#!/usr/bin/env python3
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

try:
    import numpy as np
except Exception:
    np = None  # optional


# ----------------------------
# Data containers (internal)
# ----------------------------

@dataclass
class SchemaTag:
    name: str
    version: str

@dataclass
class Model:
    method: str
    basis: Optional[str] = None
    aux_basis: Optional[str] = None

@dataclass
class SCFKeywords:
    maxiter: Optional[int] = None
    tolerance: Optional[float] = None

@dataclass
class FragmentationCutoffs:
    dimer: Optional[float] = None
    trimer: Optional[float] = None
    # easy to extend later: tetramer, etc.

@dataclass
class FragmentationKeywords:
    method: str = "MBE"
    allow_overlapping_fragments: bool = False
    level: int = 2
    embedding: str = "none"
    cutoff_method: str = "distance"
    distance_metric: str = "min"
    cutoffs: Optional[FragmentationCutoffs] = None

@dataclass
class Keywords:
    scf: Optional[SCFKeywords] = None
    fragmentation: Optional[FragmentationKeywords] = None
    raw: Optional[Dict[str, Any]] = None  # keep everything as-is too

@dataclass
class Molecule:
    symbols: List[str]                          # length N
    geometry: Any                               # Nx3 (list-of-lists or numpy array)
    connectivity: List[Tuple[int, int, int]]    # (i,j,order), 0-based
    fragments: List[List[int]]                  # 0-based atom indices
    fragment_charges: List[int]
    fragment_multiplicities: List[int]
    source_xyz: Optional[str] = None            # if read from file

@dataclass
class ProcessedInput:
    schema: SchemaTag
    molecules: List[Molecule]
    model: Model
    keywords: Keywords
    driver: str
    raw: Dict[str, Any]


# ----------------------------
# Utility / validation helpers
# ----------------------------

def _die(msg: str) -> None:
    raise ValueError(msg)

def _as_schema_tag(x: Any) -> SchemaTag:
    """
    Accept either:
      "mqc-frag-1.0"
    or:
      {"name": "mqc-frag", "version": "1.0"}
    """
    if isinstance(x, str):
        if "-" in x:
            name, version = x.rsplit("-", 1)
            return SchemaTag(name=name, version=version)
        return SchemaTag(name=x, version="unknown")

    if isinstance(x, dict):
        name = x.get("name")
        version = x.get("version")
        if not isinstance(name, str) or not isinstance(version, str):
            _die(f"schema object must contain string fields name/version; got: {x}")
        return SchemaTag(name=name, version=version)

    _die(f"schema must be string or object; got type {type(x)}")

def _check_symbols(symbols: Any) -> List[str]:
    if not isinstance(symbols, list) or not all(isinstance(s, str) for s in symbols):
        _die("molecules[*].symbols must be a list of strings")
    if len(symbols) == 0:
        _die("molecules[*].symbols must be non-empty")
    return symbols

def _reshape_geometry(geom: Any, natom: int) -> Any:
    """
    Accept flat list length 3*N, return Nx3.
    If numpy exists, returns np.ndarray (float64). Otherwise list-of-lists.
    """
    if not isinstance(geom, list):
        _die("molecules[*].geometry must be a flat list of numbers (length 3*N)")

    if len(geom) != 3 * natom:
        _die(f"geometry length mismatch: expected {3*natom}, got {len(geom)}")

    # basic numeric validation
    for k, v in enumerate(geom[:min(len(geom), 12)]):
        if not isinstance(v, (int, float)):
            _die(f"geometry contains non-numeric entry at position {k}: {v}")

    if np is not None:
        return np.asarray(geom, dtype=float).reshape((natom, 3))

    return [[float(geom[3*i+0]), float(geom[3*i+1]), float(geom[3*i+2])]
            for i in range(natom)]

def _check_index_0based(i: Any, natom: int, *, context: str) -> int:
    if not isinstance(i, int):
        _die(f"{context} must be an integer; got {type(i)}")
    if i < 0 or i >= natom:
        _die(f"{context} out of range: {i} for natom={natom} (0-based expected)")
    return i

def _check_fragments(frags: Any, natom: int) -> List[List[int]]:
    if not isinstance(frags, list):
        _die("molecules[*].fragments must be a list of fragment atom-index lists")

    out: List[List[int]] = []
    for fi, frag in enumerate(frags):
        if not isinstance(frag, list) or len(frag) == 0:
            _die(f"fragment {fi} must be a non-empty list of atom indices")

        frag_checked: List[int] = []
        for aj, a in enumerate(frag):
            frag_checked.append(_check_index_0based(a, natom, context=f"fragment[{fi}][{aj}]"))

        out.append(frag_checked)

    return out

def _check_parallel_fragment_arrays(nfrag: int, charges: Any, mults: Any) -> Tuple[List[int], List[int]]:
    if not isinstance(charges, list) or not all(isinstance(x, int) for x in charges):
        _die("molecules[*].fragment_charges must be a list of integers")
    if not isinstance(mults, list) or not all(isinstance(x, int) for x in mults):
        _die("molecules[*].fragment_multiplicities must be a list of integers")

    if len(charges) != nfrag:
        _die(f"fragment_charges length mismatch: expected {nfrag}, got {len(charges)}")
    if len(mults) != nfrag:
        _die(f"fragment_multiplicities length mismatch: expected {nfrag}, got {len(mults)}")

    for i, m in enumerate(mults):
        if m <= 0:
            _die(f"fragment multiplicity must be >= 1; fragment {i} has {m}")

    return charges, mults

def _check_connectivity(conn: Any, natom: int) -> List[Tuple[int, int, int]]:
    if conn is None:
        return []

    if not isinstance(conn, list):
        _die("molecules[*].connectivity must be a list of [i,j,order] triples")

    bonds: List[Tuple[int, int, int]] = []
    for bi, triple in enumerate(conn):
        if not isinstance(triple, list) or len(triple) != 3:
            _die(f"connectivity[{bi}] must be a list [i,j,order]")
        i, j, order = triple

        i0 = _check_index_0based(i, natom, context=f"connectivity[{bi}].i")
        j0 = _check_index_0based(j, natom, context=f"connectivity[{bi}].j")

        if not isinstance(order, int) or order <= 0:
            _die(f"connectivity[{bi}].order must be a positive int; got {order}")
        if i0 == j0:
            _die(f"connectivity[{bi}] invalid self-bond ({i0},{j0})")

        bonds.append((i0, j0, order))

    return bonds


# ----------------------------
# XYZ reader (optional molecule source)
# ----------------------------

def read_xyz_file(xyz_path: Union[str, Path]) -> Tuple[List[str], Any]:
    """
    Standard XYZ:
      line1: natom
      line2: comment
      line3+: sym x y z
    Returns (symbols, geometry Nx3).
    """
    xyz_path = Path(xyz_path)
    text = xyz_path.read_text(encoding="utf-8").strip().splitlines()
    if len(text) < 3:
        _die(f"XYZ file too short: {xyz_path}")

    try:
        natom = int(text[0].strip().split()[0])
    except Exception as e:
        _die(f"XYZ first line must start with natom int: {xyz_path} ({e})")

    if natom <= 0:
        _die(f"XYZ natom must be > 0: {xyz_path}")

    lines = text[2:]  # skip natom + comment
    if len(lines) < natom:
        _die(f"XYZ expected {natom} atom lines, got {len(lines)}: {xyz_path}")

    symbols: List[str] = []
    coords: List[float] = []

    for i in range(natom):
        parts = lines[i].split()
        if len(parts) < 4:
            _die(f"XYZ atom line {i+3} malformed in {xyz_path}: '{lines[i]}'")
        sym = parts[0]
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except Exception as e:
            _die(f"XYZ atom line {i+3} has non-numeric coords in {xyz_path}: '{lines[i]}' ({e})")
        symbols.append(sym)
        coords.extend([x, y, z])

    geom = _reshape_geometry(coords, natom)
    return symbols, geom


# ----------------------------
# Keywords parsing
# ----------------------------

def _parse_scf_keywords(d: Any) -> SCFKeywords:
    if not isinstance(d, dict):
        _die("keywords.scf must be an object")

    maxiter = d.get("maxiter")
    tol = d.get("tolerance")

    if maxiter is not None and (not isinstance(maxiter, int) or maxiter <= 0):
        _die(f"keywords.scf.maxiter must be positive int; got {maxiter}")
    if tol is not None and (not isinstance(tol, (int, float)) or tol <= 0.0):
        _die(f"keywords.scf.tolerance must be positive number; got {tol}")

    return SCFKeywords(maxiter=maxiter, tolerance=float(tol) if tol is not None else None)

def _parse_fragmentation_keywords(d: Any) -> FragmentationKeywords:
    if not isinstance(d, dict):
        _die("keywords.fragmentation must be an object")

    method = d.get("method", "MBE")
    if not isinstance(method, str) or not method:
        _die("keywords.fragmentation.method must be a non-empty string")

    allow_ov = d.get("allow_overlapping_fragments", False)
    if not isinstance(allow_ov, bool):
        _die("keywords.fragmentation.allow_overlapping_fragments must be boolean")

    level = d.get("level", 2)
    if not isinstance(level, int) or level <= 0:
        _die("keywords.fragmentation.level must be a positive int")

    embedding = d.get("embedding", "none")
    if not isinstance(embedding, str):
        _die("keywords.fragmentation.embedding must be a string")

    cutoff_method = d.get("cutoff_method", "distance")
    if not isinstance(cutoff_method, str):
        _die("keywords.fragmentation.cutoff_method must be a string")

    distance_metric = d.get("distance_metric", "min")
    if not isinstance(distance_metric, str):
        _die("keywords.fragmentation.distance_metric must be a string")

    cutoffs_obj = d.get("cutoffs")
    cutoffs = None
    if cutoffs_obj is not None:
        if not isinstance(cutoffs_obj, dict):
            _die("keywords.fragmentation.cutoffs must be an object")
        dimer = cutoffs_obj.get("dimer")
        trimer = cutoffs_obj.get("trimer")

        def _pos_float(x: Any, field: str) -> Optional[float]:
            if x is None:
                return None
            if not isinstance(x, (int, float)) or x <= 0.0:
                _die(f"keywords.fragmentation.cutoffs.{field} must be positive number; got {x}")
            return float(x)

        cutoffs = FragmentationCutoffs(
            dimer=_pos_float(dimer, "dimer"),
            trimer=_pos_float(trimer, "trimer"),
        )

    return FragmentationKeywords(
        method=method,
        allow_overlapping_fragments=allow_ov,
        level=level,
        embedding=embedding,
        cutoff_method=cutoff_method,
        distance_metric=distance_metric,
        cutoffs=cutoffs,
    )


# ----------------------------
# Main reader
# ----------------------------

def load_mqc_json(path: Union[str, Path]) -> ProcessedInput:
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        _die("Top-level JSON must be an object")

    schema = _as_schema_tag(data.get("schema", "mqc-frag-unknown"))

    # model
    model_d = data.get("model", {})
    if not isinstance(model_d, dict):
        _die("model must be an object")
    method = model_d.get("method")
    if not isinstance(method, str) or not method:
        _die("model.method must be a non-empty string")

    model = Model(
        method=method,
        basis=model_d.get("basis"),
        aux_basis=model_d.get("aux_basis"),
    )

    # keywords (scf + fragmentation at same level)
    kw_d = data.get("keywords", {}) or {}
    if not isinstance(kw_d, dict):
        _die("keywords must be an object")

    scf_kw = _parse_scf_keywords(kw_d["scf"]) if "scf" in kw_d else None
    frag_kw = _parse_fragmentation_keywords(kw_d["fragmentation"]) if "fragmentation" in kw_d else None
    keywords = Keywords(scf=scf_kw, fragmentation=frag_kw, raw=kw_d)

    # driver
    driver = data.get("driver", "Energy")
    if not isinstance(driver, str):
        _die("driver must be a string")

    # molecules
    mols = data.get("molecules")
    if not isinstance(mols, list) or len(mols) == 0:
        _die("molecules must be a non-empty list")

    molecules: List[Molecule] = []
    for mi, m in enumerate(mols):
        if not isinstance(m, dict):
            _die(f"molecules[{mi}] must be an object")

        # Two ways to specify structure:
        #   A) symbols + geometry
        #   B) xyz: "file.xyz"  (if coords/elements not specified)
        xyz_ref = m.get("xyz", None)

        has_symbols = "symbols" in m
        has_geom = "geometry" in m

        if xyz_ref is not None:
            if has_symbols or has_geom:
                _die(f"molecules[{mi}] provides xyz AND symbols/geometry. Choose one.")
            if not isinstance(xyz_ref, str) or not xyz_ref:
                _die(f"molecules[{mi}].xyz must be a non-empty string path")

            xyz_path = (path.parent / xyz_ref).resolve()
            symbols, geom = read_xyz_file(xyz_path)
            source_xyz = str(xyz_ref)
        else:
            if not (has_symbols and has_geom):
                _die(f"molecules[{mi}] must provide either (symbols+geometry) OR xyz")
            symbols = _check_symbols(m.get("symbols"))
            natom = len(symbols)
            geom = _reshape_geometry(m.get("geometry"), natom)
            source_xyz = None

        natom = len(symbols)

        # connectivity + fragments use 0-based indexing
        conn = _check_connectivity(m.get("connectivity", []), natom)
        frags = _check_fragments(m.get("fragments", []), natom)
        nfrag = len(frags)

        charges, mults = _check_parallel_fragment_arrays(
            nfrag,
            m.get("fragment_charges", [0] * nfrag),
            m.get("fragment_multiplicities", [1] * nfrag),
        )

        molecules.append(Molecule(
            symbols=symbols,
            geometry=geom,
            connectivity=conn,
            fragments=frags,
            fragment_charges=charges,
            fragment_multiplicities=mults,
            source_xyz=source_xyz,
        ))

    return ProcessedInput(
        schema=schema,
        molecules=molecules,
        model=model,
        keywords=keywords,
        driver=driver,
        raw=data,
    )


# ----------------------------
# Optional summary / debug dump
# ----------------------------

def summarize(inp: ProcessedInput) -> None:
    print(f"Schema: {inp.schema.name}-{inp.schema.version}")
    print(f"Driver: {inp.driver}")
    print(f"Model : method={inp.model.method} basis={inp.model.basis} aux={inp.model.aux_basis}")

    if inp.keywords.scf:
        print(f"SCF   : maxiter={inp.keywords.scf.maxiter} tol={inp.keywords.scf.tolerance}")

    if inp.keywords.fragmentation:
        fk = inp.keywords.fragmentation
        print("FRAG  : "
              f"method={fk.method} level={fk.level} overlap={fk.allow_overlapping_fragments} "
              f"embedding={fk.embedding} cutoff_method={fk.cutoff_method} metric={fk.distance_metric}")
        if fk.cutoffs:
            print(f"        cutoffs: dimer={fk.cutoffs.dimer} trimer={fk.cutoffs.trimer}")

    for i, mol in enumerate(inp.molecules):
        nat = len(mol.symbols)
        nfrag = len(mol.fragments)
        nbond = len(mol.connectivity)
        src = f" (xyz='{mol.source_xyz}')" if mol.source_xyz else ""
        print(f"Molecule[{i}]: natom={nat} nfrag={nfrag} nbond={nbond}{src}")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Read/validate MQC fragmentation JSON (0-based indices).")
    p.add_argument("json_file")
    args = p.parse_args()
    inp = load_mqc_json(args.json_file)
    summarize(inp)
