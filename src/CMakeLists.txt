target_sources(
  ${main_lib}
  PRIVATE # Core modules
          core/mqc_elements.f90
          core/mqc_geometry.f90
          # Basis set modules
          basis/mqc_cgto.f90
          basis/mqc_basis_reader.f90
          # I/O modules
          io/mqc_basis_file_reader.f90
          io/mqc_xyz_reader.f90
          io/mqc_input_parser.f90
          # Parallel modules
          parallel/mqc_mpi_tags.f90
          # Chemistry modules
          chemistry/mqc_chemistry_algorithms.f90
          chemistry/mqc_mbe.f90
          chemistry/mqc_physical_fragment.f90
          # CLI modules
          cli/mqc_cli_parser.f90
          cli/mqc_driver.f90)
