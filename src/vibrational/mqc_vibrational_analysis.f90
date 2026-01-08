!! Vibrational frequency analysis from Hessian matrix
module mqc_vibrational_analysis
   !! Computes vibrational frequencies from the mass-weighted Hessian matrix.
   !! Uses LAPACK eigenvalue decomposition via pic-blas interfaces.
   use pic_types, only: dp
   use pic_lapack_interfaces, only: pic_syev, pic_gesvd
   use pic_logger, only:
   use mqc_elements, only: element_mass
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: compute_vibrational_frequencies
   public :: mass_weight_hessian
   public :: project_translation_rotation

   ! Conversion factor from atomic units (Hartree/Bohr²/amu) to cm⁻¹
   ! Derived from fundamental constants:
   !   sqrt(Hartree/(Bohr²·amu)) → s⁻¹ → cm⁻¹
   real(dp), parameter :: AU_TO_CM1 = 2.642461e7_dp

contains

   subroutine compute_vibrational_frequencies(hessian, element_numbers, frequencies, &
                                              eigenvalues_out, eigenvectors, &
                                              coordinates, project_trans_rot, &
                                              projected_hessian_out)
      !! Compute vibrational frequencies from the Hessian matrix.
      !!
      !! Algorithm:
      !! 1. Mass-weight the Hessian: H_mw = M^{-1/2} * H * M^{-1/2}
      !! 2. Optionally project out translation/rotation modes
      !! 3. Diagonalize H_mw to get eigenvalues
      !! 4. Convert eigenvalues to frequencies in cm⁻¹
      !!
      !! Negative eigenvalues produce negative frequencies (imaginary modes,
      !! indicating transition states or saddle points).
      real(dp), intent(in) :: hessian(:, :)
         !! Hessian matrix in Hartree/Bohr² (3*N x 3*N)
      integer, intent(in) :: element_numbers(:)
         !! Atomic numbers for each atom (N atoms)
      real(dp), allocatable, intent(out) :: frequencies(:)
         !! Vibrational frequencies in cm⁻¹ (3*N modes, or 3*N-6 if projected)
      real(dp), allocatable, intent(out), optional :: eigenvalues_out(:)
         !! Raw eigenvalues from diagonalization (Hartree/Bohr²/amu)
      real(dp), allocatable, intent(out), optional :: eigenvectors(:, :)
         !! Normal mode eigenvectors (3*N x 3*N), columns are modes
      real(dp), intent(in), optional :: coordinates(:, :)
         !! Atomic coordinates in Bohr (3, N) - required for projection
      logical, intent(in), optional :: project_trans_rot
         !! If true, project out translation/rotation modes (requires coordinates)
      real(dp), allocatable, intent(out), optional :: projected_hessian_out(:, :)
         !! Mass-weighted Hessian after trans/rot projection (before diagonalization)

      real(dp), allocatable :: mw_hessian(:, :)
      real(dp), allocatable :: eigenvalues(:)
      integer :: n_coords, info, i
      logical :: do_projection

      n_coords = size(hessian, 1)

      ! Check if projection is requested
      do_projection = .false.
      if (present(project_trans_rot)) then
         if (project_trans_rot) then
            if (.not. present(coordinates)) then
               ! Cannot project without coordinates - fall back to no projection
               call logger%warning("Missing coordinates, not projecting out tran/rot motions")
               do_projection = .false.
            else
               do_projection = .true.
            end if
         end if
      end if

      ! Mass-weight the Hessian
      call mass_weight_hessian(hessian, element_numbers, mw_hessian)

      ! Optionally project out translation/rotation modes
      if (do_projection) then
         call project_translation_rotation(mw_hessian, coordinates, element_numbers)
      end if

      ! Return projected Hessian if requested (before diagonalization destroys it)
      if (present(projected_hessian_out)) then
         allocate (projected_hessian_out(n_coords, n_coords))
         projected_hessian_out = mw_hessian
      end if

      ! Allocate eigenvalue storage
      allocate (eigenvalues(n_coords))

      ! Diagonalize the mass-weighted Hessian
      ! pic_syev overwrites mw_hessian with eigenvectors (if jobz='V', default)
      call pic_syev(mw_hessian, eigenvalues, info=info)

      if (info /= 0) then
         ! Eigenvalue decomposition failed
         call logger%error("Eigenvalue decomposition in vibrational frequencies failed")
         allocate (frequencies(n_coords))
         frequencies = 0.0_dp
         return
      end if

      ! Convert eigenvalues to frequencies in cm⁻¹
      allocate (frequencies(n_coords))
      do i = 1, n_coords
         if (eigenvalues(i) >= 0.0_dp) then
            ! Real frequency
            frequencies(i) = sqrt(eigenvalues(i)*AU_TO_CM1)
         else
            ! Imaginary frequency (negative eigenvalue) - report as negative
            frequencies(i) = -sqrt(abs(eigenvalues(i))*AU_TO_CM1)
         end if
      end do

      ! Return eigenvalues if requested
      if (present(eigenvalues_out)) then
         allocate (eigenvalues_out(n_coords))
         eigenvalues_out = eigenvalues
      end if

      ! Return eigenvectors if requested
      if (present(eigenvectors)) then
         allocate (eigenvectors(n_coords, n_coords))
         eigenvectors = mw_hessian
      end if

      deallocate (eigenvalues, mw_hessian)

   end subroutine compute_vibrational_frequencies

   subroutine mass_weight_hessian(hessian, element_numbers, mw_hessian)
      !! Apply mass weighting to Hessian matrix.
      !!
      !! H_mw(i,j) = H(i,j) / sqrt(m_i * m_j)
      !!
      !! where m_i is the mass of the atom corresponding to coordinate i.
      !! Each atom contributes 3 coordinates (x, y, z).
      real(dp), intent(in) :: hessian(:, :)
         !! Input Hessian in Hartree/Bohr² (3*N x 3*N)
      integer, intent(in) :: element_numbers(:)
         !! Atomic numbers for each atom (N atoms)
      real(dp), allocatable, intent(out) :: mw_hessian(:, :)
         !! Mass-weighted Hessian (3*N x 3*N)

      real(dp), allocatable :: inv_sqrt_mass(:)
      integer :: n_atoms, n_coords, iatom, icoord, i, j
      real(dp) :: mass

      n_atoms = size(element_numbers)
      n_coords = 3*n_atoms

      ! Build inverse square root mass vector (each mass repeated 3x for x,y,z)
      allocate (inv_sqrt_mass(n_coords))
      do iatom = 1, n_atoms
         mass = element_mass(element_numbers(iatom))
         do icoord = 1, 3
            inv_sqrt_mass(3*(iatom - 1) + icoord) = 1.0_dp/sqrt(mass)
         end do
      end do

      ! Apply mass weighting: H_mw(i,j) = H(i,j) * inv_sqrt_mass(i) * inv_sqrt_mass(j)
      allocate (mw_hessian(n_coords, n_coords))
      do j = 1, n_coords
         do i = 1, n_coords
            mw_hessian(i, j) = hessian(i, j)*inv_sqrt_mass(i)*inv_sqrt_mass(j)
         end do
      end do

      deallocate (inv_sqrt_mass)

   end subroutine mass_weight_hessian

   subroutine project_translation_rotation(mw_hessian, coordinates, element_numbers)
      !! Project out translation and rotation modes from mass-weighted Hessian.
      !!
      !! Builds 6 vectors (3 translation + 3 rotation) in mass-weighted coordinates,
      !! orthonormalizes them using SVD, then projects them out:
      !!   H_proj = (I - D @ D^T) @ H @ (I - D @ D^T)
      !!
      !! This sets the 6 translation/rotation eigenvalues to exactly zero.
      real(dp), intent(inout) :: mw_hessian(:, :)
         !! Mass-weighted Hessian (modified in place)
      real(dp), intent(in) :: coordinates(:, :)
         !! Atomic coordinates in Bohr (3, N)
      integer, intent(in) :: element_numbers(:)
         !! Atomic numbers for each atom (N atoms)

      real(dp), allocatable :: D(:, :)        ! Translation/rotation vectors (3N, 6)
      real(dp), allocatable :: com(:)         ! Center of mass
      real(dp), allocatable :: r(:, :)        ! Coordinates relative to COM
      real(dp), allocatable :: sqrt_mass(:)   ! sqrt(mass) for each atom
      real(dp), allocatable :: S(:)           ! Singular values
      real(dp), allocatable :: U(:, :)        ! Left singular vectors
      real(dp), allocatable :: VT(:, :)       ! Right singular vectors (transposed)
      real(dp), allocatable :: D_orth(:, :)   ! Orthonormalized D vectors
      real(dp), allocatable :: proj(:, :)     ! Projector matrix
      real(dp), allocatable :: temp(:, :)     ! Temporary matrix
      real(dp) :: total_mass, mass, norm
      integer :: n_atoms, n_coords, iatom, i, j, k, n_modes, info
      integer :: idx

      n_atoms = size(element_numbers)
      n_coords = 3*n_atoms

      ! Allocate arrays
      allocate (D(n_coords, 6))
      allocate (com(3))
      allocate (r(3, n_atoms))
      allocate (sqrt_mass(n_atoms))

      ! Compute sqrt(mass) for each atom and total mass
      total_mass = 0.0_dp
      do iatom = 1, n_atoms
         mass = element_mass(element_numbers(iatom))
         sqrt_mass(iatom) = sqrt(mass)
         total_mass = total_mass + mass
      end do

      ! Compute center of mass
      com = 0.0_dp
      do iatom = 1, n_atoms
         mass = element_mass(element_numbers(iatom))
         com(:) = com(:) + mass*coordinates(:, iatom)
      end do
      com = com/total_mass

      ! Compute coordinates relative to center of mass
      do iatom = 1, n_atoms
         r(:, iatom) = coordinates(:, iatom) - com(:)
      end do

      ! Initialize D to zero
      D = 0.0_dp

      ! Build translation vectors (mass-weighted)
      ! D_trans_k: displacement along axis k, weighted by sqrt(mass)
      do iatom = 1, n_atoms
         idx = 3*(iatom - 1)
         ! Translation along x
         D(idx + 1, 1) = sqrt_mass(iatom)
         ! Translation along y
         D(idx + 2, 2) = sqrt_mass(iatom)
         ! Translation along z
         D(idx + 3, 3) = sqrt_mass(iatom)
      end do

      ! Build rotation vectors (mass-weighted)
      ! D_rot_k: rotation around axis k, proportional to r × e_k, weighted by sqrt(mass)
      do iatom = 1, n_atoms
         idx = 3*(iatom - 1)
         ! Rotation around x-axis: r × e_x = (0, r_z, -r_y)
         D(idx + 2, 4) = sqrt_mass(iatom)*r(3, iatom)
         D(idx + 3, 4) = -sqrt_mass(iatom)*r(2, iatom)
         ! Rotation around y-axis: r × e_y = (-r_z, 0, r_x)
         D(idx + 1, 5) = -sqrt_mass(iatom)*r(3, iatom)
         D(idx + 3, 5) = sqrt_mass(iatom)*r(1, iatom)
         ! Rotation around z-axis: r × e_z = (r_y, -r_x, 0)
         D(idx + 1, 6) = sqrt_mass(iatom)*r(2, iatom)
         D(idx + 2, 6) = -sqrt_mass(iatom)*r(1, iatom)
      end do

      ! Normalize each column of D
      do k = 1, 6
         norm = sqrt(sum(D(:, k)**2))
         if (norm > 1.0e-10_dp) then
            D(:, k) = D(:, k)/norm
         end if
      end do

      ! Orthonormalize D using SVD: D = U @ S @ VT
      ! The orthonormal basis is given by the columns of U corresponding to non-zero singular values
      allocate (S(6))
      allocate (U(n_coords, 6))
      allocate (VT(6, 6))

      ! pic_gesvd(A, S, U, VT, info) - A is input, U and VT are separate outputs
      call pic_gesvd(D, S, U, VT, info=info)

      ! Count non-zero singular values (determines number of modes to project)
      n_modes = 0
      do k = 1, 6
         if (S(k) > 1.0e-10_dp) n_modes = n_modes + 1
      end do

      ! Build orthonormalized D matrix from U (columns with non-zero singular values)
      allocate (D_orth(n_coords, n_modes))
      j = 0
      do k = 1, 6
         if (S(k) > 1.0e-10_dp) then
            j = j + 1
            D_orth(:, j) = U(:, k)
         end if
      end do

      ! Build projector: P = I - D_orth @ D_orth^T
      allocate (proj(n_coords, n_coords))
      proj = 0.0_dp
      do i = 1, n_coords
         proj(i, i) = 1.0_dp
      end do

      ! Subtract D_orth @ D_orth^T
      do i = 1, n_coords
         do j = 1, n_coords
            do k = 1, n_modes
               proj(i, j) = proj(i, j) - D_orth(i, k)*D_orth(j, k)
            end do
         end do
      end do

      ! Apply projection: H_proj = P @ H @ P
      allocate (temp(n_coords, n_coords))

      ! temp = H @ P
      do i = 1, n_coords
         do j = 1, n_coords
            temp(i, j) = 0.0_dp
            do k = 1, n_coords
               temp(i, j) = temp(i, j) + mw_hessian(i, k)*proj(k, j)
            end do
         end do
      end do

      ! H_proj = P @ temp
      do i = 1, n_coords
         do j = 1, n_coords
            mw_hessian(i, j) = 0.0_dp
            do k = 1, n_coords
               mw_hessian(i, j) = mw_hessian(i, j) + proj(i, k)*temp(k, j)
            end do
         end do
      end do

      ! Cleanup
      deallocate (D, com, r, sqrt_mass, S, U, VT, D_orth, proj, temp)

   end subroutine project_translation_rotation

end module mqc_vibrational_analysis
