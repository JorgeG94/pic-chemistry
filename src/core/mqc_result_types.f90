!! Quantum chemistry calculation result containers
module mqc_result_types
   !! Defines data structures for storing and managing results from
   !! quantum chemistry calculations including energies, gradients, and properties.
   use pic_types, only: dp
   use pic_mpi_lib, only: comm_t, isend, irecv, send, recv, wait, request_t, MPI_Status
   implicit none
   private

   public :: energy_t              !! Energy components type
   public :: calculation_result_t  !! Main result container type
   public :: result_send, result_isend  !! Send result over MPI
   public :: result_recv, result_irecv  !! Receive result over MPI

   type :: energy_t
      !! Container for quantum chemistry energy components
      !!
      !! Stores energy contributions from different levels of theory.
      !! Total energy is computed as: scf + mp2_correction + ccsd_correction + ccsd_t_correction
      real(dp) :: scf = 0.0_dp               !! SCF/HF reference energy (Hartree)
      real(dp) :: mp2_correction = 0.0_dp    !! MP2 correlation correction (Hartree)
      real(dp) :: ccsd_correction = 0.0_dp   !! CCSD correlation correction (Hartree)
      real(dp) :: ccsd_t_correction = 0.0_dp !! CCSD(T) triples correction (Hartree)
   contains
      procedure :: total => energy_total     !! Compute total energy from components
      procedure :: reset => energy_reset     !! Reset all components to zero
   end type energy_t

   type :: calculation_result_t
      !! Container for quantum chemistry calculation results
      !!
      !! Stores computed quantities from QC calculations with flags
      !! indicating which properties have been computed.
      type(energy_t) :: energy                  !! Energy components (Hartree)
      real(dp), allocatable :: gradient(:, :)   !! Energy gradient (3, natoms) (Hartree/Bohr)
      real(dp), allocatable :: hessian(:, :)    !! Energy hessian (future implementation)
      real(dp), allocatable :: dipole(:)        !! Dipole moment vector (3) (Debye)

      ! Computation status flags
      logical :: has_energy = .false.    !! Energy has been computed
      logical :: has_gradient = .false.  !! Gradient has been computed
      logical :: has_hessian = .false.   !! Hessian has been computed
      logical :: has_dipole = .false.    !! Dipole moment has been computed
   contains
      procedure :: destroy => result_destroy  !! Clean up allocated memory
      procedure :: reset => result_reset      !! Reset all values and flags
   end type calculation_result_t

contains

   pure function energy_total(this) result(total)
      !! Compute total energy from all components
      class(energy_t), intent(in) :: this
      real(dp) :: total

      total = this%scf + this%mp2_correction + this%ccsd_correction + this%ccsd_t_correction
   end function energy_total

   subroutine energy_reset(this)
      !! Reset all energy components to zero
      class(energy_t), intent(inout) :: this
      this%scf = 0.0_dp
      this%mp2_correction = 0.0_dp
      this%ccsd_correction = 0.0_dp
      this%ccsd_t_correction = 0.0_dp
   end subroutine energy_reset

   subroutine result_destroy(this)
      !! Clean up allocated memory in calculation_result_t
      class(calculation_result_t), intent(inout) :: this
      if (allocated(this%gradient)) deallocate (this%gradient)
      if (allocated(this%hessian)) deallocate (this%hessian)
      if (allocated(this%dipole)) deallocate (this%dipole)
      call this%reset()
   end subroutine result_destroy

   subroutine result_reset(this)
      !! Reset all values and flags in calculation_result_t
      class(calculation_result_t), intent(inout) :: this
      call this%energy%reset()
      this%has_energy = .false.
      this%has_gradient = .false.
      this%has_hessian = .false.
      this%has_dipole = .false.
   end subroutine result_reset

   subroutine result_send(result, comm, dest, tag)
      !! Send calculation result over MPI (blocking)
      !! Sends energy components and conditionally sends gradient based on has_gradient flag
      type(calculation_result_t), intent(in) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: dest, tag

      ! Send energy components
      call send(comm, result%energy%scf, dest, tag)
      call send(comm, result%energy%mp2_correction, dest, tag)
      call send(comm, result%energy%ccsd_correction, dest, tag)
      call send(comm, result%energy%ccsd_t_correction, dest, tag)

      ! Send gradient flag and data if present
      call send(comm, result%has_gradient, dest, tag)
      if (result%has_gradient) then
         call send(comm, result%gradient, dest, tag)
      end if
   end subroutine result_send

   subroutine result_isend(result, comm, dest, tag, req)
      !! Send calculation result over MPI (non-blocking)
      !! Sends SCF energy (non-blocking) and other components (blocking)
      type(calculation_result_t), intent(in) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: dest, tag
      type(request_t), intent(out) :: req

      ! Send SCF energy (non-blocking)
      call isend(comm, result%energy%scf, dest, tag, req)

      ! Send other energy components (blocking to avoid needing multiple request handles)
      call send(comm, result%energy%mp2_correction, dest, tag)
      call send(comm, result%energy%ccsd_correction, dest, tag)
      call send(comm, result%energy%ccsd_t_correction, dest, tag)

      ! Send gradient flag and data (blocking to avoid needing multiple request handles)
      call send(comm, result%has_gradient, dest, tag)
      if (result%has_gradient) then
         call send(comm, result%gradient, dest, tag)
      end if
   end subroutine result_isend

   subroutine result_recv(result, comm, source, tag, status)
      !! Receive calculation result over MPI (blocking)
      !! Receives energy components and conditionally receives gradient based on flag
      type(calculation_result_t), intent(inout) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: source, tag
      type(MPI_Status), intent(out) :: status

      ! Receive energy components
      call recv(comm, result%energy%scf, source, tag, status)
      call recv(comm, result%energy%mp2_correction, source, tag, status)
      call recv(comm, result%energy%ccsd_correction, source, tag, status)
      call recv(comm, result%energy%ccsd_t_correction, source, tag, status)
      result%has_energy = .true.

      ! Receive gradient flag and data if present
      call recv(comm, result%has_gradient, source, tag, status)
      if (result%has_gradient) then
         ! Receive allocatable gradient array (MPI lib handles allocation)
         call recv(comm, result%gradient, source, tag, status)
      end if
   end subroutine result_recv

   subroutine result_irecv(result, comm, source, tag, req)
      !! Receive calculation result over MPI (non-blocking)
      !! Receives SCF energy (non-blocking) and other components (blocking)
      type(calculation_result_t), intent(inout) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: source, tag
      type(request_t), intent(out) :: req
      type(MPI_Status) :: status

      ! Receive SCF energy (non-blocking)
      call irecv(comm, result%energy%scf, source, tag, req)

      ! Receive other energy components (blocking to avoid needing multiple request handles)
      call recv(comm, result%energy%mp2_correction, source, tag, status)
      call recv(comm, result%energy%ccsd_correction, source, tag, status)
      call recv(comm, result%energy%ccsd_t_correction, source, tag, status)
      result%has_energy = .true.

      ! Receive gradient flag and data (blocking to avoid needing multiple request handles)
      call recv(comm, result%has_gradient, source, tag, status)
      if (result%has_gradient) then
         ! Receive allocatable gradient array (MPI lib handles allocation)
         call recv(comm, result%gradient, source, tag, status)
      end if
   end subroutine result_irecv

end module mqc_result_types
