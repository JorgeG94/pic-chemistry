
program main
   !use mctc_io_math
   use pic_types, only: int32
   use mctc_env, only : wp, error_type, fatal_error
   !use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
   !   & test_failed
   use mctc_io, only : structure_type, new
   use tblite_context_type, only : context_type
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, eeq_guess
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint

   use pic_types, only: int64
   use pic_timer, only: timer_type
   implicit none
   type(error_type), allocatable :: error
   type(timer_type) :: my_timer
   type(structure_type) :: mol
   real(wp), allocatable :: xyz(:, :)
   integer, allocatable :: num(:)
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp) :: energy
   type(context_type) :: ctx
   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = 10*sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp

   num = [6, 1, 1, 1, 1]
   xyz = reshape([ &
     &  0.00000000000000_wp, -0.00000000000000_wp,  0.00000000000000_wp, &
     & -1.19220800552211_wp,  1.19220800552211_wp,  1.19220800552211_wp, &
     &  1.19220800552211_wp, -1.19220800552211_wp,  1.19220800552211_wp, &
     & -1.19220800552211_wp, -1.19220800552211_wp, -1.19220800552211_wp, &
     &  1.19220800552211_wp,  1.19220800552211_wp, -1.19220800552211_wp],&
     & [3, size(num)])

   call new(mol, num, xyz, charge=0.0_wp, uhf=0)
   call new_gfn1_calculator(calc, mol, error)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)

   energy = 0.0_wp
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)
   print *, "ENERGY IS ", energy

print *, "hello wolrd"

end program main
