!! Adapter module to convert between new mqc_config_t and existing structures
!! This allows using the new .mqc file format with the existing backend
module mqc_config_adapter
   !! Provides conversion utilities between new mqc_config_t and legacy structures
   !! to support new input format without modifying existing backend code
   use pic_types, only: dp, int32
   use mqc_config_parser, only: mqc_config_t
   use mqc_input_parser, only: input_config_t
   use mqc_physical_fragment, only: system_geometry_t, to_bohr
   use mqc_elements, only: element_symbol_to_number
   implicit none
   private

   public :: config_to_legacy, geometry_to_system
   public :: config_to_system_geometry

contains

   subroutine config_to_legacy(mqc_config, legacy_config)
      !! Convert new mqc_config_t to legacy input_config_t
      !! Maps only the fields needed by the driver
      type(mqc_config_t), intent(in) :: mqc_config
      type(input_config_t), intent(out) :: legacy_config

      ! Copy method and calc_type (already integers)
      legacy_config%method = mqc_config%method
      legacy_config%calc_type = mqc_config%calc_type

      ! Set fragmentation level
      ! For unfragmented calculations (nfrag=0), nlevel must be 0
      if (mqc_config%nfrag == 0) then
         legacy_config%nlevel = 0
      else
         legacy_config%nlevel = mqc_config%frag_level
      end if

      ! Copy log level
      if (allocated(mqc_config%log_level)) then
         legacy_config%log_level = mqc_config%log_level
      else
         legacy_config%log_level = "info"
      end if

      ! Note: geom_file and monomer_file are not set
      ! because geometry comes from mqc_config directly

   end subroutine config_to_legacy

   subroutine config_to_system_geometry(mqc_config, sys_geom, stat, errmsg)
      !! Convert mqc_config_t geometry to system_geometry_t
      !! For unfragmented calculations (nfrag=0), treats entire system as single unit
      !! For fragmented calculations, currently assumes monomer-based fragmentation
      type(mqc_config_t), intent(in) :: mqc_config
      type(system_geometry_t), intent(out) :: sys_geom
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      integer :: i
      logical :: use_angstrom

      stat = 0

      ! Check if geometry is loaded
      if (mqc_config%geometry%natoms == 0) then
         stat = 1
         errmsg = "No geometry loaded in mqc_config"
         return
      end if

      ! Determine units
      use_angstrom = .true.
      if (allocated(mqc_config%units)) then
         if (trim(mqc_config%units) == 'bohr') then
            use_angstrom = .false.
         end if
      end if

      if (mqc_config%nfrag == 0) then
         ! Unfragmented calculation: entire system is one "monomer"
         call geometry_to_system_unfragmented(mqc_config%geometry, sys_geom, use_angstrom)
      else
         ! Fragmented calculation with explicit fragments
         ! For now, we'll try to infer monomer size from first fragment
         ! TODO: This is a temporary solution - need better fragmentation support
         call geometry_to_system_fragmented(mqc_config, sys_geom, use_angstrom, stat, errmsg)
         if (stat /= 0) return
      end if

   end subroutine config_to_system_geometry

   subroutine geometry_to_system_unfragmented(geom, sys_geom, use_angstrom)
      !! Convert geometry to system_geometry_t for unfragmented calculation
      !! Treats entire system as a single monomer
      use mqc_geometry, only: geometry_type

      type(geometry_type), intent(in) :: geom
      type(system_geometry_t), intent(out) :: sys_geom
      logical, intent(in) :: use_angstrom

      integer :: i

      ! For unfragmented: n_monomers=1, atoms_per_monomer=natoms
      sys_geom%n_monomers = 1
      sys_geom%atoms_per_monomer = geom%natoms
      sys_geom%total_atoms = geom%natoms

      allocate (sys_geom%element_numbers(sys_geom%total_atoms))
      allocate (sys_geom%coordinates(3, sys_geom%total_atoms))

      ! Convert element symbols to atomic numbers
      do i = 1, sys_geom%total_atoms
         sys_geom%element_numbers(i) = element_symbol_to_number(geom%elements(i))
      end do

      ! Store coordinates (convert to Bohr if needed)
      if (use_angstrom) then
         sys_geom%coordinates = to_bohr(geom%coords)
      else
         sys_geom%coordinates = geom%coords
      end if

   end subroutine geometry_to_system_unfragmented

   subroutine geometry_to_system_fragmented(mqc_config, sys_geom, use_angstrom, stat, errmsg)
      !! Convert geometry to system_geometry_t for fragmented calculation
      !! Tries to infer monomer structure from fragments
      type(mqc_config_t), intent(in) :: mqc_config
      type(system_geometry_t), intent(out) :: sys_geom
      logical, intent(in) :: use_angstrom
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      integer :: i, atoms_in_first_frag
      integer, allocatable :: fragment_sizes(:)
      logical :: all_same_size

      stat = 0

      ! Get fragment sizes
      allocate (fragment_sizes(mqc_config%nfrag))
      do i = 1, mqc_config%nfrag
         fragment_sizes(i) = size(mqc_config%fragments(i)%indices)
      end do

      ! Check if all fragments have the same size (monomer-based)
      atoms_in_first_frag = fragment_sizes(1)
      all_same_size = all(fragment_sizes == atoms_in_first_frag)

      if (.not. all_same_size) then
         stat = 1
         errmsg = "Fragmented calculations currently only support identical fragment sizes"
         deallocate (fragment_sizes)
         return
      end if

      deallocate (fragment_sizes)

      ! Set up system geometry based on identical monomers
      sys_geom%n_monomers = mqc_config%nfrag
      sys_geom%atoms_per_monomer = atoms_in_first_frag
      sys_geom%total_atoms = mqc_config%geometry%natoms

      ! Validate that total atoms matches
      if (sys_geom%n_monomers*sys_geom%atoms_per_monomer /= sys_geom%total_atoms) then
         stat = 1
         errmsg = "Fragment atoms don't match total atoms"
         return
      end if

      allocate (sys_geom%element_numbers(sys_geom%total_atoms))
      allocate (sys_geom%coordinates(3, sys_geom%total_atoms))

      ! Convert element symbols to atomic numbers
      do i = 1, sys_geom%total_atoms
         sys_geom%element_numbers(i) = element_symbol_to_number(mqc_config%geometry%elements(i))
      end do

      ! Store coordinates (convert to Bohr if needed)
      if (use_angstrom) then
         sys_geom%coordinates = to_bohr(mqc_config%geometry%coords)
      else
         sys_geom%coordinates = mqc_config%geometry%coords
      end if

   end subroutine geometry_to_system_fragmented

   subroutine geometry_to_system(geom, sys_geom, use_angstrom)
      !! Simple conversion for backward compatibility
      !! Deprecated: use config_to_system_geometry instead
      use mqc_geometry, only: geometry_type

      type(geometry_type), intent(in) :: geom
      type(system_geometry_t), intent(out) :: sys_geom
      logical, intent(in) :: use_angstrom

      call geometry_to_system_unfragmented(geom, sys_geom, use_angstrom)

   end subroutine geometry_to_system

end module mqc_config_adapter
