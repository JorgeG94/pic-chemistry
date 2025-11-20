module pic_fragment
   use pic_types, only: default_int, dp, int64
   use pic_io, only: to_char
   implicit none
   private
   public :: pic_fragment_block
   public :: count_nonzeros
   public :: new_pic_fragment

   type :: pic_fragment_block
      integer(default_int), allocatable :: indices(:)
      real(dp), allocatable :: matrix(:, :)
   contains
      procedure :: print_indices
   end type pic_fragment_block

   type :: pic_fragment_type
      integer(default_int), allocatable :: indices(:)
      integer(default_int) :: n_momonomers
      integer(default_int) :: n_rows
      integer(default_int) :: estimated_cost

   contains
      procedure :: get_estimated_cost
   end type

contains

   function get_estimated_cost(self) result(estimated_cost)
      class(pic_fragment_type), intent(in) :: self
      integer(default_int) :: estimated_cost

      estimated_cost = self%estimated_cost
   end function get_estimated_cost

   function new_pic_fragment(indices, n_rows) result(fragment)
      integer(default_int), intent(in) :: indices(:)
      integer(default_int), intent(in) :: n_rows
      type(pic_fragment_type) :: fragment

      fragment%n_rows = n_rows
      fragment%estimated_cost = n_rows*n_rows
      allocate (fragment%indices, source=indices)
   end function new_pic_fragment

   function count_nonzeros(row) result(count)
      integer(default_int), intent(in) :: row(:)
      integer(default_int) :: i
      integer(default_int) :: count
      count = 0
      do i = 1, size(row)
         if (row(i) /= 0) count = count + 1
      end do
   end function count_nonzeros

   subroutine print_indices(self)
      class(pic_fragment_block), intent(in) :: self

      print *, to_char(self%indices)

   end subroutine print_indices

end module pic_fragment
