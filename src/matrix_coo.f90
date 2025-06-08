module matrix_class
use iso_Fortran_env, only: real64
implicit none

type :: matrix
  private
  integer :: n, m
  integer, allocatable :: row(:), col(:)
  real(real64), allocatable :: val(:)
end type

interface read;   module procedure read_matrix;  end interface
interface size;   module procedure size_matrix;  end interface
interface set;    module procedure set_sca;      end interface
interface lu;     module procedure lu_matrix;    end interface
interface print;  module procedure print_matrix; end interface
interface get
  module procedure get_sca
  module procedure get_row
  module procedure get_col
end interface
private :: read_matrix, size_matrix, get_sca, get_col, get_row, set_sca, lu_matrix, print_matrix

contains
subroutine dealloc(this)
type(matrix) :: this
integer :: i, res
character(255) :: cad

if (allocated(this%row)) deallocate(this%row, stat = res, errmsg = cad)
if (res /= 0) error stop '(matrix/dealloc) Unable to deallocate variable row: '//trim(cad)
if (allocated(this%col)) deallocate(this%col, stat = res, errmsg = cad)
if (res /= 0) error stop '(matrix/dealloc) Unable to deallocate variable col: '//trim(cad)
if (allocated(this%val)) deallocate(this%val, stat = res, errmsg = cad)
if (res /= 0) error stop '(matrix/dealloc) Unable to deallocate variable val: '//trim(cad)
end subroutine

subroutine alloc(this, row, col, val)
type(matrix)             :: this
integer,      intent(in) :: row(:), col(:)
real(real64), intent(in) :: val(:)
character(255) ::  cad

if (size(row,1) /= size(col,1) .or. size(row,1) /= size(val,1)) &
  error stop '(matrix/alloc) The three input arrays have different size.'
if (size(row,1) <= 0) then
  write(cad, *) size(row,1)
  error stop '(matrix/alloc) The input arrays have an invalid size: '//trim(adjustl(cad))
end if
!alojamiento
if (allocated(this%row) .or. allocated(this%col) .or. allocated(this%val)) call dealloc(this)
allocate(this%row(size(row,1)), this%col(size(row,1)), this%val(size(row,1)))
!valores iniciales
this%row = row
this%col = col
this%val = val
this%n = maxval(row)
this%m = maxval(col)
end subroutine

subroutine read_matrix(this, filename)
type(matrix)             :: this
character(*), intent(in) :: filename
character(255) :: cad
integer :: iu, ios, n, m, nval, i
integer, allocatable :: row(:), col(:)
real(real64), allocatable :: val(:)

open(newunit=iu, file = filename, iostat = ios, position='rewind')
if (ios /= 0) error stop '(matrix/read_matrix) Unable to open '//trim(filename)
read(iu,*) cad
read(iu,*) this%n, this%m, nval
print*, 'Read matrix ', trim(filename), ':', this%n, 'rows', this%m, 'cols', nval, 'values'
if (allocated(this%row) .or. allocated(this%col) .or. allocated(this%val)) call dealloc(this)
allocate(this%row(nval), this%col(nval), this%val(nval))
do i = 1, nval
  read(iu,*) this%row(i), this%col(i), this%val(i)
end do
close(iu)
end subroutine

function size_matrix(this, d) result(res)
type(matrix), intent(in) :: this
integer,      intent(in) :: d
integer :: res
character(255) :: cad

select case (d)
case(1)
  res = this%n
case(2)
  res = this%m
case default
  write(cad, *) d
  error stop '(matrix/size_matrix) The requested dimension is invalid: '//trim(adjustl(cad))
end select
end function

function get_sca(this, i, j) result(val)
type(matrix), intent(in) :: this
integer,      intent(in) :: i, j
real(real64)             :: val
integer :: k
character(255) :: cad

if (.not. allocated(this%row) .or. .not. allocated(this%col) .or. .not. allocated(this%val)) &
  error stop '(matrix/get_sca) Matrix is not allocated.'
if (size(this%row,1) < 0) then 
  write(cad, *) size(this%row, 1)
  error stop '(matrix/get_sca) The matrix arrays have an invalid size: '//trim(adjustl(cad))
end if
do k = 1, size(this%row,1)
  if (this%row(k) == i .and. this%col(k) == j) then
    val = this%val(k)
    return
  end if
end do
val = 0._real64
end function

function get_col(this, i, j) result(val)
type(matrix), intent(in)  :: this
integer,      intent(in)  :: i(:), j
real(real64)              :: val(size(i,1))
integer :: k

if (size(i,1) < 0) error stop '(matrix/get_col) Insufficient number of indices.'
do k = 1, size(i,1)
  val(k) = get(this, i(k), j)
end do
end function

function get_row(this, i, j) result(val)
type(matrix), intent(in)  :: this
integer,      intent(in)  :: i, j(:)
real(real64)              :: val(size(j,1))
integer :: k

if (size(j,1) < 0) error stop '(matrix/get_row) Insufficient number of indices.'
do k = 1, size(j,1)
  val(k) = get(this, i, j(k))
end do
end function

subroutine set_sca(this, i, j, val)
type(matrix), intent(inout) :: this
integer,      intent(in)    :: i, j
real(real64), intent(in)    :: val
integer ::  k, s
integer, allocatable :: trow(:), tcol(:)
real(real64), allocatable :: tval(:)

if (.not. allocated(this%row) .or. .not. allocated(this%col) .or. .not. allocated(this%val)) then
  call alloc(this, [i], [j], [val])
else
  s = size(this%row,1)
  do k = 1, s
    if (this%row(k) == i .and. this%col(k) == j) then
      this%val(k) = val
      return
    end if
  end do
  allocate(trow(s+1), tcol(s+1), tval(s+1))
  trow(1:s) =  this%row; trow(s+1) = i
  tcol(1:s) =  this%col; tcol(s+1) = j
  tval(1:s) =  this%val; tval(s+1) = val
  call move_alloc(trow, this%row)
  call move_alloc(tcol, this%col)
  call move_alloc(tval, this%val)
  this%n = max(this%n, i)
  this%m = max(this%m, j)
end if
end subroutine

subroutine lu_matrix (a, l, u)
type(matrix), intent(in) :: a
type(matrix), intent(out) :: l, u
integer :: i, j, k

do j = 1, size(a,2)
  do i = 1, j
    call set(u, i, j, get(a, i, j) - dot_product(get(l, i, [(k, k=1,i-1)]), get(u, [(k, k=1,i-1)], j)))
  end do
  call set(l, j, j, 1._real64)
  do i = j+1, size(a,2)
    call set(l, i, j, (get(a, i, j) - dot_product(get(l, i, [(k, k=1,j-1)]), get(u, [(k, k=1,j-1)], j))) / get(u, j, j))
  end do
end do
end subroutine

subroutine print_matrix(this)
type(matrix), intent(in)  :: this
integer :: i, j
do i = 1, this%n
  print*, get(this, i, [(j, j = 1, this%m)])
end do
end subroutine
end module

