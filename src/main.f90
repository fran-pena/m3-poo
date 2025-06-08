program main
use iso_Fortran_env, only: real64
use matrix_class
implicit none
type(matrix) :: a, l, u
character(260) :: filename
integer :: iu

print *, 'Programa para testear los modulos para matrices de tipo matrix.'
print*, ' '
print*, 'Nome del fichero .mtx:'
read(*,'(A)') filename
print*, trim(filename)
print*, 'Lectura de la matriz...'
call read(a, filename)
print*, 'Factorizacion LU...'
call lu(a, l, u)
print*, 'Preparacion del fichero que genera la gráfica...'
open(newunit=iu, file = 'output.m', position='rewind')
write(iu,*) 'h = figure(''visible'',''off'');'
call spy(iu, a, 'A', 1)
call spy(iu, L, 'L', 2)
call spy(iu, u, 'U', 3)
write(iu,*) 'saveas(h,''grafica.png'')'
close(iu)

contains

subroutine spy(iu, M, name, iplot)
integer, intent(in) :: iu
type(matrix), intent(in) :: M
character(*), intent(in) :: name
integer, intent(in) :: iplot
integer :: i, j

write(iu,*) 'subplot(1, 3,', iplot, ')'
write(iu,*) name//' = [', ((get(M, i, j), j = 1, size(M,2)),';', i = 1, size(M,1)), '];'
write(iu,*) 'spy(', name, ')'
write(iu,*) 'title(''', name, ''')'
end subroutine

end program