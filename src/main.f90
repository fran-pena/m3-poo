program main
use iso_Fortran_env, only: real64
use matrix_class
implicit none
type(matrix) :: a, l, u
character(260) :: filename

print *, 'Programa para testear los modulos para matrices de tipo matrix.'
print*, ' '
print*, 'Nome del fichero .mtx:'
read(*,'(A)') filename
print*, trim(filename)
print*, 'Lectura de la matriz...'
call read(a, filename)
print*, 'Factorizacion LU...'
call lu(a, l, u)
print*, 'Matriz A:'
call print(a)
print*, 'Matriz L:'
call print(l)
print*, 'Matriz U:'
call print(u)
end program