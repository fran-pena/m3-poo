# Un ejemplo de programación orientada a objetos: matrices dispersas

## Formatos de matrices dispersas

Las matrices `sparse` de Matlab se guardan en formato _compressed sparse row_ (CSR) [^1], mediante tres arreglos unidimensionales:
- los valores no nulos, 
- las longitudes de las filas y 
- los índices de las columnas.

Los valores se pueden indicar en formato _coordinate list_ (COO) [^2],[^3], es decir, indicando _(fila, columna, valor)_.

El formato _Matrix Market_ [^4] fue desarrollado (.mtx) por el National Institute of Standards and Technology (NIST) para facilitar el intercambio de matrices de prueba para computación numérica, especialmente en álgebra lineal dispersa. Es un formato basado en COO:
1. La primera línea comienza con `%%MatrixMarket` seguido de varios campos :

```text
%%MatrixMarket matrix <format> <field> <symmetry>
```
donde:
- `<format>` = coordinate (dispersa) o array (densa)
- `<field>` = real, complex, integer o pattern
- `<symmetry>` = general, symmetric, skew-symmetric, hermitian

Nuestro ejemplo estará adaptado al caso `%%MatrixMarket matrix coordinate real symmetric`. En particular usaremos el fichero [bcsstk04.mtx](https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc1/bcsstk04.html) indicado en un analisis estructural ejecutado por John Lewis (Boeing Computer Services) [^5].

Para uso intensivo de _Matrix Market_, se recomienda usar la libreria proporcionada en
[https://math.nist.gov/MatrixMarket/mmio/f/mmio.f](https://math.nist.gov/MatrixMarket/mmio/f/mmio.f).

2. La siguiente línea contiene el número de filas, columnas y (si es dispersa) número de entradas no nulas.

3. Las siguientes filas representan una entrada en formato COO: fila columna valor

Por ejemplo:

```text
%%MatrixMarket matrix coordinate real symmetric
132 132 1890
1 1  1.9960268182200e+03
2 1  5.6751562001600e+02
3 1  7.7541907594400e+02
...
```

Para usar otros formatos de matriz, se recomienda revisar [^6].

**OBSERVACIÓN:**
- Ni el formato COO ni la implementación de la clase matrix es la mejor para matrices dispersas, pero se usa por su interes pedagógico.
- La idea a transmitir es que hay que usar las librerías especializadas.

[^1] J. R. Gilbert, C. Moler, and R. Schreiber, “Sparse matrices in MATLAB: Design and implementation,” _SIAM Journal on Matrix Analysis and Applications_, vol. 13, no. 1, pp. 333–356, 1992, [doi: 10.1137/0613024](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=79bb9997f4a9f80f5a7e34868fb458d553abf01c).

[^2] MathWorks, "Constructing sparse matrices," _MathWorks_, 2021. [https://es.mathworks.com/help/matlab/math/constructing-sparse-matrices.html](https://es.mathworks.com/help/matlab/math/constructing-sparse-matrices.html).

[^3] Wikipedia, "Sparse matrix: Coordinate list (COO)," _Wikipedia_, 2025. [https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)](https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)).

[^4] NIST, "Matrix Market Exchange Formats", National Institute of Standards and Technology, 2022. [https://math.nist.gov/MatrixMarket/formats.html#mtx](https://math.nist.gov/MatrixMarket/formats.html#mtx).

[^5] NIST, "BCSSTK04: Stiffness matrix from structural engineering problem", Matrix Market, 2022. [https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc1/bcsstk04.html](https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc1/bcsstk04.html).

[^6] J. Burkardt, _SPARSEKIT: Sparse matrix software in FORTRAN77_. Department of Scientific Computing, Florida State University. [https://people.sc.fsu.edu/~jburkardt/f77_src/sparsekit/sparsekit.html](https://people.sc.fsu.edu/~jburkardt/f77_src/sparsekit/sparsekit.html).

### Clase `matrix`

Crearemos un módulo (clase) para matrices dispersas en formato COO.



Al ser un módulo largo, pegamos trozos en este Jupyter Notebook, lo que permite la coloración pero lo hace más rígido (no conocerá los cambios hechos en el código).

```fortran
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
interface get
  module procedure get_sca
  module procedure get_row
  module procedure get_col
end interface
private :: read_matrix, size_matrix, get_sca, get_col, get_row, set_sca, lu_matrix

contains
```

Se declara un tipo de dato `matrix` que guardará las variables (atributos de objeto):
  - `n`, número de filas,
  - `m`, número de columnas,
  - `row`, array de índices de fila
  - `row`, array de índices de columna
  - `val`, array de valores no nulos

Para aquellos procedimientos que vayan a tener otro nombre público, como `read_matrix` (bajo el nombre `read`):
  - el nombre público se indica en su interfaz y
  - se privatiza `matrix_read` para que no entre en conflicto con entidades declaradas en otros módulos.

Empezamos los métodos con cómo desalojar un objeto de clase `matrix`: si sus atributos están alojados se desalojan.

```fortran
subroutine dealloc(this)
type(matrix) :: this
integer :: i, res
character(260) :: cad

if (allocated(this%row)) deallocate(this%row, stat = res, errmsg = cad)
if (res /= 0) error stop '(matrix/dealloc) Unable to deallocate variable row: '//trim(cad)
if (allocated(this%col)) deallocate(this%col, stat = res, errmsg = cad)
if (res /= 0) error stop '(matrix/dealloc) Unable to deallocate variable col: '//trim(cad)
if (allocated(this%val)) deallocate(this%val, stat = res, errmsg = cad)
if (res /= 0) error stop '(matrix/dealloc) Unable to deallocate variable val: '//trim(cad)
end subroutine
```
- Se advierte de cualquier problema surgido, mediante el estado `stat` y el mensaje de error `errmsg`. 
- Se toman cadenas de 260 porque es la longitud máxima de una ruta en Windows API [https://learn.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation](https://learn.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation).
- `trim` sirve para no imprimir los blancos a la derecha.

Ahora, cómo alojar un objeto: 
```fortran
subroutine alloc(this, row, col, val)
type(matrix)             :: this
integer,      intent(in) :: row(:), col(:)
real(real64), intent(in) :: val(:)
character(260) ::  cad

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
```
- Se le deben pasar tres arrays de igual tamaño.
- Si alguno de los atributos está alojado, el objeto se desaloja primero.
- Se asignan los atributos del objeto.
- `error stop` comunica una salida por fallo (desde Fortran 2008). 


Seguidamente, escribimos la función para leer ficheros .mtx del tipo:
  ```text
  %%MatrixMarket matrix coordinate real symmetric
  132 132 1890
  1 1  1.9960268182200e+03
  2 1  5.6751562001600e+02
  3 1  7.7541907594400e+02
  ...
  ```
```fortran
subroutine read_matrix(this, filename)
type(matrix)             :: this
character(*), intent(in) :: filename
character(260) :: cad
integer :: iu, ios, nval, i

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
```
- Abre el fichero en una nueva unidad con `newunit`.
- Advierte de un error mediante `iostat`.
- Si el objeto está alojado, se desaloja previamente.
- Para uso intensivo de _Matrix Market_, se recomienda usar la libreria proporcionada en [https://math.nist.gov/MatrixMarket/mmio/f/mmio.f](https://math.nist.gov/MatrixMarket/mmio/f/mmio.f).



Ahora extendemos la función `size` a nuestra clase:

```fortran
function size_matrix(this, d) result(res)
type(matrix), intent(in) :: this
integer,      intent(in) :: d
integer :: res
character(260) :: cad

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
```  
- `size(objeto, 1)` y `size(objeto, 2)` devuelven el número de filas y columnas, respectivamente.
- Cualquier otra dimensión devuelve error.
- El texto del error mezcla cadenas y números; en Fortran debemos transformar el número mediente una _escritura interna_ en una cadena.
- La cadena debe desplazarse a la izquierda con `adjustl` antes de aliminar los blancos de la derecha con `trim`.

Empezamos con las funciones que sirven para obtener valores de la matriz. Todas ellas están bajo la interfaz:
  ```fortran
  interface get 
    module procedure get_sca
    module procedure get_row
    module procedure get_col
  end interface
  ```
La primera obtiene valores escalares:
 
```fortran
function get_sca(this, i, j) result(val)
type(matrix), intent(in) :: this
integer,      intent(in) :: i, j
real(real64)             :: val
integer :: k
character(260) :: cad

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
```
- Si el objeto no está alojado, devuelve un error.
- Si está alojado pero sin contenido, devuelve un error. Basta comprobar si `row` está alojado pues sólo `alloc` y `set` pueden crearlo.
- Se busca la posición de _(i,j)_. Si no está, se devuelve cero.


Ahora la función que devuelve una sección columna:
```fortran
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
```

Ahora la que devuelve una sección fila:
```fortran
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
```

Ahora, la función para establecer valores. Se indica el caso escalar y que acoge a la interfaz `get`, que podría ser ampliada en el futuro con funciones que establezcan secciones de la matriz.

```fortran
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
```    
- Si no está alojada, se crea una matriz $i\times j$ y se inserta el (unico) valor.
- Si está alojada, 
  - Si ya existe la posición $(i,j)$, se colocal el nuevo valor.
  - Si no existe, se crean tres arrays temporales, `trow`, `tcol` y `tval`, con una posición más y luego el objeto apunta a los nuevos arrayscon `move_alloc`. Esto reduce el trabajo de _crear, copiar, destruir_.


Veamos cómo queda la factorización LU, echando mano de los métodos `get` y `set`. Recordamos la fórmula para calcular $U_{ij}$ en matrices llenas:
   ```fortran
   u(i,j) = a(i,j) - dot_product(l(i,1:i-1), u(1:i-1,j))
   ```
   que ahora se transforma en:
   ```fortran
   call set(u, i, j, get(a, i, j) - dot_product(get(l, i, [(k, k=1,i-1)]), get(u, [(k, k=1,i-1)], j)))
   ```

```fortran
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
```

Ya podemos usar la clase para trabajar con una matrix _Matrix Market_.


```powershell
!pygmentize -g m3-poo/src/main.f90
```

## Ejercicios
1. Usa las funciones para devolver una sección fila o columna. Como el resultado es un array de tamaño variable, se recomienda usar la cláusula `source` introducida en Fortran 2003:
  ```fortran
  real(real64), allocatable :: b(:)`
  type(matrix) :: a
  ...
  allocate(b, source=get(a, 1, [(i, i=1,size(a,2))]))
  ```
Equivaldría en Matlab a `b = a(1,:)`.

2. Implementa una función para establecer una sección fila o columna.

3. Implementa la multiplicación matriz por vector y extiende `matmul` con ella.

## Más ejercicios a proponer
1. Implementar remonte y descenso para el tipo `matrix`.
2. Mostrar el tiempo de cálculo de la factorización en función de la dimensión; ver en el curso de Cocalc, la carpeta *fortran/met_num_mat/ANM/lu/llena/tiempo/*. Usar el paquete ogpf.


