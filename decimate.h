
#if defined(IN_FORTRAN_CODE)

interface
  function N_decimated(npts, nd) result(n) bind(C,name='N_decimated')
    import :: C_INT
    implicit none
    integer(C_INT), intent(IN), value        :: npts
    integer(C_INT), intent(IN), value        :: nd
    integer(C_INT)  :: n
  end function N_decimated
  function Decimate_2d(src, factor, dst, ni, li, nj) result(status) bind(C,name='Decimate_2d')
    import :: C_FLOAT, C_INT
    implicit none
    real(C_FLOAT), dimension(*), intent(IN)  :: src
    integer(C_INT), intent(IN), value        :: factor
    real(C_FLOAT), dimension(*), intent(OUT) :: dst
    integer(C_INT), intent(IN), value        :: ni
    integer(C_INT), intent(IN), value        :: li
    integer(C_INT), intent(IN), value        :: nj
    integer(C_INT)  :: status
  end function Decimate_2d
  function UnDecimate_2d(src, factor, dst, ni, li, nj) result(status) bind(C,name='UnDecimate_2d')
    import :: C_FLOAT, C_INT
    implicit none
    real(C_FLOAT), dimension(*), intent(IN)  :: src
    integer(C_INT), intent(IN), value        :: factor
    real(C_FLOAT), dimension(*), intent(OUT) :: dst
    integer(C_INT), intent(IN), value        :: ni
    integer(C_INT), intent(IN), value        :: li
    integer(C_INT), intent(IN), value        :: nj
    integer(C_INT)  :: status
  end function UnDecimate_2d
end interface

#define ALLOC_DECIMATED(array, nx, ny, nby) allocate( array( n_decimated(nx,nby), n_decimated(ny,nby) ) )

#else

int N_decimated(int npts, int nd);
int Decimate_2d(float *src, int factor, float *dst, int ni, int li, int nj);
int Undecimate_2d(float *src, int factor, float *dst, int ni, int li, int nj);

#endif