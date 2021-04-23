
#if defined(IN_FORTRAN_CODE)

interface
  function N_decimated(npts, nd) result(n) bind(C,name='N_decimated')
    import :: C_INT
    implicit none
    integer(C_INT), intent(IN), value        :: npts
    integer(C_INT), intent(IN), value        :: nd
    integer(C_INT)  :: n
  end function N_decimated
  function Decimate_1d(src, factor, dst, ni, li) result(status) bind(C,name='Decimate_1d')
    import :: C_FLOAT, C_INT
    implicit none
    real(C_FLOAT), dimension(*), intent(IN)  :: src
    integer(C_INT), intent(IN), value        :: factor
    real(C_FLOAT), dimension(*), intent(OUT) :: dst
    integer(C_INT), intent(IN), value        :: ni
    integer(C_INT), intent(IN), value        :: li
    integer(C_INT)  :: status
  end function Decimate_1d
  function Decimate_2d(src, factor, dst, ni, li, nj) result(status) bind(C,name='Decimate_2d')
    import :: C_FLOAT, C_INT
    implicit none
    real(C_FLOAT), dimension(li,*), intent(IN)  :: src     ! 2 dimensional array
    integer(C_INT), intent(IN), value           :: factor
    real(C_FLOAT), dimension(1,*), intent(OUT)  :: dst     ! 2 dimensional array
    integer(C_INT), intent(IN), value           :: ni
    integer(C_INT), intent(IN), value           :: li
    integer(C_INT), intent(IN), value           :: nj
    integer(C_INT)  :: status
  end function Decimate_2d
  function UnDecimate_1d(src, factor, dst, ni) result(status) bind(C,name='UnDecimate_1d')
    import :: C_FLOAT, C_INT
    implicit none
    real(C_FLOAT), dimension(*), intent(IN)  :: src
    integer(C_INT), intent(IN), value        :: factor
    real(C_FLOAT), dimension(*), intent(OUT) :: dst
    integer(C_INT), intent(IN), value        :: ni
    integer(C_INT)  :: status
  end function UnDecimate_1d
  function UnDecimate_2d(src, factor, dst, ni, li, nj) result(status) bind(C,name='UnDecimate_2d')
    import :: C_FLOAT, C_INT
    implicit none
    real(C_FLOAT), dimension(1,*), intent(IN)   :: src     ! fudged notation to indicate 2 dimensional array
    integer(C_INT), intent(IN), value           :: factor
    real(C_FLOAT), dimension(li,*), intent(OUT) :: dst     ! 2 dimensional array
    integer(C_INT), intent(IN), value           :: ni
    integer(C_INT), intent(IN), value           :: li
    integer(C_INT), intent(IN), value           :: nj
    integer(C_INT)  :: status
  end function UnDecimate_2d
  subroutine QuantizeRestore(z, ni, li, nj, nbits) bind(C,name='QuantizeRestore')
    import :: C_FLOAT, C_INT
    implicit none
    real(C_FLOAT), dimension(*), intent(INOUT)  :: z     ! array to submit to quantize/restore
    integer(C_INT), intent(IN), value           :: ni
    integer(C_INT), intent(IN), value           :: li
    integer(C_INT), intent(IN), value           :: nj
    integer(C_INT), intent(IN), value           :: nbits ! number of bits to keep after quantization
  end subroutine QuantizeRestore
end interface

#else

int N_decimated(int npts, int nd);
int Decimate_1d(float *src, int factor, float *dst, int ni, int li);
int UnDecimate_1d(float *src, int factor, float *dst, int ni);
int Decimate_2d(float *src, int factor, float *dst, int ni, int li, int nj);
int Undecimate_2d(float *src, int factor, float *dst, int ni, int li, int nj);

#endif
