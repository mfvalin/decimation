program decimate_test
  use ISO_C_BINDING
  implicit none
#if ! defined(NI)
#define NI 12
#endif
#if ! defined(NJ)
#define NJ 12
#endif
#include <decimate.hf>

  real, dimension(NI,NJ) :: src1, src2
  real, dimension(:,:), allocatable :: dst1
  integer :: i, j, status, by

  write(6,'(A,I3,A,I3,A)')" ==== original Fortran data (",NI,",",NJ,')'
  do j = 1, NJ
    do i = 1, NI
      src1(i,j) = i + j - 1.0
    enddo
    write(6,'(30F5.1)')src1(:,j)
  enddo
  write(6,*)

  do by = 2, 5
    write(6,'(A,I3,A)')" ==== decimation by",by,' ===='
    ALLOC_DECIMATED(dst1, NI, NJ, by)
    dst1 = 0.0
    status = Decimate_2d(src1, by, dst1, NI, NI, NJ)
    do j = 1, N_decimated(NJ, by)
      write(6,'(30F5.1)')dst1(:,j)
    enddo
    write(6,*)
    src2 = 0.0
    status = UnDecimate_2d(dst1, by, src2, NI, NI, NJ)
    do j = 1, NJ
      write(6,'(30F5.1)')src2(:,j)
    enddo
    write(6,*)
    deallocate(dst1)
 enddo
  
end program

