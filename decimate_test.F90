program decimate_test
  use decimate_array
  implicit none
#if ! defined(NI)
#define NI 12
#endif
#if ! defined(NJ)
#define NJ 12
#endif

  real, dimension(NI,NJ) :: src1, src2
  real, dimension(:,:), pointer :: dst1
  integer :: i, j, status, by

  write(6,'(A,I3,A,I3,A)')" ==== original Fortran data (",NI,",",NJ,')'
  write(6,*)'bi-linear function : f(i,j) = (i + 1) * (j+1)'
  do j = 1, NJ
    do i = 1, NI
!       src1(i,j) = (i + j - 1.0) * 1.0
      src1(i,j) = (i + 1) * (j+1)
    enddo
!     write(6,'(30F6.1)')src1(:,j)
  enddo
  write(6,*)

  do by = 2, 11
    write(6,'(A,I3,A)')" ==== decimation by",by,' ===='
    dst1 => decimated_array(NI, NJ, by)   ! allocate container for decimated array
    write(6,*)'   1 dimension test'
    dst1 = 0
    status = Decimate_1d(src1, by, dst1, NI, 0)
!     write(6,'(30F6.1)')dst1(:,1)

    src2 = 0.0
    status = UnDecimate_1d(dst1, by, src2, NI)
    if(any(src1(:,1) - src2(:,1) > .001 )) then
      write(6,*) 'ERROR: max difference =',abs(maxval(src2(:,1)-src1(:,1))), &
                  ', max rel error  =',maxval(abs(src2(:,1)-src1(:,1))/src1(:,1))
!       write(6,'(30F6.1)')abs(src2(:,1) - src1(:,1))
    else
      write(6,*) 'SUCCESS: max difference =',abs(maxval(src2(:,1)-src1(:,1))), &
                 ', max rel error  =',maxval(abs(src2(:,1)-src1(:,1))/src1(:,1))
    endif

    write(6,*)'   2 dimension test'
    dst1 = 0.0
    status = Decimate_2d(src1, by, dst1, NI, NI, NJ)
!     do j = 1, N_decimated(NJ, by)
!       write(6,'(30F6.1)')dst1(:,j)
!     enddo

    src2 = 0.0
    status = UnDecimate_2d(dst1, by, src2, NI, NI, NJ)
!     if(any(abs(src1 - src2) > .0001 )) then
    if(maxval(abs(src2-src1)/src1) > 1.0E-6) then
       write(6,*) 'ERROR: max difference =',abs(maxval(src2-src1)), &
                  ', max rel error  =',(maxval(abs(src2-src1)/src1))
!       do j = 1, NJ
!         write(6,'(30F6.1)')abs(src2(:,j) - src1(:,j))
!       enddo
      write(6,*)
    else
      write(6,*) 'SUCCESS: max difference =',abs(maxval(src2-src1)), &
                ', max rel error  =',maxval(abs(src2-src1)/src1)
    endif
    deallocate(dst1)
 enddo
  
end program

