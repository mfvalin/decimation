program decimate_test
  use decimate_array
  implicit none
#if ! defined(NI)
#define NI 12
#endif
#if ! defined(NJ)
#define NJ 12
#endif

  real, dimension(:,:,:), pointer :: src0
  real, dimension(:,:)  , pointer :: src1
  real, dimension(:,:)  , pointer :: src2
  real, dimension(:,:), pointer :: dst1
  real, dimension(:)  , pointer :: dst1d, dst1d2
  real, dimension(:,:), pointer :: dst2d
  integer :: i, j, status, by

  write(6,'(A,I3,A,I3,A)')" ==== original Fortran data (",NI,",",NJ,')'
  write(6,*)'bi-linear function : f(i,j) = (i + 1) * (j+1)'
  allocate( src0(NI,NJ,3) )
  allocate( src1(NI,NJ) )
  allocate( src2(NI,NJ) )
  write(6,*) 'allocated arrays'

  do j = 1, NJ
    do i = 1, NI
      src1(i,j)   = (i + 1) * (j+1)
    enddo
!     write(6,'(30F6.1)')src1(:,j)
  enddo

  ! syntax test
  src0 = 0
  dst1d => Decimate(src0(:,1,1), 2, NI)
  deallocate(dst1d)
  dst2d => Decimate(src0(:,:,1), 2, NI, NI, NJ)
  deallocate(dst2d)

  do by = 0, 11     ! testing 1 and 2 dimensional decimation/restore by factors of 0 through 11
    write(6,'(A,I3,A)')" ==== decimation by",by,' ===='
    dst1 => decimated_array(NI, NJ, by)           ! allocate container for 2D decimated array
    dst1 = 0
    write(6,*)'   1 dimensional test, along NI, then along NJ'
    ! =======================  1 D test along row  =======================
!   status = Decimate_1d(src1, by, dst1, NI, 0)   ! direct call to specific function
    dst1d => Decimate(src1(:,1), by, NI)          ! generic call, test along first row
    src2 = 0.0
!     status = UnDecimate_1d(dst1, by, src2, NI)
    status = UnDecimate(dst1d, by, src2(:,1), NI)
    if(maxval( abs(src2(:,1)-src1(:,1)) / src1(:,1) )  > 1.0E-6 ) then
      write(6,*) 'ERROR: max difference =',abs(maxval(src2(:,1)-src1(:,1))), &
                 ', max rel error  =',maxval(abs(src2(:,1)-src1(:,1))/src1(:,1))
      write(6,'(30F6.1)')src2(:,1)
      write(6,'(30F6.1)')abs(src2(:,1) - src1(:,1))
    else
      write(6,*) 'SUCCESS: max difference =',abs(maxval(src2(:,1)-src1(:,1))), &
                 ', max rel error  =',maxval(abs(src2(:,1)-src1(:,1))/src1(:,1))
    endif
    deallocate(dst1d)
    ! =======================  1 D test along column  ====================
    dst1d2 => Decimate(src1(1,:), by, NJ)       ! generic call, test along first column
    src2 = 0.0
    status = UnDecimate(dst1d2, by, src2(1,:), NJ)
    if(maxval( abs(src2(1,:)-src1(1,:)) / src1(1,:) )  > 1.0E-6 ) then
      write(6,*) 'ERROR: max difference =',abs(maxval(src2(1,:)-src1(1,:))), &
                 ', max rel error  =',maxval(abs(src2(1,:)-src1(1,:))/src1(1,:))
!       write(6,'(30F6.1)')src2(1,:)
!       write(6,'(30F6.1)')abs(src2(1,:) - src1(1,:))
    else
      write(6,*) 'SUCCESS: max difference =',abs(maxval(src2(1,:)-src1(1,:))), &
                 ', max rel error  =',maxval(abs(src2(1,:)-src1(1,:))/src1(1,:))
    endif
    deallocate(dst1d2)
    ! ===============================  2 D test  ===========================
    write(6,*)'   2 dimensional test'
    dst1 = 0.0
!   status = Decimate_2d(src1, by, dst1, NI, NI, NJ)   ! direct call to specific function
    dst2d => Decimate(src1, by, NI, NI, NJ)
!     do j = 1, N_decimated(NJ, by)
!       write(6,'(30F6.1)')dst1(:,j)
!     enddo

    src2 = 0.0
!   status = UnDecimate_2d(dst1, by, src2, NI, NI, NJ)   ! direct call to specific function
    status = UnDecimate(dst2d, by, src2, NI, NI, NJ)
!     if(any(abs(src1 - src2) > .0001 )) then
    if(maxval(abs(src2-src1)/src1) > 1.0E-6) then
       write(6,*) 'ERROR: max difference =',abs(maxval(src2-src1)), &
                  ', max rel error  =',(maxval(abs(src2-src1)/src1))
      do j = 1, NJ
        write(6,'(30F6.1)')src2(:,j)
!         write(6,'(30F6.1)')abs(src2(:,j) - src1(:,j))
      enddo
      write(6,*)
    else
      write(6,*) 'SUCCESS: max difference =',abs(maxval(src2-src1)), &
                ', max rel error  =',maxval(abs(src2-src1)/src1)
    endif
    deallocate(dst2d)
    deallocate(dst1)
 enddo

 contains

 subroutine ulp_diff_1(f1, f2, ni)
 end subroutine ulp_diff_1
  
end program

