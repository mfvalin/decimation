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

  write(6,'(A,I4,A,I4,A)')" ==== original Fortran data (",NI,",",NJ,')'
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

  do by = 2, 8     ! testing 1 and 2 dimensional decimation/restore by factors of 2 through 8
    write(6,'(A,I3,A)')" ==== decimation by",by,' ===='
    dst1 => decimated_array(NI, NJ, by)           ! allocate container large enough for 2D decimated array
    dst1 = 0
    write(6,*)'   1 dimensional test, along NI, then along NJ'
    ! =======================  1 D test along row  =======================
    dst1d => Decimate(src1(:,1), by, NI)           ! generic call, with auto allocate
    status = Decimate(src1(:,1), by, dst1d, NI, 0) ! alternate form for generic call
    status = Decimate_1d(src1,   by, dst1d, NI, 0) ! explicit call to specific function
    src2 = 0.0

    status = UnDecimate(   dst1d, by, src2(:,1), NI) ! generic call
    status = UnDecimate_1d(dst1d, by, src2(:,1), NI) ! explicit call to specific function
    if(maxval( abs(src2(:,1)-src1(:,1)) / src1(:,1) )  > 1.0E-6 ) then
      write(6,1) 'ERROR: max difference =',abs(maxval(src2(:,1)-src1(:,1))), &
                 ', max rel error  =',maxval(abs(src2(:,1)-src1(:,1))/src1(:,1)), &
                 ', maxulp =',ulp_diff_1(src2(:,1),src1(:,1),NI)
      write(6,'(30F6.1)')src2(:,1)
      write(6,'(30F6.1)')abs(src2(:,1) - src1(:,1))
    else
      write(6,1) 'SUCCESS: max difference =',abs(maxval(src2(:,1)-src1(:,1))), &
                 ', max rel error  =',maxval(abs(src2(:,1)-src1(:,1))/src1(:,1)), &
                 ', maxulp =',ulp_diff_1(src2(:,1),src1(:,1),NI)
    endif
    deallocate(dst1d)
    ! =======================  1 D test along column  ====================
    dst1d2 => Decimate(src1(1,:), by, NJ)                ! generic call, with auto allocate
    status  = Decimate(src1(1,:),    by, dst1d2, NJ, 0)  ! alternate form for generic call
    status  = Decimate_1d(src1(1,:), by, dst1d2, NJ, 0)  ! explicit call to specific function
    src2 = 0.0

    status = UnDecimate(   dst1d2, by, src2(1,:), NJ)    ! generic call
    status = UnDecimate_1d(dst1d2, by, src2(1,:), NJ)    ! explicit call to specific function
    if(maxval( abs(src2(1,:)-src1(1,:)) / src1(1,:) )  > 1.0E-6 ) then
      write(6,1) 'ERROR: max difference =',abs(maxval(src2(1,:)-src1(1,:))), &
                 ', max rel error  =',maxval(abs(src2(1,:)-src1(1,:))/src1(1,:)), &
                 ', maxulp =',ulp_diff_1(src2(1,:),src1(1,:),NJ)
!       write(6,'(30F6.1)')src2(1,:)
!       write(6,'(30F6.1)')abs(src2(1,:) - src1(1,:))
    else
      write(6,1) 'SUCCESS: max difference =',abs(maxval(src2(1,:)-src1(1,:))), &
                 ', max rel error  =',maxval(abs(src2(1,:)-src1(1,:))/src1(1,:)), &
                 ', maxulp =',ulp_diff_1(src2(1,:),src1(1,:),NJ)
    endif
    deallocate(dst1d2)
    ! ===============================  2 D test  ===========================
    write(6,*)'   2 dimensional test'
    dst1 = 0.0
    dst2d => Decimate(   src1, by, NI, NI, NJ)          ! generic call, with auto allocate
    status = Decimate(   src1, by, dst2d, NI, NI, NJ)   ! alternate form for generic call
    status = Decimate_2d(src1, by, dst2d, NI, NI, NJ)   ! explicit call to specific function
!     do j = 1, N_decimated(NJ, by)
!       write(6,'(30F6.1)')dst1(:,j)
!     enddo

    src2 = 0.0
    status = UnDecimate(   dst2d, by, src2, NI, NI, NJ)   ! generic call
    status = UnDecimate_2d(dst2d, by, src2, NI, NI, NJ)   ! explicit call to specific function
    if(maxval(abs(src2-src1)/src1) > 1.0E-6) then
       write(6,*) 'ERROR: max difference =',abs(maxval(src2-src1)), &
                  ', max rel error  =',(maxval(abs(src2-src1)/src1)), &
                  ', maxulp =',ulp_diff_2(src2,src1,NI,NI,NJ)
      do j = 1, NJ
        write(6,'(30F6.1)')src2(:,j)
!         write(6,'(30F6.1)')abs(src2(:,j) - src1(:,j))
      enddo
      write(6,*)
    else
      write(6,1) 'SUCCESS: max difference =',abs(maxval(src2-src1)), &
                 ', max rel error  =',maxval(abs(src2-src1)/src1), &
                 ', maxulp =',ulp_diff_2(src2,src1,NI,NI,NJ)
    endif
    deallocate(dst2d)
    deallocate(dst1)
 enddo
1 format(1X,A,G12.4,A,G12.4,A,I6)
 contains

 function ulp_diff_1(f1, f2, ni) result(ulp)  ! max distance in units of last place for 1 D arrays
   implicit none
   integer, intent(IN), value :: ni
   real, dimension(ni) :: f1, f2
   integer :: ulp
   integer :: i
   ulp = 0
   do i=1,ni
     ulp = max(ulp, abs(transfer(f1(i),1)-transfer(f2(i),1)))
   enddo
 end function ulp_diff_1

 function ulp_diff_2(f1, f2, ni, li, nj) result(ulp)  ! max distance in units of last place for 2 D arrays
   implicit none
   integer, intent(IN), value :: ni, li, nj
   real, dimension(li,nj) :: f1, f2
   integer :: ulp
   integer :: i, j
   ulp = 0
   do j=1,nj
   do i=1,ni
     ulp = max(ulp, abs(transfer(f1(i,j),1)-transfer(f2(i,j),1)))
   enddo
   enddo
 end function ulp_diff_2
  
end program

