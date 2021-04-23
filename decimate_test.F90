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
  real, dimension(:)  , pointer :: dst1d, dst1d2
  real, dimension(:,:), pointer :: dst2d
  integer :: i, j, status, by, nbits

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
    ! =======================  1 D test along row  =======================
    write(6,*)'   1 dimensional test, along NI (row 1)'
    dst1d => Decimate(src1(:,1), by, NI)           ! generic call, with auto allocate
    call undecimate_row

    status = Decimate(src1(:,1), by, dst1d, NI, 0) ! alternate form for generic call
    call undecimate_row

    status = Decimate_1d(src1,   by, dst1d, NI, 0) ! explicit call to specific function
    call undecimate_row
    deallocate(dst1d)
    ! =======================  1 D test along column  ====================
    write(6,*)'   1 dimensional test, along NJ (column 1)'
    dst1d2 => Decimate(src1(1,:), by, NJ)                ! generic call, with auto allocate
    call undecimate_column

    status  = Decimate(src1(1,:),    by, dst1d2, NJ, 0)  ! alternate form for generic call
    call undecimate_column

    status  = Decimate_1d(src1(1,:), by, dst1d2, NJ, 0)  ! explicit call to specific function
    call undecimate_column
    deallocate(dst1d2)
    ! ===============================  2 D test  ===========================
    write(6,*)'   2 dimensional test'
    dst2d => Decimate(   src1, by, NI, NI, NJ)          ! generic call, with auto allocate
    call undecimate_array

    status = Decimate(   src1, by, dst2d, NI, NI, NJ)   ! alternate form for generic call
    call undecimate_array

    status = Decimate_2d(src1, by, dst2d, NI, NI, NJ)   ! explicit call to specific function
    call undecimate_array
    deallocate(dst2d)
  enddo

  nbits = 16
  write(6,'(A,I3,A)')" ==== quantization with",nbits,' bits ===='
  call quantize_test(nbits)
  nbits = 12
  write(6,'(A,I3,A)')" ==== quantization with",nbits,' bits ===='
  call quantize_test(nbits)
  nbits = 8
  write(6,'(A,I3,A)')" ==== quantization with",nbits,' bits ===='
  call quantize_test(nbits)

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

  subroutine verify_col
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
1 format(1X,A,G12.4,A,G12.4,A,I6)
  end subroutine verify_col

  subroutine verify_row
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
1 format(1X,A,G12.4,A,G12.4,A,I6)
  end subroutine verify_row

  subroutine verify_all
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
1 format(1X,A,G12.4,A,G12.4,A,I6)
  end subroutine verify_all

  subroutine undecimate_row
    src2 = 0.0
    status = UnDecimate(   dst1d, by, src2(:,1), NI) ! generic call
    call verify_row
    src2 = 0.0
    status = UnDecimate_1d(dst1d, by, src2(:,1), NI) ! explicit call to specific function
    call verify_row
  end subroutine undecimate_row

  subroutine undecimate_column
    src2 = 0.0
    status = UnDecimate(   dst1d2, by, src2(1,:), NJ)    ! generic call
    call verify_col
    src2 = 0.0
    status = UnDecimate_1d(dst1d2, by, src2(1,:), NJ)    ! explicit call to specific function
    call verify_col
  end subroutine undecimate_column

  subroutine undecimate_array
    src2 = 0.0
    status = UnDecimate(   dst2d, by, src2, NI, NI, NJ)   ! generic call
    call verify_all
    src2 = 0.0
    status = UnDecimate_2d(dst2d, by, src2, NI, NI, NJ)   ! explicit call to specific function
    call verify_all
  end subroutine undecimate_array
  
end program

subroutine quantize_test(nbits)
  use ISO_C_BINDING
  use decimate_array
  implicit none
  real(kind=4), dimension(NI+2,NJ) :: z
  integer, intent(IN) :: nbits
  integer :: i, j

  do j = 1, NJ
  do i = 1, NI
    z(i,j) = (i - NI*.5)*(i - NI*.5) + (j - NJ*.5)*(j - NJ*.5)
    z(i,j) = sqrt(z(i,j))
  enddo
  enddo
  call QuantizeRestore(z, NI, NI+2, NJ, nbits)
  write(6,*)
  call QuantizeRestore(z, NI, NI+2, NJ, nbits)
end subroutine quantize_test


