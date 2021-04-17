! Copyright (C) 2021  Environnement Canada
!
! This is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
!
! This software is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! Author:
!     M. Valin,   Recherche en Prevision Numerique, 2020/2021
!
module decimate_array
  use ISO_C_BINDING
  implicit none
#include <decimate.hf>

  interface decimate
    procedure decimate1
    procedure Decimate_1d
    procedure decimate2
    procedure Decimate_2d
  end interface

  interface undecimate
    procedure UnDecimate_1d    ! (src, factor, dst, ni)
    procedure UnDecimate_2d    ! (src, factor, dst, ni, li, nj)
  end interface
 contains
 
  function decimated_array(ni, nj, by) result (p)  ! allocate container for decimated result
    implicit none
    integer, intent(IN), value :: ni, nj, by
    real(kind=4), dimension(:,:), pointer :: p
    allocate( p(n_decimated(ni,by), n_decimated(nj,by) ) )
    return
  end function decimated_array

  function decimate1(what, by, ni) result(d)
    implicit none
    real(kind=4), dimension(*), intent(IN) :: what
    integer, intent(IN), value :: ni, by
    real(kind=4), dimension(:), pointer :: d

    integer :: status

    allocate( d(n_decimated(ni,by)) )
    status = Decimate_1d(what, by, d, ni, 0)
  end function decimate1

  function decimate2(what, by, ni, li, nj) result(d)
    implicit none
    real(kind=4), dimension(li,*), intent(IN) :: what
    integer, intent(IN), value :: ni, li, nj, by
    real(kind=4), dimension(:,:), pointer :: d

    integer :: status

    allocate( d(n_decimated(ni,by), n_decimated(nj,by)) )
    status = Decimate_2d(what, by, d, ni, li, nj)
  end function decimate2

end module
