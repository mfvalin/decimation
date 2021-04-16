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
 contains
 
  function decimated_array(ni, nj, by) result (p)
    implicit none
    integer, intent(IN), value :: ni, nj, by
    real(kind=4), dimension(:,:), pointer :: p
    allocate( p(n_decimated(ni,by), n_decimated(nj,by) ) )
    return
  end function decimated_array

end module
