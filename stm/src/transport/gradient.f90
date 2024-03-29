!<license>
!    Copyright (C) 1996, 1997, 1998, 2001, 2007, 2009 State of California,
!    Department of Water Resources.
!    This file is part of DSM2.
!
!    The Delta Simulation Model 2 (DSM2) is free software: 
!    you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    DSM2 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with DSM2.  If not, see <http://www.gnu.org/licenses>.
!</license>

!> Module containing routines for calculating differences and limiters
!>@ingroup transport
module gradient
contains

!> Calculate the undivided lo, hi, and centered differences
subroutine difference(grad_lo,grad_hi,grad_center,vals,ncell,nvar)

use stm_precision

implicit none

!---- args
integer,intent(in)  :: ncell                          !< Number of cells
integer,intent(in)  :: nvar                           !< Number of variables
real(stm_real),intent(in)  :: vals(ncell,nvar)        !< Data to be differenced
real(stm_real),intent(out) :: grad_lo(ncell,nvar)     !< Difference on lo side, LARGEREAL in first index
real(stm_real),intent(out) :: grad_hi(ncell,nvar)     !< Difference on hi side (n+1) minus (n) LARGEREAL for last index
real(stm_real),intent(out) :: grad_center(ncell,nvar) !< Dentered diff, LARGEREAL for undefined boundary cells

!----local
integer :: ivar

do ivar = 1, nvar
  grad_center(2:(ncell-1),ivar) = (vals(3:ncell,ivar) - vals(1:(ncell-2),ivar))/two
  grad_center(1,ivar)=LARGEREAL
  grad_center(ncell,ivar)=LARGEREAL
  grad_hi(1:(ncell-1),ivar) = (vals(2:ncell,ivar) - vals(1:(ncell-1),ivar))
  grad_hi(ncell,ivar)=LARGEREAL
  grad_lo(2:ncell,ivar)=grad_hi(1:(ncell-1),ivar)
  grad_lo(1,ivar)=LARGEREAL
end do

return
end subroutine


!> Apply a flux limiter (van Leer) given one-sided and centered differences
subroutine limiter(grad_lim,grad_lo,grad_hi,grad_center,ncell,nvar)

use stm_precision
implicit none

!--- args
integer,intent(in)  :: ncell                         !< Number of cells
integer,intent(in)  :: nvar                          !< Number of variables
real(stm_real),intent(in) :: grad_lo(ncell,nvar)     !< Difference on lo side, LARGEREAL in first index
real(stm_real),intent(in) :: grad_hi(ncell,nvar)     !< Difference on hi side (n+1) minus (n) LARGEREAL for last index
real(stm_real),intent(in) :: grad_center(ncell,nvar) !< Centered difference, LARGEREAL for undefined boundary cells 
real(stm_real),intent(out) :: grad_lim(ncell,nvar)   !< Limited difference

!---locals
real(stm_real) :: delta_limit(ncell,nvar) ! Intermediate quantity
real(stm_real) :: sign                           
integer        :: ivar, icell             ! Counting variables

do ivar = 1,nvar
    do icell = 1,ncell
        delta_limit(icell,ivar) = two*min(abs(grad_lo(icell,ivar)), &
                                          abs(grad_hi(icell,ivar)) )
                                          
        if (grad_center(icell,ivar) < zero)then
            sign = minus
        else
            sign = one
        end if
        grad_lim(icell,ivar) = min(abs(grad_center(icell,ivar)), &
                                 abs(delta_limit(icell,ivar)))*sign
    end do
end do
where (grad_lo*grad_hi < zero)
    grad_lim = zero
end where

! Boundary values are not defined
do ivar = 1,nvar
  grad_lim(1,ivar)    = LARGEREAL   !todo: is this really what we want? 
  grad_lim(ncell,ivar)= LARGEREAL   !todo: is this really what we want? 
end do

return
end subroutine

end module

