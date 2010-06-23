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

!> Routines containing error metrics for assessing convergence or accuracy
!>@ingroup test
module error_metric

contains

!> Calculate the L-1, L-2 and L-inf error norms for calculated values and reference solution
subroutine error_norm(norm_1,norm_2,norm_inf,vals,reference,ncell,dx)

use stm_precision

implicit none

integer, intent(in) :: ncell                     !< Number of cells
real(stm_real), intent(out) :: norm_1            !< L-1   error norm
real(stm_real), intent(out) :: norm_2            !< L-2   error norm
real(stm_real), intent(out) :: norm_inf          !< L-inf error norm

real(stm_real), intent(in) :: vals(ncell)        !< Calculated values
real(stm_real), intent(in) :: reference(ncell)   !< Reference or 'other' values
real(stm_real), intent(in) :: dx                 !< Spatial step !todo: do we use this????

!------ locals                                    
integer :: icell                                
! todo: should we also report which_cell ?
integer :: which_cell                           !< The cell in which L-inf occurs
real(stm_real) :: err
real(stm_real) :: sq_error
real(stm_real) :: abs_error
!> initial value of all norms is zero
norm_1=zero
norm_2=zero
norm_inf=zero
!> sum up the L-1 and L-2, and fid the largest error in the domain (L-inf) 
do icell=1,ncell
   err = vals(icell) - reference(icell)
   abs_error = abs(err)
   sq_error = err*err
   if (abs_error > norm_inf) then
       norm_inf = abs_error
       which_cell = icell
   end if
   norm_1 = norm_1 + abs_error
   norm_2 = norm_2 + sq_error
end do
norm_1 = norm_1/dble(ncell)
norm_2 = norm_2/dble(ncell)
!todo: remove
!print*,which_cell,"/",ncell
return
end subroutine

end module