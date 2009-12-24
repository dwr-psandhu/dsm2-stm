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

!> Sample flow fields for tests
!>@ingroup test
module example_hydro_data

contains

!> Constant uniform flow
!> todo: needs to satisfy the hydro_data_if interface, has to change
subroutine constant_uniform(flow,flow_lo,flow_hi,ncell,time,q_const)
use stm_precision
implicit none
!--- args
integer,intent(in)  :: ncell  !< Number of cells

real(STM_REAL),intent(out) :: flow(ncell)    !< cell-centered flow at time
real(STM_REAL),intent(out) :: flow_lo(ncell) !< flow on lo side of cells at time
real(STM_REAL),intent(out) :: flow_hi(ncell) !< flow on hi side at time
real(STM_REAL), intent(in) :: time           !< time
real(STM_REAL),intent(in) :: q_const         !< constant flow to be used


!-------
flow = q_const
flow_lo=flow
flow_hi=flow

return
end subroutine

end module