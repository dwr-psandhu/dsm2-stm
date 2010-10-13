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

!> Testing the subroutine which provides the hydrodynamics of the tidal flow
!>@ingroup test
module test_tidal_flow_provider

use stm_precision
use test_advection_tidal_experience
use hydro_data

contains

! the subroutine tests the flow to satisfy the continuity of mass for water 
! dQ/dx + dA/dt = 0 
! here unit width is assumed Q=AU*1
subroutine test_tidal_hydro_provider

implicit none

integer,parameter :: ncell = 64           !< number of cells
integer,parameter :: nstep = 128 
real(stm_real),parameter :: start_time = zero
real(stm_real),parameter :: sec_per_hr = 60.d0*60.d0          !< Convert factor of hour to second 
real(stm_real),parameter :: total_time
 = 12.4d0*sec_per_hr     !< M2 tidal period 
real(stm_real),parameter :: domain_length = 128000.0d0 
  
real(stm_real):: time            !< time of request
real(stm_real):: dx              !< spatial step 
real(stm_real):: dt              !< time step 
real(stm_real):: flow(ncell)     !< cell centered flow
real(stm_real):: flow_lo(ncell)  !< lo face flow
real(stm_real):: flow_hi(ncell)  !< hi face flow
real(stm_real):: area(ncell)     !< cell center area
real(stm_real):: area_lo(ncell)  !< area lo face
real(stm_real):: area_hi(ncell)  !< area hi face
procedure(hydro_data_if),pointer :: tidal_hydro          !< The pointer points to tidal flow data

!--- local
integer :: itime
integer :: icell
real(stm_real):: flow_new(ncell)     !< cell centered flow at time n+1
real(stm_real):: flow_old(ncell)     !< cell centered flow at time n
real(stm_real):: area_new(ncell)     !< cell centered area at time n+1
real(stm_real):: area_old(ncell)     !< cell centered area at time n
real(stm_real):: mass_difference(ncell)

tidal_hydro=> tidal_flow_cell_average

time = start_time 
dt = total_time/nstep
dx = domain_length/ncell

do itime=1,nstep

call tidal_hydro(flow_old,&
                 flow_lo, &
                 flow_hi, &
                 area_old,    &
                 area_lo, &
                 area_hi, &
                 ncell,   &
                 time,    &
                 dx,      &                  
                 dt)
           
time = time + dt
           
call tidal_hydro(flow_new,&
                 flow_lo, &
                 flow_hi, &
                 area_new,&
                 area_lo, &
                 area_hi, &
                 ncell,   &
                 time,    &
                 dx,      &                  
                 dt)
 
 do icell=1,ncell-1               
    mass_difference(ncell) = (flow_old(icell+1)-flow_old(icell))/dx  - (area_new(icell)-area_old(icell))/dt
 end do 
 
! todo: remove
! print *, mass_difference         
! pause          
           


end do



end subroutine

end module