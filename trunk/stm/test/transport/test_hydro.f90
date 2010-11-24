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
module test_tidal_flow

use stm_precision
use test_advection_tidal_experience
use hydro_data
use fruit
use error_metric


contains

! the subroutine tests the flow to satisfy the continuity of mass for water 
! dQ/dx + dA/dt = 0 
! here unit width is assumed Q=AU*1
subroutine test_tidal_hydro

implicit none

! it must be equal to test_advection_tidal_experience
integer,parameter :: ncell = 32          !< number of cells
integer,parameter :: nstep = 64
  
real(stm_real):: time            !< time of request
real(stm_real):: dx              !< spatial step 
real(stm_real):: dt              !< time step 
real(stm_real):: flow(ncell)     !< cell centered flow
real(stm_real):: flow_lo(ncell)  !< lo face flow
real(stm_real):: flow_hi(ncell)  !< hi face flow
real(stm_real):: area(ncell)     !< cell center area
real(stm_real):: area_lo(ncell)  !< area lo face
real(stm_real):: area_hi(ncell)  !< area hi face
real(stm_real):: domain_length
procedure(hydro_data_if),pointer :: tidal_hydro          !< The pointer points to tidal flow data

!--- local
integer :: itime
integer :: icell
real(stm_real):: flow_new(ncell)     !< cell centered flow at time n+1
real(stm_real):: flow_hi_half(ncell)  !< high side flow at time n+1/2
real(stm_real):: flow_lo_half(ncell)  !< low side flow at time n+1/2
real(stm_real):: area_new(ncell)     !< cell centered area at time n+1
real(stm_real):: area_old(ncell)     !< cell centered area at time n
real(stm_real):: mass_difference(ncell)
real(stm_real):: max_mass_diff(nstep)
real(stm_real):: l1_mass_diff(nstep)
real(stm_real):: l2_mass_diff(nstep)
integer       :: which_cell
real(stm_real):: all_zero(ncell)


tidal_hydro=> tidal_flow_modified

time = start_time
dt = total_time/nstep
dx = domain_length/ncell

all_zero = zero

call tidal_hydro(flow_new,&
                 flow_lo_half, &
                 flow_hi_half, &
                 area_old,    &
                 area_lo, &
                 area_hi, &
                 ncell,   &
                 time,    &
                 dx,      &                  
                 dt)

do itime=1,nstep
  time = time + dt        
  call tidal_hydro(flow_new,&
                   flow_lo_half, &
                   flow_hi_half, &
                   area_new,&
                   area_lo, &
                   area_hi, &
                   ncell,   &
                   time,    &
                   dx,      &
                   dt)
          
     mass_difference = (flow_hi_half-flow_lo_half)*dt  + (area_new-area_old)*dx
  
   call error_norm(l1_mass_diff(itime),   &
                   l2_mass_diff(itime),   &
                   max_mass_diff(itime),  &
                   which_cell,            &
                   mass_difference,       &
                   all_zero,              &
                   ncell,     &
                   dx)  
  
   area_old = area_new
end do

call assert_true (maxval(max_mass_diff) < weak_eps ,'water mass balance error in tidal flow generator')

return
end subroutine

subroutine test_zoppou_flow()

use test_zoppou_advection_dispersion

implicit none

integer :: ncell
integer :: nstep
real(stm_real) :: start_t
real(stm_real) :: total_time
real(stm_real) :: domain_length
procedure(hydro_data_if),pointer :: hydrodynamics => null()
character (LEN=64) :: hydro_label

!todo: these numbers are hardwired 
ncell = 256
nstep = 128
start_t= zero
total_time = 2000.d0
domain_length = 2048.0d0
hydro_label = 'zoppou_flow'

hydrodynamics => zoppou_flow

call test_flow_continuity(ncell,        &
                          nstep,        &
                          start_t,   &
                          total_time,   &
                          domain_length,&
                          hydrodynamics,&
                          hydro_label)



end subroutine





!>the subroutine tests the flow to satisfy the continuity of mass for water 
!> dQ/dx + dA/dt = 0 
subroutine test_flow_continuity(ncell,        &
                                nstep,        &
                                start_time,   &
                                total_time,   &
                                domain_length,&
                                hydrodynamics,&
                                label)


implicit none

! it must be equal to test_advection_tidal_experience
integer,intent(in) :: ncell      !< number of cells
integer,intent(in) :: nstep      !< number of time 

real(stm_real),intent(in):: total_time             !< total time
real(stm_real),intent(in):: start_time             !< start time
real(stm_real),intent(in):: domain_length          !< domain_length
character(LEN=*),intent(in) :: label               !< label for test
procedure(hydro_data_if), pointer, intent(in) :: hydrodynamics 

!--- local
integer :: itime
integer :: icell
integer :: which_cell
real(stm_real):: time           !< time of request
real(stm_real):: flow(ncell)     !< cell centered flow
real(stm_real):: flow_lo(ncell)  !< lo face flow
real(stm_real):: flow_hi(ncell)  !< hi face flow
real(stm_real):: area(ncell)     !< cell center area
real(stm_real):: area_lo(ncell)  !< area lo face
real(stm_real):: area_hi(ncell)  !< area hi face
real(stm_real):: flow_new(ncell)     !< cell centered flow at time n+1
real(stm_real):: flow_hi_half(ncell)  !< high side flow at time n+1/2
real(stm_real):: flow_lo_half(ncell)  !< low side flow at time n+1/2
real(stm_real):: area_new(ncell)     !< cell centered area at time n+1
real(stm_real):: area_old(ncell)     !< cell centered area at time n
real(stm_real):: mass_difference(ncell)
real(stm_real):: max_mass_diff(nstep)
real(stm_real):: l1_mass_diff(nstep)
real(stm_real):: l2_mass_diff(nstep)
real(stm_real):: all_zero(ncell)
real(stm_real):: dt
real(stm_real):: dx

time = start_time
dt = total_time/dble(nstep)
dx = domain_length/dble(ncell)

all_zero = zero

call hydrodynamics(flow_new,&
                   flow_lo_half, &
                   flow_hi_half, &
                   area_old,    &
                   area_lo, &
                   area_hi, &
                   ncell,   &
                   time,    &
                   dx,      &                  
                   dt)

do itime=1,nstep
  time = time + dt        
  call hydrodynamics(flow_new,&
                   flow_lo_half, &
                   flow_hi_half, &
                   area_new,&
                   area_lo, &
                   area_hi, &
                   ncell,   &
                   time,    &
                   dx,      &
                   dt)
          
     mass_difference = (flow_hi_half-flow_lo_half)*dt  + (area_new-area_old)*dx
  
   call error_norm(l1_mass_diff(itime),   &
                   l2_mass_diff(itime),   &
                   max_mass_diff(itime),  &
                   which_cell,            &
                   mass_difference,       &
                   all_zero,              &
                   ncell,     &
                   dx)  
  
   area_old = area_new
end do

call assert_true (maxval(max_mass_diff) < weak_eps ,'water mass balance error in flow generator'//label)

return
end subroutine

end module