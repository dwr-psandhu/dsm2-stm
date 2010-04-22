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

!> Test of transport in uniform flow
!>@ingroup test
module test_uniform_flow
use stm_precision
integer :: istep = 0
integer, parameter  :: nstep_base = 80
integer, parameter  :: nx_base = 256
real(stm_real), parameter :: total_time = 6400.D0

contains

subroutine uniform_flow(flow,    &
                        flow_lo, &
                        flow_hi, &
                        area,    &
                        area_lo, &
                        area_hi, &
                        ncell,   &
                        time,    &
                        dx,      &                        
                        dt)
    use stm_precision
    implicit none
    integer, intent(in) :: ncell                   !< number of cells
    real(stm_real), intent(in) :: time             !< time of request
    real(stm_real), intent(in) :: dx               !< spatial step
    real(stm_real), intent(in) :: dt               !< time step 
    real(stm_real), intent(out) :: flow(ncell)     !< cell and time centered flow
    real(stm_real), intent(out) :: flow_lo(ncell)  !< lo face flow, time centered
    real(stm_real), intent(out) :: flow_hi(ncell)  !< hi face flow, time centered
    real(stm_real), intent(out) :: area(ncell)     !< cell center area, old time
    real(stm_real), intent(out) :: area_lo(ncell)  !< area lo face, time centered
    real(stm_real), intent(out) :: area_hi(ncell)  !< area hi face, time centered

    
    !> local
    real(stm_real), parameter :: constant_flow = 1.D2
    real(stm_real), parameter :: constant_area = 1.D2 
    real(stm_real) :: fine_initial_condition(nx_base,nvar)  !< initial condition at finest resolution
    real(stm_real) :: fine_solution(nx_base,nvar)           !< reference solution at finest resolution


    if (time <= total_time/two) then
      flow = constant_flow
    else
      flow = minus * constant_flow
    end if
    flow_hi = flow
    flow_lo = flow
    area = constant_area
    area_lo = constant_area
    area_hi = constant_area
    return
end subroutine

!> Subroutine that runs a small advective simulation
subroutine test_uniform_flow_advection()
use test_single_channel_advection
use hydro_data
procedure(hydro_data_if),pointer :: uniform_hydro
integer, parameter  :: nstep_base = 40 
integer, parameter  :: nx_base = 256
real(stm_real) :: domain_length = 51200.d0
character(LEN=12),parameter :: label = "uniform flow"

call fill_gaussian(fine_initial_condition,...)

uniform_hydro=> uniform_flow
call test_round_trip(label,         &
                     uniform_hydro, &
                     domain_length, &
                     total_time,    &
                     fine_initial_condition, &
                     fine_solution, &                     
                     nstep_base,    &
                     nx_base)

end subroutine

end module