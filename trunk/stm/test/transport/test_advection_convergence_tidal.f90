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

!> Testing advection of mass which is subjected to a tidal boundary
!>@ingroup test
module test_advection_tidal
use stm_precision
!----- module variables
! todo: make the names more meaningful
real(stm_real),parameter :: origin = zero 
real(stm_real),parameter :: domain_length = 409600.0d0
real(stm_real),parameter :: amp = half
real(stm_real),parameter :: gravity = 9.80d0
real(stm_real),parameter :: depth =16.0d0
real(stm_real),parameter :: omega= 0.506708d0 ! hr
real(stm_real),parameter :: big_a = 0.4188704167062d0
real(stm_real),parameter :: big_b = 0.040465522644d0
real(stm_real),parameter :: k_0 = 10.066637844459d0
real(stm_real),parameter :: start_time = zero 
real(stm_real),parameter :: end_time = 124d0 ! 12.4 is one exact M2 cycle of tide in hours

contains

subroutine tidal_flow(flow,    &
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
    real(stm_real), intent(in) :: time             !< time of request "old time"
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
    integer :: icell

    do icell = 1,ncell
      flow_lo =  area_lo(icell)*big_a*sin(big_b*(domain_length - (dble(icell)*dx-dx)))*sin(omega*time)
      flow_hi = area_hi(icell)*big_a*sin(big_b*(domain_length - (dble(icell)*dx)))*sin(omega*time)
      flow = area(icell)*big_a*sin(big_b*(domain_length - (dble(icell)*dx-dx/two)))*sin(omega*time)
    end do

    area = constant_area
    area_lo = constant_area
    area_hi = constant_area
    return
end subroutine

subroutine test_tidal_advection_convergence
use test_single_channel_advection
use hydro_data
procedure(hydro_data_if),pointer :: tidal_hydro
integer, parameter  :: nstep_base = 40 
integer, parameter  :: nx_base    = 256
real(stm_real), parameter :: total_time = 124.D0
real(stm_real), parameter :: domain_length = 124.D0
character(LEN=10),parameter :: label = "tidal flow"
tidal_hydro=> tidal_flow
call test_round_trip(label,         &
                     tidal_hydro,   &
                     domain_length, &
                     total_time,    &                   
                     nstep_base,    &
                     nx_base)

end subroutine

end module


