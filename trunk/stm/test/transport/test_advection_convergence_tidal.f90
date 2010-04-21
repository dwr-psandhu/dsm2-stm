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
real(stm_real),parameter :: domain_length = 40960.0d0  ! meter
real(stm_real),parameter :: amplitude = half ! meter
real(stm_real),parameter :: gravity = 9.80d0 ! m/s^2
real(stm_real),parameter :: depth =16.0d0    ! meter

real(stm_real),parameter :: start_time = zero 
real(stm_real),parameter :: sec_per_hr = 60.d0*60.d0 
real(stm_real),parameter :: m2_period = 12.4d0*sec_per_hr
real(stm_real),parameter :: end_time = ten*m2_period
real(stm_real),parameter :: freq=two*pi/m2_period


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
    real(stm_real) :: big_b 
    real(stm_real) :: big_a 
    integer :: icell



big_b = freq/sqrt(gravity*depth)
big_a = amplitude* sqrt(gravity*depth)/depth/cos(big_b*domain_length)




    do icell = 1,ncell  ! width is equal to  1 meter 
      area(icell)    = depth + amplitude * cos(big_b*(domain_length-(dble(icell)-half)*dx))/cos(big_b*domain_length)*cos(freq*time)  
      area_lo(icell) = depth + amplitude * cos(big_b*(domain_length-(dble(icell-1)*dx)))    /cos(big_b*domain_length)*cos(freq*time)  
      area_hi(icell) = depth + amplitude * cos(big_b*(domain_length-(dble(icell)*dx)))      /cos(big_b*domain_length)*cos(freq*time)  
    
      flow(icell)    = area(icell)   *big_a*sin(big_b*(domain_length - (dble(icell)*dx-dx/two)))*sin(freq*time)
      flow_lo(icell) = area_lo(icell)*big_a*sin(big_b*(domain_length - (dble(icell)*dx-dx    )))*sin(freq*time)
      flow_hi(icell) = area_hi(icell)*big_a*sin(big_b*(domain_length - (dble(icell)*dx       )))*sin(freq*time)
    end do  ! icell


    
    return
end subroutine tidal_flow

subroutine test_tidal_advection_convergence
    use test_single_channel_advection
    use hydro_data
procedure(hydro_data_if),pointer :: tidal_hydro
integer, parameter  :: nstep_base = 40 
integer, parameter  :: nx_base    = 256
real(stm_real), parameter :: total_time = ten*m2_period
real(stm_real), parameter :: domain_length = 40960.0d0


character(LEN=10),parameter :: label = "tidal flow"
tidal_hydro=> tidal_flow


call test_round_trip(label,         &
                     tidal_hydro,   &
                     domain_length, &
                     total_time,    &                   
                     nstep_base,    &
                     nx_base)

return
end subroutine

end module


