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

!> Test of mass transport convergence in uniform flow
!>@ingroup test
module test_linear_decay_no_flow

use stm_precision
use example_initial_conditions

integer :: istep = 0
integer, parameter  :: nstep_base = 40
integer, parameter  :: nx_base = 256
real(stm_real), parameter :: total_time = 128.d0
real(stm_real), parameter :: rate_1 = 0.01d0
real(stm_real), parameter :: rate_2 = 0.01d0

contains
!> Subroutine that runs a decay in the absence of diffusion
subroutine test_linear_decay_convergence(verbose)

use test_single_channel_advection
use hydro_data
use source_sink

implicit none
procedure(hydro_data_if),pointer :: no_flow_hydro
!procedure(source_if)    ,pointer :: compute_source

integer, parameter  :: nstep_base = 40 
integer, parameter  :: nx_base = 256
integer, parameter  :: nconc = 2
logical  :: verbose

real(stm_real),parameter :: domain_length = 51200.d0
real(stm_real),parameter :: origin =zero
real(stm_real),parameter :: peak_initial = sixteen
real(stm_real) :: fine_initial_condition(nx_base,nconc) !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)          !< reference solution at finest resolution
real(stm_real) :: xposition (nx_base)
real(stm_real) :: x_start = domain_length/four
real(stm_real) :: x_end = domain_length*three/four
real(stm_real) :: other_domain_value = zero
real(stm_real) :: delta_x = domain_length/nx_base

integer :: icell
character(LEN=*),parameter :: label = "linear decay no flow"

no_flow_hydro => no_flow
compute_source => linear_decay_source

!> Subroutine which generates fine initial values and reference values to compare with 
!> and feed the covvergence test subroutine.
do icell =1,nx_base
   xposition(icell)= dble(icell-1)*delta_x + delta_x/two
end do


fine_initial_condition(:,1) = ten
fine_initial_condition(:,2) = fine_initial_condition(:,1)


fine_solution(:,1) = fine_initial_condition(:,1) *exp(- rate_1*total_time)
fine_solution(:,2) = fine_initial_condition(:,2) *exp(- rate_2*total_time)


call test_advection_convergence(label,       &
                     no_flow_hydro,          &
                     domain_length,          &
                     total_time,             &
                     fine_initial_condition, &
                     fine_solution,          &            
                     nstep_base,             &
                     nx_base,                &
                     nconc,                  &
                     verbose)

end subroutine
!=========================
!>generat constant area and constant flow foreward and backward
subroutine no_flow(flow,    &
                   flow_lo, &
                   flow_hi, &
                   area,    &
                   area_lo, &
                   area_hi, &
                   ncell,   &
                   time,    &
                   dx,      &                        
                   dt)
                   
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

!--- local
real(stm_real), parameter :: constant_flow = zero 
real(stm_real), parameter :: constant_area = 27.d1 

if (time <= (total_time/two)) then
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
!========================================
subroutine linear_decay_source(source, & 
                               conc,   &
                               area,   &
                               flow,   &
                               ncell,  &
                               nvar,   &
                               time)
                                     

 use  primitive_variable_conversion
 implicit none
 
 !--- args
integer,intent(in)  :: ncell                      !< Number of cells
integer,intent(in)  :: nvar                       !< Number of variables
real(stm_real),intent(inout) :: source(ncell,nvar)!< cell centered source 
real(stm_real),intent(in)  :: conc(ncell,nvar)    !< Concentration
real(stm_real),intent(in)  :: area(ncell)         !< area at source     
real(stm_real),intent(in)  :: flow(ncell)         !< flow at source location
real(stm_real),intent(in)  :: time                !< time 
!--- local just for test
real(stm_real) :: mass (ncell,nvar)

! source must be in primitive variable 
call prim2cons(mass,conc,area,ncell,nvar)
source(:,1) = -rate_1*mass(:,1)
source(:,2) = -rate_2*mass(:,2) 
 
return
end subroutine 

end module