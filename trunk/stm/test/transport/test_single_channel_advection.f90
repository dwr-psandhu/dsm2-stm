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

!> Test advection in a single channel for a flow regime that
!> goes back and forth and ends up in the same spot periodically
!>@ingroup test
module test_single_channel_advection

contains




subroutine initial_fine_solution(fine_initial_condition, &
                                 fine_solution,          &
                                 nx_base,                &
                                 nconc                   &
                                    )




use example_initial_conditions
use stm_precision

implicit none

integer :: nconc
integer :: nx_base
real(stm_real) :: fine_initial_condition(nx_base,nconc)!< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)         !< reference solution at finest resolution
real(stm_real) :: ic_center
real(stm_real) :: ic_gaussian_sd
real(stm_real),parameter :: origin = zero
real(stm_real) :: doamin_length  



!call fill_gaussian(conc(:,1),nx,origin,dx, &
!                   three*fourth*domain_length,ic_gaussian_sd)
!call fill_gaussian(conc(:,2),nx,origin,dx, &
!                   one*fourth*domain_length,ic_gaussian_sd)
!call prim2cons( mass_prev,conc,area,nx,nconc)
!mass = mass_prev
!allocate(reference(ncell))  ! reference copy of initial state
!
!! todo: Here you coarsen the provided fine reference solution    
!reference = conc(:,2)

return
end subroutine




!> Subroutine that tests advection convergence of flow that 
!> makes a "round trip". In other words, it sloshes back and forth
!> and ends up in the same spot. 
subroutine test_round_trip(label,        &
                           hydro,        &
                           domain_length,&                          
                           total_time,   &
                           fine_initial_condition, &
                           fine_solution, &                            
                           nstep_base,   &
                           nx_base,      &
                           nconc)
    use hydro_data
    use boundary_advection_module
    use stm_precision
    use state_variables
    use primitive_variable_conversion
    use advection
    use example_initial_conditions
    use example_hydro_data
    use example_sources
    use error_metric
    use fruit
    use logging

implicit none

!--- Problem variables
procedure(hydro_data_if), pointer :: hydro
integer, intent(in) :: nconc
character(LEN=*) :: label
integer, intent(in) :: nstep_base
integer, intent(in) :: nx_base
real(stm_real), intent(in) :: fine_initial_condition(nx_base,nconc)  !< initial condition at finest resolution
real(stm_real), intent(in) :: fine_solution(nx_base,nconc)  !< reference solution at finest resolution
real(stm_real), intent(in) :: total_time
real(stm_real), intent(in) :: domain_length 

!----local
integer, parameter :: nrefine = 3
integer, parameter :: coarsen_factor = 2      ! coarsening factor used for convergence test
real(stm_real), parameter :: cfl = 0.8 
real(stm_real), parameter :: origin = zero   ! meters
real(stm_real), parameter :: constant_flow = 1.D2
real(stm_real), parameter :: constant_area = 1.D2

integer :: itime = 0
integer :: icell ! debug only -- remove later
integer :: icoarse = 0
integer :: nstep
integer :: nx
integer :: coarsening

character(LEN=64) filename

logical, parameter :: limit_slope = .false.

real(stm_real), allocatable :: reference(:)
real(stm_real), allocatable :: x_center(:)
real(stm_real) :: dt              ! seconds
real(stm_real) :: dx              ! meters
real(stm_real) :: ic_center
real(stm_real) :: ic_gaussian_sd
real(stm_real) :: vel
real(stm_real) :: time
real(stm_real) :: norm_error(3,nrefine)

!todo: this is really "no flux"
boundary_advection=>neumann_advective_flux
ic_center = three*fourth*domain_length
ic_gaussian_sd = domain_length/sixteen

! coarsening factor in convergence test
do icoarse = 1,nrefine
    coarsening = coarsen_factor**(icoarse-1)
    nx = nx_base/(coarsening)
    nstep = nstep_base/(coarsening)
    call allocate_state(nx,nvar)
    allocate(x_center(nx))
    dx = origin + domain_length/dble(nx)
    dt = total_time/dble(nstep)


    do icell = 1,nx
        x_center(icell) = (dble(icell)-half)*dx
    end do

    time = zero
    call hydro(flow,    &
               flow_lo, &
               flow_hi, &
               area,    &
               area_lo, &
               area_hi, &
               ncell,   &
               time,    &
               dx,      &               
               dt)
    area_prev = area
    
        ! todo: Here you coarsen the provided ic, move fill gaussian to calling routine
    call fill_gaussian(conc(:,1),nx,origin,dx, &
                       three*fourth*domain_length,ic_gaussian_sd)
    call fill_gaussian(conc(:,2),nx,origin,dx, &
                       one*fourth*domain_length,ic_gaussian_sd)
    call prim2cons( mass_prev,conc,area,nx,nconc)
    mass = mass_prev
    allocate(reference(ncell))  ! reference copy of initial state
    
    ! todo: Here you coarsen the provided fine reference solution    
    reference = conc(:,2)
    
    ! forwards
    do itime = 1,nstep
       time = time + dt
       call hydro(flow,    &
                  flow_lo, &
                  flow_hi, &
                  area,    &
                  area_lo, &
                  area_hi, &
                  ncell,   &
                  time,    &
                  dx,      &                  
                  dt)
      ! call transport using no_source as the callback 
      call advect(mass,     &
                  mass_prev,&  
                  flow,     &
                  flow_lo,  &
                  flow_hi,  &
                  area,     &
                  area_prev,&
                  area_lo,  &
                  area_hi,  &
                  ncell,    &
                  nvar,     &
                  time,     &
                  dt,       &
                  dx,       &
                  limit_slope)

      mass_prev = mass
      area_prev = area
      call cons2prim(conc,mass,area,nx,nconc) 
    end do
    
    write(filename, "(a\i3\'.txt')"), "uniform_gaussian_start_", ncell 
    call printout(reference,x_center,filename)
    write(filename, "(a\i3\'.txt')"), "uniform_gaussian_end_", ncell 
    call printout(conc(:,2),x_center,filename)
    ! test error norm over part of domain
    call error_norm(norm_error(1,icoarse), &
                    norm_error(2,icoarse), &
                    norm_error(3,icoarse), &
                    conc(:,2),reference,ncell,dx)

    deallocate(reference)
    deallocate(x_center)
    call deallocate_state
end do

! todo: four?
call assert_true(norm_error(1,2)/norm_error(1,1) > four,"L-1 second order convergence on " // trim(label))
call assert_true(norm_error(2,2)/norm_error(2,1) > four,"L-2 second order convergence on " // trim(label))
! This is known not to pass for second order convergence if limiter is on
call assert_true(norm_error(3,2)/norm_error(3,1) > four,"L-inf second order convergence on " // trim(label))

!todo:
!print *,label
!print *, 'L-inf = ', norm_error(3,2)/norm_error(3,1), 'L-2 = ',norm_error(2,2)/norm_error(2,1),'L-1 = ',norm_error(1,2)/norm_error(1,1)
!print *, 'dt',dt,'dx',dx, ' CFL = ' , dt/dx
!print *, '========'
return
end subroutine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Coarsen a solution at a fine level of resolution
subroutine coarsen(coarse_data,fine_data,ncell_fine,ncell_coarse, nvar)

use stm_precision
use error_handling

implicit none
!---arg
integer,intent(in) :: ncell_coarse
integer,intent(in) :: ncell_fine
integer,intent(in) :: nvar
real(stm_real), intent(in) :: fine_data(ncell_fine,nvar)
real(stm_real), intent(out):: coarse_data(ncell_coarse,nvar)

!---locals
real(stm_real) :: coarsen_factor
integer :: ivar
integer :: icell
integer :: icoarse

if ( mod(ncell_fine , ncell_coarse) /= 0) then

    call stm_fatal("Coarsening factor is not an integer!")
   
else

coarsen_factor = ncell_fine/ncell_coarse

    do ivar=1,nvar
        do icell=1,ncell_coarse
            coarse_data(icell,ivar) = zero
            icoarse = 0
            do while (icoarse < coarsen_factor) 
              coarse_data(icell,ivar) = coarse_data(icell,ivar)+ fine_data(icell*coarsen_factor-icoarse,ivar)
              icoarse= icoarse + 1   
            end do
            coarse_data(icell,ivar)= coarse_data(icell,ivar)/dble(coarsen_factor)
        end do
    end do
    
end if

!< test that coarsen_factor is correct multiple if not call stm_fatal
!< coarsen using averaging
!< don't forget a unit test. should cover cases where coarsen_factor does not
!< work, should test that first, middle last value are good in ivar = 1 and ivar =nvar
return
end subroutine






end module