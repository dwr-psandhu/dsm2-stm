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

!> examples of advection
!>@ingroup example
module example_advection

    use stm_precision

    contains
    !> Subroutine that runs a uniform flow advection
    subroutine example_advection_uniform_flow

        use state_variables
        use primitive_variable_conversion
        use advection
        use example_initial_conditions
        use logging
        !use m_vstring

        !use example_hydro_data
        !use example_sources
        !use error_metric
        !use fruit
        implicit none

        !--- Problem variables

        integer, parameter  :: nstep  = 40
        integer, parameter  :: nx = 20
        real(STM_REAL), parameter :: cfl = 0.8

        integer, parameter  :: nconc = 2
        real(STM_REAL), parameter :: origin = zero        ! meters
        real(STM_REAL), parameter :: domain_length = 100  ! meters
        real(STM_REAL) :: dt              ! seconds
        real(STM_REAL) :: dx              ! meters
        real(STM_REAL), parameter :: ic_center      = three*fourth*domain_length
        real(STM_REAL), parameter :: ic_gaussian_sd = domain_length/sixteen
        real(STM_REAL), parameter :: constant_flow = 100
        real(STM_REAL), parameter :: constant_area = 100
        real(STM_REAL) :: vel
        real(STM_REAL) :: time
        integer :: itime = 0

        !------ local
        real(STM_REAL), allocatable :: reference(:)
        real(STM_REAL), allocatable :: x(:), x_shifted(:)
        real(STM_REAL), allocatable :: conc_exact(:)
        character(LEN=64) :: filename
        integer :: i
        integer :: funit, funit_ref   


        call allocate_state(nx,nconc)
        
        area = constant_area
        area_prev = area
        area_lo = area     ! todo: used?
        area_hi = area     ! todo: used? remove from advect?

        flow = constant_flow 
        flow_hi = flow
        flow_lo = flow
        vel = constant_flow/constant_area
        dx = domain_length/dble(nx)
        dt = cfl*dx/vel   
        
        allocate(x(nx),x_shifted(nx))
        do i = 1,nx                
            x(i) = dx*(dble(i)-half)+origin
        end do        
        
        call fill_rectangular(conc(:,1),x,nx,10D0,40D0,5D0,0D0)
        call fill_rectangular(conc(:,2),x,nx,10D0,40D0,5D0,0D0)
        
        call prim2cons( mass_prev,conc,area,nx,nconc)
        mass = mass_prev
        allocate(reference(ncell))  ! reference copy of initial state
        reference = conc(:,2)

        time = zero

        write(filename, "(a)"), "uniform_rectangular" 

        !call printout(conc(:,2),x,filename)

        funit     = 11
        funit_ref = 12
        open(unit = funit,     file = trim(filename)//".dat",     status = 'REPLACE')         
        open(unit = funit_ref, file = trim(filename)//"_ref.dat", status = 'REPLACE')        
        call printout(conc(:,2),x,time,11)
        call printout(conc(:,2),x,time,12)
        

        ! forwards
        do itime = 1,nstep
            time = itime * dt
           
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
                          dx)

              mass_prev = mass
              call cons2prim(conc,mass,area,nx,nconc) 
              
              call printout(conc(:,2),x,time,funit)

              call x_shift(x,x_shifted,nx,vel,time)
              
              call printout(reference,x_shifted,time,funit_ref)

        end do
        close(funit)
        close(funit_ref)
        deallocate(reference)
        call deallocate_state
        
    end subroutine example_advection_uniform_flow

    !< simple shift function for uniform advection plot
    subroutine x_shift(x_ref,x_shifted,nloc,vel,time)

        implicit none

        real(STM_REAL), intent(out)    :: x_shifted(nloc)  !< shifted x coordinate
        real(STM_REAL), intent(in)     :: x_ref(nloc)      !< reference x coordinate        
        real(STM_REAL), intent(in)     :: vel          !< constant velocity
        real(STM_REAL), intent(in)     :: time         !< solution at this time
        integer,        intent(in)     :: nloc         !< size of array

        !--local 
        integer :: i
        
        do i = 1, nloc
            x_shifted(i) = x_ref(i) + vel*time
        end do
        
        
    end subroutine  
    
end module example_advection