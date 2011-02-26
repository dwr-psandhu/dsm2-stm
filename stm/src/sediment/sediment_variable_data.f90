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

!> Hydrodynamics interface to provide depth and flow velocity to be fed by user
!>@ingroup sediment
module suspended_sediment_variable
      !> Generic interface for fetching hydrodynamic data (
    abstract  interface
       !> Get hydrodynamic data for sediment module.
       !> This data might be calculated from a function or provided by another module
        subroutine spatiotemporal_data_sediment_if(velocity,  &
                                                   depth,     &
                                                   ncell,     &
                                                   time,      &
                                                   dx,        &
                                                   dt)
        use stm_precision
        implicit none
        integer, intent(in) :: ncell                    !< Number of cells (in)
        real(stm_real), intent(in)  :: time             !< Time of request (in)
        real(stm_real), intent(in)  :: dx               !< Spatial step (in)
        real(stm_real), intent(in)  :: dt               !< Time step  (in)
        real(stm_real), intent(out) :: velocity(ncell)  !< Cell and time centered velocity (out)
        real(stm_real), intent(out) :: depth(ncell)     !< Cell center depth (out)
       ! todo: the signature of this interface may be subjected to change 
        
        end subroutine
      end interface
      
 !> This pointer should be set by the driver or client code to specify the 
 !> depth and velocity for the sediment source sink routine
 procedure(spatiotemporal_data_sediment_if),pointer :: fill_spatiotemporal_sediment_data  => null()
 
 !> Generic interface for fetching Manning's n and width (spacially variable)inputs
 abstract interface 
       !> Get spacial varible data (Manning's n and  for sediment module.
       !> This data might be calculated from a function or provided by a subroutine
       subroutine spatial_data_sediment_if(ncell,     &
                                           time,      &
                                           dx,        &
                                           dt,        &
                                           manning,   &
                                           width)
       
            use stm_precision
            implicit none
            
            integer, intent(in) :: ncell                     !< Number of cells (in)
            real(stm_real), intent(in)  :: time              !< Time of request (in) (do we need this for future extension?)
            real(stm_real), intent(in)  :: dx                !< Spatial step : do we need this?  (in)
            real(stm_real), intent(in)  :: dt                !< do we need this ? todo:
            real(stm_real), intent(out) :: manning(ncell)    !< Cell centerd Manning's n (intent out)
            real(stm_real), intent(out) :: width(ncell)      !< Cell center width (intent out)
                                                             !< todo: any other argument we need here? 
  
     end subroutine
 
 end interface
 !> This pointer should be set by the driver or client code to specify the 
 !> Manning's n and channel width for the sediment source sink routine.
  procedure(spatial_data_sediment_if),pointer :: fill_spatial_sediment_data  => null()
 
 !> Generic interface for feeding the constants to the sediment module
 abstract interface 
       !> Get spacial varible data Gravitational acceleration. von Karman's , Kinematic viscosity of water
       !> Water density, Sediment density, Floc density, Particle diameter ,Representative cohesive mixture diameter  
       !> Critical shear stress for full deposition, Wet bulk density of the deposit, Critical shear stress for partial deposition 
       !> Critical shear stress for surface erosion, Dry bulk density, Floc strength 

       !> This data might be provided by a subroutine
       subroutine constant_data_sediment_if(g_accelration,     &
                                            kapa,              &
                                            nu,                &
                                            rho_w,             &
                                            rho_s,             &
                                            rho_floc,          &
                                            diameter,          &
                                            cohesive_diameter, &
                                            ta_d_full,         &
                                            rho_wet_bulk,      &
                                            ta_d_particle,     &
                                            ta_c_surf_erosion, &
                                            rho_dry_bulk,      &
                                            ta_floc,           &
                                            time)
                                    
        use stm_precision
        implicit none
        
        ! todo: NOTE: UNIT is important here
        ! todo: Eli how we can set a whistleblower for wrong unit(s)[I mean SI and British] here?
        ! todo: 
        
    real(stm_real), intent(out) ::  g_accelration       !< Gravitational acceleration(out)
    real(stm_real), intent(out) ::  kapa                !< von Karman's constant(out)
    real(stm_real), intent(out) ::  nu                  !< Kinematic viscosity of water(out)
    real(stm_real), intent(out) ::  rho_w               !< Water density(out)
    real(stm_real), intent(out) ::  rho_s               !< Sediment density(out)
    real(stm_real), intent(out) ::  rho_floc            !< Floc density(out)
    real(stm_real), intent(out) ::  diameter            !< Particle diameter(out)
    real(stm_real), intent(out) ::  cohesive_diameter   !< Representative cohesive mixture diameter(out) 
    real(stm_real), intent(out) ::  ta_d_full           !< Critical shear stress for full deposition (out)
    real(stm_real), intent(out) ::  rho_wet_bulk        !< Wet bulk density of the deposit (out)
    real(stm_real), intent(out) ::  ta_d_particle       !< Critical shear stress for partial deposition (out)
    real(stm_real), intent(out) ::  ta_c_surf_erosion   !< Critical shear stress for surface erosion (out)
    real(stm_real), intent(out) ::  rho_dry_bulk        !< Dry bulk density (out)
    real(stm_real), intent(out) ::  ta_floc             !< Floc strength (out)
    real(stm_real), intent(in)  ::  time                !< Time todo: do we need more parameters here? 

                                
     end subroutine
 
 end interface
 !> This pointer should be set by the driver or client code to specify the 
 !> Manning's n and channel width for the sediment source sink routine.
  procedure(constant_data_sediment_if),pointer :: fill_constant_sediment_data  => null()
 
 !> todo: fill here 
 contains
 
 !> Example spatial and time variables that prints an error and bails
 subroutine example_spatiotemporal_data_sediment(velocity,  &
                                                 depth,     &
                                                 ncell,     &
                                                 time,      &
                                                 dx,        &
                                                 dt)
     use stm_precision
     use error_handling
     implicit none
     
        integer, intent(in) :: ncell                    !< Number of cells (in)
        real(stm_real), intent(in)  :: time             !< Time of request (in)
        real(stm_real), intent(in)  :: dx               !< Spatial step (in)
        real(stm_real), intent(in)  :: dt               !< Time step  (in)
        real(stm_real), intent(out) :: velocity(ncell)  !< Cell and time centered velocity (out)
        real(stm_real), intent(out) :: depth(ncell)     !< Cell center depth (out)
        ! just to avoid warning
        velocity = LARGEREAL
        depth = minus* LARGEREAL
     
        call stm_fatal('ERROR IN SPATIOTEMPORAL DATA !')
 
 end subroutine 
 
  !> Example spatial variables that prints an error and bails
 subroutine example_spatial_data_sediment(ncell,     &
                                          time,      &
                                          dx,        &
                                          dt,        &
                                          manning,   &
                                          width)
 
 
 
    use stm_precision
    use error_handling
    implicit none
               
            integer, intent(in) :: ncell                     !< Number of cells (in)
            real(stm_real), intent(in)  :: time              !< Time of request (in) (do we need this for future extension?)
            real(stm_real), intent(in)  :: dx                !< Spatial step : do we need this?  (in)
            real(stm_real), intent(in)  :: dt                !< do we need this ? todo:
            real(stm_real), intent(out) :: manning(ncell)    !< Cell centerd Manning's n (intent out)
            real(stm_real), intent(out) :: width(ncell)      !< Cell center width (intent out)
                                                             !< todo: any other argument we need here? 
       ! just to avoid warning
        manning = minus*LARGEREAL
        width = LARGEREAL
     
        call stm_fatal('ERROR IN SPATIAL DATA !')
 
 
 end subroutine
 
 !> Example constants that prints an error and bails
 subroutine example_constant_data_sedimet(g_accelration,     &
                                          kapa,              &
                                          nu,                &
                                          rho_w,             &
                                          rho_s,             &
                                          rho_floc,          &
                                          diameter,          &
                                          cohesive_diameter, &
                                          ta_d_full,         &
                                          rho_wet_bulk,      &
                                          ta_d_particle,     &
                                          ta_c_surf_erosion, &
                                          rho_dry_bulk,      &
                                          ta_floc,           &
                                          time)
 
     use stm_precision
     use error_handling
     implicit none
    real(stm_real), intent(out) ::  g_accelration       !< Gravitational acceleration(out)
    real(stm_real), intent(out) ::  kapa                !< von Karman's constant(out)
    real(stm_real), intent(out) ::  nu                  !< Kinematic viscosity of water(out)
    real(stm_real), intent(out) ::  rho_w               !< Water density(out)
    real(stm_real), intent(out) ::  rho_s               !< Sediment density(out)
    real(stm_real), intent(out) ::  rho_floc            !< Floc density(out)
    real(stm_real), intent(out) ::  diameter            !< Particle diameter(out)
    real(stm_real), intent(out) ::  cohesive_diameter   !< Representative cohesive mixture diameter(out) 
    real(stm_real), intent(out) ::  ta_d_full           !< Critical shear stress for full deposition (out)
    real(stm_real), intent(out) ::  rho_wet_bulk        !< Wet bulk density of the deposit (out)
    real(stm_real), intent(out) ::  ta_d_particle       !< Critical shear stress for partial deposition (out)
    real(stm_real), intent(out) ::  ta_c_surf_erosion   !< Critical shear stress for surface erosion (out)
    real(stm_real), intent(out) ::  rho_dry_bulk        !< Dry bulk density (out)
    real(stm_real), intent(out) ::  ta_floc             !< Floc strength (out)
    real(stm_real), intent(in)  ::  time                !< Time todo: do we need more parameters here? 

! just to avoid warning
        g_accelration    = LARGEREAL  
        kapa             = LARGEREAL  
        nu               = LARGEREAL 
        rho_w            = LARGEREAL 
        rho_s            = LARGEREAL 
        rho_floc         = LARGEREAL 
        diameter         = LARGEREAL 
        cohesive_diameter= LARGEREAL 
        ta_d_full        = LARGEREAL 
        rho_wet_bulk     = LARGEREAL 
        ta_d_particle    = LARGEREAL 
        ta_c_surf_erosion= LARGEREAL 
        rho_dry_bulk     = LARGEREAL 
        ta_floc          = LARGEREAL 
                    
   
        call stm_fatal('ERROR IN CONSTANT PARAMETERS OF SEDIMENT!')
   
 end subroutine  
  
end module