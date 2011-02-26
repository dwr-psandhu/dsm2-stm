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
       subroutine spatial_data_sediment_if
                                    !< Number of cells (in)
                                    !< Time of request (in) (do we need this for future extension?)
                                    !< Spatial step : do we need this?  (in)
                                    !< Cell centerd Manning's n (intent out)
                                    !< Cell center width (intent out)
                                    !< any other argument we need here? 
                                        
     
  
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
       subroutine constant_data_sediment_if
                                    
                                 !< Gravitational acceleration(out)
                                 !< von Karman's constant(out)
                                 !< Kinematic viscosity of water(out)
                                 !< Water density(out)
                                 !< Sediment density(out)
                                 !< Floc density(out)
                                 !< Particle diameter(out)
                                 !< Representative cohesive mixture diameter(out) 
                                 !< Critical shear stress for full deposition (out)
                                 !< Wet bulk density of the deposit (out)
                                 !< Critical shear stress for partial deposition (out)
                                 !< Critical shear stress for surface erosion (out)
                                 !< Dry bulk density (out)
                                 !< Floc strength (out)

                           
     end subroutine
 
 end interface
 !> This pointer should be set by the driver or client code to specify the 
 !> Manning's n and channel width for the sediment source sink routine.
  procedure(constant_data_sediment_if),pointer :: fill_constant_sediment_data  => null()
 
 
 contains
  
  
  
  
  
end module