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

!> Defines the input variables for the sediment transport sources.
!> This replaces the notion of "rough I/O" suggested by Jamie Anderson originally.
!> The module defines the basic parameters for the sediment transport sources based
!> on the table agreed on 2/11/2011.
!>@ingroup sediment
module sediment_variables
 use stm_precision
 
 integer :: ncell                                 !< Number of computational cells
 integer :: nvar                                  !< Number of variables
 !> Sediment constants   
 real(stm_real), save :: gravity                  !< Acceleration of gravity; it must be in SI units (constant)
 real(stm_real), save :: water_density            !< Water density
 real(stm_real), save :: sediment_density         !< Sediment density
 real(stm_real), save :: gravity                  !< Gravitational acceleration in SI unit
 real(stm_real), save :: water_density            !< Water density
 real(stm_real), save :: sediment_density         !< Sediment density
 real(stm_real), save :: kapa                     !< von Karman's constant
 real(stm_real), save :: kinematic_viscosity      !< Kinematic viscosity of water    
 real(stm_real), save :: floc_density             !< Floc density
 real(stm_real), save :: diameter                 !< Particle diameter
 real(stm_real), save :: cohesive_diameter        !< Representative cohesive mixture diameter     
 real(stm_real), save :: ta_d_full                !< Critical shear stress for full deposition    
 real(stm_real), save :: rho_wet_bulk             !< Wet bulk density of the deposit
 real(stm_real), save :: ta_d_particle            !< Critical shear stress for partial deposition
 real(stm_real), save :: ta_c_surf_erosion        !< Critical shear stress for surface erosion 
 real(stm_real), save :: rho_dry_bulk             !< Dry bulk density 
 real(stm_real), save :: ta_floc                  !< Floc strength 
 !> Spatial sediment parameters
 real(stm_real), save, allocatable :: manning_n(:)!< Manning's n 
 real(stm_real), save, allocatable :: width(:)    !< Channel width 
    
 contains
      
!> Set sediment constants and calculated sediment constants
 subroutine set_sediment_constant()
 use error_handling
 implicit none
 
 real(stm_real) :: gravity                     !< Gravitational acceleration in SI unit
 real(stm_real) :: water_density               !< Water density
 real(stm_real) :: sediment_density            !< Sediment density
 real(stm_real) :: kapa                        !< von Karman's constant
 real(stm_real) :: kinematic_viscosity         !< Kinematic viscosity of water    
 real(stm_real) :: floc_density                !< Floc density
 real(stm_real) :: diameter                    !< Particle diameter
 real(stm_real) :: cohesive_diameter           !< Representative cohesive mixture diameter       
 real(stm_real) :: ta_d_full                   !< Critical shear stress for full deposition    
 real(stm_real) :: rho_wet_bulk                !< Wet bulk density of the deposit
 real(stm_real) :: ta_d_particle               !< Critical shear stress for partial deposition
 real(stm_real) :: ta_c_surf_erosion           !< Critical shear stress for surface erosion 
 real(stm_real) :: rho_dry_bulk                !< Dry bulk density 
 real(stm_real) :: ta_floc                     !< Floc strength 
  
 ! todo: QUESTION this is a derived variable based on constants so should it be here at all?
 !real(stm_real) :: specific_submerged_gravity 
 !specific_submerged_gravity = LARGEREAL
 
 gravity             = LARGEREAL
 water_density       = LARGEREAL
 sediment_density    = LARGEREAL
 kapa                = LARGEREAL
 kinematic_viscosity = LARGEREAL
 floc_density        = LARGEREAL
 diameter            = LARGEREAL
 cohesive_diameter   = LARGEREAL
 ta_d_full           = LARGEREAL
 rho_wet_bulk        = LARGEREAL
 ta_d_particle       = LARGEREAL
 ta_c_surf_erosion   = LARGEREAL
 rho_dry_bulk        = LARGEREAL
 ta_floc             = LARGEREAL             
       
 return
 end subroutine
    
  !> Allocate spatial sediment parameters 
 subroutine allocate_sediment_spatial_parameters(ncell)
     use error_handling
     implicit none
     integer,intent(in):: ncell    !<Number of cells
                
     allocate(manning_n(ncell))
     allocate(width(ncell))
     
     manning_n = LARGEREAL
     width = LARGEREAL
     
     return
 end subroutine
 
 !> Deallocate the sediment static variable
 !> and reset ncell and nvar to zero.
 subroutine deallocate_sediment_static()
     implicit none

     ncell = 0
   
     deallocate(manning_n)
   
     return
 end subroutine


!! todo: Considering structural point of view I think it has to move to suspended sediment varible 
!!> Calculate shear stress
!    subroutine calculate_shear_stress(shear_stress,    &
!                                      velocity,        &
!                                      ncell,           &
!                                      density)
!        implicit none
!        real(stm_real), intent(out) :: shear_stress(ncell)        
!        integer, intent(in)         :: ncell
!        real(stm_real), intent(in)  :: velocity(ncell)
!        real(stm_real), intent(in)  :: density
!        
!        return  
!    end subroutine


end module



