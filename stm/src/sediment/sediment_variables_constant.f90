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
 
 ! todo: I am not sure about the correct place of d(nvar) and 
 ! todo: if we need d_cohesive(nvar) or not 
 
 integer :: ncell                                 !< Number of computational cells
 integer :: nvar                                  !< Number of variables non cohesive
 !integer :: nvar_cohesive                        !< Number of variables cohesive
 !> Sediment constants   
 real(stm_real), save :: gravity                  !< Acceleration of gravity; it must be in SI units (constant)
 real(stm_real), save :: water_density            !< Water density
 real(stm_real), save :: sediment_density         !< Sediment density
 real(stm_real), save :: kappa                    !< von Karman's constant
 real(stm_real), save :: kinematic_viscosity      !< Kinematic viscosity of water    
 real(stm_real), save :: floc_density             !< Floc density
 real(stm_real), save, allocatable :: diameter(:) !< Particle diameter
 real(stm_real), save :: cohesive_diameter        !< Representative cohesive mixture diameter     
 real(stm_real), save :: crit_stress_full_dep     !< Critical shear stress for full deposition    
 real(stm_real), save :: density_wet_bulk         !< Wet bulk density of the deposit
 real(stm_real), save :: crit_stress_partial_dep  !< Critical shear stress for partial deposition
 real(stm_real), save :: crit_stress_surf_erosion !< Critical shear stress for surface erosion 
 real(stm_real), save :: density_dry_bulk         !< Dry bulk density 
 real(stm_real), save :: ta_floc                  !< Floc strength 
 !> Spatial sediment parameters
 real(stm_real), save, allocatable :: manning_n(:)!< Manning's n 
 real(stm_real), save, allocatable :: width(:)    !< Channel width 
    
 contains
      
!> Set sediment constants 
 subroutine set_sediment_constants
 use error_handling
 implicit none
    
 real(stm_real) :: gravity                     !< Gravitational acceleration in SI unit
 real(stm_real) :: water_density               !< Water density
 real(stm_real) :: sediment_density            !< Sediment density
 real(stm_real) :: kappa                       !< von Karman's constant
 real(stm_real) :: kinematic_viscosity         !< Kinematic viscosity of water    
 real(stm_real) :: floc_density                !< Floc density
 real(stm_real) :: cohesive_diameter           !< Representative cohesive mixture diameter       
 real(stm_real) :: crit_stress_full_dep        !< Critical shear stress for full deposition    
 real(stm_real) :: density_wet_bulk            !< Wet bulk density of the deposit
 real(stm_real) :: crit_stress_partial_dep     !< Critical shear stress for partial deposition
 real(stm_real) :: crit_stress_surf_erosion    !< Critical shear stress for surface erosion 
 real(stm_real) :: density_dry_bulk            !< Dry bulk density 
 real(stm_real) :: ta_floc                     !< Floc strength 
  
 ! todo: QUESTION this is a derived variable based on constants so should it be here at all?
 !real(stm_real) :: specific_submerged_gravity 
 !specific_submerged_gravity = LARGEREAL

 
 gravity                   = LARGEREAL
 water_density             = LARGEREAL
 sediment_density          = LARGEREAL
 kappa                     = LARGEREAL
 kinematic_viscosity       = LARGEREAL
 floc_density              = LARGEREAL
 cohesive_diameter         = LARGEREAL
 crit_stress_full_dep      = LARGEREAL
 density_wet_bulk          = LARGEREAL
 crit_stress_partial_dep   = LARGEREAL
 crit_stress_surf_erosion  = LARGEREAL
 density_dry_bulk          = LARGEREAL
 ta_floc                   = LARGEREAL             
       
 return
 end subroutine
    
  !> Allocate spatial sediment parameters 
 subroutine allocate_sediment_parameters(ncell,nvar)
     use error_handling
     implicit none
     integer,intent(in):: ncell    !< Number of cells
     integer,intent(in):: nvar     !< Number of sediment classes
            
     allocate(manning_n(ncell))
     allocate(width(ncell))
     allocate(diameter(nvar))
         
     manning_n = LARGEREAL
     width     = LARGEREAL
     diameter  = LARGEREAL
     
     return
 end subroutine
 
 !> Deallocate the sediment static variable
 !> and reset ncell and nvar to zero.
 subroutine deallocate_sediment_static()
     implicit none

     ncell = 0
     nvar  = 0 
   
     deallocate(manning_n)
     deallocate(width)
     deallocate(diameter)
   
     return
 end subroutine

end module



