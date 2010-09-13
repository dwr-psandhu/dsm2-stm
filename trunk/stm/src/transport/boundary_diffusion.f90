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

!> boundary diffusive flux interface to be fulfilled by driver or application
!>@ingroup transport
module boundary_diffusion

 !> Calculate boundary diffusion flux
 abstract interface
   !> Generic interface for boundary diffusion that should be fulfilled by client programs
   subroutine boundary_diffusive_flux_if(diffusive_flux_lo, &
                                             diffusive_flux_hi, &
                                             conc,              &
                                             area_lo,           &
                                             area_hi,           &
                                             disp_coef_lo,      &  
                                             disp_coef_hi,      &
                                             ncell,             &
                                             nvar,              &
                                             time,              &
                                             dx,                &
                                             dt)
        
         
        
         use stm_precision
         implicit none
         !--- args
         integer, intent(in)  :: ncell                                    !< Number of cells
         integer, intent(in)  :: nvar                                     !< Number of variables
         real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)   !< Face flux, lo side
         real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)   !< Face flux, hi side
         real(stm_real), intent (in)   :: area_lo         (ncell)         !< Low side area centered at time
         real(stm_real), intent (in)   :: area_hi         (ncell)         !< High side area centered at time
         real(stm_real), intent (in)   ::  time                           !< Time
         real(stm_real), intent (in)   ::  conc(ncell,nvar)               !< Concentration 
         real(stm_real), intent (in)   :: disp_coef_lo (ncell)            !< Low side constituent dispersion coef.
         real(stm_real), intent (in)   :: disp_coef_hi (ncell)            !< High side constituent dispersion coef.
         real(stm_real), intent (in)   :: dt
         real(stm_real), intent (in)   :: dx
         
   end subroutine 
 end interface

 !> This pointer should be set by the driver or client code to specify the 
 !> treatment at the boundaries
 procedure(boundary_diffusive_flux_if),pointer :: boundary_diffusion_flux  => null()
!=========================================================
 abstract interface
       !> Generic interface for calculating BC of matrix that should be fulfilled by
       !> the driver or the client programs
       subroutine boundary_diffusive_matrix_if(center_diag ,           &
                                               up_diag,                &     
                                               down_diag,              &
                                               right_hand_side,        &
                                               conc,                   &
                                               explicit_diffuse_op,    &
                                               area,                   &
                                               area_lo,                &
                                               area_hi,                &          
                                               disp_coef_lo,           &
                                               disp_coef_hi,           &
                                               theta_stm,              &
                                               ncell,                  &
                                               time,                   & 
                                               nvar,                   & 
                                               dx,                     &
                                               dt)
                                                      
                                                      
        
        use stm_precision
        implicit none
        !--- args
                                       
        integer, intent (in) :: ncell                                               !< Number of cells
        integer, intent (in) :: nvar                                                !< Number of variables

        real(stm_real),intent (inout):: down_diag(ncell,nvar)                       !< Values of the coefficients below diagonal in matrix
        real(stm_real),intent (inout):: center_diag(ncell,nvar)                     !< Values of the coefficients at the diagonal in matrix
        real(stm_real),intent (inout):: up_diag(ncell,nvar)                         !< Values of the coefficients above the diagonal in matrix
        real(stm_real),intent (inout):: right_hand_side(ncell,nvar)                 !< Values of the coefficients of right  hand side vector
        real(stm_real), intent (in)  :: conc(ncell,nvar)                            !< concentration
        real(stm_real), intent (in)  :: explicit_diffuse_op(ncell,nvar)              
        real(stm_real), intent (in)  :: area (ncell)                                !< Cell centered area at new time 
        real(stm_real), intent (in)  :: area_lo(ncell)                              !< Low side area at new time
        real(stm_real), intent (in)  :: area_hi(ncell)                              !< High side area at new time 
        real(stm_real), intent (in)  :: disp_coef_lo (ncell)                        !< Low side constituent dispersion coef. at new time
        real(stm_real), intent (in)  :: disp_coef_hi (ncell)                        !< High side constituent dispersion coef. at new time
        real(stm_real), intent (in)  :: time                                        !< Current time
        real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
        real(stm_real), intent (in)  :: dx                                          !< Spatial step  
        real(stm_real), intent (in)  :: dt                                          !< Time step     

       end subroutine 
 end interface

 !> This pointer should be set by the driver or client code to specify the 
 !> treatment at the first and last row of coefficient matrix
 ! todo: check the "boundary_diffusion_matrix"
 procedure(boundary_diffusive_matrix_if),pointer :: boundary_diffusion_matrix  => null()

contains
 


subroutine set_diffusion_boundary(bc_flux,bc_matrix)
use error_handling
implicit none
procedure(boundary_diffusive_flux_if),pointer :: bc_flux      !< Diffusion flux routine
procedure(boundary_diffusive_matrix_if),pointer :: bc_matrix  !< Diffusion matrix routine
if ( (associated(bc_flux) .and. associated(bc_matrix)) .ne. &
     (associated(bc_flux) .or.  associated(bc_matrix)) )then
    call stm_fatal("Boundary diffusive flux and boundary diffusive matrix not consistently assigned (one null, other not)")
end if

boundary_diffusion_flux => bc_flux
boundary_diffusion_matrix => bc_matrix

return
end subroutine
 
 
 
 !> Example diffusive flux that prints an error and bails
 subroutine no_diffusion_flux(diffusive_flux_lo, &
                                         diffusive_flux_hi, &
                                         conc,              &
                                         area_lo,           &
                                         area_hi,           &
                                         disp_coef_lo,      &  
                                         disp_coef_hi,      &
                                         ncell,             &
                                         nvar,              &
                                         time,              &
                                         dx,                &
                                         dt)
     
     
     
     use stm_precision
     use error_handling
     implicit none
     !--- args
     integer, intent(in)  :: ncell                                    !< Number of cells
     integer, intent(in)  :: nvar                                     !< Number of variables
     real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)   !< Face flux, lo side
     real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)   !< Face flux, hi side
     real(stm_real), intent (in)   :: area_lo         (ncell)         !< Low side area centered at time
     real(stm_real), intent (in)   :: area_hi         (ncell)         !< High side area centered at time
     real(stm_real), intent (in)   ::  time                           !< Time
     real(stm_real), intent (in)   ::  conc(ncell,nvar)               !< Concentration 
     real(stm_real), intent (in)   :: disp_coef_lo (ncell)            !< Low side constituent dispersion coef.
     real(stm_real), intent (in)   :: disp_coef_hi (ncell)            !< High side constituent dispersion coef.
     real(stm_real), intent (in)   :: dt
     real(stm_real), intent (in)   :: dx
    
     call stm_fatal("boundary not implemented")
     
     return
 end subroutine 
 
!> Diffusion boundary condition that enforces zero constituent flux 
 subroutine neumann_zero_diffusive_flux(diffusive_flux_lo, &
                                           diffusive_flux_hi, &
                                           conc,              &
                                           area_lo,           &
                                           area_hi,           &
                                           disp_coef_lo,      &  
                                           disp_coef_hi,      &
                                           ncell,             &
                                           nvar,              &
                                           time,              &
                                           dx,                &
                                           dt)
    use stm_precision
    implicit none
    !--- args
    integer, intent(in)  :: ncell                                   !< number of cells
    integer, intent(in)  :: nvar                                    !< number of variables
    real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)  !< face flux, lo side
    real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)  !< face flux, hi side
    real(stm_real), intent (in)   :: area_lo(ncell)        !< Low side area centered at time
    real(stm_real), intent (in)   :: area_hi(ncell)        !< High side area centered at time
    real(stm_real), intent (in)   :: time                           !< time
    real(stm_real), intent (in)   :: conc(ncell,nvar)               !< concentration 
    real(stm_real), intent (in)   :: disp_coef_lo(ncell)      !< Low side constituent dispersion coef.
    real(stm_real), intent (in)   :: disp_coef_hi(ncell)      !< High side constituent dispersion coef.
    real(stm_real), intent (in)   :: dt
    real(stm_real), intent (in)   :: dx
    
    diffusive_flux_lo(1,:) = zero
    diffusive_flux_hi(ncell,:) = zero
    
       
    return
 end subroutine
 
 
!> Example diffusive flux that imposes sinusoidal time dependent Neumann boundary flux at
!> both ends of the channel.
 subroutine neumann_sin_diffusive_flux(diffusive_flux_lo, &
                                       diffusive_flux_hi, &
                                       conc,              &
                                       area_lo,           &
                                       area_hi,           &
                                       disp_coef_lo,      &  
                                       disp_coef_hi,      &
                                       ncell,             &
                                       nvar,              &
                                       time,              &
                                       dx,                &
                                       dt)
   use stm_precision
   implicit none
   !--- args
   integer, intent(in)  :: ncell                                   !< number of cells
   integer, intent(in)  :: nvar                                    !< number of variables
   real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)  !< face flux, lo side
   real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)  !< face flux, hi side
   real(stm_real), intent (in)   :: area_lo(ncell)        !< Low side area centered at time
   real(stm_real), intent (in)   :: area_hi(ncell)        !< High side area centered at time
   real(stm_real), intent (in)   ::  time                          !< time
   real(stm_real), intent (in)   ::  conc(ncell,nvar)              !< concentration 
   real(stm_real), intent (in)   :: disp_coef_lo(ncell)      !< Low side constituent dispersion coef.
   real(stm_real), intent (in)   :: disp_coef_hi(ncell)      !< High side constituent dispersion coef.
   real(stm_real), intent (in)   :: dt
   real(stm_real), intent (in)   :: dx
    
   diffusive_flux_lo(1,:) = two*cos(pi*time/three)               !Just for test 
   diffusive_flux_hi(ncell,:) = five*sin (pi*time/seven)
       
   return
 end subroutine
 

 
!===========================================
!> No-diffusion implementation for use as a pointer when diffusion is off.
subroutine no_diffusion_matrix(center_diag ,           &
                               up_diag,                &     
                               down_diag,              &
                               right_hand_side,        &
                               conc,                   &
                               explicit_diffuse_op,    &
                               area,                   &
                               area_lo,                &
                               area_hi,                &          
                               disp_coef_lo,           &
                               disp_coef_hi,           &
                               theta_stm,              &
                               ncell,                  &
                               time,                   & 
                               nvar,                   & 
                               dx,                     &
                               dt)
                                                      
                                                      
        
    use stm_precision
    use error_handling
    implicit none
    !--- args
                                   
    integer, intent (in) :: ncell                                               !< Number of cells
    integer, intent (in) :: nvar                                                !< Number of variables

    real(stm_real),intent (inout):: down_diag(ncell,nvar)                       !< Values of the coefficients below diagonal in matrix
    real(stm_real),intent (inout):: center_diag(ncell,nvar)                     !< Values of the coefficients at the diagonal in matrix
    real(stm_real),intent (inout):: up_diag(ncell,nvar)                         !< Values of the coefficients above the diagonal in matrix
    real(stm_real),intent (inout):: right_hand_side(ncell,nvar)                 !< Values of the coefficients of right  hand side vector
    real(stm_real), intent (in)  :: conc(ncell,nvar)                            !< concentration
    real(stm_real), intent (in)  :: explicit_diffuse_op(ncell,nvar)              
    real(stm_real), intent (in)  :: area (ncell)                                !< Cell centered area at new time 
    real(stm_real), intent (in)  :: area_lo(ncell)                              !< Low side area at new time
    real(stm_real), intent (in)  :: area_hi(ncell)                              !< High side area at new time 
    real(stm_real), intent (in)  :: disp_coef_lo (ncell)                        !< Low side constituent dispersion coef. at new time
    real(stm_real), intent (in)  :: disp_coef_hi (ncell)                        !< High side constituent dispersion coef. at new time
    real(stm_real), intent (in)  :: time                                        !< Current time
    real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
    real(stm_real), intent (in)  :: dx                                          !< Spatial step  
    real(stm_real), intent (in)  :: dt                                          !< Time step     
    
    call stm_fatal("boundary not implemented!")
    return
end subroutine 
 

 
 end module
