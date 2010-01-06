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

!> Module orchestrating the diffusion scheme. The main
!> routine in the module is diffuse().
!> Explicit and implicit diffusion operators are included here.
!>@ingroup transport
module diffusion

contains

! ///////////////////////////////////////////////////////////////////

!> Calculates the diffusive portion of the constituent transport.
!> It contains an explicit version of the diffusion operator and a general (involving all
!> potential cases) diffusion operator as well, with a coefficient theta_stm for 
!> selecting the level of implicitness. (theta_stm=0.5 is Crank Nicolson.).
!> The matrix is solved via a tri-diagonal solver.  !
!> The algoritm looks like this:
!>   - This creates the diffusive fluxes sends them for modification for boundaries and then differences the fluxes to get the operator d/dx(Ad/dx). 
!>         - Calculate interior and boundary fluxes 
!>   - Construct right hand side with Nuemann boundary condition imposed
!>   - Construct diffusion coefficeint matrix with Nuemann boundary condition imposed
!>   - Solve the system
subroutine diffuse(conc,             &
                  conc_prev,         &
                  mass,              &
                  mass_prev,         &
                  area,              &
                  area_prev,         &
                  area_lo,           &
                  area_hi,           &
                  area_lo_prev,      &
                  area_hi_prev,      &
                  disp_coef_lo,      &  
                  disp_coef_hi,      &
                  disp_coef_lo_prev, &  
                  disp_coef_hi_prev, &
                  ncell,             &
                  nvar,              &
                  time,              &
                  theta_stm,         &
                  dt,                &
                  dx                 )

use stm_precision
use primitive_variable_conversion 

implicit none

! ---- args

integer, intent (in) :: ncell                                !< Number of cells
integer, intent (in) :: nvar                                 !< Number of variables

real(stm_real), intent (out) :: conc(ncell,nvar)             !< Concentration at new time
real(stm_real), intent (out) :: mass(ncell,nvar)             !< Mass (A*C) at new time
real(stm_real), intent (in) :: mass_prev(ncell,nvar)         !< Mass (A*C) at old time
real(stm_real), intent (in) :: conc_prev(ncell,nvar)         !< Concentration at old time
real(stm_real), intent (in) :: area (ncell)                  !< Cell-centered area at new time
real(stm_real), intent (in) :: area_prev (ncell)             !< Cell-centered area at old time
real(stm_real), intent (in) :: area_lo (ncell)               !< Low side area centered in time
real(stm_real), intent (in) :: area_hi (ncell)               !< High side area centered in time 
real(stm_real), intent (in) :: area_lo_prev (ncell)          !< Low side area centered at old time
real(stm_real), intent (in) :: area_hi_prev (ncell)          !< High side area centered at old time 
real(stm_real), intent (in) :: disp_coef_lo (ncell,nvar)     !< Low side constituent dispersion coef. at new time
real(stm_real), intent (in) :: disp_coef_hi (ncell,nvar)     !< High side constituent dispersion coef. at new time
real(stm_real), intent (in) :: disp_coef_lo_prev(ncell,nvar) !< Low side constituent dispersion coef. at old time
real(stm_real), intent (in) :: disp_coef_hi_prev(ncell,nvar) !< High side constituent dispersion coef. at old time
real(stm_real), intent (in) :: time                          !< Current time
real(stm_real), intent (in) :: theta_stm                     !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
real(stm_real), intent (in) :: dt                            !< Time step   
real(stm_real), intent (in) :: dx                            !< Spacial step 

! ---- locals

real(stm_real) :: explicit_diffuse_op(ncell,nvar)
real(stm_real) :: down_diag(ncell)                            !< Values of the coefficients below diagonal in matrix
real(stm_real) :: center_diag(ncell)                          !< Values of the coefficients at the diagonal in matrix
real(stm_real) :: up_diag(ncell)                              !< Values of the coefficients above the diagonal in matrix
real(stm_real) :: diffusive_flux_boundary(ncell,nvar)         !< Explicit diffusive boundary flux
real(stm_real) :: diffusive_flux_interior(ncell,nvar)         !< Explicit diffusive interior flux
real(stm_real) :: right_hand_side(ncell,nvar)         

!todo: remove this
mass=0
conc=0

! This routine gives the effects of diffusion fluxes on each cell
! for a single time step (ie, explicit). This is needed for the advection step.
! It is also probably part of the right hand side of the implicit diffusion solver 
! matrix calculation. 

! Explicit diffusion operator construction. This creates the diffusive fluxes
! sends them for modification for boundaries and then differences the fluxes
! to get the operator d/dx(Ad/dx). 
call explicit_diffusion_operator(explicit_diffuse_op,     &
                                        conc_prev,        &
                                        area_lo_prev,     &
                                        area_hi_prev,     &
                                        disp_coef_lo_prev,&  
                                        disp_coef_hi_prev,&
                                        ncell,            &
                                        nvar,             &
                                        time,             &
                                        dx)
                                                                  


call construct_right_hand_side( right_hand_side,   & 
                                  explicit_diffuse_op,   & 
                                  area_prev,             &
                                  area_lo_prev,          &
                                  area_hi_prev,          &
                                  disp_coef_lo_prev,     &
                                  disp_coef_hi_prev,     &
                                  conc_prev,             &
                                  theta_stm,             &
                                  ncell,                 &
                                  time,                  &
                                  nvar,                  &  
                                  dx,                    &
                                  dt)
                                        
! Construct the matrix for the diffusion solver
! without boundary condition modification or structure on interior of domain
call construct_diffusion_matrix( center_diag ,      &
                                  up_diag,          &     
                                  down_diag,        &
                                  area,             &
                                  area_lo,          &
                                  area_hi,          &
                                  conc,             &
                                  conc_prev,        &
                                  disp_coef_lo,     &
                                  disp_coef_hi,     &
                                  theta_stm,        &
                                  ncell,            &
                                  time,             & 
                                  nvar,             & 
                                  dx,               &
                                  dt)

! this function will add boundary conditions to the matrix


call solve ( center_diag ,          &
                  up_diag,          &     
                  down_diag,        &
                  right_hand_side,  &
                  conc,             &
                  ncell,            &
                  nvar)


return
end subroutine diffuse


!> Calculate the explicit diffusion operator.
!> Explicit diffusion operator construction. This creates the diffusive fluxes
!> sends them for modification for boundaries and then differences the fluxes
!> to get the operator d/dx(Ad/dx). 
subroutine explicit_diffusion_operator (explicit_diffuse_op,  &
                                            conc_prev,        &
                                            area_lo_prev,     &
                                            area_hi_prev,     &
                                            disp_coef_lo_prev,&  
                                            disp_coef_hi_prev,&
                                            ncell,            &
                                            nvar,             &
                                            time,             &
                                            dx)
                                                                                          
use stm_precision

implicit none 

!----args

integer, intent (in) :: ncell !< Number of cells
integer, intent (in) :: nvar  !< Number of variables

real(stm_real), intent (out) :: explicit_diffuse_op(ncell,nvar)             !< Explicit diffusion operator
real(stm_real), intent (in)  :: conc_prev(ncell,nvar)                       !< Concentration at old time
real(stm_real), intent (in)  :: area_lo_prev (ncell)                        !< Low side area at old time
real(stm_real), intent (in)  :: area_hi_prev (ncell)                        !< High side area at old time 
real(stm_real), intent (in)  :: disp_coef_lo_prev (ncell,nvar)              !< Low side constituent dispersion coef. at old time
real(stm_real), intent (in)  :: disp_coef_hi_prev (ncell,nvar)              !< High side constituent dispersion coef. at old time
real(stm_real), intent (in)  :: time                                        !< Current time
real(stm_real), intent (in)  :: dx                                          !< Spacial step  

!--- locals
integer :: ivar
integer :: icell
real(stm_real):: diffusive_flux_interior (ncell,nvar) 
real(stm_real):: diffusive_flux_boundary (ncell,nvar) 

!todo: remove this
explicit_diffuse_op =0

call interior_diffusive_flux ( diffusive_flux_interior,       &
                                            conc_prev,        &
                                            area_lo_prev,     &
                                            area_hi_prev,     &
                                            disp_coef_lo_prev,&  
                                            disp_coef_hi_prev,&
                                            ncell,            &
                                            nvar,             &
                                            time,             &
                                            dx)
                                                        
call boundary_diffusive_flux ( diffusive_flux_boundary,       &
                                            conc_prev,        &
                                            area_lo_prev,     &
                                            area_hi_prev,     &
                                            disp_coef_lo_prev,&  
                                            disp_coef_hi_prev,&
                                            ncell,            &
                                            nvar,             &
                                            time,             &
                                            dx)

        
!do ivar = 1, nvar 
!    do icell = 2,ncell-1
!    
!       
!        explicit_diffusion_term(icell,ivar) =  (area_hi(icell,ivar)*disp_coef_hi(icell,ivar)* conc_prev(icell+1,ivar) &
!                                              -  area_hi(icell,ivar)*disp_coef_hi(icell,ivar)*conc_prev(icell,ivar) &
!                                              -  area_lo(icell,ivar)*disp_coef_lo(icell,ivar)*conc_prev(icell,ivar)   &
!                                              +  area_lo(icell,ivar)*disp_coef_lo(icell,ivar)*conc_prev(icell-1,ivar) )/dx/dx
!                
!               
!           
!    end do
!
!
!  ! ----  boundary conditions on on lo side of domain
!  ! todo: make dtbydx2 just dx and separate into diffuse_flux_hi and diffuse_flux_lo
!  ! then call boundary_diffusive_flux
!  ! then difference the fluxes one more time
!        explicit_diffusion_term(1,ivar) = (area_hi(1,ivar)*disp_coef_hi(1,ivar)* conc_prev(2,ivar) &
!                                            -  area_hi(1,ivar)*disp_coef_hi(1,ivar)*conc_prev(1,ivar) &
!                                            -  area_lo(1,ivar)*disp_coef_lo(1,ivar)*conc_prev(1,ivar) &
!                                            +  area_lo(1,ivar)*disp_coef_lo(1,ivar)*(conc_prev(2,ivar)- two * dx * bc_start_flux ) )/dx/dx
!        
!   
!        
!        
!  ! ----  boundary conditions on on hi side of domain
!        explicit_diffusion_term(ncell,ivar) =  (area_hi(ncell,ivar)*disp_coef_hi(ncell,ivar)* (conc_prev(ncell-1,ivar)+ two*dx*bc_end_flux) &
!                                                -  area_hi(ncell,ivar)*disp_coef_hi(ncell,ivar)*conc_prev(ncell,ivar)  &
!                                                -  area_lo(ncell,ivar)*disp_coef_lo(ncell,ivar)*conc_prev(ncell,ivar) &
!                                                +  area_lo(ncell,ivar)*disp_coef_lo(ncell,ivar)*conc_prev(icell-1,ivar) )/dx/dx
!        
!        
!    end do                          
                                  
         
return
end subroutine explicit_diffusion_operator 


subroutine interior_diffusive_flux (diffusive_flux_interior,  &
                                            conc_prev,        &
                                            area_lo_prev,     &
                                            area_hi_prev,     &
                                            disp_coef_lo_prev,&  
                                            disp_coef_hi_prev,&
                                            ncell,            &
                                            nvar,             &
                                            time,             &
                                            dx)

use stm_precision
! --- args
                                            
integer, intent (in) :: ncell !< Number of cells
integer, intent (in) :: nvar  !< Number of variables

real(stm_real), intent (out) :: diffusive_flux_interior(ncell,nvar)         !< Explicit diffusive flux
real(stm_real), intent (in)  :: conc_prev(ncell,nvar)                       !< Concentration at old time
real(stm_real), intent (in)  :: area_lo_prev (ncell)                        !< Low side area at old time
real(stm_real), intent (in)  :: area_hi_prev (ncell)                        !< High side area at old time 
real(stm_real), intent (in)  :: disp_coef_lo_prev (ncell,nvar)              !< Low side constituent dispersion coef. at old time
real(stm_real), intent (in)  :: disp_coef_hi_prev (ncell,nvar)              !< High side constituent dispersion coef. at old time
real(stm_real), intent (in)  :: time                                        !< Current time
real(stm_real), intent (in)  :: dx                                          !< Spatial step   

!todo: remove this
diffusive_flux_interior=0
                                            
return
end subroutine interior_diffusive_flux
                                                      
subroutine boundary_diffusive_flux(diffusive_flux_boundary,       &
                                                conc_prev,        &
                                                area_lo_prev,     &
                                                area_hi_prev,     &
                                                disp_coef_lo_prev,&  
                                                disp_coef_hi_prev,&
                                                ncell,            &
                                                nvar,             &
                                                time,             &
                                                dx)
use stm_precision

 ! --- args
                                            
integer, intent (in) :: ncell !< Number of cells
integer, intent (in) :: nvar  !< Number of variables

real(stm_real), intent (out) :: diffusive_flux_boundary(ncell,nvar)         !< Explicit diffusive boundary flux
real(stm_real), intent (in)  :: conc_prev(ncell,nvar)                       !< Concentration at old time
real(stm_real), intent (in)  :: area_lo_prev (ncell)                        !< Low side area at old time
real(stm_real), intent (in)  :: area_hi_prev (ncell)                        !< High side area at old time 
real(stm_real), intent (in)  :: disp_coef_lo_prev (ncell,nvar)              !< Low side constituent dispersion coef. at old time
real(stm_real), intent (in)  :: disp_coef_hi_prev (ncell,nvar)              !< High side constituent dispersion coef. at old time
real(stm_real), intent (in)  :: time                                        !< Current time
real(stm_real), intent (in)  :: dx                                          !< Spatial step   
                                                                                           
!todo: remove this
diffusive_flux_boundary = 0                                               
                                                
return
end subroutine boundary_diffusive_flux



!/////////////////////////////////////////////////
!> Construct the right hand side vector from previous step,
!> and impose Neumann boundary condition on it.
pure subroutine construct_right_hand_side( right_hand_side,   & 
                                  explicit_diffuse_op,   & 
                                  area_prev,             &
                                  area_lo_prev,          &
                                  area_hi_prev,          &
                                  disp_coef_lo_prev,     &
                                  disp_coef_hi_prev,     &
                                  conc_prev,             &
                                  theta_stm,             &
                                  ncell,                 &
                                  time,                  &
                                  nvar,                  &  
                                  dx,                    &
                                  dt)
use stm_precision   
  ! ---args                                
integer, intent (in) :: ncell !< Number of cells
integer, intent (in) :: nvar  !< Number of variables

real(stm_real), intent (out) :: right_hand_side(ncell,nvar)                 !< The right hand side vector
real(stm_real), intent (in)  :: explicit_diffuse_op (ncell,nvar)            !< Explicit diffusion operator
real(stm_real), intent (in)  :: area_prev (ncell)                           !< Cell centered area at old time 
real(stm_real), intent (in)  :: conc_prev(ncell,nvar)                       !< Concentration at old time
real(stm_real), intent (in)  :: area_lo_prev (ncell)                        !< Low side area at old time
real(stm_real), intent (in)  :: area_hi_prev (ncell)                        !< High side area at old time 
real(stm_real), intent (in)  :: disp_coef_lo_prev (ncell,nvar)              !< Low side constituent dispersion coef. at old time
real(stm_real), intent (in)  :: disp_coef_hi_prev (ncell,nvar)              !< High side constituent dispersion coef. at old time
real(stm_real), intent (in)  :: time                                        !< Current time
real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
real(stm_real), intent (in)  :: dx                                          !< Spatial step  
real(stm_real), intent (in)  :: dt                                          !< Time step                                   
                                  
!todo: remove this
right_hand_side=0                                  
                                  
return

end subroutine construct_right_hand_side

!/////////////////////////////////////////////////

subroutine construct_diffusion_matrix( center_diag ,      &
                                  up_diag,          &     
                                  down_diag,        &
                                  area,             &
                                  area_lo,          &
                                  area_hi,          &
                                  conc,             &
                                  conc_prev,        &
                                  disp_coef_lo,     &
                                  disp_coef_hi,     &
                                  theta_stm,        &
                                  ncell,            &
                                  time,             & 
                                  nvar,             & 
                                  dx,               &
                                  dt)

use stm_precision
                                  
 ! ---args                                
integer, intent (in) :: ncell !< Number of cells
integer, intent (in) :: nvar  !< Number of variables

real(stm_real),intent (out)  :: down_diag(ncell)                            !< Values of the coefficients below diagonal in matrix
real(stm_real),intent (out)  :: center_diag(ncell)                          !< Values of the coefficients at the diagonal in matrix
real(stm_real),intent (out)  :: up_diag(ncell)                              !< Values of the coefficients above the diagonal in matrix
real(stm_real), intent (in)  :: area (ncell)                           !< Cell centered area at new time 
real(stm_real), intent (in)  :: conc(ncell,nvar)                            !< Concentration at new time
real(stm_real), intent (in)  :: conc_prev(ncell,nvar)                       !< Concentration at old time
real(stm_real), intent (in)  :: area_lo(ncell)                   !< Low side area at new time
real(stm_real), intent (in)  :: area_hi(ncell)                   !< High side area at new time 
real(stm_real), intent (in)  :: disp_coef_lo (ncell,nvar)                   !< Low side constituent dispersion coef. at new time
real(stm_real), intent (in)  :: disp_coef_hi (ncell,nvar)                   !< High side constituent dispersion coef. at new time
real(stm_real), intent (in)  :: time                                        !< Current time
real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
real(stm_real), intent (in)  :: dx                                          !< Spatial step  
real(stm_real), intent (in)  :: dt                                          !< Time step                                   
                                  
!---local                                  
real(stm_real) :: explicit_diffuse_op (ncell,nvar)            !< Explicit diffusion operator
                                  
                                  
!todo: remove this                                  
center_diag=0
up_diag=0    
down_diag=0                                 
                                  
                                  
                                  
                                  
                                  
                                  
return
end subroutine construct_diffusion_matrix
!/////////////////////////////////////////////////
!> Solve the system of linear equations
subroutine solve ( center_diag ,           &
                          up_diag,              &     
                          down_diag,            &
                          right_hand_side,      &
                          conc,                 &
                          ncell,                &
                          nvar)

use matrix_solver
                                                        
! ----- args
use stm_precision
integer,intent (in) :: ncell                          !< Number of volumes
integer, intent (in):: nvar                           !< Number of variables 

real(stm_real),intent (in)  :: down_diag(ncell)       !< Values of the coefficients below diagonal in matrix
real(stm_real),intent (in)  :: center_diag(ncell)     !< Values of the coefficients at the diagonal in matrix
real(stm_real),intent (in)  :: up_diag(ncell)         !< Values of the coefficients above the diagonal in matrix
real(stm_real),intent (in)  :: right_hand_side(ncell) !< Values of the right hand side vector
real(stm_real),intent (out) :: conc(ncell)            !< Values of the computed solution


call tridi_solver ( center_diag ,               &
                          up_diag,              &     
                          down_diag,            &
                          right_hand_side,      &
                          conc,                 &
                          ncell,                &
                          nvar)


return
end subroutine solve

end module diffusion    