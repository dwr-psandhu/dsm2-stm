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
module new_diffusion

use stm_type

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
subroutine New_diffuse(conc,             &
                  conc_prev,         &
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
                  dx,               &                 
                  bc_single_canal)

use stm_precision
use primitive_variable_conversion 


implicit none

! ---- args

integer, intent (in) :: ncell                                !< Number of cells
integer, intent (in) :: nvar                                 !< Number of variables

real(stm_real), intent (out) :: conc(ncell,nvar)             !< Concentration at new time
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
type (bc_type), intent (in) :: bc_single_canal               !< boundary condition of single channel -Type up- Value Up-Type Down- Value down

! ---- locals

real(stm_real) :: explicit_diffuse_op(ncell,nvar)
real(stm_real) :: down_diag(ncell,nvar)                            !< Values of the coefficients below diagonal in matrix
real(stm_real) :: center_diag(ncell,nvar)                          !< Values of the coefficients at the diagonal in matrix
real(stm_real) :: up_diag(ncell,nvar)                              !< Values of the coefficients above the diagonal in matrix
real(stm_real) :: diffusive_flux_boundary(ncell,nvar)         !< Explicit diffusive boundary flux
real(stm_real) :: diffusive_flux_interior(ncell,nvar)         !< Explicit diffusive interior flux
real(stm_real) :: right_hand_side(ncell,nvar)         



! This routine gives the effects of diffusion fluxes on each cell
! for a single time step (ie, explicit). This is needed for the advection step.
! It is also probably part of the right hand side of the implicit diffusion solver 
! matrix calculation. 

! Explicit diffusion operator construction. This creates the diffusive fluxes
! sends them for modification for boundaries and then differences the fluxes
! to get the operator d/dx(Ad/dx). 
call explicit_diffusion_operator(explicit_diffuse_op,  &
                                            conc_prev,        &
                                            area_lo_prev,     &
                                            area_hi_prev,     &
                                            disp_coef_lo_prev,&  
                                            disp_coef_hi_prev,&
                                            ncell,            &
                                            nvar,             &
                                            time,             &
                                            dx,                 &
                                            bc_single_canal)
                                                                  

! need to change this to use just diffusive_flux_hi/lo
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
                                  ! todo: check these 
!                                  conc,             &
!                                  conc_prev,        &
                                  disp_coef_lo,     &
                                  disp_coef_hi,     &
                                  theta_stm,        &
                                  ncell,            &
                                  time,             & 
                                  nvar,             & 
                                  dx,               &
                                  dt)
call
! todo: this function will add boundary conditions to the matrix
!call boundary_diffusion_matrix(center_diag...

call solve ( center_diag ,          &
                  up_diag,          &     
                  down_diag,        &
                  right_hand_side,  &
                  conc,             &
                  ncell,            &
                  nvar)



return
end subroutine New_diffuse


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
                                            dx,                 &
                                            bc_single_canal)
                                                                                          
use stm_precision
use stm_type
use boundary_diffusion

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
type (bc_type), intent (in) :: bc_single_canal               !< boundary condition of single channel -Type up- Value Up-Type Down- Value down
!--- locals
integer :: ivar
integer :: icell

real(stm_real):: diffusive_flux_lo(ncell,nvar)
real(stm_real):: diffusive_flux_hi(ncell,nvar)
real(stm_real):: diffusive_flux_lo_prev(ncell,nvar)
real(stm_real):: diffusive_flux_hi_prev(ncell,nvar)



do ivar = 1,nvar
    do icell = 2,ncell
        diffusive_flux_lo(icell,ivar) = (area_lo_prev(icell)*disp_coef_lo_prev(icell,ivar)* (conc_prev(icell,ivar)- conc_prev(icell-1,ivar)))/dx                       
    end do
end do 
diffusive_flux_hi(1:ncell-1,:) =  diffusive_flux_lo(2:ncell,:)  
! todo: 
diffusive_flux_hi(ncell,:)=bc_single_canal.hi_value
diffusive_flux_lo(1,:)=bc_single_canal.lo_value
  
  
  
  
  
                                  

   explicit_diffuse_op (1,:) = (diffusive_flux_hi(1,:) - diffusive_flux_lo(1,:) )/dx
   explicit_diffuse_op (ncell,:) = (diffusive_flux_hi(ncell,:) - diffusive_flux_lo(ncell,:) )/dx
         
    do ivar = 1,nvar
    ! to do : 1, ncell
        do icell = 1,ncell 
         
         explicit_diffuse_op (icell,ivar) = (diffusive_flux_hi (icell,ivar) - diffusive_flux_lo (icell,ivar))/dx
        
        end do
    end do
    
    
return
end subroutine explicit_diffusion_operator 


                                                      




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

  
  !---- locals
   integer :: ivar
   integer :: icell
   
 do ivar = 1,nvar 
    do icell = 1, ncell
    
         right_hand_side(icell,ivar) = area_prev(icell)*conc_prev(icell,ivar) &
                                       + (1-theta_stm)*dt * explicit_diffuse_op (icell,ivar) 
            
    end do      
    
 end do
 ! to do :todo : take care of BC of right hand side
! right_hand_side(1,ivar) = right_hand_side(1,ivar) - dt*theta_stm* diffusive_flux_boundary_lo(ivar)  /dx 
! right_hand_side(ncell,ivar) = right_hand_side(ncell,ivar) + dt*theta_stm* diffusive_flux_boundary_hi(ivar)  /dx
return
end subroutine construct_right_hand_side

!/////////////////////////////////////////////////

subroutine construct_diffusion_matrix( center_diag ,      &
                                  up_diag,          &     
                                  down_diag,        &
                                  area,             &
                                  area_lo,          &
                                  area_hi,          &
!                                  conc,             &
!                                  conc_prev,        &
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

real(stm_real),intent (out)  :: down_diag(ncell,nvar)                            !< Values of the coefficients below diagonal in matrix
real(stm_real),intent (out)  :: center_diag(ncell,nvar)                          !< Values of the coefficients at the diagonal in matrix
real(stm_real),intent (out)  :: up_diag(ncell,nvar)                              !< Values of the coefficients above the diagonal in matrix
real(stm_real), intent (in)  :: area (ncell)                                !< Cell centered area at new time 
!real(stm_real), intent (in)  :: conc(ncell,nvar)                            !< Concentration at new time
!real(stm_real), intent (in)  :: conc_prev(ncell,nvar)                       !< Concentration at old time
real(stm_real), intent (in)  :: area_lo(ncell)                              !< Low side area at new time
real(stm_real), intent (in)  :: area_hi(ncell)                              !< High side area at new time 
real(stm_real), intent (in)  :: disp_coef_lo (ncell,nvar)                   !< Low side constituent dispersion coef. at new time
real(stm_real), intent (in)  :: disp_coef_hi (ncell,nvar)                   !< High side constituent dispersion coef. at new time
real(stm_real), intent (in)  :: time                                        !< Current time
real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
real(stm_real), intent (in)  :: dx                                          !< Spatial step  
real(stm_real), intent (in)  :: dt                                          !< Time step                                   
                                  
!---local                                  
real(stm_real) :: explicit_diffuse_op (ncell,nvar)                          !< Explicit diffusion operator
real(stm_real) :: d_star

integer :: icell
integer :: ivar

d_star = dt/dx/dx  
    do ivar = 1,nvar 
        do icell = 2,ncell-1
         down_diag(icell,ivar) = - theta_stm*d_star*area_lo(icell)*disp_coef_lo(icell,ivar) 
         center_diag(icell,ivar) = area(icell) + theta_stm*d_star*(area_hi(icell)*disp_coef_hi(icell,ivar) + area_lo(icell)*disp_coef_lo(icell,ivar))
         up_diag(icell,ivar) = - theta_stm*d_star*area_hi(icell)*disp_coef_hi(icell,ivar)             
        end do                            
        center_diag(1,ivar) = area(1) + theta_stm*d_star*(area_hi(1)*disp_coef_hi(1,ivar) )
        up_diag(1,ivar) = - theta_stm*d_star*(area_hi(1)*disp_coef_hi(1,ivar)+ area_lo(1)*disp_coef_lo(1,ivar) )        
        down_diag(ncell,ivar)  =   - theta_stm*d_star*(area_hi(ncell)*disp_coef_hi(ncell,ivar)+ area_lo(ncell)*disp_coef_lo(ncell,ivar) )
        center_diag(ncell,ivar) = area(ncell) + theta_stm*d_star*(area_lo(ncell)*disp_coef_lo(ncell,ivar) )    
    end do   
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

! --- local
integer :: ivar
do ivar = 1 ,nvar
    call tridi_solver ( center_diag ,               &
                              up_diag,              &     
                              down_diag,            &
                              right_hand_side,      &
                              conc,                 &
                              ncell)
end do
return
end subroutine solve

end module new_diffusion    