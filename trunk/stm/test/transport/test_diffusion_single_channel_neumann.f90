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

!> todo: write tests for explicit operator and implicit diffusion
!>@ingroup test
module test_diffusion_single_channel_neumann

use fruit
use stm_precision
use primitive_variable_conversion
use diffusion
use boundary_diffusion
use boundary_diffusion_matrix_module
use example_initial_conditions

contains





subroutine test_new_diffusion_calc

integer,parameter :: ncell = 1001                              !< Number of cells
integer,parameter :: nvar = 1                                  !< Number of variables

real(stm_real) :: conc(ncell,nvar)              !< Concentration at new time
real(stm_real) :: mass(ncell,nvar)              !< Mass (A*C) at new time
real(stm_real) :: mass_prev(ncell,nvar)         !< Mass (A*C) at old time
real(stm_real) :: conc_prev(ncell,nvar)         !< Concentration at old time
real(stm_real) :: area (ncell)                  !< Cell-centered area at new time
real(stm_real) :: area_prev (ncell)             !< Cell-centered area at old time
real(stm_real) :: area_lo (ncell)               !< Low side area centered in time
real(stm_real) :: area_hi (ncell)               !< High side area centered in time 
real(stm_real) :: area_lo_prev (ncell)          !< Low side area centered at old time
real(stm_real) :: area_hi_prev (ncell)          !< High side area centered at old time 
real(stm_real) :: disp_coef_lo (ncell,nvar)     !< Low side constituent dispersion coef. at new time
real(stm_real) :: disp_coef_hi (ncell,nvar)     !< High side constituent dispersion coef. at new time
real(stm_real) :: disp_coef_lo_prev(ncell,nvar) !< Low side constituent dispersion coef. at old time
real(stm_real) :: disp_coef_hi_prev(ncell,nvar) !< High side constituent dispersion coef. at old time
real(stm_real) :: time                          !< Current time
real(stm_real) :: theta_stm                     !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
real(stm_real) :: dt                            !< Time step   
real(stm_real) :: dx                            !< Spacial step 


!--- locals
integer :: iivar
integer :: jjvar
real(stm_real) :: xpos(ncell)
! todo: remove this
real(stm_real) :: dummy_higher
real(stm_real) :: dummy_lower
procedure(boundary_diffusive_matrix_if),pointer :: boundary_diffusion_matrix  
procedure(boundary_diffusive_flux_if),pointer :: boundary_diffusion_flux  

boundary_diffusion_matrix  => single_channel_neumann_matrix
boundary_diffusion_flux  => channel_neumann_gaussian_diffusive_flux
 
! ---- these will remain same in the process
time = LARGEREAL
dt = 0.001d0
dx = 0.05d0 * (1000d0/(ncell-1))
theta_stam = 0.6d0

area (:)= 1.0d0                 
area_prev (:) = 1.0d0            
area_lo (:)= 1.0d0              
area_hi (:)= 1.0d0              
area_lo_prev (:)= 1.0d0         
area_hi_prev (:)= 1.0d0         
disp_coef_lo (:,:) = 0.1d0   
disp_coef_hi (:,:) = 0.1d0   
disp_coef_lo_prev(:,:) = 0.1d0 
disp_coef_hi_prev(:,:) = 0.1d0 





!---- t initial is t=0 sec 

call fill_gaussian(conc_prev,ncell, -25.0d0 ,dx,zero,one)
!do iivar = 1, ncell
! xpos(iivar) = -25.0d0 + dx* (iivar-1)
! conc_prev (iivar, nvar) = exp(-(xpos(iivar)**2.0d0)/(4.0d0*disp_coef_lo_prev(iivar,nvar)))
!end do



call prim2cons(mass_prev,conc_prev,area,ncell,nvar)



!! add use boundary_diffusion
!! write neumann_diffusion_bc_flux and neumann_diffusion_bc_matrix

  

timemarch: do jjvar = 1,1000


     call diffuse(conc,             &
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
                  dx)


    conc_prev(:,nvar) = conc(:,nvar)
    

          
    call prim2cons(mass_prev,conc_prev,area,ncell,nvar)
   
    
end do timemarch

open (3,file="results.txt")

do ivar=1,ncell
        if  ((-5.0d0 - 1d-8 < xpos(ivar)) .and. (xpos(ivar) < 5d0 + 1d-8 )) then
        write (3,*) xpos(ivar),conc(ivar,1)
    end if
end do

continue
!!!
!!!!! todo: remove these
!!!!print *, xpos(500) ,conc_prev(500,nvar)
!!!!print *, xpos(501) ,conc_prev(501,nvar)
!!!!print *, xpos(502) ,conc_prev(502,nvar)
!!!!print *, epsilon   
!!!!pause
!!!   
!!!  ! ---check symmetry in solutions
!!!  call assertEquals(conc_prev((ncell-1)/2 +1 + 5,nvar),conc_prev((ncell-1)/2 + 1 - 5,nvar),1d-9,"Diffusion solution is not symmetric! theta=0.6")
!!!!  call assertEquals(conc_prev((ncell-1)/2 +1 + 50,nvar),conc_prev((ncell-1)/2 + 1 - 50,nvar),1d-9,"Diffusion solution is not symmetric! theta=0.6")
!!!!  call assertEquals(conc_prev((ncell-1)/2 +1 + 250,nvar),conc_prev((ncell-1)/2 + 1 - 250,nvar),1d-9,"Diffusion solution is not symmetric!theta=0.6")
!!!!  
!!!!    !----- check with exact solution 
!!!!  call assertEquals(conc_prev(501 ,nvar),0.707106781d0,1d-8,"Diffusion solution is not same as exact! theta=0.6")
!!!!  call assertEquals(conc_prev(551,nvar),0.148217633d0,1d-8,"Diffusion solution is not same as exact! theta=0.6")
!!!!  call assertEquals(conc_prev(601,nvar),0.001365037d0,1d-8,"Diffusion solution is not same as exact! theta=0.6")
!!!  
!!!  
!!!  !--- test for theta = 1
!!!  theta_stm = 1d0
!!!  
!!!  !---- t initial is t=1 sec 
!!!
!!!do iivar = 1, ncell
!!! xpos(iivar) = (iivar -1 - (ncell-1)/2)*dx
!!! conc_prev (iivar, nvar) = exp(-(xpos(iivar)**2)/4.0d0/disp_coef_lo_prev(iivar,nvar))
!!!
!!!end do
!!!
!!!call prim2cons(mass_prev,conc_prev,area,ncell,nvar)
!!!
!!!! !todo: remove these 
!!!!print *, xpos(4) ,conc_prev(4,nvar)
!!!!print *, xpos(501) ,conc_prev(501,nvar)
!!!!print *, xpos(998) ,conc_prev(998,nvar)   
!!!!pause
!!!
!!!!---- march
!!!
!!!timemarch1: do jjvar = 1,1000
!!!
!!!   !xmarch: do iivar = 1,ncell ! do I need this? I dont think so 
!!!
!!!     call diffuse(conc,             &
!!!                  conc_prev,         &
!!!                  mass,              &
!!!                  mass_prev,         &
!!!                  area,              &
!!!                  area_prev,         &
!!!                  area_lo,           &
!!!                  area_hi,           &
!!!                  area_lo_prev,      &
!!!                  area_hi_prev,      &
!!!                  diffusive_flux_boundary_lo, &
!!!                  diffusive_flux_boundary_hi, &
!!!                  disp_coef_lo,      &  
!!!                  disp_coef_hi,      &
!!!                  disp_coef_lo_prev, &  
!!!                  disp_coef_hi_prev, &
!!!                  ncell,             &
!!!                  nvar,              &
!!!                  time,              &
!!!                  theta_stm,         &
!!!                  dt,                &
!!!                  dx                 )
!!!
!!!    conc_prev(:,nvar) = conc(:,nvar)
!!!    
!!!
!!!          
!!!    call prim2cons(mass_prev,conc_prev,area,ncell,nvar)
!!!   
!!!    
!!!end do timemarch1
!!!
!!!!! todo: remove these
!!!!print *, xpos(500) ,conc_prev(500,nvar)
!!!!print *, xpos(501) ,conc_prev(501,nvar)
!!!!print *, xpos(502) ,conc_prev(502,nvar)
!!!!print *, epsilon   
!!!!pause
!!!   
!!!  ! ---check symmetry in solutions
!!!  call assertEquals(conc_prev((ncell-1)/2 +1 + 5,nvar),conc_prev((ncell-1)/2 + 1 - 5,nvar),1d-9,"Diffusion solution is not symmetric! theta=1")
!!!!  call assertEquals(conc_prev((ncell-1)/2 +1 + 50,nvar),conc_prev((ncell-1)/2 + 1 - 50,nvar),1d-9,"Diffusion solution is not symmetric! theta=1")
!!!!  call assertEquals(conc_prev((ncell-1)/2 +1 + 250,nvar),conc_prev((ncell-1)/2 + 1 - 250,nvar),1d-9,"Diffusion solution is not symmetric!theta=1")
!!!!  
!!!!    !----- check with exact solution 
!!!!  call assertEquals(conc_prev(501 ,nvar),0.707106781d0,1d-8,"Diffusion solution is not same as exact! theta=1")
!!!!  call assertEquals(conc_prev(551,nvar),0.148217633d0,1d-8,"Diffusion solution is not same as exact! theta=1")
!!!!  call assertEquals(conc_prev(601,nvar),0.001365037d0,1d-8,"Diffusion solution is not same as exact! theta=1")
!!!  
!!!  
!!!  !--- test for theta = 0.0
!!!  theta_stm = 0.0d0
!!!  
!!!  !---- t initial is t=1 sec 
!!!
!!!do iivar = 1, ncell
!!! xpos(iivar) = (iivar -1 - (ncell-1)/2)*dx
!!! conc_prev (iivar, nvar) = exp(-(xpos(iivar)**2)/4.0d0/disp_coef_lo_prev(iivar,nvar))
!!!
!!!end do
!!!
!!!call prim2cons(mass_prev,conc_prev,area,ncell,nvar)
!!!
!!!! !todo: remove these 
!!!!print *, xpos(4) ,conc_prev(4,nvar)
!!!!print *, xpos(501) ,conc_prev(501,nvar)
!!!!print *, xpos(998) ,conc_prev(998,nvar)   
!!!!pause
!!!
!!!!---- march
!!!
!!!timemarch2: do jjvar = 1,1000
!!!
!!!   !xmarch: do iivar = 1,ncell ! do I need this? I dont think so 
!!!
!!!     call diffuse(conc,             &
!!!                  conc_prev,         &
!!!                  mass,              &
!!!                  mass_prev,         &
!!!                  area,              &
!!!                  area_prev,         &
!!!                  area_lo,           &
!!!                  area_hi,           &
!!!                  area_lo_prev,      &
!!!                  area_hi_prev,      &
!!!                  diffusive_flux_boundary_lo, &
!!!                  diffusive_flux_boundary_hi, &
!!!                  disp_coef_lo,      &  
!!!                  disp_coef_hi,      &
!!!                  disp_coef_lo_prev, &  
!!!                  disp_coef_hi_prev, &
!!!                  ncell,             &
!!!                  nvar,              &
!!!                  time,              &
!!!                  theta_stm,         &
!!!                  dt,                &
!!!                  dx                 )
!!!
!!!    conc_prev(:,nvar) = conc(:,nvar)
!!!    
!!!
!!!          
!!!    call prim2cons(mass_prev,conc_prev,area,ncell,nvar)
!!!   
!!!    
!!!end do timemarch2
!!!
!!!!! todo: remove these
!!!!print *, xpos(500) ,conc_prev(500,nvar)
!!!!print *, xpos(501) ,conc_prev(501,nvar)
!!!!print *, xpos(502) ,conc_prev(502,nvar)
!!!!print *, epsilon   
!!!!pause
!!!   
!!!  ! ---check symmetry in solutions
!!!  call assertEquals(conc_prev((ncell-1)/2 +1 + 5,nvar),conc_prev((ncell-1)/2 + 1 - 5,nvar),1d-9,"Diffusion solution is not symmetric! theta=0")
!!!!  call assertEquals(conc_prev((ncell-1)/2 +1 + 50,nvar),conc_prev((ncell-1)/2 + 1 - 50,nvar),1d-9,"Diffusion solution is not symmetric! theta=0")
!!!!  call assertEquals(conc_prev((ncell-1)/2 +1 + 250,nvar),conc_prev((ncell-1)/2 + 1 - 250,nvar),1d-9,"Diffusion solution is not symmetric!theta=0")
!!!!  
!!!!    !----- check with exact solution 
!!!!  call assertEquals(conc_prev(501 ,nvar),0.707106781d0,1d-8,"Diffusion solution is not same as exact! theta=0")
!!!!  call assertEquals(conc_prev(551,nvar),0.148217633d0,1d-8,"Diffusion solution is not same as exact! theta=0")
!!!!  call assertEquals(conc_prev(601,nvar),0.001365037d0,1d-8,"Diffusion solution is not same as exact! theta=0")
!!!! 
!!!!todo remove
!!!

  return
end subroutine test_new_diffusion_calc


! may be fine to put this in diffusion
subroutine channel_neumann_gaussian_diffusive_flux (diffusive_flux_lo, &
                                                     diffusive_flux_hi, &
                                                     conc,              &
                                                     ncell,             &
                                                     nvar,              &
                                                     time)
                                                     
 use stm_precision
 implicit none
 !--- args
 integer,intent(in)  :: ncell                                   !< number of cells
 integer,intent(in)  :: nvar                                    !< number of variables
 real(stm_real),intent(inout) :: diffusive_flux_lo(ncell,nvar)  !< face flux, lo side
 real(stm_real),intent(inout) :: diffusive_flux_hi(ncell,nvar)  !< face flux, hi side
 real(stm_real),intent(in)  :: time                             !< time
 real(stm_real),intent(in)  :: conc(ncell,nvar)                 !< concentration 
       
                                                     
                                                     
                                                     
                                                     
                                                     
                                                     
! must follow interface
return
end subroutine



 
subroutine single_channel_neumann_matrix (center_diag ,      &
                                                      up_diag,          &     
                                                      down_diag,        &
                                                      area,             &
                                                      area_lo,          &
                                                      area_hi,          &          
                                                      disp_coef_lo,     &
                                                      disp_coef_hi,     &
                                                      theta_stm,        &
                                                      ncell,            &
                                                      time,             & 
                                                      nvar,             & 
                                                      dx,               &
                                                      dt)
                                                      
                                              
 use stm_precision
 implicit none
 !--- args

                                       
        integer, intent (in) :: ncell                                                 !< Number of cells
        integer, intent (in) :: nvar                                                  !< Number of variables

        real(stm_real),intent (inout)  :: down_diag(ncell,nvar)                            !< Values of the coefficients below diagonal in matrix
        real(stm_real),intent (inout)  :: center_diag(ncell,nvar)                          !< Values of the coefficients at the diagonal in matrix
        real(stm_real),intent (inout)  :: up_diag(ncell,nvar)                              !< Values of the coefficients above the diagonal in matrix
        real(stm_real), intent (in)  :: area (ncell)                                !< Cell centered area at new time 
        real(stm_real), intent (in)  :: area_lo(ncell)                              !< Low side area at new time
        real(stm_real), intent (in)  :: area_hi(ncell)                              !< High side area at new time 
        real(stm_real), intent (in)  :: disp_coef_lo (ncell,nvar)                   !< Low side constituent dispersion coef. at new time
        real(stm_real), intent (in)  :: disp_coef_hi (ncell,nvar)                   !< High side constituent dispersion coef. at new time
        real(stm_real), intent (in)  :: time                                        !< Current time
        real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
        real(stm_real), intent (in)  :: dx                                          !< Spatial step  
        real(stm_real), intent (in)  :: dt       


! must follow interface
return
end subroutine




end module
