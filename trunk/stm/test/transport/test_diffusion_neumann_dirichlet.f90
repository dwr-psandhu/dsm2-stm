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

!> todo: write tests for diffusion subroutine in case of Neumann BC at the left and Dirichlet BC at the right side 
!>@ingroup test
module test_diffusion_neumann_dirichlet

use diffusion
use fruit
use stm_precision
use primitive_variable_conversion

contains

subroutine test_diffusion_n_d
  use diffusion
  
  implicit none
  
integer,parameter :: ncell = 9                      !< Number of cells
integer,parameter :: nvar = 1                       !< Number of variables



real(stm_real) :: conc(ncell,nvar)                  !< Concentration at new time
real(stm_real) :: mass(ncell,nvar)                  !< Mass (A*C) at new time
real(stm_real) :: mass_prev(ncell,nvar)             !< Mass (A*C) at old time
real(stm_real) :: conc_prev(ncell,nvar)             !< Concentration at old time
real(stm_real) :: area (ncell)                      !< Cell-centered area at new time
real(stm_real) :: area_prev (ncell)                 !< Cell-centered area at old time
real(stm_real) :: area_lo (ncell)                   !< Low side area centered in time
real(stm_real) :: area_hi (ncell)                   !< High side area centered in time 
real(stm_real) :: area_lo_prev (ncell)              !< Low side area centered at old time
real(stm_real) :: area_hi_prev (ncell)              !< High side area centered at old time 
real(stm_real) :: disp_coef_lo (ncell,nvar)         !< Low side constituent dispersion coef. at new time
real(stm_real) :: disp_coef_hi (ncell,nvar)         !< High side constituent dispersion coef. at new time
real(stm_real) :: disp_coef_lo_prev(ncell,nvar)     !< Low side constituent dispersion coef. at old time
real(stm_real) :: disp_coef_hi_prev(ncell,nvar)     !< High side constituent dispersion coef. at old time
real(stm_real) :: time                              !< Current time
real(stm_real) :: theta_stm                         !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
real(stm_real) :: dt                                !< Time step   
real(stm_real) :: dx                                !< Spacial step 
real(stm_real) :: diffusive_flux_boundary_lo(nvar)  !< Neumann BC on low side    
real(stm_real) :: diffusive_flux_boundary_hi (nvar) !< Neumann BC on high side

!---- locals
integer :: jvar
integer :: ivar
real(stm_real) :: xpos


 !--- Small numbers on center diag large numbers on up and down diag
  
time = zero
dt = 0.001d0
dx = 0.9d0 /ncell
theta_stm = 0.5d0

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

! this analitical solution only valid for disp_coef =constant
! todo: these must remove
!diffusive_flux_boundary_lo(nvar) = two - two * pi* sin(0.05d0*pi)
!diffusive_flux_boundary_hi(nvar) = two - two * pi* sin(0.5d0*pi)
! 
 !---initial condition
  do ivar=1,ncell
  
    xpos = 0.1d0+ (ivar-half)*dx
    conc_prev(ivar,nvar) = two*xpos +two*two*cos(pi*xpos/two) 
!    print* ,"IC-----",ivar,xpos,conc_prev(ivar,nvar)
  
  end do
  
!  print *, "C at B lo", two*0.05d0 +two*two*cos(pi*0.05d0/two)
!  print *, " slope at B lo",(conc_prev(1,1) -   two*0.05d0 +two*two*cos(pi*0.05d0/two))/dx
!  pause
!  
  
 call prim2cons(mass_prev,conc_prev,area,ncell,nvar)
 
do jvar=1,1 
 
    call diffuse(conc,               &
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
                  dx                 )

   
   time = (jvar)*dt
   conc_prev(:,nvar) = conc(:,nvar) 
           
call prim2cons(mass_prev,conc_prev,area,ncell,nvar)

   diffusive_flux_boundary_lo(nvar) = two -two*pi*sin(0.05d0*pi)*exp(-0.1d0*time*pi*pi/4) 
   diffusive_flux_boundary_hi (nvar) =  two - two *pi*sin(0.5d0*pi)*exp(-0.1d0*time*pi*pi/4)

!print *,diffusive_flux_boundary_lo



end do



do ivar=1,ncell
  
    xpos = 0.1d0+ (ivar-half)*dx
    conc_prev(ivar,nvar) = two*xpos +two*two*cos(pi*xpos/two)*exp(-0.1d0*time*pi*pi/4) 
      
end do

!todo fill this part
!conc_exact

!  call assertEquals (
!  call assertEquals (
!  call assertEquals (


return
end subroutine test_diffusion_n_d

end module