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

!> Test of convergence of error norms in diffusion subroutine 
!>@ingroup test

module test_diffusion_norms

use diffusion 
use fruit
use stm_precision
use error_metric
use primitive_variable_conversion
use state_variables
use example_initial_conditions

contains

subroutine test_diffusion_error_norms
implicit none
integer :: ncell                                      !< Number of cells
integer,parameter :: nvar = 1                         !< Number of variables

real(stm_real),allocatable :: conc(:,:)               !< Concentration at new time
real(stm_real),allocatable :: mass(:,:)               !< Mass (A*C) at new time
real(stm_real),allocatable :: mass_prev(:,:)          !< Mass (A*C) at old time
real(stm_real),allocatable :: conc_prev(:,:)          !< Concentration at old time
real(stm_real),allocatable :: area (:)                !< Cell-centered area at new time
real(stm_real),allocatable :: area_prev (:)           !< Cell-centered area at old time
real(stm_real),allocatable :: area_lo (:)             !< Low side area centered in time
real(stm_real),allocatable :: area_hi (:)             !< High side area centered in time 
real(stm_real),allocatable :: area_lo_prev (:)        !< Low side area centered at old time
real(stm_real),allocatable :: area_hi_prev (:)        !< High side area centered at old time 
real(stm_real),allocatable :: disp_coef_lo (:,:)      !< Low side constituent dispersion coef. at new time
real(stm_real),allocatable :: disp_coef_hi (:,:)      !< High side constituent dispersion coef. at new time
real(stm_real),allocatable :: disp_coef_lo_prev(:,:)  !< Low side constituent dispersion coef. at old time
real(stm_real),allocatable :: disp_coef_hi_prev(:,:)  !< High side constituent dispersion coef. at old time
real(stm_real) :: time                                !< Current time
real(stm_real) :: theta_stm                           !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
real(stm_real) :: dt                                  !< Time step   
real(stm_real) :: dx                                  !< Spacial step 
real(stm_real) :: diffusive_flux_boundary_lo(nvar)    !< Neumann BC on low side    
real(stm_real) :: diffusive_flux_boundary_hi (nvar)   !< Neumann BC on high side


!--- locals
integer :: icell
integer :: itime
integer :: kkvar
integer,parameter :: n_err_try = 4
real(stm_real),allocatable :: xpos(:)
real(stm_real),allocatable :: conc_exact(:,:)
real(stm_real) :: norm_1(n_err_try)
real(stm_real) :: norm_2(n_err_try)
real(stm_real) :: norm_inf(n_err_try)

 

 do kkvar=1,n_err_try
 
dx = 0.2d0 / dble(2**(kkvar-1))
ncell = 250 * (2**(kkvar-1))+ 1

allocate (xpos(ncell)) 
allocate (conc_exact(ncell,nvar))
allocate (conc(ncell,nvar))             
allocate (mass(ncell,nvar))            
allocate (mass_prev(ncell,nvar))        
allocate (conc_prev(ncell,nvar))       
allocate (area (ncell))      
allocate (area_prev (ncell))     
allocate (area_lo (ncell))    
allocate (area_hi (ncell))   
allocate (area_lo_prev (ncell))  
allocate (area_hi_prev (ncell)) 
allocate (disp_coef_lo (ncell,nvar))
allocate (disp_coef_hi (ncell,nvar))
allocate (disp_coef_lo_prev(ncell,nvar))
allocate (disp_coef_hi_prev(ncell,nvar))

 
 
time = LARGEREAL
dt = 0.001d0
 
theta_stm = 0.7d0 

area (:)= 1.0d0                 
area_prev (:) = 1.0d0            
area_lo (:)= 1.0d0              
area_hi (:)= 1.0d0              
area_lo_prev (:)= 1.0d0         
area_hi_prev (:)= 1.0d0         
disp_coef_lo (:,:) = 0.5d0   
disp_coef_hi (:,:) = 0.5d0   
disp_coef_lo_prev(:,:) = 0.5d0 
disp_coef_hi_prev(:,:) = 0.5d0 
diffusive_flux_boundary_lo(nvar) = zero      
diffusive_flux_boundary_hi (nvar) = zero 
 
  
!---- t initial is t=1 sec 
! todo : re write with fill gaussian 
do icell = 1, ncell
 xpos(icell) = -25.0d0 + dx* (icell-1)
 conc_prev (icell, :) = exp(-(xpos(icell)**2.0d0)/4.0d0/disp_coef_lo_prev(icell,nvar))
end do

call prim2cons(mass_prev,conc_prev,area,ncell,nvar)

!---- march

timemarch: do itime = 1,1000

     call diffuse(conc,              &
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

    conc_prev(:,nvar) = conc(:,nvar)
        
    call prim2cons(mass_prev,conc_prev,area,ncell,nvar)
    
end do timemarch


!call fill_gaussian(vals,nloc,origin,dx,mean,sd)
!todo : check above works or not it is not exact gaussian 
do icell = 1, ncell
     xpos(icell) = -25.0d0 + dx* (icell-1)
     conc_exact (icell, nvar) = sqrt(0.5d0)* exp(-(xpos(icell)**2.0d0)/4.0d0/disp_coef_lo_prev(icell,nvar))
end do

call error_norm(norm_1(kkvar),norm_2(kkvar),norm_inf(kkvar),conc,conc_exact,ncell,dx)



deallocate(conc, conc_prev, mass,mass_prev)
deallocate(area, area_prev, area_lo, area_hi)

deallocate (area_lo_prev ,      &
            area_hi_prev ,      & 
            disp_coef_lo ,      &
            disp_coef_hi ,      &
            disp_coef_lo_prev,  &
            disp_coef_hi_prev,  &
            xpos,               &
            conc_exact)


end do 

return
end subroutine test_diffusion_error_norms

end module 
