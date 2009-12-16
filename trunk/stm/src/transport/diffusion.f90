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

use stm_precision
use primitive_variable_conversion

contains

! ///////////////////////////////////////////////////////////////////

!> Calculates the diffusive portion of the constituent transport.
!> It contains an explicit version of the diffusion operator and a general (involving all
!> potential cases) diffusion operator as well, with a coefficient theta_stm for 
!> selecting the level of implicitness. (theta_stm=0.5 is Crank Nicolson.).
!> The matrix is solved via a tri-diagonal solver.  
subroutine diffuse(conc,     &
                  conc_prev,&
                  mass,     &
                  mass_prev, &
                  area,     &
                  area_prev,&
                  area_lo,  &
                  area_hi,  &
                  ks_lo,    &
                  ks_hi,    &
                  ncell,    &
                  nvar,     &
                  time,     &
                  theta_stm,&
                  dt,       &
                  dx,       &
                  end_flx_bc_prev, &
                  strt_flx_bc_prev, &
                  explicit_diffusion_term )

use primitive_variable_conversion ! does it needed here 

implicit none

! ---- args

integer, intent (in) :: ncell !< Number of cells
integer, intent (in) :: nvar  !< Number of variables

real(stm_real), intent (out) :: conc(ncell,nvar)     !< Concentration at new time
real(stm_real), intent (out) :: explicit_diffusion_term(ncell,nvar)     !< Explicit diffusion predictor
real(stm_real), intent (out) :: mass(ncell,nvar)     !< mass (A*C) at new time
real(stm_real), intent (in) :: mass_prev(ncell,nvar) !< mass (A*C) at old time
real(stm_real), intent (in) :: conc_prev(ncell,nvar) !< Concentration at old time
real(stm_real), intent (in) :: area (ncell,nvar)     !< Cell-centered area at new time
real(stm_real), intent (in) :: area_prev (ncell,nvar)!< Cell-centered area at old time
real(stm_real), intent (in) :: area_lo (ncell,nvar)  !< Low side area centered intime
real(stm_real), intent (in) :: area_hi (ncell,nvar)  !< High side area centered in time 
real(stm_real), intent (in) :: ks_lo (ncell,nvar)    !< Low side sediment dispersion coef. centered in time
real(stm_real), intent (in) :: ks_hi (ncell,nvar)    !< High side sediment dispersion coef. centered in time
real(stm_real), intent (in) :: time                  !< Current time
real(stm_real), intent (in) :: theta_stm             !< Expilicitness coefficient 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
real(stm_real), intent (in) :: dt                    !< Time step   
real(stm_real), intent (in) :: dx                    !< Spacial step 
real(stm_real), intent (in) :: end_flx_bc_prev         !< flux at right hand side old time
real(stm_real), intent (in) :: strt_flx_bc_prev         !< flux at left hand side old time

! TODO: remove this
 conc = 0
 mass = 0


! This routine gives the effects of diffusion fluxes on each cell
! for a single time step (ie, explicit). This is needed for the advection step.
! It is also probably part of the right hand side of the implicit diffusion solver 
! matrix calculation. 



call explicit_diffusion_operator (explicit_diffusion_term, &
                                                        conc,             &
                                                        conc_prev,        &
                                                        mass,             &
                                                        mass_prev,        & 
                                                        area,             &      
                                                        area_lo,          &
                                                        area_hi,          &
                                                        ks_lo,            &
                                                        ks_hi,            &
                                                        ncell,            &
                                                        nvar,             &
                                                        time,             &
                                                        strt_flx_bc_prev, &   
                                                        end_flx_bc_prev,  &
                                                        dt,               &
                                                        dx)
! 
call construct_diffusion_matrix()

call construct_right_hand_side()

call apply_diffusion_boundary()

call solve


return
end subroutine diffuse



!> Calculate the explicit diffusion term.
subroutine explicit_diffusion_operator (explicit_diffusion_term, &
                                        conc,             &
                                        conc_prev,        &
                                        mass,             &
                                        mass_prev,        & 
                                        area,             &      
                                        area_lo,          &
                                        area_hi,          &
                                        ks_lo,            &
                                        ks_hi,            &
                                        ncell,            &
                                        nvar,             &
                                        time,             &
                                        strt_flx_bc_prev, &   
                                        end_flx_bc_prev,  &
                                        dt,               &
                                        dx)
                                                                                          
use stm_precision

implicit none 

!----args

integer, intent (in) :: ncell !< Number of cells
integer, intent (in) :: nvar  !< Number of variables

real(stm_real), intent (out) :: explicit_diffusion_term(ncell,nvar)     !< Explicit diffusion predictor
real(stm_real), intent (out) :: conc(ncell,nvar)      !< Concentration at new time
real(stm_real), intent (in)  :: conc_prev(ncell,nvar)  !< Concentration at old time
real(STM_REAL), intent (out) :: mass(ncell,nvar)      !< mass at new time
real(STM_REAL), intent (in)  :: mass_prev(ncell,nvar) !< mass at old time
real(stm_real), intent (in)  :: area (ncell,nvar)      !< Cell-centered area at new time
real(stm_real), intent (in)  :: area_lo (ncell,nvar)   !< Low side area at old time
real(stm_real), intent (in)  :: area_hi (ncell,nvar)   !< High side area  
real(stm_real), intent (in)  :: ks_lo (ncell,nvar)     !< Low side sediment dispersion coef. at old time
real(stm_real), intent (in)  :: ks_hi (ncell,nvar)     !< High side sediment dispersion coef. at old time
real(stm_real), intent (in)  :: time                   !< Current time
real(stm_real), intent (in)  :: end_flx_bc_prev        !< flux at right hand side old time
real(stm_real), intent (in)  :: strt_flx_bc_prev       !< flux at left hand side old time  
real(stm_real), intent (in)  :: dt                     !< Time step   
real(stm_real), intent (in)  :: dx                     !< Spacial step  

!--- locals
integer :: ivar
integer :: icell
real(stm_real) :: dtbydx2


dtbydx2 = dt/dx/dx

do ivar = 1, nvar 
    do icell = 2,ncell-1
    
        if ( ks_lo(icell,ivar)*dtbydx2> half .or. ks_hi(icell,ivar)*dtbydx2 > half) then 
            print *, 'Diffusion Unstable!'
            pause
        end if 
        
        explicit_diffusion_term(icell,ivar) = dtbydx2 * (area_hi(icell,ivar)*ks_hi(icell,ivar)* conc_prev(icell+1,ivar) &
                                              -  area_hi(icell,ivar)*ks_hi(icell,ivar)*conc_prev(icell,ivar) &
                                              -  area_lo(icell,ivar)*ks_lo(icell,ivar)*conc_prev(icell,ivar)   &
                                              +  area_lo(icell,ivar)*ks_lo(icell,ivar)*conc_prev(icell-1,ivar) )
                
        mass(icell,ivar) = mass_prev(icell,ivar) + explicit_diffusion_term(icell,ivar)
        
           
    end do


  ! ---- rhs boundary
        explicit_diffusion_term(1,ivar) = dtbydx2 *(area_hi(1,ivar)*ks_hi(1,ivar)* conc_prev(2,ivar) &
                                            -  area_hi(1,ivar)*ks_hi(1,ivar)*conc_prev(1,ivar) &
                                            -  area_lo(1,ivar)*ks_lo(1,ivar)*conc_prev(1,ivar) &
                                            +  area_lo(1,ivar)*ks_lo(1,ivar)*(conc_prev(2,ivar)- two * dx * strt_flx_bc_prev ) )
        
        mass(1,ivar) = mass_prev(1,ivar) + explicit_diffusion_term(1,ivar)
        
        
  ! --- lhs boundary
        explicit_diffusion_term(ncell,ivar) =  dtbydx2 * (area_hi(ncell,ivar)*ks_hi(ncell,ivar)* (conc_prev(ncell-1,ivar)+ two*dx*end_flx_bc_prev) &
                                                -  area_hi(ncell,ivar)*ks_hi(ncell,ivar)*conc_prev(ncell,ivar)  &
                                                -  area_lo(ncell,ivar)*ks_lo(ncell,ivar)*conc_prev(ncell,ivar) &
                                                +  area_lo(ncell,ivar)*ks_lo(ncell,ivar)*conc_prev(icell-1,ivar) )
               
        mass(ncell,ivar) = mass_prev(ncell,ivar) + explicit_diffusion_term(ncell,ivar)
        
    end do                          
                                  
         call cons2prim(conc,mass,area,ncell,nvar)                                
return
end subroutine explicit_diffusion_operator 


end module diffusion

!////////////////////////////


                                 

