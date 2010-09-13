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

!> boundary conditions for a single channel
!>@ingroup transport
module single_channel_boundary

use boundary_diffusion
use boundary_advection

!> Advectove fluxes on lo and hi side of single channel.
!> Each is expected to operate on only one end
procedure(boundary_advective_flux_if),  pointer :: advective_flux_lo  => null()
procedure(boundary_advective_flux_if),  pointer :: advective_flux_hi  => null()


!> Diffusive fluxes on lo and hi side of single channel.
!> Each is expected to operate on only one end
procedure(boundary_diffusive_flux_if),  pointer :: diffusion_flux_lo  => null()
procedure(boundary_diffusive_flux_if),  pointer :: diffusion_flux_hi  => null()


!> Diffusive matrix on lo and hi side of single channel.
!> Each is expected to operate on only one end
procedure(boundary_diffusive_matrix_if),pointer :: diffusion_matrix_lo  => null()
procedure(boundary_diffusive_matrix_if),pointer :: diffusion_matrix_hi  => null()
 
!> Provide data for boundary
abstract interface
!> Generic interface for boundary diffusion that should be fulfilled by client programs
subroutine boundary_data_if(bc_data,           &
                            xloc,              &
                            conc,              &
                            ncell,             &
                            nvar,              &
                            origin,            &
                            time,              &
                            dx,                &
                            dt)
    use stm_precision
    implicit none
    !--- args
    integer, intent(in)  :: ncell                                    !< Number of cells
    integer, intent(in)  :: nvar                                     !< Number of variables
    real(stm_real), intent(out)   :: bc_data(nvar)                   !< concentration or gradient data
    real(stm_real), intent(in)    :: xloc                            !< location where data is requested
    real(stm_real), intent (in)   :: time                            !< Time
    real(stm_real), intent (in)   :: origin                          !< Space origin
    real(stm_real), intent (in)   :: conc(ncell,nvar)                !< Concentration 
    real(stm_real), intent (in)   :: dt
    real(stm_real), intent (in)   :: dx
    
    end subroutine 
 end interface

!> User settable function that provides the data (concentration or derivative) used for the BC
procedure(boundary_data_if),pointer :: advection_data_lo  => null()
procedure(boundary_data_if),pointer :: advection_data_hi  => null()
procedure(boundary_data_if),pointer :: diffusion_data_lo  => null()
procedure(boundary_data_if),pointer :: diffusion_data_hi  => null()

contains

subroutine set_single_channel_boundary(advect_bc_lo, advect_data_lo, &
                                       advect_bc_hi, advect_data_hi, &
                                       diffuse_bc_lo, diffuse_data_lo, &
                                       diffuse_bc_hi, diffuse_data_hi)
   use error_handling
   implicit none

    procedure(boundary_advective_flux_if),  pointer :: advect_bc_lo
    procedure(boundary_advective_flux_if),  pointer :: advect_bc_hi
    procedure(boundary_diffusive_flux_if),  pointer :: diffuse_bc_lo
    procedure(boundary_diffusive_flux_if),  pointer :: diffuse_bc_hi

    procedure(boundary_data_if),pointer :: advect_data_lo
    procedure(boundary_data_if),pointer :: advect_data_hi
    procedure(boundary_data_if),pointer :: diffuse_data_lo
    procedure(boundary_data_if),pointer :: diffuse_data_hi

    advective_flux_lo => advect_bc_lo
    advective_flux_hi => advect_bc_hi
    diffusion_flux_lo  => diffuse_bc_lo
    diffusion_flux_hi  => diffuse_bc_hi

    advection_data_lo => advect_data_lo
    advection_data_hi => advect_data_hi
    diffusion_data_lo => diffuse_data_lo
    diffusion_data_hi => diffuse_data_hi

    if(associated(diffusion_flux_lo,dirichlet_diffusive_flux_lo))then
       diffusion_matrix_lo => dirichlet_diffusive_matrix_lo
    elseif(associated(diffusion_flux_lo,neumann_diffusive_flux_lo))then
       diffusion_matrix_lo => neumann_diffusive_matrix_lo
    else
       call stm_fatal("Unable to infer diffusion lo-side boundary condition (for matrix)")
    end if
    if(associated(diffusion_flux_hi,dirichlet_diffusive_flux_hi))then
       diffusion_matrix_hi => dirichlet_diffusive_matrix_hi
    elseif(associated(diffusion_flux_hi,neumann_diffusive_flux_hi))then
       diffusion_matrix_hi => neumann_diffusive_matrix_hi
    else
       call stm_fatal("Unable to infer diffusion hi-side boundary condition (for matrix)")
    end if


return
end subroutine

!> Advective flux boundary condition for a single channel, delagates to an implementation at each end
subroutine single_channel_boundary_advective_flux(flux_lo,    &
                                                  flux_hi,    &
                                                  conc_lo,    &
                                                  conc_hi,    &
                                                  flow_lo,    &
                                                  flow_hi,    &
                                                  ncell,      &
                                                  nvar,       &
                                                  time,       &
                                                  dt,         &
                                                  dx)
 
 use stm_precision
  
 implicit none
  !--- args          
 integer,intent(in)  :: ncell  !< Number of cells
 integer,intent(in)  :: nvar   !< Number of variables
 ! todo: check the intents
 real(stm_real),intent(inout) :: flux_lo(ncell,nvar)     !< flux on lo side of cell, time centered
 real(stm_real),intent(inout) :: flux_hi(ncell,nvar)     !< flux on hi side of cell, time centered
 real(stm_real),intent(in)    :: flow_lo(ncell)          !< flow on lo side of cells centered in time
 real(stm_real),intent(in)    :: flow_hi(ncell)          !< flow on hi side of cells centered in time
 real(stm_real),intent(in)    :: conc_lo(ncell,nvar)     !< concentration extrapolated to lo face
 real(stm_real),intent(in)    :: conc_hi(ncell,nvar)     !< concentration extrapolated to hi face
 real(stm_real),intent(in)    :: time                    !< Current time
 real(stm_real),intent(in)    :: dx                      !< Spatial step  
 real(stm_real),intent(in)    :: dt                      !< Time step     
 call advective_flux_lo(flux_lo,    &
                        flux_hi,    &
                        conc_lo,    &
                        conc_hi,    &
                        flow_lo,    &
                        flow_hi,    &
                        ncell,      &
                        nvar,       &
                        time,       &
                        dt,         &
                        dx)
 call advective_flux_hi(flux_lo,    &
                        flux_hi,    &
                        conc_lo,    &
                        conc_hi,    &
                        flow_lo,    &
                        flow_hi,    &
                        ncell,      &
                        nvar,       &
                        time,       &
                        dt,         &
                        dx)                        

return
end subroutine single_channel_boundary_advective_flux


!> Diffusive flux boundary condition for a single channel, delagates to an implementation at each end
subroutine single_channel_boundary_diffusive_flux(diffusive_flux_lo, &
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
    real(stm_real), intent (in)   :: area_lo(ncell)         !< Low side area centered at time
    real(stm_real), intent (in)   :: area_hi(ncell)         !< High side area centered at time
    real(stm_real), intent (in)   ::  time                           !< Time
    real(stm_real), intent (in)   ::  conc(ncell,nvar)               !< Concentration 
    real(stm_real), intent (in)   :: disp_coef_lo(ncell)       !< Low side constituent dispersion coef.
    real(stm_real), intent (in)   :: disp_coef_hi(ncell)       !< High side constituent dispersion coef.
    real(stm_real), intent (in)   :: dt
    real(stm_real), intent (in)   :: dx
    call diffusion_flux_lo(diffusive_flux_lo, &
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

    call diffusion_flux_hi(diffusive_flux_lo, &
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

return
end subroutine


!> Matrix boundary condition, delagates to an implementation at each end
subroutine single_channel_boundary_diffusive_matrix(center_diag ,           &
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
 real(stm_real), intent (in)  :: disp_coef_lo (ncell)                   !< Low side constituent dispersion coef. at new time
 real(stm_real), intent (in)  :: disp_coef_hi (ncell)                   !< High side constituent dispersion coef. at new time
 real(stm_real), intent (in)  :: time                                        !< Current time
 real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
 real(stm_real), intent (in)  :: dx                                          !< Spatial step  
 real(stm_real), intent (in)  :: dt                                          !< Time step     

 call diffusion_matrix_lo(center_diag ,           &
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
                          
 call diffusion_matrix_hi(center_diag ,           &
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
return
end subroutine 






!> Example advective flux that imposes dirichlet boundaries on lo side
subroutine dirichlet_advective_flux_lo(flux_lo,    &
                                       flux_hi,    &
                                       conc_lo,    &
                                       conc_hi,    &
                                       flow_lo,    &
                                       flow_hi,    &
                                       ncell,      &
                                       nvar,       &
                                       time,       &
                                       dt,         &
                                       dx)
     
     use stm_precision
     implicit none
      !--- args          
     integer,intent(in)  :: ncell  !< Number of cells
     integer,intent(in)  :: nvar   !< Number of variables
     ! todo: check the intents
     real(stm_real),intent(inout) :: flux_lo(ncell,nvar)     !< flux on lo side of cell at time 
     real(stm_real),intent(inout) :: flux_hi(ncell,nvar)     !< flux on hi side of cell at time
     real(stm_real),intent(in)    :: flow_lo(ncell)          !< flow on lo side of cells centered at time
     real(stm_real),intent(in)    :: flow_hi(ncell)          !< flow on hi side of cells centered at time
     real(stm_real),intent(in)    :: conc_lo(ncell,nvar)     !< concentration extrapolated to lo face
     real(stm_real),intent(in)    :: conc_hi(ncell,nvar)     !< concentration extrapolated to hi face
     real(stm_real), intent (in)  :: time                    !< Current time
     real(stm_real), intent (in)  :: dx                      !< Spatial step  
     real(stm_real), intent (in)  :: dt                      !< Time step     
    !---local
    real(stm_real) :: bc_data(nvar)
    real(stm_real) :: xloc
    real(stm_real) :: origin = zero !todo: HARDWIRE

    xloc = origin
    call advection_data_lo( bc_data,           &
                            xloc,              &
                            conc_lo,           &
                            ncell,             &
                            nvar,              &
                            origin,            &
                            time,              &
                            dx,                &
                            dt)
   if (flow_lo(1) .ge. zero)then
     flux_lo(1,:)=bc_data*flow_lo(1)
   end if
   return
 end subroutine

!> Example advective flux that imposes dirichlet boundaries on hi side
subroutine dirichlet_advective_flux_hi(flux_lo,    &
                                       flux_hi,    &
                                       conc_lo,    &
                                       conc_hi,    &
                                       flow_lo,    &
                                       flow_hi,    &
                                       ncell,      &
                                       nvar,       &
                                       time,       &
                                       dt,         &
                                       dx)
     
    use stm_precision
    implicit none
   !--- args          
    integer,intent(in)  :: ncell  !< Number of cells
    integer,intent(in)  :: nvar   !< Number of variables
    ! todo: check the intents
    real(stm_real),intent(inout) :: flux_lo(ncell,nvar)     !< flux on lo side of cell, time centered
    real(stm_real),intent(inout) :: flux_hi(ncell,nvar)     !< flux on hi side of cell, time centered
    real(stm_real),intent(in)    :: flow_lo(ncell)          !< flow on lo side of cells centered in time
    real(stm_real),intent(in)    :: flow_hi(ncell)          !< flow on hi side of cells centered in time
    real(stm_real),intent(in)    :: conc_lo(ncell,nvar)     !< concentration extrapolated to lo face
    real(stm_real),intent(in)    :: conc_hi(ncell,nvar)     !< concentration extrapolated to hi face
    real(stm_real), intent (in)  :: time                    !< Current time
    real(stm_real), intent (in)  :: dx                      !< Spatial step  
    real(stm_real), intent (in)  :: dt                      !< Time step     
    !---local
    real(stm_real) :: bc_data(nvar)
    real(stm_real) :: xloc
    real(stm_real) :: origin = zero !todo: HARDWIRE

    xloc = origin
    call advection_data_hi( bc_data,           &
                            xloc,              &
                            conc_hi,           &
                            ncell,             &
                            nvar,              &
                            origin,            &
                            time,              &
                            dx,                &
                            dt)
    if (flow_hi(ncell) .le. zero) then
       flux_hi(ncell,:)=bc_data*flow_hi(ncell)
    end if
    return
 end subroutine


!> dirichlet boundary condition that sets only the low side boundary
subroutine dirichlet_diffusive_flux_lo(diffusive_flux_lo, &
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
    real(stm_real), intent (in)   :: disp_coef_lo (ncell)       !< Low side constituent dispersion coef.
    real(stm_real), intent (in)   :: disp_coef_hi (ncell)       !< High side constituent dispersion coef.
    real(stm_real), intent (in)   :: dt
    real(stm_real), intent (in)   :: dx
    !---local
    real(stm_real) :: bc_data(nvar)
    real(stm_real) :: xloc
    real(stm_real) :: origin = zero !todo: HARDWIRE

    xloc = origin
    call diffusion_data_lo( bc_data,           &
                            xloc,              &
                            conc,              &
                            ncell,             &
                            nvar,              &
                            origin,            &
                            time,              &
                            dx,                &
                            dt)
    
   diffusive_flux_lo(1,:)=-two*area_lo(1)*disp_coef_lo(1)*(conc(1,:)-bc_data(:))/dx
   
   
return
end subroutine

!> dirichlet boundary condition that sets only the hi side boundary
subroutine dirichlet_diffusive_flux_hi(diffusive_flux_lo, &
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
    real(stm_real), intent (in)   :: area_lo         (ncell)         !< Low side area
    real(stm_real), intent (in)   :: area_hi         (ncell)         !< High side area
    real(stm_real), intent (in)   :: time                           !< Time
    real(stm_real), intent (in)   :: conc(ncell,nvar)               !< Concentration 
    real(stm_real), intent (in)   :: disp_coef_lo (ncell)       !< Low side constituent dispersion coef.
    real(stm_real), intent (in)   :: disp_coef_hi (ncell)       !< High side constituent dispersion coef.
    real(stm_real), intent (in)   :: dt
    real(stm_real), intent (in)   :: dx
    real(stm_real) :: bc_data(nvar)
    real(stm_real) :: xloc
    real(stm_real) :: origin = zero !todo: HARDWIRE

    xloc = origin + dx*dble(ncell)
    call diffusion_data_hi( bc_data,           &
                            xloc,              &
                            conc,              &
                            ncell,             &
                            nvar,              &
                            origin,            &
                            time,              &
                            dx,                &
                            dt)
    
   diffusive_flux_hi(ncell,:)=-two*area_hi(ncell)*disp_coef_hi(ncell)*(bc_data(:)-conc(ncell,:))/dx


return
end subroutine



!> Matrix boundary condition for dirichlet, only operates on lo end
subroutine dirichlet_diffusive_matrix_lo(center_diag ,           &
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
 real(stm_real), intent (in)  :: disp_coef_lo(ncell)                   !< Low side constituent dispersion coef. at new time
 real(stm_real), intent (in)  :: disp_coef_hi(ncell)                   !< High side constituent dispersion coef. at new time
 real(stm_real), intent (in)  :: time                                        !< Current time
 real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
 real(stm_real), intent (in)  :: dx                                          !< Spatial step  
 real(stm_real), intent (in)  :: dt                                          !< Time step    
  !---local
 real(stm_real) :: dt_by_dxsq
 real(stm_real) :: xloc
 real(stm_real) :: bc_data(nvar)
 real(stm_real) :: origin = zero !todo: HARDWIRE

 xloc = origin
 call diffusion_data_lo( bc_data,           &
                         xloc,              &
                         conc,              &
                         ncell,             &
                         nvar,              &
                         origin,            &
                         time,              &
                         dx,                &
                         dt)
 
 
 dt_by_dxsq = dt/(dx*dx)
 xloc = origin
! todo: one part of center diag is based on old time and other part new time
!       is this really true?
 center_diag(1,:)= center_diag(1,:)+theta_stm*dt_by_dxsq*(area_lo(1)*disp_coef_lo(1))                  
 right_hand_side(1,:) = right_hand_side(1,:)&
             + two*theta_stm*dt_by_dxsq*(area_lo(1)*disp_coef_lo(1))*bc_data



return
end subroutine


!> Matrix boundary condition for dirichlet, only operates on hi end
subroutine dirichlet_diffusive_matrix_hi(center_diag ,           &
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
 real(stm_real), intent (in)  :: disp_coef_lo (ncell)                   !< Low side constituent dispersion coef. at new time
 real(stm_real), intent (in)  :: disp_coef_hi (ncell)                   !< High side constituent dispersion coef. at new time
 real(stm_real), intent (in)  :: time                                        !< Current (new) time
 real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
 real(stm_real), intent (in)  :: dx                                          !< Spatial step  
 real(stm_real), intent (in)  :: dt                                          !< Time step    
  !---local
 real(stm_real) :: dt_by_dxsq
 real(stm_real) :: xloc
 real(stm_real) :: bc_data(nvar)
 real(stm_real) :: origin = zero !todo: HARDWIRE

 dt_by_dxsq = dt/(dx*dx)
 xloc = origin + dx*dble(ncell)
 call diffusion_data_hi( bc_data,           &
                         xloc,              &
                         conc,              &
                         ncell,             &
                         nvar,              &
                         origin,            &
                         time,              &
                         dx,                &
                         dt)
 
 ! todo: one part of center diag is based on old time and other part new time?
 ! todo: is this really true?
 center_diag(ncell,:)= center_diag(ncell,:)&
                    +  theta_stm*dt_by_dxsq*(area_hi(ncell)*disp_coef_hi(ncell))
 right_hand_side(ncell,:) = right_hand_side(ncell,:)&
                    + two*theta_stm*dt_by_dxsq*(area_hi(ncell)*disp_coef_hi(ncell))*bc_data

return
end subroutine


!=============================================================================================



 !> Example advective flux that imposes Neumann boundaries on lo side of channel
 subroutine neumann_advective_flux_lo(flux_lo,    &
                                      flux_hi,    &
                                      conc_lo,    &
                                      conc_hi,    &
                                      flow_lo,    &
                                      flow_hi,    &
                                      ncell,      &
                                      nvar,       &
                                      time,       &
                                      dt,         &
                                      dx)
     
   use stm_precision
   use error_handling
   implicit none
   !--- args          
   integer,intent(in)  :: ncell  !< Number of cells
   integer,intent(in)  :: nvar   !< Number of variables
   ! todo: check the intents
   real(stm_real),intent(inout) :: flux_lo(ncell,nvar)     !< flux on lo side of cell, time centered
   real(stm_real),intent(inout) :: flux_hi(ncell,nvar)     !< flux on hi side of cell, time centered
   real(stm_real),intent(in)    :: flow_lo(ncell)          !< flow on lo side of cells centered in time
   real(stm_real),intent(in)    :: flow_hi(ncell)          !< flow on hi side of cells centered in time
   real(stm_real),intent(in)    :: conc_lo(ncell,nvar)     !< concentration extrapolated to lo face
   real(stm_real),intent(in)    :: conc_hi(ncell,nvar)     !< concentration extrapolated to hi face
   real(stm_real), intent (in)  :: time                    !< Current time
   real(stm_real), intent (in)  :: dx                      !< Spatial step  
   real(stm_real), intent (in)  :: dt                      !< Time step     
   call stm_fatal("Single channel advection neumann boundary not implemented")
   
   return
 end subroutine

 !> Example advective flux that imposes Neumann boundaries on hi side of channel
 subroutine neumann_advective_flux_hi(flux_lo,    &
                                      flux_hi,    &
                                      conc_lo,    &
                                      conc_hi,    &
                                      flow_lo,    &
                                      flow_hi,    &
                                      ncell,      &
                                      nvar,       &
                                      time,       &
                                      dt,         &
                                      dx)
     
   use stm_precision
   use error_handling

   implicit none
   !--- args          
   integer,intent(in)  :: ncell  !< Number of cells
   integer,intent(in)  :: nvar   !< Number of variables
   ! todo: check the intents
   real(stm_real),intent(inout) :: flux_lo(ncell,nvar)     !< flux on lo side of cell, time centered
   real(stm_real),intent(inout) :: flux_hi(ncell,nvar)     !< flux on hi side of cell, time centered
   real(stm_real),intent(in)    :: flow_lo(ncell)          !< flow on lo side of cells centered in time
   real(stm_real),intent(in)    :: flow_hi(ncell)          !< flow on hi side of cells centered in time
   real(stm_real),intent(in)    :: conc_lo(ncell,nvar)     !< concentration extrapolated to lo face
   real(stm_real),intent(in)    :: conc_hi(ncell,nvar)     !< concentration extrapolated to hi face
   real(stm_real), intent (in)  :: time                    !< Current time
   real(stm_real), intent (in)  :: dx                      !< Spatial step  
   real(stm_real), intent (in)  :: dt                      !< Time step     

   call stm_fatal("Single channel advection neumann boundary not implemented")
      
   return
 end subroutine


!> neumann boundary condition that sets only the lo side boundary
subroutine neumann_diffusive_flux_lo(diffusive_flux_lo, &
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
    real(stm_real), intent (in)   :: disp_coef_lo (ncell)       !< Low side constituent dispersion coef.
    real(stm_real), intent (in)   :: disp_coef_hi (ncell)       !< High side constituent dispersion coef.
    real(stm_real), intent (in)   :: dt
    real(stm_real), intent (in)   :: dx
    !---local
    real(stm_real) :: dt_by_dxsq
    real(stm_real) :: xloc
    real(stm_real) :: bc_data(nvar)
    real(stm_real) :: origin = zero !todo: HARDWIRE

    xloc = origin
    call diffusion_data_lo( bc_data,           &
                         xloc,              &
                         conc,              &
                         ncell,             &
                         nvar,              &
                         origin,            &
                         time,              &
                         dx,                &
                         dt)    
    
        
    diffusive_flux_lo(1,:) = -area_lo(1)*disp_coef_lo(1)*bc_data
    
return
end subroutine

!> neumann boundary condition that sets only the hi side boundary
subroutine neumann_diffusive_flux_hi(diffusive_flux_lo, &
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
    real(stm_real), intent (in)   :: disp_coef_lo (ncell)       !< Low side constituent dispersion coef.
    real(stm_real), intent (in)   :: disp_coef_hi (ncell)       !< High side constituent dispersion coef.
    real(stm_real), intent (in)   :: dt
    real(stm_real), intent (in)   :: dx

  !---local
 real(stm_real) :: dt_by_dxsq
 real(stm_real) :: xloc
 real(stm_real) :: bc_data(nvar)
 real(stm_real) :: origin = zero !todo: HARDWIRE

 dt_by_dxsq = dt/(dx*dx)
 xloc = origin + dx*dble(ncell)
 call diffusion_data_hi( bc_data,           &
                         xloc,              &
                         conc,              &
                         ncell,             &
                         nvar,              &
                         origin,            &
                         time,              &
                         dx,                &
                         dt)
    
  diffusive_flux_hi(ncell,:) = -area_hi(ncell)*disp_coef_hi(ncell)*bc_data


return
end subroutine


!> Matrix boundary condition for neumann, only operates on lo end
subroutine neumann_diffusive_matrix_lo(center_diag ,           &
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
 real(stm_real), intent (in)  :: disp_coef_lo (ncell)                   !< Low side constituent dispersion coef. at new time
 real(stm_real), intent (in)  :: disp_coef_hi (ncell)                   !< High side constituent dispersion coef. at new time
 real(stm_real), intent (in)  :: time                                        !< Current time
 real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
 real(stm_real), intent (in)  :: dx                                          !< Spatial step  
 real(stm_real), intent (in)  :: dt                                          !< Time step    

  !---local
 real(stm_real) :: dt_by_dxsq
 real(stm_real) :: xloc
 real(stm_real) :: bc_data(nvar)
 real(stm_real) :: origin = zero !todo: HARDWIRE
 real(stm_real) :: flux(nvar)
 dt_by_dxsq = dt/(dx*dx) 
 xloc = origin
 call diffusion_data_lo( bc_data,           &
                         xloc,              &
                         conc,              &
                         ncell,             &
                         nvar,              &
                         origin,            &
                         time,              &
                         dx,                &
                         dt)
 
 
     ! todo: there may be an issue of area vs area previous    
    
    flux(:) = -area_lo(1)*disp_coef_lo(1)*bc_data
    
    center_diag(1,:)= area(1)+ theta_stm*dt_by_dxsq* area_hi(1)*disp_coef_hi(1)  
    right_hand_side(1,:) = right_hand_side(1,:)+ theta_stm*(dt/dx)*flux(:)
     

return
end subroutine


!> Matrix boundary condition for neumann, only operates on hi end
subroutine neumann_diffusive_matrix_hi(center_diag ,           &
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
 real(stm_real), intent (in)  :: disp_coef_lo (ncell)                   !< Low side constituent dispersion coef. at new time
 real(stm_real), intent (in)  :: disp_coef_hi (ncell)                   !< High side constituent dispersion coef. at new time
 real(stm_real), intent (in)  :: time                                        !< Current time
 real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
 real(stm_real), intent (in)  :: dx                                          !< Spatial step  
 real(stm_real), intent (in)  :: dt                                          !< Time step    

  !---local
 real(stm_real) :: dt_by_dxsq
 real(stm_real) :: xloc
 real(stm_real) :: bc_data(nvar)
 real(stm_real) :: origin = zero !todo: HARDWIRE
 real(stm_real) :: flux(nvar)
 dt_by_dxsq = dt/(dx*dx) 
 xloc = origin + dx*dble(ncell)
 call diffusion_data_hi( bc_data,           &
                         xloc,              &
                         conc,              &
                         ncell,             &
                         nvar,              &
                         origin,            &
                         time,              &
                         dx,                &
                         dt)
 
 
     ! todo: there may be an issue of area vs area previous    
    flux(:) = -area_hi(ncell)*disp_coef_hi(ncell)*bc_data
    
    
    center_diag(ncell,:)= area(ncell)+ theta_stm*dt_by_dxsq* area_lo(ncell)*disp_coef_lo(1)
    right_hand_side(ncell,:)= right_hand_side(ncell,:) - theta_stm*(dt/dx)*flux(:)

return
end subroutine

end module