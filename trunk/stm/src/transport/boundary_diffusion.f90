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
                                             dt)
        
         
        
         use stm_precision
         implicit none
         !--- args
         integer, intent(in)  :: ncell                                    !< Number of cells
         integer, intent(in)  :: nvar                                     !< Number of variables
         real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)   !< Face flux, lo side
         real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)   !< Face flux, hi side
         real(stm_real), intent (in)   :: area_lo         (ncell)         !< Low side area centered at old time
         real(stm_real), intent (in)   :: area_hi         (ncell)         !< High side area centered at old time
         real(stm_real), intent (in)   ::  time                           !< Time
         real(stm_real), intent (in)   ::  conc(ncell,nvar)               !< Concentration 
         real(stm_real), intent (in)   :: disp_coef_lo (ncell,nvar)       !< Low side constituent dispersion coef.
         real(stm_real), intent (in)   :: disp_coef_hi (ncell,nvar)       !< High side constituent dispersion coef.
         real(stm_real), intent (in)   :: dt
         
       end subroutine boundary_diffusive_flux_if
 end interface

 !> This pointer should be set by the driver or client code to specify the 
 !> treatment at the boundaries
 procedure(boundary_diffusive_flux_if),pointer :: boundary_diffusion_flux  => null()

!=========================================================

 abstract interface
       !> Generic interface for calculating BC of matrix that should be fulfilled by
       !> the driver or the client programs
       subroutine boundary_diffusive_matrix_if( center_diag ,           &
                                                      up_diag,          &     
                                                      down_diag,        &
                                                      right_hand_side,  &
                                                      explicit_diffuse_op, &
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
                                       
        integer, intent (in) :: ncell                                               !< Number of cells
        integer, intent (in) :: nvar                                                !< Number of variables

        real(stm_real),intent (inout):: down_diag(ncell,nvar)                       !< Values of the coefficients below diagonal in matrix
        real(stm_real),intent (inout):: center_diag(ncell,nvar)                     !< Values of the coefficients at the diagonal in matrix
        real(stm_real),intent (inout):: up_diag(ncell,nvar)                         !< Values of the coefficients above the diagonal in matrix
        real(stm_real),intent (inout):: right_hand_side(ncell,nvar)                 !< Values of the coefficients of right  hand side vector
        real(stm_real), intent (in)  :: explicit_diffuse_op(ncell,nvar)              
        real(stm_real), intent (in)  :: area (ncell)                                !< Cell centered area at new time 
        real(stm_real), intent (in)  :: area_lo(ncell)                              !< Low side area at new time
        real(stm_real), intent (in)  :: area_hi(ncell)                              !< High side area at new time 
        real(stm_real), intent (in)  :: disp_coef_lo (ncell,nvar)                   !< Low side constituent dispersion coef. at new time
        real(stm_real), intent (in)  :: disp_coef_hi (ncell,nvar)                   !< High side constituent dispersion coef. at new time
        real(stm_real), intent (in)  :: time                                        !< Current time
        real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
        real(stm_real), intent (in)  :: dx                                          !< Spatial step  
        real(stm_real), intent (in)  :: dt                                          !< Time step     

       end subroutine boundary_diffusive_matrix_if
 end interface

 !> This pointer should be set by the driver or client code to specify the 
 !> treatment at the first and last row of coefficient matrix
 ! todo: check the "boundary_diffusion_impose"
 procedure(boundary_diffusive_matrix_if),pointer :: boundary_diffusion_impose  => null()

 contains
 
 !> Example diffusive flux that prints an error and bails
 subroutine uninitialized_diffusive_flux(diffusive_flux_lo, &
                                         diffusive_flux_hi, &
                                         conc,              &
                                         area_lo,           &
                                         area_hi,           &
                                         disp_coef_lo,      &  
                                         disp_coef_hi,      &
                                         ncell,             &
                                         nvar,              &
                                         time)
                                         
     use stm_precision 
     use error_handling
 
         implicit none
         !--- args
         integer, intent(in)  :: ncell                                   !< number of cells
         integer, intent(in)  :: nvar                                    !< number of variables
         real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)  !< face flux, lo side
         real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)  !< face flux, hi side
         real(stm_real), intent (in)   :: area_lo         (ncell)        !< Low side area centered at old time
         real(stm_real), intent (in)   :: area_hi         (ncell)        !< High side area centered at old time
         real(stm_real), intent (in)   ::  time                          !< time
         real(stm_real), intent (in)   ::  conc(ncell,nvar)              !< concentration 
         real(stm_real), intent (in)   :: disp_coef_lo (ncell,nvar)      !< Low side constituent dispersion coef.
         real(stm_real), intent (in)   :: disp_coef_hi (ncell,nvar)      !< High side constituent dispersion coef.
    
     call stm_fatal("boundary not implemented")
     
     return
 end subroutine 
 
 !> Example diffusive flux that imposes Neumann boundaries with zero flux at
 !> both ends of the channel.
 subroutine neumann_no_flow_diffusive_flux(diffusive_flux_lo,  &
                                             diffusive_flux_hi, &
                                             conc,              &
                                             area_lo,           &
                                             area_hi,           &
                                             disp_coef_lo,      &  
                                             disp_coef_hi,      &
                                             ncell,             &
                                             nvar,              &
                                             time,              &
                                             dt)
    use stm_precision
         implicit none
         !--- args
         integer, intent(in)  :: ncell                                   !< number of cells
         integer, intent(in)  :: nvar                                    !< number of variables
         real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)  !< face flux, lo side
         real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)  !< face flux, hi side
         real(stm_real), intent (in)   :: area_lo         (ncell)        !< Low side area centered at old time
         real(stm_real), intent (in)   :: area_hi         (ncell)        !< High side area centered at old time
         real(stm_real), intent (in)   ::  time                          !< time
         real(stm_real), intent (in)   ::  conc(ncell,nvar)              !< concentration 
         real(stm_real), intent (in)   :: disp_coef_lo (ncell,nvar)      !< Low side constituent dispersion coef.
         real(stm_real), intent (in)   :: disp_coef_hi (ncell,nvar)      !< High side constituent dispersion coef.
         real(stm_real), intent (in)   :: dt
     
     diffusive_flux_lo(1,:) = zero
     diffusive_flux_hi(ncell,:) = zero
     
    return
 end subroutine
 
!> Example diffusive flux that imposes sinusoidal time dependent Neumann boundary flux at
!> both ends of the channel.
 subroutine neumann_sin_diffusive_flux(diffusive_flux_lo,   &
                                             diffusive_flux_hi, &
                                             conc,              &
                                             area_lo,           &
                                             area_hi,           &
                                             disp_coef_lo,      &  
                                             disp_coef_hi,      &
                                             ncell,             &
                                             nvar,              &
                                             time,              &
                                             dt)
    use stm_precision
         implicit none
         !--- args
         integer, intent(in)  :: ncell                                   !< number of cells
         integer, intent(in)  :: nvar                                    !< number of variables
         real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)  !< face flux, lo side
         real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)  !< face flux, hi side
         real(stm_real), intent (in)   :: area_lo         (ncell)        !< Low side area centered at old time
         real(stm_real), intent (in)   :: area_hi         (ncell)        !< High side area centered at old time
         real(stm_real), intent (in)   ::  time                          !< time
         real(stm_real), intent (in)   ::  conc(ncell,nvar)              !< concentration 
         real(stm_real), intent (in)   :: disp_coef_lo (ncell,nvar)      !< Low side constituent dispersion coef.
         real(stm_real), intent (in)   :: disp_coef_hi (ncell,nvar)      !< High side constituent dispersion coef.
         real(stm_real), intent (in)   :: dt
          
     diffusive_flux_lo(1,:) = two * cos( pi* time / three)               !Just for test 
     diffusive_flux_hi(ncell,:) = five * sin (pi*time / seven)
        
     return
 end subroutine
 
subroutine n_d_test_diffusive_flux(diffusive_flux_lo,   &
                                     diffusive_flux_hi, &
                                     conc,              &
                                     area_lo,           &
                                     area_hi,           &
                                     disp_coef_lo,      &  
                                     disp_coef_hi,      &
                                     ncell,             &
                                     nvar,              &
                                     time,              &
                                     dt)
    use stm_precision
         implicit none
         !--- args
         integer, intent(in)  :: ncell                                   !< number of cells
         integer, intent(in)  :: nvar                                    !< number of variables
         real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)  !< face flux, lo side
         real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)  !< face flux, hi side
         real(stm_real), intent (in)   :: area_lo         (ncell)        !< Low side area centered at old time
         real(stm_real), intent (in)   :: area_hi         (ncell)        !< High side area centered at old time
         real(stm_real), intent (in)   :: time                           !< time
         real(stm_real), intent (in)   :: conc(ncell,nvar)               !< concentration 
         real(stm_real), intent (in)   :: disp_coef_lo (ncell,nvar)      !< Low side constituent dispersion coef.
         real(stm_real), intent (in)   :: disp_coef_hi (ncell,nvar)      !< High side constituent dispersion coef.
         real(stm_real), intent (in)   :: dt
    !--local
    real(stm_real)::xstart = 0.1d0
    real(stm_real)::xend = one
    real(stm_real) :: dx
    dx = (xend-xstart)/ncell
    
    ! todo : time /time -dt
    ! dt is hard wired 
    !two do what about negative???
    diffusive_flux_lo(1,:) = -100.0d0*half *(two-two*pi *sin(pi*xstart/two)*exp(-half*pi*pi*(time-dt)/four))
    diffusive_flux_hi(ncell,:) = -two*100.0d0*half*(two-conc(ncell,:))/dx
        ! todo : A and Ks are hard wired here (A=100.0 and Ks=three)
     return
 end subroutine
 
 
 !===========================================
  !> Example matrix that prints an error and bails
 subroutine uninitialized_diffusive_bc_matrix(center_diag,      &
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
     use error_handling
     
     implicit none
         !--- args
                                       
        integer, intent (in) :: ncell                                               !< Number of cells
        integer, intent (in) :: nvar                                                !< Number of variables

        real(stm_real),intent (inout):: down_diag(ncell,nvar)                       !< Values of the coefficients below diagonal in matrix
        real(stm_real),intent (inout):: center_diag(ncell,nvar)                     !< Values of the coefficients at the diagonal in matrix
        real(stm_real),intent (inout):: up_diag(ncell,nvar)                         !< Values of the coefficients above the diagonal in matrix
        real(stm_real), intent (in)  :: area (ncell)                                !< Cell centered area at new time 
        real(stm_real), intent (in)  :: area_lo(ncell)                              !< Low side area at new time
        real(stm_real), intent (in)  :: area_hi(ncell)                              !< High side area at new time 
        real(stm_real), intent (in)  :: disp_coef_lo (ncell,nvar)                   !< Low side constituent dispersion coef. at new time
        real(stm_real), intent (in)  :: disp_coef_hi (ncell,nvar)                   !< High side constituent dispersion coef. at new time
        real(stm_real), intent (in)  :: time                                        !< Current time
        real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
        real(stm_real), intent (in)  :: dx                                          !< Spatial step  
        real(stm_real), intent (in)  :: dt                                          !< Time step     
                                                 
      
     call stm_fatal("boundary not implemented!")
     
     return
 end subroutine 
 
 !> Example diffusive flux that imposes Neumann boundaries with zero flux at
 !> both ends of the channel.
 !todo: make sure this is generic for all neumann bc
 subroutine neumann_diffusion_matrix( center_diag ,       &
                                      up_diag,            &     
                                      down_diag,          &
                                      right_hand_side,    & 
                                      explicit_diffuse_op,&
                                      area,               &
                                      area_lo,            &
                                      area_hi,            &          
                                      disp_coef_lo,       &
                                      disp_coef_hi,       &
                                      theta_stm,          &
                                      ncell,              &
                                      time,               & 
                                      nvar,               & 
                                      dx,                 &
                                      dt)
     use stm_precision
     implicit none
         !--- args
                                       
        integer, intent (in) :: ncell                                               !< Number of cells
        integer, intent (in) :: nvar                                                !< Number of variables

        real(stm_real),intent (inout):: down_diag(ncell,nvar)                       !< Values of the coefficients below diagonal in matrix
        real(stm_real),intent (inout):: center_diag(ncell,nvar)                     !< Values of the coefficients at the diagonal in matrix
        real(stm_real),intent (inout):: up_diag(ncell,nvar)                         !< Values of the coefficients above the diagonal in matrix
        real(stm_real),intent (inout):: right_hand_side(ncell,nvar)                 !< Values of the coefficients of the right hand side
        real(stm_real), intent (in)  :: explicit_diffuse_op(ncell,nvar)
        real(stm_real), intent (in)  :: area (ncell)                                !< Cell centered area at new time 
        real(stm_real), intent (in)  :: area_lo(ncell)                              !< Low side area at new time
        real(stm_real), intent (in)  :: area_hi(ncell)                              !< High side area at new time 
        real(stm_real), intent (in)  :: disp_coef_lo (ncell,nvar)                   !< Low side constituent dispersion coef. at new time
        real(stm_real), intent (in)  :: disp_coef_hi (ncell,nvar)                   !< High side constituent dispersion coef. at new time
        real(stm_real), intent (in)  :: time                                        !< Current time
        real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
        real(stm_real), intent (in)  :: dx                                          !< Spatial step  
        real(stm_real), intent (in)  :: dt                                          !< Time step     
      
        !---local
        real(stm_real) :: d_star 
        d_star = dt/(dx*dx)  
      
     ! todo: add types of other BC 
          
     center_diag(1,nvar)= area(1) + theta_stm*d_star*(area_hi(1)*disp_coef_hi(1,nvar) )
     center_diag(ncell,nvar)= area(ncell) + theta_stm*d_star*(area_lo(ncell)*disp_coef_lo(ncell,nvar))
     
     ! todo: implement and test
     return
 end subroutine
 
  subroutine n_d_test_diffusion_matrix(center_diag ,       &
                                       up_diag,            &     
                                       down_diag,          &
                                       right_hand_side,    & 
                                       explicit_diffuse_op,&
                                       area,               &
                                       area_lo,            &
                                       area_hi,            &          
                                       disp_coef_lo,       &
                                       disp_coef_hi,       &
                                       theta_stm,          &
                                       ncell,              &
                                       time,               & 
                                       nvar,               & 
                                       dx,                 &
                                       dt)
     use stm_precision
     implicit none
         !--- args
                                       
        integer, intent (in) :: ncell                                               !< Number of cells
        integer, intent (in) :: nvar                                                !< Number of variables

        real(stm_real),intent (inout):: down_diag(ncell,nvar)                       !< Values of the coefficients below diagonal in matrix
        real(stm_real),intent (inout):: center_diag(ncell,nvar)                     !< Values of the coefficients at the diagonal in matrix
        real(stm_real),intent (inout):: up_diag(ncell,nvar)                         !< Values of the coefficients above the diagonal in matrix
        real(stm_real),intent (inout):: right_hand_side(ncell,nvar)                 !< Values of the coefficients of right hand side vector
        real(stm_real), intent (in)  :: explicit_diffuse_op(ncell,nvar) 
        real(stm_real), intent (in)  :: area (ncell)                                !< Cell centered area at new time 
        real(stm_real), intent (in)  :: area_lo(ncell)                              !< Low side area at new time
        real(stm_real), intent (in)  :: area_hi(ncell)                              !< High side area at new time 
        real(stm_real), intent (in)  :: disp_coef_lo (ncell,nvar)                   !< Low side constituent dispersion coef. at new time
        real(stm_real), intent (in)  :: disp_coef_hi (ncell,nvar)                   !< High side constituent dispersion coef. at new time
        real(stm_real), intent (in)  :: time                                        !< Current time
        real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
        real(stm_real), intent (in)  :: dx                                          !< Spatial step  
        real(stm_real), intent (in)  :: dt                                          !< Time step     
      
        !---local
        real(stm_real) :: d_star 
        d_star = dt/(dx*dx)  
      
     ! todo: error is here
          ! todo : KS=3 hard wired
     center_diag(1,:)= area(1) + theta_stm*d_star*(area_hi(1)*disp_coef_hi(1,:)) 
     up_diag(1,:) = - theta_stm*d_star*(area_hi(1)*disp_coef_hi(1,:))
     right_hand_side(1,:) = right_hand_side(1,:)&
                                - (- theta_stm*(dt/dx)*area_hi(1)*disp_coef_hi(1,:)*(two-two*pi*sin(0.05d0*pi)*exp(-half*pi*pi*time/four) ))
     
     
     
     center_diag(ncell,:)=  area(ncell) + theta_stm*d_star*(area_hi(ncell)*disp_coef_hi(ncell,:)+two*area_lo(ncell)*disp_coef_lo(ncell,:))
     down_diag(ncell,:) = - theta_stm*d_star*(area_lo(ncell)*disp_coef_lo(ncell,:))
     right_hand_side(ncell,:) = right_hand_side(ncell,:)&
                                 - two*theta_stm*d_star*area_hi(:)*disp_coef_hi(ncell,:)*two ! two is c on the boundary
     
     return
 end subroutine
 
end module
