 
 !<!license>
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

!> Initial conditions for examples
!> It can fill an array with 1D Gaussin function, rectangular, and triangular,
!> and also it can fill discontinuities. In case of multi-constituents tehy have to initialize separately     
!>@ingroup example
module gaussian_init_boundary_condition
use stm_precision   
private gaussian_cdf

contains

!> Gaussian cdf (integrated Guassian pdf) from -inf to x0
real(stm_real) function gaussian_cdf(x0,mean,sd)        
    implicit none   
    real(stm_real), intent(in) :: x0   !< End of integration
    real(stm_real), intent(in) :: mean !< Mean
    real(stm_real), intent(in) :: sd   !< Sigma (standard deviation)
    gaussian_cdf = half + half*erf((x0-mean)/(sqrt(two)*sd))
    return
end function


!> Fill array with 1D gaussian shape
!> This routine expects a 1D array, so multi-constituents
!> have to be initialized separately
!> fill_guassian(OUTPUT,num_cell,Left side of domain,dx,Center,sigma,a)
!> f(x) = a*exp(-(x-b)^2/(2c^2)) [c is sigma]
subroutine fill_gaussian(vals,nloc,origin,dx,mean,sd,scale) 
    use stm_precision
    implicit none
    integer, intent(in) :: nloc                   !< Number of cells (size of array) 
    real(stm_real), intent(out) :: vals(nloc)     !< Values to be filled
    real(stm_real), intent(in)  :: origin         !< Origin (lo side of channel)
    real(stm_real), intent(in)  :: dx             !< dx
    real(stm_real), intent(in)  :: mean           !< Center of the gaussian shape
    real(stm_real), intent(in)  :: sd             !< Standard deviation (Sigma)
    real(stm_real), intent(in), optional :: scale !< scale
    !-----locals
    real(stm_real) :: xlo                         !< Low side position
    real(stm_real) :: xhi                         !< High side position
    integer        :: iloc                        !< Cell counter
    real(stm_real) :: actual_scale                !< Scale 
   
    if (present(scale))then
        actual_scale = scale*sqrt(two*pi*sd*sd)
    else
        actual_scale = one*sqrt(two*pi*sd*sd)
    end if
    
    do iloc = 1,nloc
       xlo = origin + dble(iloc - 1)*dx
       xhi = origin + dble(iloc)*dx
      ! todo: need to populate using cell averages
       vals(iloc) =  (gaussian_cdf(xhi,mean,sd) & 
                    - gaussian_cdf(xlo,mean,sd))
    end do
    vals = vals*(actual_scale/dx)
    return
end subroutine

!> Compute a point value of a gaussian function 
!> f(x) = a*exp(-(x-b)^2/(2c^2)) , [c is sigma]   
!> where a is a scale factor, b is the center/mean and 
!> is a distance/standard deviation
subroutine gaussian(val,xposition,center,sd,scale)
    ! df(x)/dx = -a*2*(x-b)/(2c^2)*exp(-(x-b)^2/(2c^2))
    use stm_precision
    implicit none
    real(stm_real), intent(out) :: val            !< value to be produced
    real(stm_real), intent(in)  :: xposition      !< X
    real(stm_real), intent(in)  :: center         !< Center of gaussian shape
    real(stm_real), intent(in)  :: sd             !< Standard deviation (Sigma)
    real(stm_real), intent(in)  :: scale          !< scale
    !---locals
   
    val = scale*exp(-(xposition-center)**2/(two*sd*sd)) 
   
    return
end subroutine


!> Compute the slope of gaussian shape, base on sigma, center 
!> of shape, scale factor and distance from the center
!> derivative_gaussian(OUTPUT or df/dx,x,center or miu,sigma,scale or a)  
!> f(x) = a*exp(-(x-b)^2/(2c^2)) , [c is sigma]   
!> df(x)/dx = -a*2*(x-b)/(2c^2)*exp(-(x-b)^2/(2c^2))
subroutine derivative_gaussian(val,xposition,center,sd,scale)
    
    use stm_precision
    implicit none
    real(stm_real), intent(out) :: val            !< value to be produced
    real(stm_real), intent(in)  :: xposition      !< X
    real(stm_real), intent(in)  :: center         !< center of gaussian shape
    real(stm_real), intent(in)  :: sd             !< Standard deviation (Sigma)
    real(stm_real), intent(in)  :: scale          !< Scale
    !---locals
   
    val = -(scale*(xposition - center)/(sd*sd))*exp(-(xposition-center)**2/(two*sd*sd)) 
   
    return
end subroutine

!> Example diffusive flux that imposes Neumann boundaries with zero flux at
!> both ends of the channel.
!todo: make sure this is generic for all Neumann bc
subroutine neumann_diffusion_matrix(center_diag ,       &
                                     up_diag,            &     
                                     down_diag,          &
                                     right_hand_side,    & 
                                     conc,               &
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
    real(stm_real), intent (in)  :: conc(ncell,nvar)                            !< Concentration
    real(stm_real), intent (in)  :: explicit_diffuse_op(ncell,nvar)             !< Explicit diffusive operator
    real(stm_real), intent (in)  :: area (ncell)                                !< Cell centered area at new time 
    real(stm_real), intent (in)  :: area_lo(ncell)                              !< Low side area at new time
    real(stm_real), intent (in)  :: area_hi(ncell)                              !< High side area at new time 
    real(stm_real), intent (in)  :: disp_coef_lo (ncell)                        !< Low side constituent dispersion coef. at new time
    real(stm_real), intent (in)  :: disp_coef_hi (ncell)                        !< High side constituent dispersion coef. at new time
    real(stm_real), intent (in)  :: time                                        !< Current time
    real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
    real(stm_real), intent (in)  :: dx                                          !< Spatial step  
    real(stm_real), intent (in)  :: dt                                          !< Time step     
    !---local
    real(stm_real) :: dt_by_dxsq 
    real(stm_real) :: flux_start(nvar)
    real(stm_real) :: flux_end(nvar)
         
     ! todo: these must control for area or area_prev 
         
    real(stm_real) :: xend =51200d0/two
    dt_by_dxsq = dt/(dx*dx) 
    
    flux_start(:) = -area_lo(1)*disp_coef_lo(1)*(1/sqrt(four*pi*disp_coef_lo(1)))* &
                             (xend/two/disp_coef_lo(1)/time)*exp(-(xend**2)/four/disp_coef_lo(1)/time)
    flux_end(:) = -area_hi(ncell)*disp_coef_hi(ncell)*(1/sqrt(four*pi*disp_coef_hi(ncell)))* &
                                   (-xend/two/disp_coef_hi(1)/time)*exp(-(xend**2)/four/disp_coef_lo(1)/time)
    
        
    center_diag(1,:)= area(1)+ theta_stm*dt_by_dxsq* area_hi(1)*disp_coef_hi(1)  
    right_hand_side(1,:) = right_hand_side(1,:) &
                               + theta_stm*(dt/dx)*flux_start(:)
    
    center_diag(ncell,:)= area(ncell)+ theta_stm*dt_by_dxsq* area_lo(ncell)*disp_coef_lo(1)
    right_hand_side(ncell,:)= right_hand_side(ncell,:) &
                                  - theta_stm*(dt/dx)*flux_end(:)
    return
end subroutine
 
!> Example diffusion matrix values imposes Dirichlet and Nuemann boundaries at
!> the ends of the channel
subroutine n_d_test_diffusion_matrix(center_diag ,       &
                                       up_diag,            &     
                                       down_diag,          &
                                       right_hand_side,    & 
                                       conc,               &
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
     real(stm_real), intent (in)  :: conc(ncell,nvar)                            !< Concentration 
     real(stm_real), intent (in)  :: explicit_diffuse_op(ncell,nvar)             !< Explicit diffusive operator
     real(stm_real), intent (in)  :: area (ncell)                                !< Cell centered area at new time 
     real(stm_real), intent (in)  :: area_lo(ncell)                              !< Low side area at new time
     real(stm_real), intent (in)  :: area_hi(ncell)                              !< High side area at new time 
     real(stm_real), intent (in)  :: disp_coef_lo(ncell)                         !< Low side constituent dispersion coef. at new time
     real(stm_real), intent (in)  :: disp_coef_hi(ncell)                         !< High side constituent dispersion coef. at new time
     real(stm_real), intent (in)  :: time                                        !< Current time
     real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
     real(stm_real), intent (in)  :: dx                                          !< Spatial step  
     real(stm_real), intent (in)  :: dt                                          !< Time step     
      
        !---local
 
     real(stm_real) :: dt_by_dxsq
     real(stm_real) :: xstart
     real(stm_real) :: xend  
     real(stm_real) :: flux_start(nvar)
     real(stm_real) :: flux_end (nvar)
    
    dt_by_dxsq = dt/(dx*dx)
    xstart = 0.1d0
    xend = one 
    ! area is new?!? 
    !todo: call flux
    !todo: which is dirchlet?
    flux_start(:) = - area_lo(1)*disp_coef_lo(1)*(two-two*pi*sin(pi*xstart/two)*exp(-disp_coef_lo(1)*pi*pi*time/four)) 
    flux_end(:) = - area_hi(ncell)*disp_coef_hi(ncell)*(two-two*pi*sin(pi*xend/two)*exp(-disp_coef_hi(ncell)*pi*pi*time/four))
         
    center_diag(1,:)= area(1)+ theta_stm*dt_by_dxsq* area_hi(1)*disp_coef_hi(1)  
    right_hand_side(1,:) = right_hand_side(1,:) &
                                + theta_stm*(dt/dx)*flux_start(:)
    center_diag(ncell,:)= area(ncell)+ theta_stm*dt_by_dxsq* area_lo(ncell)*disp_coef_lo(1)
    right_hand_side(ncell,:)= right_hand_side(ncell,:) &
                                   - theta_stm*(dt/dx)*flux_end(:)
      
    return
end subroutine
 
!> Example diffusion matrix values imposes Dirichlet boundaries at
!> both ends of the channel. 
subroutine dirichlet_test_diffusion_matrix(center_diag ,       &
                                           up_diag,            &     
                                           down_diag,          &
                                           right_hand_side,    &
                                           conc,               & 
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
   real(stm_real), intent (in)  :: conc(ncell,nvar)                            !< Concentration 
   real(stm_real), intent (in)  :: explicit_diffuse_op(ncell,nvar)             !< Explicit diffusive operator
   real(stm_real), intent (in)  :: area (ncell)                                !< Cell centered area at new time 
   real(stm_real), intent (in)  :: area_lo(ncell)                              !< Low side area at new time
   real(stm_real), intent (in)  :: area_hi(ncell)                              !< High side area at new time 
   real(stm_real), intent (in)  :: disp_coef_lo(ncell)                         !< Low side constituent dispersion coef. at new time
   real(stm_real), intent (in)  :: disp_coef_hi(ncell)                         !< High side constituent dispersion coef. at new time
   real(stm_real), intent (in)  :: time                                        !< Current time
   real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
   real(stm_real), intent (in)  :: dx                                          !< Spatial step  
   real(stm_real), intent (in)  :: dt                                          !< Time step     
   
   !---local
   real(stm_real) :: dt_by_dxsq
   real(stm_real) :: xstart
   real(stm_real) :: xend  
   real(stm_real) :: conc_start(nvar)
   real(stm_real) :: conc_end(nvar)
  
   dt_by_dxsq = dt/(dx*dx)
   xstart = 0.1d0
   xend = one 
   ! here time is new time and area and Ks for updating rhs are for time stap n+1
   conc_end = two
   conc_start = two*xstart + four*cos(pi*xstart/two)*exp(-disp_coef_lo(1)*time*pi*pi/four)
   ! todo: one part of center diag is based on old time and other part new time
   center_diag(1,:)=  center_diag(1,:) &
                         + theta_stm*dt_by_dxsq*(area_lo(1)*disp_coef_lo(1))                  
   right_hand_side(1,:) = right_hand_side(1,:)&
               + two * theta_stm*dt_by_dxsq*(area_lo(1)*disp_coef_lo(1))*conc_start
     
   center_diag(ncell,:)= center_diag(ncell,:)&
                          +  theta_stm*dt_by_dxsq*(area_hi(ncell)*disp_coef_hi(ncell))
   right_hand_side(ncell,:) = right_hand_side(ncell,:)&
              + two * theta_stm*dt_by_dxsq*(area_hi(ncell)*disp_coef_hi(ncell))*conc_end

   
   return
 end subroutine
!> Example diffusive flux that imposes Neumann boundary at one end and Dirichlet boundary at
!> the other end of the channel.
subroutine n_d_test_diffusive_flux(diffusive_flux_lo, &
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
    integer, intent(in)  :: ncell                                   !< Number of cells
    integer, intent(in)  :: nvar                                    !< Number of variables
    real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)  !< Face flux, lo side
    real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)  !< Face flux, hi side
    real(stm_real), intent (in)   :: area_lo(ncell)                 !< Low side area centered at time
    real(stm_real), intent (in)   :: area_hi(ncell)                 !< High side area centered at time
    real(stm_real), intent (in)   :: time                           !< Time
    real(stm_real), intent (in)   :: conc(ncell,nvar)               !< Concentration 
    real(stm_real), intent (in)   :: disp_coef_lo(ncell)            !< Low side constituent dispersion coef.
    real(stm_real), intent (in)   :: disp_coef_hi(ncell)            !< High side constituent dispersion coef.
    real(stm_real), intent (in)   :: dt                             !< Spatial step
    real(stm_real), intent (in)   :: dx                             !< Time step   
    !--local
    real(stm_real) :: xstart = 0.1d0
    real(stm_real) :: xend = one
              
    diffusive_flux_lo(1,:) = -area_lo(1)*disp_coef_lo(1)*(two - two*pi*sin(pi*xstart/two)*exp(-disp_coef_lo(1)*pi*pi*time/four))
    diffusive_flux_hi(ncell,:) = -area_hi(ncell)*disp_coef_hi(ncell)*(two - two*pi*sin(pi*xend/two)*exp(-disp_coef_hi(ncell)*pi*pi*time/four))
    
    return
 end subroutine

!> Example diffusive flux that imposes Dirichlet boundaries at
!> both ends of the channel.  
subroutine dirichlet_test_diffusive_flux(diffusive_flux_lo, &
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
    integer, intent(in)  :: ncell                                   !< Number of cells
    integer, intent(in)  :: nvar                                    !< Number of variables
    real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)  !< Face flux, lo side
    real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)  !< Face flux, hi side
    real(stm_real), intent (in)   :: area_lo(ncell)                 !< Low side area centered at time
    real(stm_real), intent (in)   :: area_hi(ncell)                 !< High side area centered at time
    real(stm_real), intent (in)   :: time                           !< Time
    real(stm_real), intent (in)   :: conc(ncell,nvar)               !< Concentration 
    real(stm_real), intent (in)   :: disp_coef_lo (ncell)           !< Low side constituent dispersion coef.
    real(stm_real), intent (in)   :: disp_coef_hi (ncell)           !< High side constituent dispersion coef.
    real(stm_real), intent (in)   :: dt                             !< Time step  
    real(stm_real), intent (in)   :: dx                             !< Spatial step
    !--local
    
    real(stm_real) :: conc_start(nvar)
    real(stm_real) :: conc_end(nvar) 
    real(stm_real) :: xstart = 0.1d0
   
    conc_end(:) = two
    conc_start(:) = two*xstart + four*cos(pi*xstart/two)*exp(-disp_coef_lo(1)*pi*pi*time/four)
   
   ! todo: check convergence for second order boundary fitting 
   ! todo: this area also must be area_prev  
   diffusive_flux_lo(1,:)=-two*area_lo(1)*disp_coef_lo(1)*(conc(1,:)-conc_start(:))/dx
    
   diffusive_flux_hi(ncell,:)=-two*area_hi(ncell)*disp_coef_hi(ncell)*(conc_end(:)-conc(ncell,:))/dx
        
   return
 end subroutine


 !> Example diffusive flux that imposes Neumann boundaries with zero flux at
 !> both ends of the channel.
 subroutine neumann_gaussian_diffusive_flux(diffusive_flux_lo, &
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
    integer, intent(in)  :: ncell                                   !< Number of cells
    integer, intent(in)  :: nvar                                    !< Number of variables
    real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)  !< Face flux, lo side
    real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)  !< Face flux, hi side
    real(stm_real), intent (in)   :: area_lo(ncell)                 !< Low side area centered at time
    real(stm_real), intent (in)   :: area_hi(ncell)                 !< High side area centered at time
    real(stm_real), intent (in)   :: time                           !< Time
    real(stm_real), intent (in)   :: conc(ncell,nvar)               !< Concentration 
    real(stm_real), intent (in)   :: disp_coef_lo(ncell)            !< Low side constituent dispersion coef.
    real(stm_real), intent (in)   :: disp_coef_hi(ncell)            !< High side constituent dispersion coef.
    real(stm_real), intent (in)   :: dt                             !< Time step  
    real(stm_real), intent (in)   :: dx                             !< Spatial step
    real(stm_real) :: xend  
    xend = 51200d0/two
    diffusive_flux_lo(1,:) = -area_lo(1)*disp_coef_lo(1)*(1/sqrt(four*pi*disp_coef_lo(1)))* &
                             (xend/(two*disp_coef_lo(1)*time))*exp(-(xend**2)/four/disp_coef_lo(1)/time)
    diffusive_flux_hi(ncell,:) = -area_hi(ncell)*disp_coef_hi(ncell)*(1/sqrt(four*pi*disp_coef_hi(ncell)))* &
                             (-xend/two/disp_coef_hi(1)/time)*exp(-(xend**2)/four/disp_coef_lo(1)/time)
       
    return
 end subroutine
 
end module
