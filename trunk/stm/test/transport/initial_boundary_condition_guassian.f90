 
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
    subroutine fill_gaussian(vals,nloc,origin,dx,mean,sd,scale)
    ! fill_guassian(OUTPUT,num_cell,Left side of domain,dx,Center,sigma,a)
    ! f(x) = a*exp(-(x-b)^2/(2c^2)) [c is sigma]   
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
    real(stm_real) :: xlo
    real(stm_real) :: xhi
    integer        :: iloc
    real(stm_real) :: actual_scale
   
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
!> Compute the slope of gaussian shape, base on sigma, center 
!> of shape, scale factor and distance from the center
subroutine derivative_gaussian(val,xposition,center,sd,scale)
    ! derivative_gaussian(OUTPUT or df/dx,x,center or miu,sigma,scale or a)  
    ! f(x) = a*exp(-(x-miu)^2/(2c^2)) , [c is sigma]   
    ! df(x)/dx = -a*2*(x-miu)/(2c^2)*exp(-(x-miu)^2/(2c^2))
    use stm_precision
    implicit none
    real(stm_real), intent(out) :: val            !< value to be produced
    real(stm_real), intent(in)  :: xposition      !< X
    real(stm_real), intent(in)  :: center         !< center of gaussian shape
    real(stm_real), intent(in)  :: sd             !< Standard deviation (Sigma)
    real(stm_real), intent(in)  :: scale !< scale
    !---locals
   
    val = -(scale*(xposition - center)/(sd*sd))*exp(-(xposition-center)**2/(two*sd*sd)) 
   
    return
end subroutine


end module
