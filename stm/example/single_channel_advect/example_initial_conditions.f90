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
!>@ingroup example

module example_initial_conditions
    
    use stm_precision
    private gaussian_cdf
    
    contains
    
    !> Gaussian cdf (integrated Guassian pdf) from -inf to x0
    real(stm_real) function gaussian_cdf(x0,mean,sd)
        use stm_precision
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
        !-----------
        if ( present(scale))then
            actual_scale = scale*sqrt(two*pi*sd*sd)
        else
            actual_scale = one*sqrt(two*pi*sd*sd)
        end if
        
        do iloc = 1,nloc
           xlo = origin + dble(iloc - 1)*dx
           xhi = origin + dble(iloc)*dx
          ! need to populate using cell averages
           vals(iloc) =  (gaussian_cdf(xhi,mean,sd) & 
                        - gaussian_cdf(xlo,mean,sd))
        end do
        vals = vals*(actual_scale/dx)
        return
    end subroutine
    
    
    !> Fill array with rectanglar shape
    !> This routine expects a 1D array, so multi-constituents
    !> have to be initialized separately
    subroutine fill_rectangular(array,x,nloc,xlo,xhi,fill,fill_else)

        implicit none
        integer,        intent(in)  :: nloc         !< Size of array
        real(stm_real), intent(out) :: array(nloc)  !< Array to be filled
        real(stm_real), intent(in)  :: x(nloc)      !< Cell-centered x coordinate
        real(stm_real), intent(in)  :: xlo          !< Lo side boundary of fill
        real(stm_real), intent(in)  :: xhi          !< Hi side boundary of fill
        real(stm_real), intent(in)  :: fill         !< Filled value between xlo and xhi
        real(stm_real), intent(in)  :: fill_else    !< Filled value if not between xlo and xhi
                
        array = fill_else
        
        where (x > xlo .and. x< xhi)
            array = fill
        end where
        
        return
    end subroutine
    
        
    !> Fill array with symetric triangular shape (isosceles triangle)
    !> This routine expects a 1D array, so multi-constituents
    !> have to be initialized separately
    subroutine fill_triangular(array,xpos,nloc,xlo,xhi,vertex_hight,fill_else)

        implicit none
        integer,        intent(in)  :: nloc         !< Size of array
        real(stm_real), intent(out) :: array(nloc)  !< Array to be filled
        real(stm_real), intent(in)  :: xpos(nloc)   !< Cell-centered x coordinate
        real(stm_real), intent(in)  :: xlo          !< Lo side boundary of fill
        real(stm_real), intent(in)  :: xhi          !< Hi side boundary of fill
        real(stm_real), intent(in)  :: vertex_hight !< Filled value of the toppest point
        real(stm_real), intent(in)  :: fill_else    !< Filled value if not between xlo and xhi
        !---loc
        real(stm_real) :: xcenter 
        real(stm_real) :: lenght
        xcenter = half*(xhi-xlo)+xlo
        length = xhi-xlo
                
        array = fill_else
        !todo: check if you want this or something else
        where (xpos > xlo .and. (xpos < (half*length)))
            array = array + (xpos - xlo)*vertex_hight/(half*length)
        elsewhere (xpos < xhi .and. (xpos > (half*length)))
            array = array + vertex_hight - (xpos - xcenter)*vertex_hight/(half*length)
        end where
        
        
        return
    end subroutine     
    
    !> Initialize the concentration fields with a step function
    subroutine fill_discontinuity(vals,nloc,origin,dx,x0,value_lo,value_hi)

    implicit none
    integer, intent(in) ::  nloc                !< Size of array
    real(stm_real), intent(out) :: vals(nloc)   !< Values to be filled
    real(stm_real), intent(in)  :: origin       !< Low side of channel 
    real(stm_real), intent(in)  :: dx           !< dx
    real(stm_real), intent(in)  :: x0           !< Location of discontinuity todo: ??
    real(stm_real), intent(in)  :: value_lo     !< Value in the low side of discontinuity todo: ??
    real(stm_real), intent(in)  :: value_hi     !< Value in the high side of discontinuity todo: ??
    !---locals
    real(stm_real) :: fraction_lo               !discountinuity distance percent from lo side
    real(stm_real) :: fraction_hi               !discountinuity distance percent from hi side
    real(stm_real) :: xlo
    real(stm_real) :: xhi
    integer :: iloc
    
    !----------------------------------
    do iloc = 1,nloc
       xlo = origin + dble(iloc - 1)*dx
       xhi = origin + dble(iloc)*dx
       fraction_lo = (x0 - xlo)/dx
       ! tend to the usual cases where cell is entirely on lo/hi side of discontinuity
       fraction_lo = max(fraction_lo,zero)
       fraction_lo = min(fraction_lo,one)        
       fraction_hi = one - fraction_lo      
      ! need to populate using cell averages
       vals(iloc) =  (fraction_lo*value_lo + fraction_hi*value_hi)
    end do
    
    return
    end subroutine

end module
