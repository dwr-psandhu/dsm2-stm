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
        
        real(stm_real), intent(in) :: x0   !< end of integration
        real(stm_real), intent(in) :: mean !< mean
        real(stm_real), intent(in) :: sd   !< sigma/standard deviation
        gaussian_cdf = half + half*erf((x0-mean)/(sqrt(two)*sd))
        return
    end function
    
    
    !> Fill array with 1D gaussian shape
    !> This routine expects a 1D array, so multi-constituents
    !> have to be initialized separately
    
    subroutine fill_gaussian(vals,nloc,origin,dx,mean,sd)

        implicit none
        integer, intent(in) :: nloc                !< number of cells (size of array) 
        real(stm_real), intent(out) :: vals(nloc)  !< values to be filled
        real(stm_real), intent(in)  :: origin      !< origin (lo side of channel)
        real(stm_real), intent(in)  :: dx          !< dx
        real(stm_real), intent(in)  :: mean        !< center of the gaussian shape
        real(stm_real), intent(in)  :: sd          !< length of gaussian shape
        !-----locals
        real(stm_real) :: xlo
        real(stm_real) :: xhi
        integer        :: iloc
        !-----------
        do iloc = 1,nloc
           xlo = origin + dble(iloc - 1)*dx
           xhi = origin + dble(iloc)*dx
          ! need to populate using cell averages
           vals(iloc) =  (gaussian_cdf(xhi,mean,sd) & 
                        - gaussian_cdf(xlo,mean,sd))
           vals(iloc)=vals(iloc)*sqrt(two*pi*sd*sd)/dx   !todo: move out of loop and make pi a constant instead of acos(zero)
        end do
        return
    end subroutine
    
    
    !> Fill array with rectanglar shape
    !> This routine expects a 1D array, so multi-constituents
    !> have to be initialized separately
    subroutine fill_rectangular(array,x,nloc,xlo,xhi,fill,fill_else)

        implicit none
        real(stm_real), intent(out) :: array(nloc)  !< array to be filled
        real(stm_real), intent(in)  :: x(nloc)      !< cell-centered x coordinate
        integer,        intent(in)  :: nloc         !< size of array
        real(stm_real), intent(in)  :: xlo          !< lo side boundary of fill
        real(stm_real), intent(in)  :: xhi          !< hi side boundary of fill
        real(stm_real), intent(in)  :: fill         !< filled value between xlo and xhi
        real(stm_real), intent(in)  :: fill_else    !< filled value if not between xlo and xhi
                
        array = fill_else
        
        where (x > xlo .and. x< xhi)
            array = fill
        end where
        
        return
    end subroutine    

    subroutine fill_triangular(array,x,nloc,xlo,xhi,fill,fill_else)

        implicit none
        real(stm_real), intent(out) :: array(nloc)  !< array to be filled
        real(stm_real), intent(in)  :: x(nloc)      !< cell-centered x coordinate
        integer,        intent(in)  :: nloc         !< size of array
        real(stm_real), intent(in)  :: xlo          !< lo side boundary of fill
        real(stm_real), intent(in)  :: xhi          !< hi side boundary of fill
        real(stm_real), intent(in)  :: fill         !< filled value between xlo and xhi
        real(stm_real), intent(in)  :: fill_else    !< filled value if not between xlo and xhi
                
        array = fill_else
        
        where (x > xlo .and. x< xhi)
            array = fill
        end where
        
        return
    end subroutine     
    
    !> Initialize the concentration fields with a step function
    subroutine fill_discontinuity(vals,nloc,origin,dx,x0,value_lo,value_hi)

    implicit none
    integer, intent(in) ::  nloc                !< size of array
    real(stm_real), intent(out) :: vals(nloc)   !< values to be filled
    real(stm_real), intent(in)  :: origin       !< low side of channel 
    real(stm_real), intent(in)  :: dx           !< dx
    real(stm_real), intent(in)  :: x0           !< location of discontinuity todo: ??
    real(stm_real), intent(in)  :: value_lo     !< value in the low side of discontinuity todo: ??
    real(stm_real), intent(in)  :: value_hi     !< value in the high side of discontinuity todo: ??
    !---locals
    real(stm_real) :: fraction_lo
    real(stm_real) :: fraction_hi
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
