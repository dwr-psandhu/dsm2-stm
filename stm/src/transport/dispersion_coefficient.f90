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

!> dispersion coefficient interface to be fulfilled by driver or application
!>@ingroup transport

module dispersion_coefficient

interface
    subroutine diffusion_coef_if(disp_coef,            &
                                 disp_coef_lo,         &
                                 disp_coef_hi,         &
                                 flow,                 &
                                 flow_lo,              &
                                 flow_hi,              &
                                 time,                 &
                                 dx,                   &
                                 dt,                   &
                                 origin,               &
                                 ncell,                &
                                 nvar,                 &
                                 const_disp_coef)      
                                  
    use stm_precision    
    implicit none
    
    integer,intent(in)  :: ncell                         !< Number of cells
    integer,intent(in)  :: nvar                          !< Number of variables   
    real(stm_real),intent(in) :: time                    !< Current time
    real(stm_real),intent(in) :: dx                      !< Spatial step  
    real(stm_real),intent(in) :: dt                      !< Time step 
    real(stm_real),intent(in) :: origin                  !< Left side of the channel
    real(stm_real),intent(in) :: flow_lo(ncell)          !< flow on lo side of cells centered in time
    real(stm_real),intent(in) :: flow_hi(ncell)          !< flow on hi side of cells centered in time       
    real(stm_real),intent(in) :: flow(ncell)             !< flow on center of cells 
    real(stm_real),intent(out):: disp_coef(ncell,nvar)   !< center constituent dispersion coef. at new time
    real(stm_real),intent(out):: disp_coef_lo(ncell,nvar)!< Low side constituent dispersion coef. at new time
    real(stm_real),intent(out):: disp_coef_hi(ncell,nvar)!< High side constituent dispersion coef. at new time
    ! todo: should it be here     
    real(stm_real),intent(in),optional :: const_disp_coef(nvar)   !< Constant value of dispersion coef. 
           
    end subroutine diffusion_coef_if
end interface

!> This pointer should be set by the driver or client code to specify the 
!> treatment at the dispersion coefficients
procedure(diffusion_coef_if),pointer :: dispersion_coef  => null()

contains

subroutine set_constant_dispersion(disp_coef_lo,         &
                                   disp_coef_hi,         &
                                   flow,                 &
                                   flow_lo,              &
                                   flow_hi,              &
                                   time,                 &
                                   dx,                   &
                                   dt,                   &
                                   origin,               &
                                   ncell,                &
                                   nvar,                 &
                                   const_disp_coef)  
     
     use stm_precision
     use error_handling
     
     implicit none
      !--- args          
    integer,intent(in)  :: ncell                         !< Number of cells
    integer,intent(in)  :: nvar                          !< Number of variables   
    real(stm_real),intent(in) :: time                    !< Current time
    real(stm_real),intent(in) :: dx                      !< Spatial step  
    real(stm_real),intent(in) :: dt                      !< Time step 
    real(stm_real),intent(in) :: origin                  !< Left side of the channel
    real(stm_real),intent(in) :: flow_lo(ncell)          !< flow on lo side of cells centered in time
    real(stm_real),intent(in) :: flow_hi(ncell)          !< flow on hi side of cells centered in time       
    real(stm_real),intent(in) :: flow(ncell)             !< flow on center of cells 
    real(stm_real),intent(out):: disp_coef_lo(ncell,nvar)!< Low side constituent dispersion coef. at new time
    real(stm_real),intent(out):: disp_coef_hi(ncell,nvar)!< High side constituent dispersion coef. at new time
    real(stm_real),intent(in) :: const_disp_coef(nvar)   !< Constant value of dispersion coef. 
    !--
    integer :: ivar
   
         do ivar=1,nvar
             disp_coef_lo(:,ivar)=const_disp_coef(ivar)
         end do
         disp_coef_hi = disp_coef_lo
        
     return
 end subroutine


end module