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

!> Hydrodynamics interface to be fulfilled by driver or application
!>@ingroup transport
module hydro_data
      !> Generic interface for fetching hydrodynamic data
      interface
       !> Fill in hydrodynamic data.
       !> This data might be calculated from a function or provided by another module
       !> Note that continuity must be satisfied between time steps. The implementation
       !> must be provided by the driver or application 
        subroutine hydro_data_if(flow,    &
                                 flow_lo, &
                                 flow_hi, &
                                 area,    &
                                 area_lo, &
                                 area_hi, &
                                 ncell,   &
                                 time,    &
                                 dx,      &
                                 dt)
        use stm_precision
        implicit none
        integer, intent(in) :: ncell                   !< Number of cells
        real(stm_real), intent(in)  :: time            !< Time of request
        real(stm_real), intent(in)  :: dx              !< Spatial step 
        real(stm_real), intent(in)  :: dt              !< Time step 
        real(stm_real), intent(out) :: flow(ncell)     !< Cell and time centered flow
        real(stm_real), intent(out) :: flow_lo(ncell)  !< Lo face flow, time centered
        real(stm_real), intent(out) :: flow_hi(ncell)  !< Hi face flow, time centered
        real(stm_real), intent(out) :: area(ncell)     !< Cell center area, old time
        real(stm_real), intent(out) :: area_lo(ncell)  !< Area lo face, time centered
        real(stm_real), intent(out) :: area_hi(ncell)  !< Area hi face, time centered
        
        end subroutine
      end interface
      
 !> This pointer should be set by the driver or client code to specify the 
 !> treatment at the boundaries
 procedure(hydro_data_if),pointer :: fill_hydro_data  => null()
  !todo: some variables need to be centered in time     
      
end module