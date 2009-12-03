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
module hydro_data_if
      !> Generic interface for fetching hydrodynamic data
      interface hydro_data
       !> Fill in hydrodynamic data.
       !> This data might be calculated from a function or provided by another module
       !> Note that continuity must be satisfied between time steps. The implementation
       !> must be provided by the driver or application 
        subroutine hydro_data_impl(flow,flow_lo,flow_hi,area,area_lo,area_hi,ncell,time,dt)
        use stm_precision
        implicit none
        integer, intent(in) :: ncell                   !< number of cells
        real(STM_REAL), intent(in) :: time             !< time of request "old time"
        real(STM_REAL), intent(in) :: dt               !< time step for 
        real(STM_REAL), intent(out) :: flow(ncell)     !< cell and time centered flow
        real(STM_REAL), intent(out) :: flow_lo(ncell)  !< lo face flow, time centered
        real(STM_REAL), intent(out) :: flow_hi(ncell)  !< hi face flow, time centered
        real(STM_REAL), intent(out) :: area(ncell)     !< cell center area, old time
        real(STM_REAL), intent(out) :: area_lo(ncell)  !< area lo face, time centered
        real(STM_REAL), intent(out) :: area_hi(ncell)  !< area hi face, time centered
        end subroutine
      end interface
end module