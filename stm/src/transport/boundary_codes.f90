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

! todo: what does this module exactly do? and it should be here or in driver
!> Module provides boundary condition flags
!> routine in the module is advection().
!>@ingroup transport

module boundary_codes
!> sets the BC type
interface set_bc_type

subroutine set_boundary_codes (boundary_flag)
 implicit none 
 
 integer,intent(out)  :: boundary_flag

end subroutine set_boundary_codes
end interface

 procedure(set_boundary_codes),pointer :: boundary_flag => null()

contains

subroutine set_boundary_neumann(boundary_flag)

    implicit none
    ! todo : inout or out
    integer,intent(inout)  :: boundary_flag
    !---local
    integer,parameter  :: neumann = 2

    return
end subroutine

subroutine set_boundary_dirichlet(boundary_flag)

    implicit none
    ! todo : inout or out
    integer,intent(inout)  :: boundary_flag
    !--local
    integer,parameter  :: dirichlet = 2

    return
end subroutine

end module boundary_codes