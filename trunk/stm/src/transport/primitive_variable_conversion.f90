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
!> Module containing routines for converting mass to concentraion, and concentration to mass
!>@ingroup transport
module primitive_variable_conversion

contains
!> Convert conservative variables (ie mass) to primitive (concentration)
pure subroutine cons2prim(conc,mass,area,nloc,nvar)

use stm_precision

implicit none
real(stm_real),intent(out) :: conc(nloc,nvar)   !< concentration (converted from mass per unit length )
real(stm_real),intent(in)  :: mass(nloc,nvar)   !< mass per unit length 
real(stm_real),intent(in)  :: area(nloc)        !< area at conversion locations
!--- args
integer,intent(in)  :: nloc                     !< Number of cells or faces
integer,intent(in)  :: nvar                     !< Number of variables
!--- locals
integer :: ivar
!-------------------
do ivar = 1,nvar
    conc(:,ivar) = mass(:,ivar)/area
end do

return
end subroutine

!> Convert  primitive (concentration) to conservative variables (ie mass)
pure subroutine prim2cons(mass,conc,area,nloc,nvar)

use stm_precision

implicit none
real(stm_real),intent(out) :: mass(nloc,nvar)  !< mass per unit length (converted from concentration)
real(stm_real),intent(in)  :: conc(nloc,nvar)  !< concentrations to convert
real(stm_real),intent(in)  :: area(nloc)       !< area at conversion locations
!--- args
integer,intent(in)  :: nloc                    !< Number of cells or faces
integer,intent(in)  :: nvar                    !< Number of variables
!--- locals
integer :: ivar
!-------------------

do ivar = 1,nvar
    mass(:,ivar) = conc(:,ivar)*area
end do

return
end subroutine

end module


