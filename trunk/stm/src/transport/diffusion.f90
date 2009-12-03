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

!> Explicit and implicit diffusion operators go here
!>@ingroup transport
module diffusion

! We need a single step routine with inputs and outputs that
! do not involve computation detail (like matrices)
! This is just stolen from the interface for advection, 
! but the real diffusion API will look similar
!subroutine diffuse(mass,     &
!                  mass_prev,&
!                  flow,     &                  
!                  flow_lo,  &
!                  flow_hi,  &
!                  area,     &
!                  area_prev,&
!                  area_lo,  &
!                  area_hi,  &
!                  ncell,    &
!                  nvar,     &
!                  time,     &
!                  dt,       &
!                  dx)

! This routine should give the effects of diffusion fluxes on each cell
! for a single time step (ie, explicit). This is needed for the advection step.
! It is also probably part of the right hand side of the implicit diffusion solver 
! matrix calculation. 
!subroutine diffusion_operator(mass,     &
!                  mass_prev,&
!                  flow,     &                  
!                  flow_lo,  &
!                  flow_hi,  &
!                  area,     &
!                  area_prev,&
!                  area_lo,  &
!                  area_hi,  &
!                  ncell,    &
!                  nvar,     &
!                  time,     &
!                  dt,       &
!                  dx)

! The rest of this should be neat, but is not important to the public use of the library


end module