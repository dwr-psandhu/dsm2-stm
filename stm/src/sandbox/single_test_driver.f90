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
!    along with DSM2. If not, see <http://www.gnu.org/licenses>.
!</license>

!> Main program unit for testing a single test in transport
!>@ingroup test
program single_test_driver

use fruit
use unit_test_suspend_sed_utility
use test_non_cohesive
use test_bed_load

implicit none
logical :: verbose = .true.

call init_fruit

call test_submerged_specific_gravity
call test_explicit_particle_reynolds_number
call test_particle_reynolds_number
call test_dimless_particle_diameter
call test_critical_shields_parameter
call test_settling_velocity
call test_shear_velocity
call test_rouse_number

call test_allocation_ratio()

! Non_cohesive sink source
call test_first_einstein_integral
call test_es_garcia_parker
!bedload
call test_volumetric_bedload_transport_rate


!!!!!!!!!!!!!!!!
call fruit_summary


pause
end program 