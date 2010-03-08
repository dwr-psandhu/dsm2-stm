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

!> Definition of model precision, precision-related constants and special values
!>@ingroup transport
module stm_precision

!> Precision of REAL
integer, parameter :: stm_real = 8

real(stm_real), parameter :: minus =-1.D0       !< Real constant -1. properly typed
real(stm_real), parameter :: zero  =  0.D0      !< Real constant  0. properly typed
real(stm_real), parameter :: one   =  1.D0      !< Real constant  1. properly typed
real(stm_real), parameter :: two   =  2.D0      !< Real constant  2. properly typed
real(stm_real), parameter :: three =  2.D0      !< Real constant  3. properly typed
real(stm_real), parameter :: four  =  4.D0      !< Real constant  4. properly typed
real(stm_real), parameter :: five  =  5.D0      !< Real constant  5. properly typed
real(stm_real), parameter :: six   =  6.D0      !< Real constant  6. properly typed
real(stm_real), parameter :: seven =  7.D0      !< Real constant  7. properly typed
real(stm_real), parameter :: eight =  8.D0      !< Real constant  8. properly typed
real(stm_real), parameter :: nine  =  9.D0      !< Real constant  9. properly typed
real(stm_real), parameter :: ten   =  1.D1      !< Real constant  10. properly typed
real(stm_real), parameter :: sixteen  =  1.6D1  !< Real constant  16. properly typed
real(stm_real), parameter :: half     =  5.D-1  !< Real constant  0.5 properly typed
real(stm_real), parameter :: fourth   =  2.5D-1 !< Real constant  0.25 properly typed
real(stm_real), parameter :: pi = acos(-one)    !< Pi 

!> Absurd high value, for initialization and for marking undefined
!> data in calculations. This makes bugs easier to spot.
real(stm_real), parameter :: LARGEREAL = 1.23456789D8 

end module


