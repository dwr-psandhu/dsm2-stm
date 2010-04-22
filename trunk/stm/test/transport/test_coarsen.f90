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

!> test coarsenig subroutine subroutine
!>@ingroup test
module test_coarsening

contains

subroutine test_coarsen
use fruit
use stm_precision
use test_single_channel_advection
implicit none
  

!---arg
integer :: ncell_coarse
integer,parameter :: ncell_fine = 6
integer,parameter :: nvar = 2
real(stm_real) :: fine_data(ncell_fine,nvar)
real(stm_real), allocatable :: coarse_data(:,:)
!---- local
real(stm_real), parameter :: tol = 1.d-15

ncell_coarse = 6
allocate (coarse_data(ncell_coarse,nvar))
fine_data(:,1) = (/1:6/)
fine_data(:,2) = (/11:16/)

call coarsen(coarse_data,fine_data,ncell_fine,ncell_coarse, nvar)
call assertEquals(coarse_data(1,1),one,tol,"error in coarsening, no refinement (1,1)")
call assertEquals(coarse_data(6,2),16.d0,tol,"error in coarsening, no refinement (6,2)")
deallocate(coarse_data)

ncell_coarse=3
allocate (coarse_data(ncell_coarse,nvar))
call coarsen(coarse_data,fine_data,ncell_fine,ncell_coarse, nvar)
call assertEquals (coarse_data(1,1),1.5d0,tol,"error in coarsening (1,1)")
call assertEquals (coarse_data(3,2),15.5d0,tol,"error in coarsening constituent (2,3)")
deallocate (coarse_data)

return
end subroutine test_coarsen

end module