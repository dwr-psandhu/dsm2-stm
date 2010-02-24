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

!> Testing of the extrapolate step.
!>@ingroup test
module test_extrapolate
use fruit

contains
!///////////////////////////////////////

!> Test the extrapolation part of the predictor step of the advection algorithm
!> There are three cells and two constituents to test symmetry in space across constituents
subroutine test_extrapolation
use stm_precision
use advection
implicit none
  integer,parameter :: nx = 3       !interior and two ends
  integer,parameter :: nconst = 2
 
  real(stm_real) :: grad(nx,nconst)
  real(stm_real) :: conc(nx,nconst)
  real(stm_real) :: conc_hi(nx,nconst)
  real(stm_real) :: conc_lo(nx,nconst)
  real(stm_real) :: source(nx,nconst)

  real(stm_real) :: flow(nx)
  real(stm_real) :: area(nx)  
  real(stm_real) :: dx
  real(stm_real) :: dt
  real(stm_real) :: time
  integer        :: ix,ic
  
  conc(1,1)=62.D0
  conc(2,1)=64.D0
  conc(3,1)=65.D0
  conc(1,2)=65.D0
  conc(2,2)=64.D0
  conc(3,2)=62.D0

  grad(1,1)=2.D0
  grad(2,1)=1.5D0
  grad(3,1)=1.D0
  grad(1,2)=-1.D0
  grad(2,2)=-1.5D0
  grad(3,2)=-2.D0
  
  source = 0.D0
  flow(1) = 8.D0
  flow(2) = zero
  flow(3) = -8.D0
  area = 16.D0
  dx = 32.D0
  dt = 48.D0
  time = zero

  call extrapolate(conc_lo,  &
                   conc_hi,  &
                   conc,     &
                   grad,     &
                   source,   &
                   flow,     &  
                   area,     &
                   nx,       &
                   nconst,   &
                   time,     &
                   dt,       &
                   dx)

  call assertEquals(conc_lo(1,1),60.25D0,"extrapolate, lo(1,1)")
  call assertEquals(conc_lo(2,1),63.25D0,"extrapolate, lo(2,1)")
  call assertEquals(conc_lo(3,1),64.875D0,"extrapolate, lo(3,1)")
  call assertEquals(conc_lo(1,1),conc_hi(3,2),"extrapolate, lo-hi(1,1)")
  call assertEquals(conc_lo(2,1),conc_hi(2,2),"extrapolate, lo-hi(2,1)")
  call assertEquals(conc_lo(3,1),conc_hi(1,2),"extrapolate, lo-hi(3,1)")

  call assertEquals(conc_hi(1,1),62.25D0,"extrapolate, hi(1,1)")
  call assertEquals(conc_hi(2,1),64.75D0,"extrapolate, hi(2,1)")
  call assertEquals(conc_hi(3,1),65.875D0,"extrapolate, hi(3,1)")
  call assertEquals(conc_hi(1,1),conc_lo(3,2),"extrapolate, hi-lo(3,2)")
  call assertEquals(conc_hi(2,1),conc_lo(2,2),"extrapolate, hi-lo(2,2)")
  call assertEquals(conc_hi(3,1),conc_lo(1,2),"extrapolate, hi-lo(1,2)")

return
end subroutine



!///////////////////////////////////////

!> Test the extrapolation part of the predictor step of the advection algorithm
!> There are three cells and two constituents to test symmetry in space across constituents
subroutine test_flux_calculation
use stm_precision
use advection
implicit none
  integer,parameter :: nx = 3       !interior and two ends
  integer,parameter :: nconst = 2
 
  real(stm_real) :: flow_lo(nx)
  real(stm_real) :: flow_hi(nx)  
  real(stm_real) :: conc(nx,nconst)
  real(stm_real) :: conc_hi(nx,nconst)
  real(stm_real) :: conc_lo(nx,nconst)
  real(stm_real) :: flux_lo(nx,nconst)
  real(stm_real) :: flux_hi(nx,nconst)
  real(stm_real) :: flow(nx)
  real(stm_real) :: area(nx)  
  real(stm_real) :: dx
  real(stm_real) :: dt
  real(stm_real) :: time
  integer        :: ix,ic
  
  conc_lo(1,1)=62.D0
  conc_lo(2,1)=64.D0
  conc_lo(3,1)=65.D0
  conc_lo(1,2)=65.D0
  conc_lo(2,2)=64.D0
  conc_lo(3,2)=62.D0
  
  conc_hi(1,1)=64.D0
  conc_hi(2,1)=66.D0
  conc_hi(3,1)=68.D0
  conc_hi(1,2)=61.D0
  conc_hi(2,2)=60.D0
  conc_hi(3,2)=66.D0
    
  flow_lo(1) = 1.D0
  flow_lo(2) = 2.D0
  flow_lo(3) = -1.D0
  flow_hi(1) = 2.0D0
  flow_hi(2) = -1.D0
  flow_hi(3) = -4.D0
  

 ! Compute upwind value of fluxes. This is a naive guess based on the extrapolated states
 ! It doesn't include any node-based sources or reservoirs or the like.
 call compute_flux(flux_lo,  &
                   flux_hi,  &
                   conc_lo,  &
                   conc_hi,  &                       
                   flow_lo,  &
                   flow_hi,  &
                   nx,       &
                   nconst    &
                   )

  call assertEquals(flux_lo(1,1),LARGEREAL,"flux computation, lo(1,1)")
  call assertEquals(flux_lo(2,1),1.28D2,"flux computation, lo(2,1)")
  call assertEquals(flux_hi(3,1),LARGEREAL,"flux computation, hi(3,1)")
  call assertEquals(flux_lo(3,1),-6.2D1,"flux computation, hi(3,1)")

  flow_lo(1) = -1.D0
  flow_hi(3) =  4.D0
  call compute_flux(flux_lo,  &
                   flux_hi,  &
                   conc_lo,  &
                   conc_hi,  &                       
                   flow_lo,  &
                   flow_hi,  &
                   nx,       &
                   nconst    &
                   )
  call assertEquals(flux_lo(1,1),-6.2D1,"flux computation, outflow, lo(1,1)")
  call assertEquals(flux_hi(3,1),2.72D2,"flux computation, outflow, hi(3,1)")

return
end subroutine




end module



