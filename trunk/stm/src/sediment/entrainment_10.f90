! copyright (c) 1996, 1997, 1998, 2001, 2007, 2009 state of california,
! department of water resources.
! this file is part of dsm2.

! the delta simulation model 2 (dsm2) is free software: 
! you can redistribute it and/or modify
! it under the terms of the gnu general public license as published by
! the free software foundation, either version 3 of the license, or
! (at your option) any later version.

! dsm2 is distributed in the hope that it will be useful,
! but without any warranty; without even the implied warranty of
! merchantability or fitness for a particular purpose.  see the
! gnu general public license for more details.

! you should have received a copy of the gnu general public license
! along with dsm2.  if not, see <http://www.gnu.org/licenses>.

! delta modeling section 
! modeling support branch
! bay-delta office
! california department of water resources
! http://baydeltaoffice.water.ca.gov/modeling/deltamodeling/index.cfm
!
! purpose:
!
! it produces the sources of sediment for advection-diffusion-sink-source equation
! partial_(a c_s)/partial_t = partial/partial_x [a k_s partial_(c_s)/partial_x]
! 
! these sub-routines are being developed by the department of civil and 
! environmental engineering at the university of california, davis and dwr
!
! Record of revisions:
!       Date             Programmers                 Description of change
!       ====             ===========                 =====================
!    11/27/09           Kaveh Zamani                 Basic Structure of code 
!    11/28/09           Kaveh Zamani                 Basic Structure, Shear Velocity, Settling Velocity
!    11/29/09           Kaveh Zamani                 Miu, Nu, Salinity, Temp sub modules
!    11/30/09           Kaveh Zamani                 Garcia_Parker Formula
!    12/01/09           Kaveh Zamani                 van Rijn Formula
!    12/01/09           Kaveh Zamani                 Smith_McLean Formula
!    12/06/09           Kaveh Zamani                 Zyserman_Fredsoe Formula & debug
!    12/07/09           Kaveh Zamani                 debug
!    12/06/09           Kaveh Zamani                 debug & some clean up
!***********************************************************************************************
program entrainment

implicit none

integer,parameter:: stm_real=8

real(stm_real) :: minus = -1.d0    !< real constant -1. properly typed
real(stm_real) :: zero  =  0.d0    !< real constant  0. properly typed
real(stm_real) :: one   =  1.d0    !< real constant  1. properly typed
real(stm_real) :: two   =  2.d0    !< real constant  2. properly typed
real(stm_real) :: half   =  5.d-1  !< real constant  0.5 properly typed
real(stm_real) :: fourth =  2.5d-1 !< real constant  0.25 properly typed

real(stm_real),parameter :: g_acl =	9.800d0 ! value for San Francisco latitude and longitude
real(stm_real) :: spcfc_g  !=2.65d0 !!!!!!!!!!!!!!!!?
real(stm_real) :: rho_s  !  =2650.0d0 !!!!!!!!!!!!

real(stm_real) :: stl_v
real(stm_real) :: diameter
real(stm_real) :: d_star
real(stm_real) :: nu
real(stm_real) :: mu
real(stm_real) :: rho_0
real(stm_real) :: temp ! = 20.0d0 !!!!!!
real(stm_real) :: salinity !  = 0.0d0 !!!!!!!!!!!!!!!
real(stm_real) :: u_star
real(stm_real) :: v_lo
real(stm_real) :: v_hi
real(stm_real) :: hyd_r
real(stm_real) :: manning_n
real(stm_real) :: re_p
real(stm_real) :: area_lo
real(stm_real) :: area_hi
real(stm_real) :: d50                            !  =256.0d-3 !!!!!!!!!!!!!!!!!!!!!!!!!!!
real(stm_real) :: e_source
real(stm_real) :: dpt
real(stm_real) :: ta_c_star
real(stm_real) :: b_prmtr
real(stm_real) :: ks_bed 



character (len=8) :: entr_form ! !"grc_prkr"  "van_rijn"  "smt_mcln"  "zsr_frds"


open(3,file='input_entrainment.dat',action='READ')

	read(3,*) 	
	read(3,*) d50	
	read(3,*) diameter	
	read(3,*) rho_s	
	read(3,*) temp
	read(3,*) salinity
	read(3,*) spcfc_g
	read(3,*) entr_form
	read(3,*)	
	
close (3)



call  knm_viscosity (nu,rho_0 , &
                               salinity,   temp)
                                
call settling_velocity (stl_v,d_star ,&  
                                         nu,diameter)

call shr_v_re_prtcl (u_star,re_p,ta_c_star ,&
                                             v_lo,v_hi,hyd_r,manning_n,g_acl,d50)
                                             
call hydro_input_source(v_lo,v_hi,hyd_r,manning_n,area_lo,area_hi,dpt,ks_bed   &
                                                                                          )
                                             
                                                  
                                
select case  (entr_form)

    case ('grc_prkr')                             
                                
        call gp_entrainment (e_source, & 
                                        u_star,re_p,stl_v) 
                                        
    case ('van_rijn')
        
        
        call  vr_entrainment (e_source, & 
                                            d50,  dpt,ks_bed, ta_c_star, d_star,rho_0,v_lo,v_hi,g_acl, rho_s,diameter)

    case ('smt_mcln')
    
        call sm_entrainment (e_source, & 
                                            d50,  dpt,ks_bed, ta_c_star, d_star,rho_0,v_lo,v_hi,g_acl, rho_s,diameter)

    case ('zsr_frds')
    
        call zf_entrainment (e_source, & 
                                            d50,  dpt,ks_bed, ta_c_star, d_star,rho_0,v_lo,v_hi,g_acl, rho_s,diameter)
                                     
end select

open (4,file='Entrainment.dat',status='unknown')

write (4,*)'Entrainment Formula is :  ', entr_form
write (4,*)  'E = ' ,e_source
write (4,*) 'Temp = ', temp
write (4,*) 'Salinity = ' , salinity
write (4,*)  'Settleing Velocity = ',stl_v
write (4,*)  'Dynamic Viscosity = ', nu

end program entrainment 



!************************************
subroutine settling_velocity (stl_v,d_star ,&  
                                                nu,diameter)

    integer,parameter:: stm_real=8

    real(stm_real),parameter :: spcfc_g =2.65d0 !!!!!!!!!!!!! ???
    real(stm_real),parameter :: g_acl =	9.800d0 ! value for San Francisco latitude and longitude

    real(stm_real) :: stl_v
    real(stm_real) :: diameter
    real(stm_real) :: d_star
    real(stm_real) :: nu


    !! we should add cohesive in future
    call knm_viscosity (nu,rho_0 , &
                                 salinity ,temp)

    d_star= diameter * (g_acl*(spcfc_g-1.0d0)/nu/nu)**(1.0d0/3.0d0)
    stl_v =8.0d0*nu*(sqrt(1+0.0139d0*(d_star**3.0d0))-1.0d0)/diameter

end subroutine settling_velocity
!************************************
subroutine knm_viscosity (nu,rho_0 ,&
                                salinity ,temp)
    integer,parameter:: stm_real=8

    real(stm_real) :: minus = -1.d0    !< real constant -1. properly typed
    real(stm_real) :: zero  =  0.d0    !< real constant  0. properly typed
    real(stm_real) :: one   =  1.d0    !< real constant  1. properly typed
    real(stm_real) :: two   =  2.d0    !< real constant  2. properly typed
    real(stm_real) :: half   =  5.d-1  !< real constant  0.5 properly typed
    real(stm_real) :: fourth =  2.5d-1 !< real constant  0.25 properly typed

    real(stm_real) :: nu
    real(stm_real) :: mu
    real(stm_real) :: rho
    real(stm_real) :: temp
    real(stm_real) :: salinity
    real(stm_real) :: rho_0

    rho_0 = 999.842594d0 + 0.06793952d0*temp - 0.00909529d0*temp**two +0.0001685*temp**3.0d0 -0.00000120083*temp**(two*two) + 0.000000006536322*temp**5.0d0
    rho_0 = rho_0 + salinity *(0.824493d0 - 0.0040899d0*temp +0.000076438d0*temp**two - 0.00000082467*temp**3.0d0 +0.0000000053785*temp**(4.0d0) )
    rho_0 = rho_0 +(-0.00572466d0 + 0.00010227d0*temp -0.0000016546*temp**two)*salinity**(1.5d0)
    rho_0 = rho_0 + 0.00048314*salinity**two
    
    if (temp<1) then
        mu = 1.79d-3
    elseif (temp<4) then
        mu = 1.65d-3
    elseif (temp<6) then
        mu = 1.51d-3
    elseif (temp<8) then
        mu = 1.41d-3
    elseif (temp<11) then
        mu = 1.31d-3
    elseif (temp<14) then
        mu = 1.22d-3
    elseif (temp<16) then
        mu = 1.14d-3
    elseif  (temp<19)then
        mu = 1.07d-3
    elseif (temp<21) then
        mu = 1.0d-3
    elseif (temp<24) then
        mu = 9.1d-4
    elseif (temp<26) then 
        mu =8.91d-4
    elseif (temp<29) then 
        mu =8.43d-4
    else
        mu =7.96d-4 
    endif

    nu= mu/rho_0


end subroutine knm_viscosity
!*********************************
subroutine shr_v_re_prtcl (u_star,re_p,ta_c_star, &
                                v_lo,v_hi,hyd_r,manning_n,g_acl,d50)
                             
    integer,parameter:: stm_real=8

    real(stm_real) :: minus = -1.d0    !< real constant -1. properly typed
    real(stm_real) :: zero  =  0.d0    !< real constant  0. properly typed
    real(stm_real) :: one   =  1.d0    !< real constant  1. properly typed
    real(stm_real) :: two   =  2.d0    !< real constant  2. properly typed
    real(stm_real) :: half   =  5.d-1  !< real constant  0.5 properly typed
    real(stm_real) :: fourth =  2.5d-1 !< real constant  0.25 properly typed

    real(stm_real) :: v_lo
    real(stm_real) :: v_hi
    real(stm_real) :: hyd_r
    real(stm_real) :: manning_n
    real(stm_real) :: u_star
    real(stm_real) :: d50
    real(stm_real) :: nu
    real(stm_real) :: re_p
    real(stm_real) :: dpt
    real(stm_real) :: g_acl 
    real(stm_real) :: ta_c_star

    call knm_viscosity (nu,rho_0, &
                                salinity , temp)
                                    
    call hydro_input_source(v_lo,v_hi,hyd_r,manning_n,area_lo,area_hi,dpt,ks_bed &
                                                                                     )
                                                                               
                                                                               
    u_star = (v_lo+v_hi)*manning_n*sqrt(g_acl)/(2*hyd_r**(1.0d0/6.0d0))
    re_p= sqrt(hyd_r*g_acl*d50)*d50/nu
    ta_c_star = 0.22d0*(re_p**(-0.6))+ 0.06d0*(10.0d0**(-7.7d0*re_p**(-0.6d0)))
    
                                                                
                          
end subroutine shr_v_re_prtcl
!*****************************************
subroutine hydro_input_source(v_lo,v_hi,hyd_r,manning_n,area_lo,area_hi,dpt,ks_bed &
                                                                                         )
    integer,parameter:: stm_real=8

    real(stm_real) :: v_lo
    real(stm_real) :: v_hi
    real(stm_real) :: hyd_r
    real(stm_real) :: manning_n
    real(stm_real) :: area_lo
    real(stm_real) :: area_hi
    real(stm_real) :: dpt
    real(stm_real) :: ks_bed


    v_lo = 1.0d0
    v_hi = 1.0d0
    hyd_r = 2.5d0
    area_lo=100d0
    area_hi=100d0
    dpt =2.0d0
    ks_bed = 0.1d0
    manning_n = 0.02d0 !!!!!!!!! SI and British
                                                           
end subroutine  hydro_input_source 
!***********************************  
subroutine gp_entrainment (e_source, & 
                                 u_star,re_p,stl_v)

    integer,parameter:: stm_real=8
                                     
    real(stm_real) :: u_star
    real(stm_real) :: stl_v
    real(stm_real) :: re_p
    real(stm_real) :: e_source

    real(stm_real),parameter :: a_gp =1.3d-7 !!!!!!! local
    real(stm_real) :: zu 


    zu = u_star*sqrt(re_p**0.6d0)/stl_v
    e_source = (a_gp*zu**5)/(1+(a_gp*zu**5)/0.3d0)

!    print *,'G_P E= ',e_source                                
!    Pause
                                 
end subroutine gp_entrainment          
!**********************************
subroutine vr_entrainment (e_source, & 
                                            d50,  dpt,ks_bed, ta_c_star, d_star,rho_0,v_lo,v_hi,g_acl, rho_s,diameter)

                                            
 integer,parameter:: stm_real=8
 
    real(stm_real) :: e_source
    real(stm_real) :: dpt
    real(stm_real) :: ta_c_star
    real(stm_real) :: b_prmtr
    real(stm_real) :: ks_bed !!!!!!!!!!!!!!!!!!! it can be 3 D50
    real(stm_real) :: d50
    real(stm_real) :: d_star
    real(stm_real) :: ta_bs
    real(stm_real) :: c_fs
    real(stm_real) :: rho_0
    real(stm_real) :: rho_s
    real(stm_real) :: v_lo
    real(stm_real) :: v_hi
    real(stm_real) :: ta_s_star
    real(stm_real) :: diameter
    
    real(stm_real) :: spc_g
    real(stm_real) :: g_acl
    
    real(stm_real),parameter :: von_krmn =0.41d0
                                               
                                           
b_prmtr = max (0.01*dpt,ks_bed)
ks_bed= 3.0d0*d50

c_fs =(1.0d0/von_krmn)*log(12*dpt/ks_bed)**(-2.0d0)   !!!!! I have qustion about formula it self
ta_bs = (1.0d0/4.0d0)*rho_0*c_fs* (v_hi + v_lo)**2.0d0
ta_s_star = ta_bs/(rho_s - rho_0)/diameter/g_acl


!print *, 'v' , v_hi,v_lo
!print*,'d50',d50, g_acl
!print *,'b',b_prmtr,'tbs',ta_bs
!print *, ta_c_star, ta_s_star ,rho_0,rho_s
!print *,d_star


e_source = 0.015d0*d50*((ta_s_star/ta_c_star - 1.0d0)**1.5d0) /b_prmtr/ (d_star**0.3d0)

!print *,'E=van_Rijn',e_source
!pause

!!!!!!!! OK I check it with hand calculator

end subroutine vr_entrainment
!**************************************************
subroutine sm_entrainment (e_source, & 
                                            d50,  dpt,ks_bed, ta_c_star, d_star,rho_0,v_lo,v_hi,g_acl, rho_s,diameter)

                                            
 integer,parameter:: stm_real=8
 
    real(stm_real) :: e_source
    real(stm_real) :: dpt
    real(stm_real) :: ta_c_star
    real(stm_real) :: b_prmtr
    real(stm_real) :: ks_bed !!!!!!!!!!!!!!!!!!! it can be 3 D50
    real(stm_real) :: d50
    real(stm_real) :: d_star
    real(stm_real) :: ta_bs
    real(stm_real) :: c_fs
    real(stm_real) :: rho_0
    real(stm_real) :: rho_s
    real(stm_real) :: v_lo
    real(stm_real) :: v_hi
    real(stm_real) :: ta_s_star
    real(stm_real) :: diameter
    
    real(stm_real) :: spc_g
    real(stm_real) :: g_acl
    
    real(stm_real),parameter :: von_krmn =0.41d0
    
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ta_bs = rho_0 * (u_star ** 2.0)
!ta_s_star = ta_bs/(rho_s - rho_0)/d50/g_acl   !!! is it d 50 ro diameter? 
!
!print *,'ta_s_star with ta_bs =ro u star 2', ta_s_star
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                               
                                           
b_prmtr = max (0.01*dpt,ks_bed)
ks_bed= 3.0d0*d50

c_fs =(1.0d0/von_krmn)*log(12*dpt/ks_bed)**(-2.0d0)   !!!!! I have qustion about formula it self
ta_bs = (1.0d0/4.0d0)*rho_0*c_fs* (v_hi + v_lo)**2.0d0
ta_s_star = ta_bs/(rho_s - rho_0)/diameter/g_acl

e_source = 0.65d0* 0.0024d0*(ta_s_star/ta_c_star - 1.0d0)/(1.0d0 + 0.0024d0*(ta_s_star/ta_c_star - 1.0d0))
!print *,'E=smith McLean',e_source
!pause

end subroutine sm_entrainment
!**************************************************
subroutine zf_entrainment (e_source, & 
                                            d50,  dpt,ks_bed, ta_c_star, d_star,rho_0,v_lo,v_hi,g_acl, rho_s,diameter)

                                            
 integer,parameter:: stm_real=8
 
    real(stm_real) :: e_source
    real(stm_real) :: dpt
    real(stm_real) :: ta_c_star
    real(stm_real) :: b_prmtr
    real(stm_real) :: ks_bed !!!!!!!!!!!!!!!!!!! it can be 3 D50
    real(stm_real) :: d50
    real(stm_real) :: d_star
    real(stm_real) :: ta_bs
    real(stm_real) :: c_fs
    real(stm_real) :: rho_0
    real(stm_real) :: rho_s
    real(stm_real) :: v_lo
    real(stm_real) :: v_hi
    real(stm_real) :: ta_s_star
    real(stm_real) :: diameter
    
    real(stm_real) :: spc_g
    real(stm_real) :: g_acl
    
    real(stm_real),parameter :: von_krmn =0.41d0
                                               
                                           
b_prmtr = max (0.01*dpt,ks_bed)
ks_bed= 3.0d0*d50

c_fs =(1.0d0/von_krmn)*log(12*dpt/ks_bed)**(-2.0d0)   !!!!! I have qustion about formula it self
ta_bs = (1.0d0/4.0d0)*rho_0*c_fs* (v_hi + v_lo)**2.0d0
ta_s_star = ta_bs/(rho_s - rho_0)/diameter/g_acl

e_source = 0.331d0*((ta_s_star-ta_c_star)**1.75d0) /(1.0d0 + 0.72d0*(ta_s_star-ta_c_star)**1.75d0)

!print *,'E=Zyserman Fredsoe',e_source
!pause

end subroutine zf_entrainment
!*******************************************


!subroutine sm_entrainment (e_source , & 
!                                        u_star, rho_0 ,rho_s,d50,g_acl,ta_c_star)
!
!integer,parameter:: stm_real=8
!
!real(stm_real) :: d50
!real(stm_real) :: rho_0
!real(stm_real) :: rho_s
!real(stm_real) :: g_acl
!real(stm_real) :: u_star
!real(stm_real) :: ta_bs
!real(stm_real) :: ta_s_star
!real(stm_real) :: ta_c_star
!real(stm_real) :: e_source
!real(stm_real) :: ks_bed 
!real(stm_real) :: dpt
!real(stm_real) :: hs
!real(stm_real) :: v_lo
!real(stm_real) :: v_hi
!real(stm_real) :: wslope
!
!real(stm_real),parameter :: von_krmn =0.41d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ta_bs = rho_0 * (u_star ** 2.0)
!ta_s_star = ta_bs/(rho_s - rho_0)/d50/g_acl   !!! is it d 50 ro diameter? 
!
!print *,'ta_s_star with ta_bs =ro u star 2', ta_s_star
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!hs = (v_lo/2.0d0 + v_hi/2.0d0)**2 /g_acl/wslope
!hs = hs* 
!
!
!
!
!
!!print *, 'v' , v_hi,v_lo
!!print*,'d50',d50, g_acl
!!print *,'b',b_prmtr,'tbs',ta_bs
!!print *, ta_c_star, ta_s_star ,rho_0,rho_s
!!print *,d_star
!
!
!e_source = 0.65d0* 0.0024d0*(ta_s_star/ta_c_star - 1.0d0)/(1.0d0 + 0.0024d0*(ta_s_star/ta_c_star - 1.0d0))
!print *,'ta_c_s',ta_c_star,'ta_s_s',ta_s_star
!print *,'e_sm =',e_source
!pause
!end subroutine sm_entrainment
!**************************************************


