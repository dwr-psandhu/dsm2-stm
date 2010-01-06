subroutine  explicit_diffusion_operator (explicit_diffusion_term, &
                                                        conc,             &
                                                        conc_prev,        &
                                                        mass,             &
                                                        mass_prev,        & 
                                                        area,             &      
                                                        area_lo,          &
                                                        area_hi,          &
                                                        ks_lo,            &
                                                        ks_hi,            &
                                                        ncell,            &
                                                        nvar,             &
                                                        time,             &
                                                        strt_flx_bc_prev, &   
                                                        end_flx_bc_prev,  &
                                                        dt,               &
                                                        dx)
                                                                                          
use stm_precision

implicit none 

!----args

integer, intent (in) :: ncell !< Number of cells
integer, intent (in) :: nvar  !< Number of variables


real(stm_real), intent (out) :: explicit_diffusion_term(ncell,nvar)     !< Explicit diffusion predictor
real(stm_real), intent (out) :: conc(ncell,nvar)     !< Concentration at new time
real(stm_real), intent (in) :: conc_prev(ncell,nvar) !< Concentration at old time
real(STM_REAL), intent(out) :: mass(ncell,nvar)     !< mass at new time
real(STM_REAL), intent(in)  :: mass_prev(ncell,nvar) !< mass at old time
real(stm_real), intent (in) :: area (ncell,nvar)     !< Cell-centered area at new time
real(stm_real), intent (in) :: area_lo (ncell,nvar)  !< Low side area at old time
real(stm_real), intent (in) :: area_hi (ncell,nvar)  !< High side area  
real(stm_real), intent (in) :: ks_lo (ncell,nvar)    !< Low side sediment dispersion coef. at old time
real(stm_real), intent (in) :: ks_hi (ncell,nvar)    !< High side sediment dispersion coef. at old time
real(stm_real), intent (in) :: time                  !< Current time
real(stm_real), intent (in) :: end_flx_bc_prev         !< flux at right hand side old time
real(stm_real), intent (in) :: strt_flx_bc_prev         !< flux at left hand side old time  
real(stm_real), intent (in) :: dt                    !< Time step   
real(stm_real), intent (in) :: dx                    !< Spacial step  

!--- locals
integer :: ivar
integer :: icell
real(stm_real) :: dtbydx2


dtbydx2 = dt/dx/dx

do ivar = 1, nvar 
    do icell = 2,ncell-1
    
        if ( ks_lo(icell,ivar)*dtbydx2> half .or. ks_hi(icell,ivar)*dtbydx2 > half) then 
            print *, 'Unsatabality in solution!'
            pause
        end if 
        
        explicit_diffusion_term(icell,ivar) = dtbydx2 * (area_hi(icell,ivar)*ks_hi(icell,ivar)* conc_prev(icell+1,ivar) &
                                              -  area_hi(icell,ivar)*ks_hi(icell,ivar)*conc_prev(icell,ivar) &
                                              -  area_lo(icell,ivar)*ks_lo(icell,ivar)*conc_prev(icell,ivar)   &
                                              +  area_lo(icell,ivar)*ks_lo(icell,ivar)*conc_prev(icell-1,ivar) )
                
        mass(icell,ivar) = mass_prev(icell,ivar) + explicit_diffusion_term(icell,ivar)
        
        !call cons2prim(conc,mass,area,ncell,nvar) !!! nloc???
           
    end do


  ! ---- rhs boundary
        explicit_diffusion_term(1,ivar) = dtbydx2 *(area_hi(1,ivar)*ks_hi(1,ivar)* conc_prev(2,ivar) &
                                            -  area_hi(1,ivar)*ks_hi(1,ivar)*conc_prev(1,ivar) &
                                            -  area_lo(1,ivar)*ks_lo(1,ivar)*conc_prev(1,ivar) &
                                            +  area_lo(1,ivar)*ks_lo(1,ivar)*(conc_prev(2,ivar)- two * dx * strt_flx_bc_prev ) )
        
        mass(1,ivar) = mass_prev(1,ivar) + explicit_diffusion_term(1,ivar)
        
        !call cons2prim(conc,mass,area,ncell,ivar) !!!! nloc
        
  ! --- lhs boundary
        explicit_diffusion_term(ncell,ivar) =  dtbydx2 * (area_hi(ncell,ivar)*ks_hi(ncell,ivar)* (conc_prev(ncell-1,ivar)+ two*dx*end_flx_bc_prev) &
                                                -  area_hi(ncell,ivar)*ks_hi(ncell,ivar)*conc_prev(ncell,ivar)  &
                                                -  area_lo(ncell,ivar)*ks_lo(ncell,ivar)*conc_prev(ncell,ivar) &
                                                +  area_lo(ncell,ivar)*ks_lo(ncell,ivar)*conc_prev(icell-1,ivar) )
               
        mass(ncell,ivar) = mass_prev(ncell,ivar) + explicit_diffusion_term(ncell,ivar)
        


    end do                          
                                  
         call cons2prim(conc,mass,area,ncell,nvar) !!!! nloc                                 
return
end subroutine explicit_diffusion_operator 