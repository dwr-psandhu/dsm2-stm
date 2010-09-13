module hydro_uniform_flow
use stm_precision

real(stm_real) :: const_flow = zero 
real(stm_real) :: const_area = zero
real(stm_real) :: reversal_time = LARGEREAL

contains 


!> Set constant flow and area that will be used in the no_flow hydro interface 
subroutine set_uniform_flow_area(flow, area, reverse_time)
use stm_precision
implicit none
real(stm_real), intent(in) :: flow     !< constant flow to be set
real(stm_real), intent(in) :: area     !< constant area to be set
real(stm_real), intent(in), optional :: reverse_time !< time at which flow will reverse direction
const_flow = flow
const_area = area
if (present(reverse_time)) then
  reversal_time = reverse_time
else
  reversal_time = LARGEREAL
end if
return 
end subroutine


!>Simple hydrodynamic interface for constant area and constant flow
subroutine uniform_flow_area(flow,    &
                             flow_lo, &
                             flow_hi, &
                             area,    &
                             area_lo, &
                             area_hi, &
                             ncell,   &
                             time,    &
                             dx,      &                        
                             dt)
use stm_precision                   
implicit none
integer, intent(in) :: ncell                   !< number of cells
real(stm_real), intent(in) :: time             !< time of request
real(stm_real), intent(in) :: dx               !< spatial step
real(stm_real), intent(in) :: dt               !< time step 
real(stm_real), intent(out) :: flow(ncell)     !< cell and time centered flow
real(stm_real), intent(out) :: flow_lo(ncell)  !< lo face flow, time centered
real(stm_real), intent(out) :: flow_hi(ncell)  !< hi face flow, time centered
real(stm_real), intent(out) :: area(ncell)     !< cell center area, old time
real(stm_real), intent(out) :: area_lo(ncell)  !< area lo face, time centered
real(stm_real), intent(out) :: area_hi(ncell)  !< area hi face, time centered

if (time <= reversal_time) then
   flow = const_flow
else
   flow = minus * const_flow
end if
    
flow_hi = flow
flow_lo = flow
area = const_area
area_lo = const_area
area_hi = const_area

return
end subroutine


end module


