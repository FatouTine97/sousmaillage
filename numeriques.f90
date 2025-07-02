module lesconstantes_numeriques
implicit none 
integer , parameter ::  dp = selected_real_kind(15, 307) 
real(dp) ,parameter ::  pi =acos(-1.0_dp)
 real(dp), parameter :: EPSILON_0 = 8.854187818d-12, MU_0 = 12.5663706144d-7
 real(dp), parameter :: fmax = 1e9
 real(dp), parameter :: c = 1.0d0 / (sqrt(EPSILON_0 * MU_0))
 integer,  parameter :: i1 = 200, j1 =0 , i2 =0, j2 = 0
end module lesconstantes_numeriques