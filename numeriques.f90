module lesconstantes_numeriques
implicit none 
integer , parameter ::  dp = selected_real_kind(15, 307) 
real(dp) ,parameter ::  pi =acos(-1.0_dp)
real(dp), parameter ::fmax = 1e9, EPSILON_0 = 8.854187818d-12, MU_0 = 12.5663706144d-7
real(dp), parameter :: c = 1.0d0 / (sqrt(EPSILON_0 * MU_0))
integer,  parameter :: i1 = 300, j1 =0 

Integer, parameter       :: Npt=5
REAL,dimension(2,Npt),parameter :: Pobs = reshape( (/2.5, 2.5,&
                                                     1.0, 2.5,&
                                                     2.8, 2.5,&
                                                     3.0, 2.5,&
                                                     3.5, 2.5/), (/2,Npt/) )
REAL(dp),dimension(Npt)  :: EPtobs
end module lesconstantes_numeriques