Program FDTD_sous_maillage
use lesconstantes_numeriques
use  lexcitation_donde
use lesfonction
use lesreference
implicit none

type(tableau) :: sm
integer :: n, Nt, m, i, j,Ntr
integer :: Nx, Ny, Nx_sm, Ny_sm, Nx_r, Ny_r
integer, parameter :: r = 3
real(dp), allocatable :: Esrc(:),Esrc_r(:)
real(dp) ::  dt, dx_r, dy_r, dt_r
real(dp) :: Lx, Ly, Lx_r, Ly_r
!integer :: nr

 real(dp), parameter :: dx = 0.01  !3e8 / (fmax * 30.0d0)
 real(dp), parameter :: dy = 0.01  ! 3e8 / (fmax * 30.0d0)


Nx = 500
Ny = 500
Nt = 1200
Ntr = 3600
Lx = Nx * dx
Ly = Ny * dy
print *, "Lx = ", Lx, "Ly = ", Ly
dt = 0.98d0 / (c * sqrt(1.0d0/(dx*dx) + 1.0d0/(dy*dy)))  ! CFL
!tmax = dt * Nt 

Nx_r = Nx*r
Ny_r = Ny*r

dx_r = dx/r  !Lx / (Nx_r-1)
dy_r = dy/r   !Ly / (Ny_r-1) dt, dx,

Lx_r = Nx_r * dx_r
Ly_r = Ny_r * dy_r


!print *, "Lx_r = ", Lx_r, "Ly_r = ", Ly_r


dt_r  = dt/r
Nx_sm  = (Nx - i1)*r
Ny_sm  = (Ny -j1)*r


! Excitation
call excitation_donde(Esrc, Nt, dt)
call excitation_donde(Esrc_r, Nt*3, dt_r)

! Initialisation des champs
call sm%initialiser_champs(Nx, Ny, Nx_sm, Ny_sm, Nx_r, Ny_r, dt)

call sm%mise_a_jour_champs(Nx, Ny, Ntr,  dx, dt, dy, Nx_sm, Ny_sm, Esrc_r)

call reference1(sm, Nx, Ny, Nt, dx, dt, dy, Esrc)
call reference2(sm, Nx_r, Ny_r, Nt*3, dx_r, dy_r, dt_r, Esrc_r)

print *, "Nombre de cartes générées :", m

! Libération mémoire
deallocate(Esrc)
deallocate(sm%Hx, sm%Ez, sm%Hy, sm%c_E, sm%c_H, sm%Ez_r1, sm%Hx_r1, sm%Hy_r1)
deallocate(sm%hx_s, sm%ez_s, sm%hy_s)

end program
