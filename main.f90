program FDTD_sous_maillage
use lesconstantes_numeriques
use  lexcitation_donde
use lesfonction
use lesreference
  use cpml_2d
!use  dispersion_tmz2d_cg_fg
implicit none

type(tableau) :: sm
integer :: n, Nt, m, i, j, Ntr
integer :: Nx, Ny, Nx_sm, Ny_sm, Nx_r, Ny_r
integer, parameter :: r = 3
real(dp), allocatable :: Esrc(:),Esrc_r(:)
real(dp) ::  dt, dx_r, dy_r, dt_r
real(dp) :: Lx, Ly, Lx_r, Ly_r
real(dp), parameter :: dx = 0.01, dy = 0.01
real :: t1, t2, t_init, t_update, t_ref, t_total   ! pour cpu_time

! -------------------------
! Paramètres numériques
! -------------------------
Nx = 500
Ny = 500
Nt = 1200
Ntr = Nt*r
Lx = Nx * dx
Ly = Ny * dy
print *, "Lx = ", Lx, "Ly = ", Ly
dt = 0.98d0 / (c * sqrt(1.0d0/(dx*dx) + 1.0d0/(dy*dy)))  ! CFL

Nx_r = Nx*r
Ny_r = Ny*r
print *, "Nx_r = ", Nx_r, "Ny_r = ", Ny_r
dx_r = dx/r  
dy_r = dy/r   

Lx_r = Nx_r * dx_r
Ly_r = Ny_r * dy_r



dt_r  = dt/r
Nx_sm  = (Nx - i1)*r
Ny_sm  = (Ny - j1)*r

allocate(Esrc_r(0:Ntr-1))

! -------------------------
! Excitation
! -------------------------
call excitation_donde(Esrc, Nt, dt)
call excitation_donde(Esrc_r, Ntr, dt_r)

! -------------------------
! Initialisation
! -------------------------

call sm%initialiser_champs(Nx, Ny, Nx_sm, Ny_sm, Nx_r, Ny_r, dt)

! -------------------------
! Mise à jour des champs
! -------------------------

call sm%mise_a_jour_champs(Nx, Ny, Ntr, dx, dt, dy, Nx_sm, Ny_sm, Esrc_r)

! -------------------------
! Références
! -------------------------

call reference1(sm, Nx, Ny, Nt, dx, dt, dy, Esrc)
call reference2(sm, Nx_r, Ny_r, Ntr, dx_r, dt_r, dy_r, Esrc_r)
! -------------------------
! Dispersion
! -------------------------
!call dispersion_numerique(dx, dy, r, dx_r, dy_r, dt, dt_r)


! -------------------------
! Libération mémoire
! -------------------------
deallocate(Esrc)
deallocate(sm%Hx, sm%Ez, sm%Hy, sm%c_E, sm%c_H, sm%Ez_r1, sm%Hx_r1, sm%Hy_r1)
deallocate(sm%hx_s, sm%ez_s, sm%hy_s)

end program FDTD_sous_maillage
