Program FDTD_sous_maillage
use lesconstantes_numeriques
use  lexcitation_donde
use lesfonction
implicit none
! Déclaration des variables
type(tableau) :: sm
integer :: Nx, Ny, Nt, Nx_sm, Ny_sm
!integer, parameter :: r = 3  ! facteur de raffinement du sous-maillage
real(dp),dimension(:), allocatable :: Esrc
real(dp), parameter :: dx = c / (fmax * 30.0d0)
real(dp), parameter :: dy = c / (fmax * 30.0d0)
real(dp), parameter :: dt = 0.98d0 / (c* sqrt(1/(dx*dx) + 1/(dy*dy)))


Nx = 500
Ny = 500
NT = 1200
 
Nx_sm = (Nx - i1) 
Ny_sm = (Ny - j1) 

print*, Nx_sm,Ny_sm


! Calcul de l'excitation
call excitation_donde(Esrc, Nt, dt)
    ! Initialiser les champs
call sm%initialiser_champs(Nx, Ny, Nx_sm, Ny_sm, dt)
call sm%mise_a_jour_champs(Nx, Ny, Nt, dx, dt, dy, Nx_sm, Ny_sm, Esrc)


    ! Libération de la mémoire
    deallocate(Esrc)
    deallocate(sm%Hx, sm%Ez, sm%Hy, sm%c_E, sm%c_H)
    deallocate(sm%hx_s, sm%ez_s, sm%hy_s)

end program