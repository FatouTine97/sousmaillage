Program FDTD_sous_maillage
use lesconstantes_numeriques
use  lexcitation_donde
use lesfonction
!use lesreference
implicit none
! Déclaration des variables
type(tableau) :: sm
integer :: n, Nt, m , i, j
integer :: Nx, Ny, Nx_sm, Ny_sm, Nx_r, Ny_r
integer, parameter :: r = 3  ! facteur de raffinement du sous-maillage
real(dp),dimension(:), allocatable :: Esrc
real(dp), parameter :: dx = c / (fmax * 30.0d0)
real(dp), parameter :: dy = c / (fmax * 30.0d0)
real(dp), parameter :: dt = 0.98d0 / (c* sqrt(1/(dx*dx) + 1/(dy*dy)))
real(dp), parameter :: dt_prime = dt /r

Nx = 500
Ny = 500
NT = 1200
 
Nx_sm = (Nx - i1) *r
Ny_sm = (Ny - j1) *r
Nx_r = Nx*r
Ny_r = Ny*r

print*, Nx_sm,Ny_sm, Nx_r, Ny_r
 

! Calcul de l'excitation
call excitation_donde(Esrc, Nt, dt)

call sm%initialiser_champs(Nx, Ny, Nx_sm, Ny_sm, Nx_r,Ny_r, dt)

  !ouverture des fichiers pour l'ecriture
open(unit=10, file="GM/Ez_t.dat",           status = 'replace', action = 'write')
open(unit=11, file="GM/Hx_t.dat",           status = 'replace', action = 'write')
open(unit=12, file="GM/Hy_t.dat",           status = 'replace', action = 'write')
open(unit=13, file="PM/carto_t.dat",        status = 'replace', action = 'write')
open(unit=14, file="GM/carto_t1.dat",       status = 'replace', action = 'write')
open(unit=15, file="PM/ez_s.dat",           status = 'replace', action = 'write')
open(unit=16, file="R1/Ez_r1.dat",          status = 'replace', action = 'write')
open(unit=17, file="R1/Hx_r1.dat",          status = 'replace', action = 'write')
open(unit=18, file="R1/Hy_r1.dat",          status = 'replace', action = 'write')
open(unit=19, file="R2/Ez_r2.dat",          status = 'replace', action = 'write')
open(unit=20, file="R2/Hx_r2.dat",          status = 'replace', action = 'write')
open(unit=21, file="R2/Hy_r2.dat",          status = 'replace', action = 'write')
m =0
do n = 0, Nt - 1
if (mod(n, 100) == 0) print *, "Itération temporelle n°", n
call sm%mise_a_jour_champs(Nx, Ny, Nt,n, dx, dt, dy, Nx_sm, Ny_sm, Esrc)
!call reference1(Nx, Ny, Nt,n, dx, dt, dy, Esrc)
!call reference2(sm, Nx_r, Ny_r,Nt, n, dx, dy, dt, Esrc)

!Ecriture des champs dans les fichiers 
write(10,*) n*dt,       sm%Ez(50,250),   sm%Ez(250,250),   sm%Ez(150,250), sm%Ez(250,150)
write(11,*) n*dt,       sm%Hx(250,250),   sm%Hx(150,100),   sm%Hx(50,100), sm%Hx(100,150)
write(12,*) n*dt,       sm%Hy(250,250),    sm%Hy(250,250),   sm%Hy(50,100), sm%Hy(100,150)
write(15,*) n*dt_prime, sm%ez_s(300,750),  sm%ez_s(450,750), sm%ez_s(500,600)
write(16,*) n*dt,       sm%Ez_r1(350,350), sm%Ez_r1(250,250), sm%Ez_r1(350,250), sm%Ez_r1(250,350)
write(17,*) n*dt,       sm%Hx_r1(100,100), sm%Hx_r1(150,100), sm%Hx_r1(50,100), sm%Hx_r1(100,150)
write(18,*) n*dt,       sm%Hy_r1(100,100), sm%Hy_r1(150,100), sm%Hy_r1(50,100), sm%Hy_r1(100,150)
write(30,*) n*dt,       sm%Ez_r2(1050,1050), sm%Ez_r2(750,750), sm%Ez_r2(1050,750), sm%Ez_r2(750,1050)
write(31,*) n*dt,       sm%Hx_r2(300,300),   sm%Hx_r2(150,300), sm%Hx_r2(50,100),    sm%Hx_r2(100,150)
 write(32,*) n*dt,      sm%Hy_r2(100,100),  sm%Hy_r2(150,100), sm%Hy_r2(50,100),  sm%Hy_r2(100,150)

   if (mod(n, 100) == 0 .and. n > 0) then
    m = m + 1

    !write(14,*) ! Ligne vide pour séparer les temps
    ! on n’écrit que la valeur Ez(i,j), séparée par un espace
    do i = 1, Nx, 2
       write(14,*) (sm%Ez(i,j), j= 0, Ny, 2)
    end do
     write(14,*)

     !valeur ez_s(i,j), séparée par un espace
     do i = 1, Nx_sm, 2
         write(13,*) (sm%ez_s(i,j), j= 1, Ny_sm, 2)
    end do
 write(13,*)
  end if   

 end do
    ! Libération de la mémoire
    deallocate(Esrc)
    deallocate(sm%Hx, sm%Ez, sm%Hy, sm%c_E, sm%c_H)
    deallocate(sm%hx_s, sm%ez_s, sm%hy_s)

end program