module lesfonction
use lesconstantes_numeriques
use lexcitation_donde
 implicit none
    type :: tableau
        real(dp), dimension(:,:), allocatable :: Hx, Ez, Hy
        real(dp), dimension(:,:), allocatable :: hx_s, ez_s, hy_s
        real(dp), dimension(:,:), allocatable :: c_E, c_H
    contains
        procedure :: initialiser_champs
        procedure :: mise_a_jour_champs
    end type tableau

contains

    ! Sous-routine pour initialiser les champs
    subroutine initialiser_champs(this, Nx, Ny, Nx_sm, Ny_sm, dt)
        use lesconstantes_numeriques
        implicit none
        class(tableau), intent(inout) :: this
        integer, intent(in) :: Nx, Ny, Nx_sm, Ny_sm
        real(dp), intent(in) :: dt

        allocate(this%hx_s(0:Nx_sm, 0:Ny_sm), this%ez_s(0:Nx_sm, 0:Ny_sm), this%hy_s(0:Nx_sm, 0:Ny_sm))
        this%hx_s = 0.0_dp
        this%ez_s = 0.0_dp
        this%hy_s = 0.0_dp

       allocate( this%Ez(0:Nx, 0:Ny), this%Hx(0:Nx, 0:Ny), this%Hy(0:Nx, 0:Ny))
        this%Ez = 0.0_dp
        this%Hy = 0.0_dp
        this%Hx = 0.0_dp

        !coefficients pour les champs
allocate(this%c_H(0:Nx, 0:Ny),  this%c_E(0:Nx, 0:Ny))
        this%c_E = dt / EPSILON_0
        this%c_H = dt / MU_0
        
    end subroutine initialiser_champs


subroutine couplage_CG_vers_FG(fd, i1, j1, r)
  !  use lesconstantes_numeriques
    implicit none
    class(tableau), intent(inout) :: fd
    integer, intent(in)           :: i1, j1, r
    integer :: j, i, i_cg, j_cg

    ! Injecte Hy CG → hy_s (bord gauche FG)
    do j = 1, size(fd%hy_s, 2)-2
        j_cg = j1 + j / r
        fd%hy_s(0,j) = 0.5_dp * (fd%Hy(i1-1,j_cg) + fd%Hy(i1,j_cg))
    end do

    ! Injecte Hx CG → hx_s (bord bas FG)
    do i = 1, size(fd%hx_s, 1)-2
        i_cg = i1 + i / r
        fd%hx_s(i,0) = 0.5_dp * (fd%Hx(i_cg,j1-1) + fd%Hx(i_cg,j1))
    end do
end subroutine couplage_CG_vers_FG


    ! Subroutine pour la mise à jour des champs
    subroutine mise_a_jour_champs(fd, Nx, Ny, Nt, dx, dt, dy, Nx_sm, Ny_sm, Esrc)
        use lesconstantes_numeriques
        implicit none
    class(tableau), intent(inout)                              :: fd
    integer, intent(in)                                        :: Nx, Ny, Nt, Nx_sm, Ny_sm
    real(dp), intent(in)                                       :: dx, dy, dt
    real(dp), intent(in)                                       :: Esrc(0:Nt-1)
    integer, parameter                                         :: r = 3
    integer                                                    :: i, j, n,  m
    real(dp)                                                   :: Hz2, Hz3, dx_sm, dy_sm, dt_prime
    real(dp)                                                   :: Haux_y, dx_prime, dy_prime, f_Hy, eaux_z, fez
    !real(dp)                                                   :: LCC, LIC1, LIC2, TCC1, TCC2, alpha 
    real(dp), dimension(0 : 1       , 0 : Ny)                  :: Ezx0_n1 , Ezx0_n2                    
    real(dp), dimension(Nx - 1 : Nx , 0 : Ny)                  :: Ezx_n1 , Ezx_n2                 
    real(dp), dimension(0 : 1       , 0 : Nx)                  :: Ezy0_n1 , Ezy0_n2                     
    real(dp), dimension(Ny - 1 : Ny , 0 : Nx)                  :: Ezy_n1 , Ezy_n2
   ! real(dp), dimension(0 : 1       , 0 : Ny)                  :: Hx0_n1, Hx0_n2                    
   ! real(dp), dimension(Nx - 1 : Nx , 0 : Ny)                  :: Hx_n1, Hx_n2
  !  real(dp), dimension(0 : 1       , 0 : Nx)                  :: Hy0_n1, Hy0_n2                    
  !  real(dp), dimension(Ny - 1 : Ny , 0 : Nx)                  :: Hy_n1, Hy_n2
    real(dp)                                                   :: coef_mur1, coef_mur2, coef_mur3
   real(dp)                                                   :: energy,energy_CG, energy_FG
   real(dp), dimension(0:Nt-1)                                :: total_energy


        dx_sm        = dx/ r
        dy_sm        = dy/ r
        dx_prime     = 0.5*dx + 0.5*dx_sm
        dy_prime     = 0.5*dy + 0.5*dy_sm
        f_Hy         = dy/dy_sm
        fez          = 1.0
        dt_prime     = dt/r

!ouverture des fichiers pour stockes les resultats
open(unit=10, file="Ez_t.dat",        status = 'replace', action = 'write')
open(unit=11, file="Hx_t.dat",        status = 'replace', action = 'write')
open(unit=12, file="Hy_t.dat",        status = 'replace', action = 'write')
open(unit=13, file="carto_t.dat",     status = 'replace', action = 'write')
 open(unit=14, file="carto_complet.dat", position='append')
open(unit=15, file="ez_s.dat",        status = 'replace', action = 'write')
open(unit=16, file="ez_obs.dat",   status = 'replace', action = 'write')
open(unit=17, file="Ez_obs.dat",      status = 'replace', action = 'write')
 m = 0
Ezx0_n1 = 0.d0 ;    Ezx0_n2 = 0.0d0 
Ezx_n1  = 0.d0 ;    Ezx_n2  = 0.0d0
Ezy0_n1 = 0.d0 ;    Ezy0_n2 = 0.0d0
Ezy_n1  = 0.d0 ;    Ezy_n2  = 0.0d0
!Hx0_n1  = 0.d0 ;    Hx0_n2  = 0.0d0
!Hx_n1   = 0.d0 ;    Hx_n2   = 0.0d0
!Hy0_n1  = 0.d0 ;    Hy0_n1   = 0.0d0
!Hy_n1   = 0.d0 ;    Hy_n2   = 0.0d0

coef_mur1 = (c * dt - dx) / (c * dt + dx)
coef_mur2 = (2.d0 * dx) / (c * dt + dx)
coef_mur3 = (c * dt)**2 / ( 2 * dx * ( c * dt + dx ) )

do n = 0, Nt - 1
if (mod(n, 100) == 0) print *, "Itération temporelle n°", n
! Condition de bord x = 0
Ezx0_n2        = Ezx0_n1                        ! temps n - 2 pas i = 0 et i = 1
Ezx0_n1(0 , :) = fd%Ez(0,:)                     ! temps n - 1 pas i = 0
Ezx0_n1(1 , :) = fd%Ez(1,:)                     ! temps n - 1 pas i = 1

            ! Condition de bord x = Nx
Ezx_n2             = Ezx_n1                     ! temps n - 2 pas i = Nx - 1 et i = Nx
Ezx_n1(Nx - 1 , :) = fd%Ez(Nx - 1 , :)          ! temps n - 1 pas i = Nx - 1
Ezx_n1(Nx , :)     = fd%Ez(Nx , :)              ! temps n - 1 pas i = Nx

            ! Condition de bord y = 0
Ezy0_n2       = Ezy0_n1                         ! temps n - 2 pas j = 0 et j = 1
Ezy0_n1(0, :) = fd%Ez(:, 0)                     ! temps n - 1 pas j = 0
Ezy0_n1(1, :) = fd%Ez(:, 1)                     ! temps n - 1 pas j = 1

            ! Condition de bord y = Ny
Ezy_n2            = Ezy_n1                      ! temps n - 2 pas j = Ny - 1 et j = Ny
Ezy_n1(Ny - 1, :) = fd%Ez(:, Ny - 1)            ! temps n - 1 pas j = Ny - 1
Ezy_n1(Ny, :)     = fd%Ez(:, Ny)                ! temps n - 1 pas j = Ny
            
do i = 1, Nx - 1
 if (i < i1 ) then
do j = 1, Ny - 1
fd%Ez(i,j) = fd%Ez(i,j) + fd%c_E(i,j) * ((fd%Hy(i,j) - fd%Hy(i-1,j))/ dx  &
                        - (fd%Hx(i,j) - fd%Hx(i,j-1)) / dy)
!end if
end do
end if
end do
 
            !! Condition aux bord (Mur absorbant) x = 0
fd%Ez(0 , :)  =  - Ezx0_n2(1, :)                                             &
                 + coef_mur1 * (fd%Ez(1, :) + Ezx0_n2(0, :))                 &         ! n2 -> le temps n - 1    n1 -> n
                  + coef_mur2 * (Ezx0_n1(0, :) + Ezx0_n1(1, :))

            ! Condition aux bord (Mur absorbant) x = Nx
fd%Ez(Nx , :) =    - Ezx_n2(Nx - 1, :)                                         &
                  + coef_mur1 * (fd%Ez(Nx - 1, :) + Ezx_n2(Nx, :))             &
                  + coef_mur2 * (Ezx_n1(Nx, :) + Ezx_n1(Nx - 1, :))

            ! Condition aux bord (Mur absorbant) y = 0
fd%Ez(0:i1 , 0)  =   - Ezy0_n2(1, :)                                             &
                + coef_mur1 * (fd%Ez(:, 1) + Ezy0_n2(0, :))                      &
                + coef_mur2 * (Ezy0_n1(0, :) + Ezy0_n1(1, :))

!fd%Ez(i1:Nx, 0) = 0.0_dp 

            ! Condition aux bord (Mur absorbant) y = Ny
fd%Ez(0:i1 , Ny) =   - Ezy_n2(Ny - 1, :)                                         &
                  + coef_mur1 * (fd%Ez(:, Ny - 1) + Ezy_n2(Ny, :))               &
                  + coef_mur2 * (Ezy_n1(Ny, :) + Ezy_n1(Ny - 1,:))

!fd%Ez(i1:Nx, Ny) = 0.0_dp ! Condition aux bord (Mur absorbant) pour le coin

             ! Source dans domaine grossier
fd%Ez(350,350) = fd%Ez(350, 350) + Esrc(n)

            ! Mise à jour Hx dans le domaine grossier       
do i = 0, Nx - 1
 if (i <= i1 ) then
do j = 0, Ny - 1
!if (j < j1 ) then
     fd%Hx(i,j) = fd%Hx(i,j) - fd%c_H(i,j) / dy * (fd%Ez(i,j+1) - fd%Ez(i,j))
    ! end if 
end do
end if
end do

            ! Mise à jour  Hy dans le domaine grossier
do i = 0, Nx - 1
 if (i < i1 ) then
do j = 0, Ny - 1
fd%Hy(i,j) = fd%Hy(i,j) + fd%c_H(i,j) / dx * (fd%Ez(i+1,j) - fd%Ez(i,j))
 
 end do
    end if
end do


! Ez raffiné
do i = 1, Nx_sm - 1
do j = 1, Ny_sm- 1   
    fd%ez_s(i,j) = fd%ez_s(i,j) + (dt_prime/ epsilon_0 ) * ((fd%hy_s(i,j) - &
    fd%hy_s(i-1,j))/dx_sm - (fd%hx_s(i,j) - fd%hx_s(i,j-1))/dy_sm)
   end do
  end do
    ! Source dans le sous-maillage
!if (mod(n, r) == 0) then
   fd%ez_s(300,750) = fd%ez_s(300,750) + Esrc(n)
!end if


 !le champs magnetique(Hx) à l'intérieur du sous-maillage
do i = 0, Nx_sm - 1
do j = 0, Ny_sm - 1
    fd%hx_s(i,j) = fd%hx_s(i,j) - dt_prime/ (mu_0 * dy_sm) * (fd%ez_s(i,j+1) - fd%ez_s(i,j))
end do
end do
 
 ! le champs magnetique(Hy) à l'intérieur du sous-maillage
 do i = 1, Nx_sm - 1
 do j = 0, Ny_sm - 1
 fd%hy_s(i,j) = fd%hy_s(i,j) + dt_prime/ (mu_0 * dx_sm) * (fd%ez_s(i+1,j) - fd%ez_s(i,j))
end do
end do


! Interpolation depuis la grille fine 
do j = 0, Ny_sm - 1
     eaux_z = compute_ez_aux(fd, 0*r+1 , j)
     end do 

do i = 1, Nx - 1
   do j = 1, Ny - 1
   fd%Hy(i1, j) = fd%Hy(i1, j) + dt / (mu_0 * dx) * (fez*eaux_z- fd%Ez(i, j))
   end do 
end do


! le champs electrique dans l'interface
do j = 0, Ny - 1
Hz2 = 0.0_dp
Hz3 = 0.0_dp
Hz2 = fd%Hy(i1-1, j)  
Hz3 = fd%Hy(i1+1, j)
Haux_y = (2.0*Hz2 /9.0)+ (1.0*Hz3 /9.0)
end do


do i =  1, Nx_sm - 1
do j =  1, Ny_sm - 1
fd%ez_s(0,j) = fd%ez_s(0,j) + (dt_prime/ epsilon_0 ) *((fd%hy_s(0,j) - (f_Hy *Haux_y))/dx_prime - &
              (fd%hx_s(i,j) - fd%hx_s(i,j-1))/dy_sm)
end do
end do 


    write(10,*) n*dt,       fd%Ez(350,250), fd%Ez(250,250), fd%Ez(150,250), fd%Ez(450,150)
    write(11,*) n*dt,       fd%Hx(250,250), fd%Hx(150,100), fd%Hx(50,100), fd%Hx(100,150)
    write(12,*) n*dt,       fd%Hy(250,250), fd%Hy(250,250), fd%Hy(50,100), fd%Hy(100,150)
    write(15,*) n*dt_prime, fd%ez_s(300,750), fd%ez_s(150,250), fd%ez_s(50,600)

   if (mod(n, 100) == 0 .and. n > 0) then
    m = m + 1

    !write(14,*) ! Ligne vide pour séparer les temps
    ! on n’écrit que la valeur Ez(i,j), séparée par un espace
    do i = 1, Nx, 2
       write(14,*) (fd%Ez(i,j), j= 0, Ny, 2)
    end do
     write(14,*)

     !valeur ez_s(i,j), séparée par un espace
     do i = 1, Nx_sm, 2
         write(13,*) (fd%ez_s(i,j), j= 1, Ny_sm, 2)
    end do
 write(13,*)
  end if   

!pour calculer l'erreur entre les deux champs
 write(16,*) n*dt_prime, fd%ez_s(0,750)
 write(17,*) n*dt_prime, fd%Ez(300,250)
 end do
        print *, "Nombre de carto en temps = ", m
        print *, "Fin de la simulation"
        close(10); close(11); close(12); close(13); close(14); close(15)
        close(16); close(17)
   
 end subroutine mise_a_jour_champs
  
 ! Subroutine pour la mise à jour des champs
    subroutine reference1(fd, Nx, Ny, Nt, dx, dt, dy, Esrc)
        use lesconstantes_numeriques
        implicit none
    class(tableau), intent(inout)                              :: fd
    integer, intent(in)                                        :: Nx, Ny, Nt
    real(dp), intent(in)                                       :: dx, dy, dt
    real(dp), intent(in)                                       :: Esrc(0:Nt-1)
    integer, parameter                                         :: r = 3
    integer                                                    :: i, j, n,  m 
    real(dp), dimension(0 : 1       , 0 : Ny)                  :: Ezx0_n1 , Ezx0_n2                    
    real(dp), dimension(Nx - 1 : Nx , 0 : Ny)                  :: Ezx_n1 , Ezx_n2                 
    real(dp), dimension(0 : 1       , 0 : Nx)                  :: Ezy0_n1 , Ezy0_n2                     
    real(dp), dimension(Ny - 1 : Ny , 0 : Nx)                  :: Ezy_n1 , Ezy_n2
    real(dp)                                                   :: coef_mur1, coef_mur2, coef_mur3
   

        dx_sm        = dx/ r
        dy_sm        = dy/ r
        dx_prime     = 0.5*dx + 0.5*dx_sm
        dy_prime     = 0.5*dy + 0.5*dy_sm
        f_Hy         = dy/dy_sm
        fez          = 1.0
        dt_prime     = dt/r

!ouverture des fichiers pour stockes les resultats
open(unit=10, file="Ez_t1.dat",        status = 'replace', action = 'write')
open(unit=11, file="Hx_t1.dat",        status = 'replace', action = 'write')
open(unit=12, file="Hy_t1.dat",        status = 'replace', action = 'write')
open(unit=13, file="carto_t1.dat",     status = 'replace', action = 'write')

 m = 0
Ezx0_n1 = 0.d0 ;    Ezx0_n2 = 0.0d0 
Ezx_n1  = 0.d0 ;    Ezx_n2  = 0.0d0
Ezy0_n1 = 0.d0 ;    Ezy0_n2 = 0.0d0
Ezy_n1  = 0.d0 ;    Ezy_n2  = 0.0d0


coef_mur1 = (c * dt - dx) / (c * dt + dx)
coef_mur2 = (2.d0 * dx) / (c * dt + dx)
coef_mur3 = (c * dt)**2 / ( 2 * dx * ( c * dt + dx ) )

do n = 0, Nt - 1
if (mod(n, 100) == 0) print *, "Itération temporelle n°", n
! Condition de bord x = 0
Ezx0_n2        = Ezx0_n1                        ! temps n - 2 pas i = 0 et i = 1
Ezx0_n1(0 , :) = fd%Ez(0,:)                     ! temps n - 1 pas i = 0
Ezx0_n1(1 , :) = fd%Ez(1,:)                     ! temps n - 1 pas i = 1

            ! Condition de bord x = Nx
Ezx_n2             = Ezx_n1                     ! temps n - 2 pas i = Nx - 1 et i = Nx
Ezx_n1(Nx - 1 , :) = fd%Ez(Nx - 1 , :)          ! temps n - 1 pas i = Nx - 1
Ezx_n1(Nx , :)     = fd%Ez(Nx , :)              ! temps n - 1 pas i = Nx

            ! Condition de bord y = 0
Ezy0_n2       = Ezy0_n1                         ! temps n - 2 pas j = 0 et j = 1
Ezy0_n1(0, :) = fd%Ez(:, 0)                     ! temps n - 1 pas j = 0
Ezy0_n1(1, :) = fd%Ez(:, 1)                     ! temps n - 1 pas j = 1

            ! Condition de bord y = Ny
Ezy_n2            = Ezy_n1                      ! temps n - 2 pas j = Ny - 1 et j = Ny
Ezy_n1(Ny - 1, :) = fd%Ez(:, Ny - 1)            ! temps n - 1 pas j = Ny - 1
Ezy_n1(Ny, :)     = fd%Ez(:, Ny)                ! temps n - 1 pas j = Ny
            
do i = 1, Nx - 1
do j = 1, Ny - 1
fd%Ez(i,j) = fd%Ez(i,j) + fd%c_E(i,j) * ((fd%Hy(i,j) - fd%Hy(i-1,j))/ dx  &
                        - (fd%Hx(i,j) - fd%Hx(i,j-1)) / dy)

end do
end do
 
            !! Condition aux bord (Mur absorbant) x = 0
fd%Ez(0 , :)  =  - Ezx0_n2(1, :)                                             &
                 + coef_mur1 * (fd%Ez(1, :) + Ezx0_n2(0, :))                 &         ! n2 -> le temps n - 1    n1 -> n
                  + coef_mur2 * (Ezx0_n1(0, :) + Ezx0_n1(1, :))

            ! Condition aux bord (Mur absorbant) x = Nx
fd%Ez(Nx , :) =    - Ezx_n2(Nx - 1, :)                                         &
                  + coef_mur1 * (fd%Ez(Nx - 1, :) + Ezx_n2(Nx, :))            &
                  + coef_mur2 * (Ezx_n1(Nx, :) + Ezx_n1(Nx - 1, :))

            ! Condition aux bord (Mur absorbant) y = 0
fd%Ez(: , 0)  =   - Ezy0_n2(1, :)                                             &
                + coef_mur1 * (fd%Ez(:, 1) + Ezy0_n2(0, :))                      &
                + coef_mur2 * (Ezy0_n1(0, :) + Ezy0_n1(1, :))

!fd%Ez(i1:Nx, 0) = 0.0_dp 

            ! Condition aux bord (Mur absorbant) y = Ny
fd%Ez(: , Ny) =   - Ezy_n2(Ny - 1, :)                                         &
                  + coef_mur1 * (fd%Ez(:, Ny - 1) + Ezy_n2(Ny, :))               &
                  + coef_mur2 * (Ezy_n1(Ny, :) + Ezy_n1(Ny - 1,:))



             ! Source dans domaine grossier
fd%Ez(40,400) = fd%Ez(400, 400) + Esrc(n)

            ! Mise à jour Hx dans le domaine grossier       
do i = 1, Nx - 1
do j = 1, Ny - 1
     fd%Hx(i,j) = fd%Hx(i,j) - fd%c_H(i,j) / dy * (fd%Ez(i,j+1) - fd%Ez(i,j))
  
end do
end do

            ! Mise à jour  Hy dans le domaine grossier
do i = 1, Nx - 1
do j = 1, Ny - 1
fd%Hy(i,j) = fd%Hy(i,j) + fd%c_H(i,j) / dx * (fd%Ez(i+1,j) - fd%Ez(i,j))
 end do
end do



    write(10,*) n*dt,       fd%Ez(350,250), fd%Ez(250,250), fd%Ez(150,250), fd%Ez(250,150)
    write(11,*) n*dt,       fd%Hx(100,100), fd%Hx(150,100), fd%Hx(50,100), fd%Hx(100,150)
    write(12,*) n*dt,       fd%Hy(100,100), fd%Hy(150,100), fd%Hy(50,100), fd%Hy(100,150)


   if (mod(n, 100) == 0 .and. n > 0) then
    m = m + 1

    write(14,*) ! Ligne vide pour séparer les temps
    ! on n’écrit que la valeur Ez(i,j), séparée par un espace
    do i = 1, Nx, 2
       write(14,*) (fd%Ez(i,j), j= 0, Ny, 2)
    end do
     write(14,*)

     !valeur ez_s(i,j), séparée par un espace
     do i = 1, Nx_sm, 2
         write(13,*) (fd%ez_s(i,j), j= 1, Ny_sm, 2)
    end do
 write(13,*)
  end if   

 end do
        print *, "Nombre de carto en temps = ", m
        print *, "Fin de la simulation"
        close(10); close(11); close(12); close(13); close(14); close(15)
        close(16); close(17)
   
 end subroutine mise_a_jour_champs


  ! Interpolation conservatrice pour Ez aux interfaces FG → CG
function compute_ez_aux(fd, i_f, j_f) result(ez_aux)
    use lesconstantes_numeriques
    implicit none
    type(tableau), intent(in) :: fd
    integer, intent(in)       :: i_f, j_f  ! Indices dans la grille fine (FG)
    real(dp)                  :: ez_aux
    real(dp) :: sum
    integer :: ii, jj
   real(dp), dimension(-1:1, -1:1) :: w

    ! Pondération conservatrice 3x3 (somme = 16)
    w = reshape([ &
        1.0/3.0_dp, 2.0/3.0_dp, 1.0_dp, 2.0/3.0_dp, 1.0/3.0_dp, &
        2.0/3.0_dp, 1.0/3.0_dp, 1.0_dp, 1.0/3.0_dp ], [3,3])
    ! Moyenne pondérée 3x3 centrée sur (i_f, j_f)
    sum = 0.0_dp
    do jj = -1, 1
        do ii = -1, 1
            sum = sum + w(ii, jj) * fd%ez_s(i_f, j_f+jj)
        end do
    end do

    ez_aux = sum
    ez_aux = sum / 3.0_dp
end function compute_ez_aux

end module lesfonction  