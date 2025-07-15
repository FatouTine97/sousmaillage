module lesreference 
use lesconstantes_numeriques
use lesfonction
implicit none
contains 

! Subroutine pour la mise à jour des champs
    subroutine reference1(fd, Nx, Ny, Nt, dx, dt, dy, Esrc)
        use lesconstantes_numeriques
        implicit none
    class(tableau), intent(inout)                              :: fd
    integer, intent(in)                                        :: Nx, Ny, Nt
    real(dp), intent(in)                                       :: dx, dy, dt
    real(dp), intent(in)                                       :: Esrc(0:Nt-1)
    integer                                                    :: i, j, n,  m 
    real(dp), dimension(0 : 1       , 0 : Ny)                  :: Ezx0_n1 , Ezx0_n2                    
    real(dp), dimension(Nx - 1 : Nx , 0 : Ny)                  :: Ezx_n1 , Ezx_n2                 
    real(dp), dimension(0 : 1       , 0 : Nx)                  :: Ezy0_n1 , Ezy0_n2                     
    real(dp), dimension(Ny - 1 : Ny , 0 : Nx)                  :: Ezy_n1 , Ezy_n2
    real(dp)                                                   :: coef_mur1, coef_mur2, coef_mur3
   

!ouverture des fichiers pour stockes les resultats
open(unit=20, file="Ez_r1.dat",        status = 'replace', action = 'write')
open(unit=21, file="Hx_r1.dat",        status = 'replace', action = 'write')
open(unit=22, file="Hy_r1.dat",        status = 'replace', action = 'write')
!open(unit=23, file="carto_t2.dat",     status = 'replace', action = 'write')
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
Ezx0_n1(0 , :) = fd%Ez_r1(0,:)                     ! temps n - 1 pas i = 0
Ezx0_n1(1 , :) = fd%Ez_r1(1,:)                     ! temps n - 1 pas i = 1

            ! Condition de bord x = Nx
Ezx_n2             = Ezx_n1                     ! temps n - 2 pas i = Nx - 1 et i = Nx
Ezx_n1(Nx - 1 , :) = fd%Ez_r1(Nx - 1 , :)          ! temps n - 1 pas i = Nx - 1
Ezx_n1(Nx , :)     = fd%Ez_r1(Nx , :)              ! temps n - 1 pas i = Nx

            ! Condition de bord y = 0
Ezy0_n2       = Ezy0_n1                         ! temps n - 2 pas j = 0 et j = 1
Ezy0_n1(0, :) = fd%Ez_r1(:, 0)                     ! temps n - 1 pas j = 0
Ezy0_n1(1, :) = fd%Ez_r1(:, 1)                     ! ; close(14)temps n - 1 pas j = 1

            ! Condition de bord y = Ny
Ezy_n2            = Ezy_n1                      ! temps n - 2 pas j = Ny - 1 et j = Ny
Ezy_n1(Ny - 1, :) = fd%Ez_r1(:, Ny - 1)            ! temps n - 1 pas j = Ny - 1
Ezy_n1(Ny, :)     = fd%Ez_r1(:, Ny)                ! temps n - 1 pas j = Ny
            
do i = 1, Nx - 1
do j = 1, Ny - 1
fd%Ez_r1(i,j) = fd%Ez_r1(i,j) + fd%c_E(i,j) * ((fd%Hy_r1(i,j) - fd%Hy_r1(i-1,j))/ dx  &
                        - (fd%Hx_r1(i,j) - fd%Hx_r1(i,j-1)) / dy)

end do
end do
 
            !! Condition aux bord (Mur absorbant) x = 0
fd%Ez_r1(0 , :)  =  - Ezx0_n2(1, :)                                             &
                 + coef_mur1 * (fd%Ez_r1(1, :) + Ezx0_n2(0, :))                 &         ! n2 -> le temps n - 1    n1 -> n
                  + coef_mur2 * (Ezx0_n1(0, :) + Ezx0_n1(1, :))

            ! Condition aux bord (Mur absorbant) x = Nx
fd%Ez_r1(Nx , :) =    - Ezx_n2(Nx - 1, :)                                        &
                  + coef_mur1 * (fd%Ez_r1(Nx - 1, :) + Ezx_n2(Nx, :))            &
                  + coef_mur2 * (Ezx_n1(Nx, :) + Ezx_n1(Nx - 1, :))

            ! Condition aux bord (Mur absorbant) y = 0
fd%Ez_r1(: , 0)  =   - Ezy0_n2(1, :)                                             &
                + coef_mur1 * (fd%Ez_r1(:, 1) + Ezy0_n2(0, :))                   &
                + coef_mur2 * (Ezy0_n1(0, :) + Ezy0_n1(1, :))



            ! Condition aux bord (Mur absorbant) y = Ny
fd%Ez_r1(: , Ny) =   - Ezy_n2(Ny - 1, :)                                         &
                  + coef_mur1 * (fd%Ez_r1(:, Ny - 1) + Ezy_n2(Ny, :))            &
                  + coef_mur2 * (Ezy_n1(Ny, :) + Ezy_n1(Ny - 1,:))

             ! Source dans domaine grossier
fd%Ez_r1(250,250) = fd%Ez_r1(750, 750) + Esrc(n)

            ! Mise à jour Hx dans le domaine grossier       
do i = 0, Nx - 1
do j = 0, Ny - 1
     fd%Hx_r1(i,j) = fd%Hx_r1(i,j) - fd%c_H(i,j) / dy * (fd%Ez_r1(i,j+1) - fd%Ez_r1(i,j))
  
end do
end do

            ! Mise à jour  Hy dans le domaine grossier
do i = 0, Nx - 1
do j = 0, Ny - 1
fd%Hy_r1(i,j) = fd%Hy_r1(i,j) + fd%c_H(i,j) / dx * (fd%Ez_r1(i+1,j) - fd%Ez_r1(i,j))
 end do
end do

    write(20,*) n*dt,       fd%Ez_r1(350,350), fd%Ez_r1(250,250), fd%Ez_r1(350,250), fd%Ez_r1(250,350)
    write(21,*) n*dt,       fd%Hx_r1(100,100), fd%Hx_r1(150,100), fd%Hx_r1(50,100), fd%Hx_r1(100,150)
    write(22,*) n*dt,       fd%Hy_r1(100,100), fd%Hy_r1(150,100), fd%Hy_r1(50,100), fd%Hy_r1(100,150)


 !  if (mod(n, 100) == 0 .and. n > 0) then
 !   m = m + 1

 !   write(14,*) ! Ligne vide pour séparer les temps
 !   ! on n’écrit que la valeur Ez(i,j), séparée par un espace
   ! do i = 1, Nx, 2
  !     write(23,*) (fd%Ez(i,j), j= 0, Ny, 2)
   ! end do
 ! !!   write(23,*)


 ! end if   

 end do
    !    print *, "Nombre de carto en temps = ", m
     !   print *, "Fin de la simulation"
        close(20); close(21); close(22); close(23)
   
 end subroutine reference1

  subroutine reference2(fd, Nx_r, Ny_r, Nt, dx, dy, dt, Esrc)
        use lesconstantes_numeriques
        implicit none
    class(tableau), intent(inout)                              :: fd
    integer, intent(in)                                        :: Nx_r, Ny_r, Nt
    real(dp), intent(in)                                       :: dx, dy, dt
    real(dp)                                                   ::dx_s, dy_s
    real(dp), intent(in)                                       :: Esrc(0:Nt-1)
    integer                                                    :: i, j, n,  m 
    integer, parameter                                         :: r = 3
    real(dp), dimension(0 : 1       , 0 : Ny_r)                   :: Ezx0_n1 , Ezx0_n2                    
    real(dp), dimension(Nx_r - 1 : Nx_r , 0 : Ny_r)                  :: Ezx_n1 , Ezx_n2                 
    real(dp), dimension(0 : 1       , 0 : Nx_r)                    :: Ezy0_n1 , Ezy0_n2                     
    real(dp), dimension(Ny_r - 1 : Ny_r, 0 : Nx_r)                 :: Ezy_n1 , Ezy_n2
    real(dp)                                                       :: coef_mur1, coef_mur2, coef_mur3
   

!ouverture des fichiers pour stockes les resultats
open(unit=30, file="Ez_r2.dat",        status = 'replace', action = 'write')
open(unit=31, file="Hx_r2.dat",        status = 'replace', action = 'write')
open(unit=32, file="Hy_r2.dat",        status = 'replace', action = 'write')
!open(unit=23, file="carto_t2.dat",     status = 'replace', action = 'write')
 m = 0
Ezx0_n1 = 0.d0 ;    Ezx0_n2 = 0.0d0 
Ezx_n1  = 0.d0 ;    Ezx_n2  = 0.0d0
Ezy0_n1 = 0.d0 ;    Ezy0_n2 = 0.0d0
Ezy_n1  = 0.d0 ;    Ezy_n2  = 0.0d0

dx_s = dx / r   
dy_s = dy / r

coef_mur1 = (c * dt - dx_s) / (c * dt + dx_s)
coef_mur2 = (2.d0 * dx_s) / (c * dt + dx_s)
coef_mur3 = (c * dt)**2 / ( 2 * dx_s * ( c * dt + dx_s) )
 
do n = 0, Nt - 1
if (mod(n, 100) == 0) print *, "Itération temporelle n°", n
! Condition de bord x = 0
Ezx0_n2        = Ezx0_n1                        ! temps n - 2 pas i = 0 et i = 1
Ezx0_n1(0 , :) = fd%Ez_r2(0,:)                     ! temps n - 1 pas i = 0
Ezx0_n1(1 , :) = fd%Ez_r2(1,:)                     ! temps n - 1 pas i = 1

            ! Condition de bord x = Nx
Ezx_n2             = Ezx_n1                     ! temps n - 2 pas i = Nx - 1 et i = Nx
Ezx_n1(Nx_r - 1 , :) = fd%Ez_r2(Nx_r - 1 , :)          ! temps n - 1 pas i = Nx - 1
Ezx_n1(Nx_r , :)     = fd%Ez_r2(Nx_r , :)              ! temps n - 1 pas i = Nx

            ! Condition de bord y = 0
Ezy0_n2       = Ezy0_n1                         ! temps n - 2 pas j = 0 et j = 1
Ezy0_n1(0, :) = fd%Ez_r2(:, 0)                     ! temps n - 1 pas j = 0
Ezy0_n1(1, :) = fd%Ez_r2(:, 1)                     ! ; close(14)temps n - 1 pas j = 1

            ! Condition de bord y = Ny
Ezy_n2            = Ezy_n1                      ! temps n - 2 pas j = Ny - 1 et j = Ny
Ezy_n1(Ny_r - 1, :) = fd%Ez_r2(:, Ny_r - 1)            ! temps n - 1 pas j = Ny - 1
Ezy_n1(Ny_r, :)     = fd%Ez_r2(:, Ny_r)                ! temps n - 1 pas j = Ny
            
do i = 1, Nx_r - 1
do j = 1, Ny_r - 1
fd%Ez_r2(i,j) = fd%Ez_r2(i,j) + fd%c_E(i,j) * ((fd%Hy_r2(i,j) - fd%Hy_r2(i-1,j))/ dx_s  &
                        - (fd%Hx_r2(i,j) - fd%Hx_r2(i,j-1)) / dy_s)

end do
end do
 
            !! Condition aux bord (Mur absorbant) x = 0
fd%Ez_r2(0 , :)  =  - Ezx0_n2(1, :)                                             &
                 + coef_mur1 * (fd%Ez_r2(1, :) + Ezx0_n2(0, :))                 &         ! n2 -> le temps n - 1    n1 -> n
                  + coef_mur2 * (Ezx0_n1(0, :) + Ezx0_n1(1, :))

            ! Condition aux bord (Mur absorbant) x = Nx
fd%Ez_r2(Nx_r , :) =    - Ezx_n2(Nx_r - 1, :)                                        &
                  + coef_mur1 * (fd%Ez_r2(Nx_r - 1, :) + Ezx_n2(Nx_r, :))            &
                  + coef_mur2 * (Ezx_n1(Nx_r, :) + Ezx_n1(Nx_r - 1, :))

            ! Condition aux bord (Mur absorbant) y = 0
fd%Ez_r2(: , 0)  =   - Ezy0_n2(1, :)                                             &
                + coef_mur1 * (fd%Ez_r2(:, 1) + Ezy0_n2(0, :))                   &
                + coef_mur2 * (Ezy0_n1(0, :) + Ezy0_n1(1, :))



            ! Condition aux bord (Mur absorbant) y = Ny
fd%Ez_r2(: , Ny_r) =   - Ezy_n2(Ny_r - 1, :)                                         &
                  + coef_mur1 * (fd%Ez_r2(:, Ny_r - 1) + Ezy_n2(Ny_r, :))            &
                  + coef_mur2 * (Ezy_n1(Ny_r, :) + Ezy_n1(Ny_r - 1,:))

             ! Source dans domaine grossier
fd%Ez_r2(750,750) = fd%Ez_r2(750, 750) + Esrc(n)

            ! Mise à jour Hx dans le domaine grossier       
do i = 0, Nx_r - 1
do j = 0, Ny_r - 1
     fd%Hx_r2(i,j) = fd%Hx_r2(i,j) - fd%c_H(i,j) / dy_s * (fd%Ez_r2(i,j+1) - fd%Ez_r2(i,j))
  
end do
end do

            ! Mise à jour  Hy dans le domaine grossier
do i = 0, Nx_r - 1
do j = 0, Ny_r - 1
fd%Hy_r2(i,j) = fd%Hy_r2(i,j) + fd%c_H(i,j) / dx_s * (fd%Ez_r2(i+1,j) - fd%Ez_r2(i,j))
 end do
end do

    write(30,*) n*dt,       fd%Ez_r2(1050,1050), fd%Ez_r2(750,750), fd%Ez_r2(1050,750), fd%Ez_r2(750,1050)
    write(31,*) n*dt,       fd%Hx_r2(300,300), fd%Hx_r2(150,300), fd%Hx_r2(50,100), fd%Hx_r2(100,150)
    write(32,*) n*dt,       fd%Hy_r2(100,100), fd%Hy_r2(150,100), fd%Hy_r2(50,100), fd%Hy_r2(100,150)


 !  if (mod(n, 100) == 0 .and. n > 0) then
 !   m = m + 1

 !   write(14,*) ! Ligne vide pour séparer les temps
 !   ! on n’écrit que la valeur Ez(i,j), séparée par un espace
   ! do i = 1, Nx, 2
  !     write(23,*) (fd%Ez(i,j), j= 0, Ny, 2)
   ! end do
 ! !!   write(23,*)


 ! end if   

 end do
  
        close(30); close(31); close(32); close(33)
   
 end subroutine reference2

end module lesreference