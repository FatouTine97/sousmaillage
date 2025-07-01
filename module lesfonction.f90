module lesfonction
use lesconstantes_numeriques
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

        allocate(this%c_H(0:Nx, 0:Ny), &
        this%c_E(0:Nx, 0:Ny), this%Ez(0:Nx, 0:Ny),&
         this%Hx(0:Nx, 0:Ny), &
        this%Hy(0:Nx, 0:Ny))

        this%Ez = 0.0_dp
        this%Hy = 0.0_dp
        this%Hx = 0.0_dp

        this%c_E = dt / EPSILON_0
        this%c_H = dt / MU_0
        
    end subroutine initialiser_champs


    ! Subroutine pour la mise à jour des champs
    subroutine mise_a_jour_champs(fd, Nx, Ny, Nt, dx, dt, dy, Nx_sm, Ny_sm, Esrc)
        use lesconstantes_numeriques
        implicit none
        class(tableau), intent(inout) :: fd
        integer, intent(in) :: Nx, Ny, Nt, Nx_sm, Ny_sm
        real(dp), intent(in) :: dx, dy, dt
        real(dp), intent(in) :: Esrc(0:Nt-1)
        integer, parameter :: r = 3
        integer :: i, j, n,  m,i_f, j_f,  i_1, k
        real(dp) :: Hz2, Hz3, dx_sm, dy_sm, dt_prime
        real(dp) :: Haux_y, dx_prime, dy_prime, f_Hy1, eaux_z, fez
        real(dp) :: LCC, LIC1, LIC2, TCC1, TCC2, alpha 

        dx_sm = dx/ r
        dy_sm = dy/ r
        dx_prime = (dx/2.0_dp) + (dx_sm/2.0_dp)
        dy_prime = dy/2.0_dp +dy_sm/2.0_dp
        f_Hy1 = dy/dy_sm
        fez = 1.0
        dt_prime= dt/r

open(unit=10, file="Ez_t.dat",    status = 'replace', action = 'write')
open(unit=11, file="Hx_t.dat",    status = 'replace', action = 'write')
open(unit=12, file="Hy_t.dat",    status = 'replace', action = 'write')
open(unit=13, file="carto_t.dat", status = 'replace', action = 'write')
open(unit=14, file="carto_t1.dat",status = 'replace', action = 'write')
open(unit=15, file="ez_s.dat",status = 'replace', action = 'write')

        m = 0
        do n = 0, Nt - 1
            if (mod(n, 100) == 0) print *, "Itération temporelle n°", n

            ! Mise à jour du champ Ez dans le domaine grossier
            ! Mise à jour du champ électrique
            do i = 1, Nx - 1
             if (i < i1 ) then
                do j = 1, Ny - 1
                    fd%Ez(i,j) = fd%Ez(i,j) + fd%c_E(i,j) * ((fd%Hy(i,j) - fd%Hy(i-1,j))/ dx - &
                    (fd%Hx(i,j) - fd%Hx(i,j-1)) / dy)
               end do
               end if
            end do

            ! Source dans domaine grossier
            fd%Ez(150, 150) = fd%Ez(150, 150) + Esrc(n)
            ! Conditions aux limites Ez
            fd%Ez(0,:) = 0.0_dp
            fd%Ez(Nx,:) = 0.0_dp
            fd%Ez(:,0) = 0.0_dp
            fd%Ez(:,Ny) = 0.0_dp

            ! Mise à jour Hx dans le domaine grossier
                
            do i = 0, Nx - 1
             if (i < i1 ) then
            do j = 0, Ny - 1
             fd%Hx(i,j) = fd%Hx(i,j) - fd%c_H(i,j) / dy * (fd%Ez(i,j+1) - fd%Ez(i,j))
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


 ! Mise à jour dans le sous-maillage avec r sous-pas de temps
         ! do k = 1, r
! Ez raffiné
                do i = 0, Nx_sm - 1
                 do j = 1, Ny_sm - 1 
                 
    fd%ez_s(i,j) = fd%ez_s(i,j) + dt/ (epsilon_0 * dx_sm) * (fd%hy_s(i,j) - &
    fd%hy_s(i-1,j)) - (dt/r) / (epsilon_0 * dy_sm) * (fd%hx_s(i,j) - fd%hx_s(i,j-1))
     end do
    end do
    ! Source dans le sous-maillage
fd%ez_s(350, Ny_sm/2) = fd%ez_s(350, Ny_sm/2) + Esrc(n)

                ! Hx raffinés
do i = 0, Nx_sm - 1
do j = 0, Ny_sm - 1
    fd%hx_s(i,j) = fd%hx_s(i,j) - dt/ (mu_0 * dy_sm) * (fd%ez_s(i,j+1) - fd%ez_s(i,j))
             end do
            end do
            
 !Hy raffinés
do i = 0, Nx_sm - 1
 do j = 0, Ny_sm - 1
 fd%hy_s(i,j) = fd%hy_s(i,j) +  dt/ (mu_0 * dx_sm) * (fd%ez_s(i+1,j) - fd%ez_s(i,j))
    end do
    end do
! le champs electrique dans l'interface
 do j_f = 0, Ny_sm - 1
  i_f = i1
 alpha = 1.0/3.0
 LCC = 1.0_dp / real(r, dp)
LIC1 = 1.0_dp - alpha
LIC2 =  alpha

  TCC1 = LIC1 * LCC
  TCC2 = LIC2 * LCC       
        
Hz2 = fd%Hy(i_f, j_f)  
Hz3 = fd%Hy(i_f, j_f+1)

 Haux_y = TCC1*Hz2  + TCC2*Hz3

fd%ez_s(i_f,j_f) = fd%ez_s(i_f,j_f) +  dt/ (epsilon_0 ) *&
        (fd%hy_s(i_f,j_f) - (f_Hy1 *Haux_y)/dx_prime -&
        (fd%hx_s(i_f,j_f) - fd%hx_s(i_f,j_f-1))/(dy_sm))
     
end do

!champs magnetiques dans l'interface
  do j = 1, Nx -1
i = i1
 eaux_z = (1.0d0 / 3.0d0) * ( &
     (1.0d0 / 3.0d0 * fd%ez_s(i, j) + 2.0d0 / 3.0d0 * fd%ez_s(i, j+1) + fd%ez_s(i, j+2) + & 
      2.0d0 / 3.0d0 * fd%ez_s(i, j+3) + 1.0d0 / 3.0d0 * fd%ez_s(i, j+4)) + &
    (1.0d0 / 3.0d0 * fd%ez_s(i, j+5) + 2.0d0 / 3.0d0 * fd%ez_s(i, j+6) + fd%ez_s(i, j+7) + & 
       2.0d0 / 3.0d0 * fd%ez_s(i, j+8) + 1.0d0 / 3.0d0 * fd%ez_s(i, j+9)) )

  ! Mise à jour de Hy à l'interface
  fd%Hy(i, j) =fd%Hy(i, j) +   (dt_prime/ (mu_0 * dx)) * (fez * eaux_z - fd%Ez(i, j))
end do 
        
!end do
    write(10,*) n*dt, fd%Ez(100,100), fd%Ez(150,100), fd%Ez(50,100), fd%Ez(1,150)
    write(11,*) n*dt, fd%Hx(100,100), fd%Hx(150,100), fd%Hx(50,100), fd%Hx(1,150)
    write(12,*) n*dt, fd%Hy(100,100), fd%Hy(150,100), fd%Hy(50,100), fd%Hy(1,150)
if (mod(n, 100) == 0 .and. n > 0) then
    m = m + 1
    do i = 1, Nx_sm, 2
         write(13,*) (fd%ez_s(i,j), j=0, Ny_sm, 2)
    end do
 write(13,*)
     do i_1 = 1, Nx, 2
    write(14,*) (fd%Ez(i_1,j), j= 0, Ny, 2)
        end do
     write(14,*)
  end if   

    write(15,*) n*dt_prime, fd%ez_s(300,250)
        end do
        print *, "Nombre de carto en temps = ", m
        print *, "Fin de la simulation"
        close(10); close(11); close(12); close(13); close(14); close(15)
    end subroutine mise_a_jour_champs
  
end module lesfonction  






