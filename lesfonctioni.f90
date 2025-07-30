module lesfonction
use lesconstantes_numeriques
use lexcitation_donde
implicit none

type :: tableau
    real(dp), dimension(:,:), allocatable :: Hx, Ez, Hy
    real(dp), dimension(:,:), allocatable :: Hx_r1, Ez_r1, Hy_r1
    real(dp), dimension(:,:), allocatable :: Ez_r2, Hx_r2, Hy_r2
    real(dp), dimension(:,:), allocatable :: hx_s, ez_s, hy_s
    real(dp), dimension(:,:), allocatable :: c_E, c_H
contains
    procedure :: initialiser_champs
    procedure :: mise_a_jour_champs
end type tableau

contains

!=======================================================================
!   INITIALISATION DES CHAMPS
!=======================================================================
subroutine initialiser_champs(this, Nx, Ny, Nx_sm, Ny_sm, Nx_r, Ny_r, dt)
    use lesconstantes_numeriques
    implicit none
    class(tableau), intent(inout) :: this
    integer, intent(in) :: Nx, Ny, Nx_sm, Ny_sm, Nx_r, Ny_r
    real(dp), intent(in) :: dt

    ! Allocation sous-maillage
    allocate(this%hx_s(0:Nx_sm, 0:Ny_sm), this%ez_s(0:Nx_sm, 0:Ny_sm), this%hy_s(0:Nx_sm, 0:Ny_sm))
    this%hx_s = 0.0_dp
    this%ez_s = 0.0_dp
    this%hy_s = 0.0_dp

    ! Allocation domaine grossier
    allocate(this%Ez(0:Nx, 0:Ny), this%Hx(0:Nx, 0:Ny), this%Hy(0:Nx, 0:Ny))
    this%Ez = 0.0_dp
    this%Hy = 0.0_dp
    this%Hx = 0.0_dp

    ! Allocation référence 1
    allocate(this%Hx_r1(0:Nx, 0:Ny), this%Ez_r1(0:Nx, 0:Ny), this%Hy_r1(0:Nx, 0:Ny)) 
    this%Ez_r1 = 0.0_dp  
    this%Hy_r1 = 0.0_dp
    this%Hx_r1 = 0.0_dp

    ! Allocation référence 2
    allocate(this%Ez_r2(0:Nx_r, 0:Ny_r), this%Hx_r2(0:Nx_r, 0:Ny_r), this%Hy_r2(0:Nx_r, 0:Ny_r))
    this%Ez_r2 = 0.0_dp
    this%Hy_r2 = 0.0_dp
    this%Hx_r2 = 0.0_dp

    ! Coefficients de mise à jour
    allocate(this%c_H(0:Nx, 0:Ny), this%c_E(0:Nx, 0:Ny))
    this%c_E = dt / EPSILON_0
    this%c_H = dt / MU_0
end subroutine initialiser_champs

!=======================================================================
!   MISE À JOUR DES CHAMPS
!=======================================================================
subroutine mise_a_jour_champs(fd, Nx, Ny, Nt, dx, dt, dy, Nx_sm, Ny_sm, Esrc)
    use lesconstantes_numeriques
    implicit none
    class(tableau), intent(inout) :: fd
    integer, intent(in) :: Nx, Ny, Nt, Nx_sm, Ny_sm
    real(dp), intent(in) :: dx, dy, dt
    real(dp), intent(in) :: Esrc(0:Nt-1)

    integer, parameter :: r = 3
    integer :: i, j, jr, ig
    real(dp) :: dx_sm, dy_sm, dt_prime
    real(dp) :: dx_prime, dy_prime, coefMur
    real(dp) :: e1, e2, e3, e4, e5, eaux_z
    real(dp) :: Hz2, Hz3
    real(dp), allocatable :: Haux_y(:)
    integer :: i0, j0, n,m

      

    ! Source au centre
    i0 = Nx / 2
    j0 = Ny / 2

    ! Facteurs
    dx_sm = dx / r
    dy_sm = dy / r
    dt_prime = dt / r
    dx_prime = 0.5_dp * dx + 0.5_dp * dx_sm
    dy_prime = 0.5_dp * dy + 0.5_dp * dy_sm

    coefMur = (c * dt - dx) / (c * dt + dx)

    allocate(Haux_y(0:Ny_sm-1))

 ! Fichiers de sortie
open(unit=10, file="observation/Ez_t.dat", status='replace', action='write')
open(unit=11, file="observation/Hx_t.dat", status='replace', action='write')
open(unit=12, file="observation/Hy_t.dat", status='replace', action='write')
open(unit=13, file="observation/carto_t.dat", status='replace', action='write')
open(unit=14, file="observation/carto_t1.dat", status='replace', action='write')
open(unit=15, file="observation/ez_s.dat", status='replace', action='write')


m = 0
do n = 0, Nt- 1
   if (mod(n, 100) == 0) print *, "Itération n°", n
!===========================================================
! 1. Mise à jour Ez dans le domaine grossier (sauf interface)
!===========================================================

do i = 1, Nx - 1
    if (i < i1) then
        do j = 1, Ny - 1
            fd%Ez(i,j) = fd%Ez(i,j) + fd%c_E(i,j) * ((fd%Hy(i,j) - fd%Hy(i-1,j))/dx - &
                                                     (fd%Hx(i,j) - fd%Hx(i,j-1))/dy)
        end do
    end if
end do


!===========================================================
! 2. Conditions aux limites absorbantes (Mur)
!===========================================================
fd%Ez(0,:) = fd%Ez(1,:) + coefMur * (fd%Ez(1,:) - fd%Ez(0,:))
fd%Ez(Nx,:) = fd%Ez(Nx-1,:) + coefMur * (fd%Ez(Nx-1,:) - fd%Ez(Nx,:))
fd%Ez(:,0) = fd%Ez(:,1) + coefMur * (fd%Ez(:,1) - fd%Ez(:,0))
fd%Ez(:,Ny) = fd%Ez(:,Ny-1) + coefMur * (fd%Ez(:,Ny-1) - fd%Ez(:,Ny))

!===========================================================
! 3. Source
!===========================================================
fd%Ez(i0,j0) = fd%Ez(i0,j0) + Esrc(n)*3
!===========================================================
! 4. Mise à jour Hx, Hy (domaine grossier)
!===========================================================
do i = 0, Nx - 1
    if (i < i1-1) then
        do j = 0, Ny - 1
            fd%Hx(i,j) = fd%Hx(i,j) - fd%c_H(i,j) / dy * (fd%Ez(i,j+1) - fd%Ez(i,j))
        end do
    end if
end do

do i = 0, Nx - 1
    if (i < i1-1) then
        do j = 0, Ny - 1
            fd%Hy(i,j) = fd%Hy(i,j) + fd%c_H(i,j) / dx * (fd%Ez(i+1,j) - fd%Ez(i,j))
        end do
    end if
end do

!===========================================================
! 5. Interpolation FG → CG (interface Ez)
!===========================================================
do j = 2, Ny_sm-3
    jr = int(j/r)
    e1 = fd%ez_s(0,j-2)
    e2 = fd%ez_s(0,j-1)
    e3 = fd%ez_s(0,j)
    e4 = fd%ez_s(0,j+1)
    e5 = fd%ez_s(0,j+2)

    eaux_z = (1.0_dp/9.0_dp)*e1 + (2.0_dp/9.0_dp)*e2 + (1.0_dp/30_dp)*e3 + (2.0_dp/9.0_dp)*e4 +&
                (1.0_dp/9.0_dp)*e5 + (1.0_dp/9.0_dp)*e1 + (2.0_dp/9.0_dp)*e2 + (1.0_dp/30_dp)*e3 + (2.0_dp/9.0_dp)*e4 +&
                (1.0_dp/9.0_dp)*e5

    fd%Hy(i1-1,jr) = fd%Hy(i1-1,jr) + (fd%c_H(i1-1,jr)/dx) * (eaux_z - fd%Ez(i1-1,jr))
end do

!===========================================================
! 6. Mise à jour sous-maillage
!===========================================================
do j = 1, Ny_sm-1
    ig = int(j/r)
    Hz2 = fd%Hy(i1-1, ig)
    Hz3 = fd%Hy(i1-1, ig+1)
    Haux_y(j) = (2.0_dp/9.0_dp)*Hz2 + (1.0_dp/9.0_dp)*Hz3

    fd%ez_s(0,j) = fd%ez_s(0,j) + (dt_prime / epsilon_0) * ((fd%hy_s(0,j) - Haux_y(j))/dx_prime - &
                                                             (fd%hx_s(0,j) - fd%hx_s(0,j-1))/dy_sm)
end do


 ! -- Mise à jour Ez_s (zone interne FG)
    do i = 1, Nx_sm-1
        do j = 1, Ny_sm-1
            fd%ez_s(i,j) = fd%ez_s(i,j) + (dt_prime/epsilon_0) * &
                           ((fd%hy_s(i,j) - fd%hy_s(i-1,j))/dx_sm - &
                            (fd%hx_s(i,j) - fd%hx_s(i,j-1))/dy_sm)
        end do
    end do

    ! -- Mise à jour Hx_s
    do i = 0, Nx_sm-1
        do j = 0, Ny_sm-1   ! éviter dépassement
            fd%hx_s(i,j) = fd%hx_s(i,j) - (dt_prime/mu_0) * &
                           (fd%ez_s(i,j+1) - fd%ez_s(i,j)) / dy_sm
        end do
    end do

    ! -- Mise à jour Hy_s
    do i = 0, Nx_sm-1  ! éviter dépassement
        do j = 0, Ny_sm-1
            fd%hy_s(i,j) = fd%hy_s(i,j) + (dt_prime/mu_0) * &
                           (fd%ez_s(i+1,j) - fd%ez_s(i,j)) / dx_sm
        end do
    end do

write(10,*)n*dt_prime,   fd%Ez(250,250),     fd%Ez(250,300),   fd%ez_s(300,1200)
      ! write(15,*) n*dt,                                     , sm%ez_s(500,900)
! write(11,*) n*dt, sm%Hx(250,250), sm%Hx(150,100), sm%Hx(50,100), sm%Hx(100,150)
!!  write(12,*) n*dt, sm%Hy(250,250), sm%Hy(150,100), sm%Hy(50,100), sm%Hy(100,150)
     
    if (mod(n/3, 100) == 0) then
      m = m + 1
      do i = 0, Nx, 4
         write(14,*) (fd%Ez(i,j), j=0, Ny, 4)
      end do
      write(14,*)

      do i = 0, Nx_sm, 4
         write(13,*) (fd%ez_s(i,j), j=0, Ny_sm, 4)
      end do
      write(13,*)
   end if

end do 

 end subroutine mise_a_jour_champs
  

 real(dp) function compute_ez_aux(fd, i_f, j_f) result(ez_aux)
    use lesconstantes_numeriques
    implicit none
    type(tableau), intent(in) :: fd
    integer, intent(in)       :: i_f, j_f
    integer :: ii, jj
    real(dp) :: sum
    real(dp), dimension(-1:1, -1:1) :: w

    ! Pondération conservatrice 3x3, somme = 1
    w = reshape([ &
        1.0_dp/16.0_dp, 2.0_dp/16.0_dp, 1.0_dp/16.0_dp, &
        2.0_dp/16.0_dp, 4.0_dp/16.0_dp, 2.0_dp/16.0_dp, &
        1.0_dp/16.0_dp, 2.0_dp/16.0_dp, 1.0_dp/16.0_dp], &
        [3,3])

       sum = 0.0_dp
    do ii = -1, 1
        do jj = -1, 1
            if (i_f >= lbound(fd%ez_s,1) .and. i_f <= ubound(fd%ez_s,1) .and. &
                j_f+jj >= lbound(fd%ez_s,2) .and. j_f+jj <= ubound(fd%ez_s,2)) then
                sum = sum + w(ii+2, jj+2) * fd%ez_s(i_f, j_f + jj)
            end if
        end do
    end do
    ez_aux = sum
end function compute_ez_aux



end module lesfonction  