module lesreference
use lesconstantes_numeriques
use lesfonction
implicit none
contains

!======================================================================
! Subroutine mise à jour des champs - Domaine grossier (r1)
!======================================================================
subroutine reference1(fd, Nx, Ny, Nt, n, dx, dt, dy, Esrc)
    use lesconstantes_numeriques
    implicit none
    class(tableau), intent(inout) :: fd
    integer, intent(in)           :: Nx, Ny, Nt
    real(dp), intent(in)          :: dx, dy, dt
    real(dp), intent(in)          :: Esrc(0:Nt-1)
    integer                       :: i, j, n, isc_x, isc_y
    real(dp)                      :: alpha_x, alpha_y
    real(dp), allocatable, save   :: Ez_prev(:,:)

    !--- Allocation Ez_prev à la première itération ---
    if (.not.allocated(Ez_prev)) allocate(Ez_prev(0:Nx,0:Ny), source=0.0_dp)

    !--- Position de la source ---
    isc_x = Nx/2
    isc_y = Ny/2

    !--- Coefficients ABC 1er ordre ---
    alpha_x = (c*dt - dx) / (c*dt + dx)
    alpha_y = (c*dt - dy) / (c*dt + dy)

    !==============================
    !   Mise à jour Ez
    !==============================
    do i = 1, Nx-1
        do j = 1, Ny-1
            fd%Ez_r1(i,j) = fd%Ez_r1(i,j) + fd%c_E(i,j) * ((fd%Hy_r1(i,j) - fd%Hy_r1(i-1,j))/dx &
                                 - (fd%Hx_r1(i,j) - fd%Hx_r1(i,j-1))/dy)
        end do
    end do

    !=== Application des conditions ABC 1er ordre ===
    ! Bord gauche (i=0)
    do j = 0, Ny
        fd%Ez_r1(0,j) = Ez_prev(1,j) + alpha_x*(fd%Ez_r1(1,j) - Ez_prev(0,j))
    end do
    ! Bord droit (i=Nx)
    do j = 0, Ny
        fd%Ez_r1(Nx,j) = Ez_prev(Nx-1,j) + alpha_x*(fd%Ez_r1(Nx-1,j) - Ez_prev(Nx,j))
    end do
    ! Bord bas (j=0)
    do i = 0, Nx
        fd%Ez_r1(i,0) = Ez_prev(i,1) + alpha_y*(fd%Ez_r1(i,1) - Ez_prev(i,0))
    end do
    ! Bord haut (j=Ny)
    do i = 0, Nx
        fd%Ez_r1(i,Ny) = Ez_prev(i,Ny-1) + alpha_y*(fd%Ez_r1(i,Ny-1) - Ez_prev(i,Ny))
    end do

    ! Sauvegarde pour l'itération suivante
    Ez_prev = fd%Ez_r1

    !=== Source dans le domaine grossier ===
    fd%Ez_r1(isc_x, isc_y) = fd%Ez_r1(isc_x, isc_y) + Esrc(n)

    !==============================
    !   Mise à jour Hx
    !==============================
    do i = 0, Nx-1
        do j = 0, Ny-1
            fd%Hx_r1(i,j) = fd%Hx_r1(i,j) - (fd%c_H(i,j)/dy) * (fd%Ez_r1(i,j+1) - fd%Ez_r1(i,j))
        end do
    end do

    !==============================
    !   Mise à jour Hy
    !==============================
    do i = 0, Nx-1
        do j = 0, Ny-1
            fd%Hy_r1(i,j) = fd%Hy_r1(i,j) + (fd%c_H(i,j)/dx) * (fd%Ez_r1(i+1,j) - fd%Ez_r1(i,j))
        end do
    end do

end subroutine reference1

!======================================================================
! Subroutine mise à jour des champs - Domaine raffiné (r2)
!======================================================================
subroutine reference2(fd, Nx_r, Ny_r, Nt, n, dx_r, dy_r, dt_r, Esrc)
    use lesconstantes_numeriques
    implicit none
    class(tableau), intent(inout) :: fd
    integer, intent(in)           :: Nx_r, Ny_r, Nt
    real(dp), intent(in)          :: dx_r, dy_r, dt_r
    real(dp), intent(in)          :: Esrc(0:Nt-1)
    integer                       :: i, j, n, ics_x, ics_y
    real(dp)                      :: alpha_x, alpha_y
    real(dp), allocatable, save   :: Ez_prev(:,:)

    if (.not.allocated(Ez_prev)) allocate(Ez_prev(0:Nx_r,0:Ny_r), source=0.0_dp)

    ! Position source
    ics_x = Nx_r/2
    ics_y = Ny_r/2

    ! Coefficients ABC 1er ordre
    alpha_x = (c*dt_r - dx_r) / (c*dt_r + dx_r)
    alpha_y = (c*dt_r - dy_r) / (c*dt_r + dy_r)

    !==============================
    !   Mise à jour Ez raffiné
    !==============================
    do i = 1, Nx_r-1
        do j = 1, Ny_r-1
            fd%Ez_r2(i,j) = fd%Ez_r2(i,j) + (dt_r/EPSILON_0) * ((fd%Hy_r2(i,j) - fd%Hy_r2(i-1,j))/dx_r &
                                     - (fd%Hx_r2(i,j) - fd%Hx_r2(i,j-1))/dy_r)
        end do
    end do

    !=== Application ABC 1er ordre ===
    do j = 0, Ny_r
        fd%Ez_r2(0,j) = Ez_prev(1,j) + alpha_x*(fd%Ez_r2(1,j) - Ez_prev(0,j))
        fd%Ez_r2(Nx_r,j) = Ez_prev(Nx_r-1,j) + alpha_x*(fd%Ez_r2(Nx_r-1,j) - Ez_prev(Nx_r,j))
    end do
    do i = 0, Nx_r
        fd%Ez_r2(i,0) = Ez_prev(i,1) + alpha_y*(fd%Ez_r2(i,1) - Ez_prev(i,0))
        fd%Ez_r2(i,Ny_r) = Ez_prev(i,Ny_r-1) + alpha_y*(fd%Ez_r2(i,Ny_r-1) - Ez_prev(i,Ny_r))
    end do

    Ez_prev = fd%Ez_r2

    !=== Source dans domaine raffiné ===
    fd%Ez_r2(ics_x, ics_y) = fd%Ez_r2(ics_x, ics_y) + Esrc(n)

    !==============================
    !   Mise à jour Hx raffiné
    !==============================
    do i = 0, Nx_r-1
        do j = 0, Ny_r-1
            fd%Hx_r2(i,j) = fd%Hx_r2(i,j) - (dt_r/(mu_0*dy_r)) * (fd%Ez_r2(i,j+1) - fd%Ez_r2(i,j))
        end do
    end do

    !==============================
    !   Mise à jour Hy raffiné
    !==============================
    do i = 0, Nx_r-1
        do j = 0, Ny_r-1
            fd%Hy_r2(i,j) = fd%Hy_r2(i,j) + (dt_r/(mu_0*dx_r)) * (fd%Ez_r2(i+1,j) - fd%Ez_r2(i,j))
        end do
    end do

end subroutine reference2

end module lesreference
