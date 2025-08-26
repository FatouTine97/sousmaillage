module lesreference
use lesconstantes_numeriques
use lesfonction
implicit none
contains

!======================================================================
! Subroutine mise à jour des champs - Domaine grossier (r1)
!======================================================================
subroutine reference1(fd, Nx, Ny, Nt, dx, dt, dy, Esrc)
    use lesconstantes_numeriques
    implicit none
    class(tableau), intent(inout) :: fd
    integer, intent(in)           :: Nx, Ny, Nt
    real(dp), intent(in)          :: dx, dy, dt
    real(dp), intent(in)          :: Esrc(0:Nt-1)
    integer                       :: i, j, n, isc_x, isc_y,m
    real(dp)                      :: alpha_x, alpha_y
    real(dp), allocatable, save   :: Ez_prev(:,:)
    real(dp)                      :: t_start, t_end

    if (.not.allocated(Ez_prev)) allocate(Ez_prev(0:Nx,0:Ny), source=0.0_dp)

    isc_x = Nx/2
    isc_y = Ny/2

    alpha_x = (c*dt - dx) / (c*dt + dx)
    alpha_y = (c*dt - dy) / (c*dt + dy)

    open(unit=16, file="observation/Ez_r1.dat", status='replace', action='write')

    call cpu_time(t_start)   ! <<< Début du chrono

    do n = 0, Nt - 1
        ! Mise à jour Ez
        do i = 1, Nx-1
            do j = 1, Ny-1
                fd%Ez_r1(i,j) = fd%Ez_r1(i,j) + fd%c_E(i,j) * ((fd%Hy_r1(i,j) - fd%Hy_r1(i-1,j))/dx &
                                     - (fd%Hx_r1(i,j) - fd%Hx_r1(i,j-1))/dy)
            end do
        end do

        ! Conditions ABC
        do j = 0, Ny
            fd%Ez_r1(0,j)   = Ez_prev(1,j)   + alpha_x*(fd%Ez_r1(1,j)   - Ez_prev(0,j))
            fd%Ez_r1(Nx,j)  = Ez_prev(Nx-1,j)+ alpha_x*(fd%Ez_r1(Nx-1,j)- Ez_prev(Nx,j))
        end do
        do i = 0, Nx
            fd%Ez_r1(i,0)   = Ez_prev(i,1)   + alpha_y*(fd%Ez_r1(i,1)   - Ez_prev(i,0))
            fd%Ez_r1(i,Ny)  = Ez_prev(i,Ny-1)+ alpha_y*(fd%Ez_r1(i,Ny-1)- Ez_prev(i,Ny))
        end do

        Ez_prev = fd%Ez_r1

        ! Source
        fd%Ez_r1(isc_x, isc_y) = fd%Ez_r1(isc_x, isc_y) + Esrc(n)

        ! Hx
        do i = 0, Nx-1
            do j = 0, Ny-1
                fd%Hx_r1(i,j) = fd%Hx_r1(i,j) - (fd%c_H(i,j)/dy) * (fd%Ez_r1(i,j+1) - fd%Ez_r1(i,j))
            end do
        end do

        ! Hy
        do i = 0, Nx-1
            do j = 0, Ny-1
                fd%Hy_r1(i,j) = fd%Hy_r1(i,j) + (fd%c_H(i,j)/dx) * (fd%Ez_r1(i+1,j) - fd%Ez_r1(i,j))
            end do
        end do

        write(16,*) n*dt, fd%Ez_r1(250,250),  fd%Ez_r1(250,300), fd%Ez_r1(400,400)
    end do

    call cpu_time(t_end)     ! <<< Fin du chrono
    print *, 'Temps CPU (reference1) :', t_end - t_start, ' secondes'

end subroutine reference1


!======================================================================
! Subroutine mise à jour des champs - Domaine raffiné (r2)
!======================================================================
subroutine reference2(fd, Nx_r, Ny_r, Nt, dx_r, dy_r, dt_r, Esrc)
    use lesconstantes_numeriques
    implicit none
    class(tableau), intent(inout) :: fd
    integer, intent(in)           :: Nx_r, Ny_r, Nt
    real(dp), intent(in)          :: dx_r, dy_r, dt_r
    real(dp), intent(in)          :: Esrc(0:Nt-1)
    integer                       :: i, j, n, ics_x, ics_y
    real(dp)                      :: alpha_x, alpha_y, S
    real(dp), allocatable, save   :: Ez_prev(:,:)
    real(dp)                      :: t_start, t_end
    real(dp)                      :: cE_r, cHx_r, cHy_r
    integer                       :: ios

    ! (Re)allocation robuste de Ez_prev
    if (.not.allocated(Ez_prev)) then
        allocate(Ez_prev(0:Nx_r,0:Ny_r), source=0.0_dp)
    else if ( size(Ez_prev,1) /= Nx_r+1 .or. size(Ez_prev,2) /= Ny_r+1 ) then
        deallocate(Ez_prev)
        allocate(Ez_prev(0:Nx_r,0:Ny_r), source=0.0_dp)
    end if

    ics_x = Nx_r/2
    ics_y = Ny_r/2

    alpha_x = (c*dt_r - dx_r) / (c*dt_r + dx_r)
    alpha_y = (c*dt_r - dy_r) / (c*dt_r + dy_r)

    ! Coeffs (vide). Si milieux, remplace comme en r1 avec fd%c_E / fd%c_H spécifiques.
    cE_r  = dt_r/EPSILON_0
    cHx_r = dt_r/(mu_0*dy_r)
    cHy_r = dt_r/(mu_0*dx_r)

    S = c*dt_r*sqrt( (1.0_dp/dx_r**2) + (1.0_dp/dy_r**2) )
    if (S >= 1.0_dp) then
        print *, 'Avertissement: CFL (r2) >= 1 (S = ', S, ') — instable en 2D'
    end if

  !  call execute_command_line('mkdir -p observation', wait=.true.)

    open(unit=19, file='observation/Ez_r2.dat', status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        print *, 'Erreur ouverture fichier Ez_r2.dat, iostat=', ios
        return
    end if

    call cpu_time(t_start)

    do n = 0, Nt - 1
        ! ----- Ez (intérieur) -----
        do i = 1, Nx_r-1
            do j = 1, Ny_r-1
                fd%Ez_r2(i,j) = fd%Ez_r2(i,j) + cE_r * ( (fd%Hy_r2(i,j) - fd%Hy_r2(i-1,j))/dx_r &
                                         - (fd%Hx_r2(i,j) - fd%Hx_r2(i,j-1))/dy_r )
            end do
        end do

        ! ----- ABC (Mur 1) -----
        do j = 0, Ny_r
            fd%Ez_r2(0,j)    = Ez_prev(1,j)       + alpha_x*(fd%Ez_r2(1,j)       - Ez_prev(0,j))
            fd%Ez_r2(Nx_r,j) = Ez_prev(Nx_r-1,j)  + alpha_x*(fd%Ez_r2(Nx_r-1,j)  - Ez_prev(Nx_r,j))
        end do
        do i = 0, Nx_r
            fd%Ez_r2(i,0)    = Ez_prev(i,1)       + alpha_y*(fd%Ez_r2(i,1)       - Ez_prev(i,0))
            fd%Ez_r2(i,Ny_r) = Ez_prev(i,Ny_r-1)  + alpha_y*(fd%Ez_r2(i,Ny_r-1)  - Ez_prev(i,Ny_r))
        end do

        Ez_prev = fd%Ez_r2

        ! ----- Source (gain 3x conservé) -----
        fd%Ez_r2(ics_x, ics_y) = fd%Ez_r2(ics_x, ics_y) + (Esrc(n)*3.0_dp)

        ! ----- Hx -----
        do i = 0, Nx_r-1
            do j = 0, Ny_r-1
                fd%Hx_r2(i,j) = fd%Hx_r2(i,j) - cHx_r*(fd%Ez_r2(i,j+1) - fd%Ez_r2(i,j))
            end do
        end do

        ! ----- Hy -----
        do i = 0, Nx_r-1
            do j = 0, Ny_r-1
                fd%Hy_r2(i,j) = fd%Hy_r2(i,j) + cHy_r*(fd%Ez_r2(i+1,j) - fd%Ez_r2(i,j))
            end do
        end do

        write(19,*) n*dt_r, fd%Ez_r2(750,750), fd%Ez_r2(750,900), fd%Ez_r2(1200,1200)
        ! flush(19)
    end do

    call cpu_time(t_end)
    print *, 'Temps CPU (reference2) :', t_end - t_start, ' secondes'

    close(19)
end subroutine reference2


end module lesreference
