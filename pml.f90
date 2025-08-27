module cpml_2d
  use lesconstantes_numeriques
  implicit none
  private
  public :: cpml_t, cpml_init, cpml_update_ez, cpml_update_h

  type cpml_t
     integer :: nx, ny, npml
     real(dp), allocatable :: kappa_x(:), sigma_x(:), alpha_x(:), a_ex(:), b_ex(:), a_hx(:), b_hx(:)
     real(dp), allocatable :: kappa_y(:), sigma_y(:), alpha_y(:), a_ey(:), b_ey(:), a_hy(:), b_hy(:)
     real(dp), allocatable :: psi_ez_x(:,:), psi_ez_y(:,:), psi_hx_y(:,:), psi_hy_x(:,:)
  end type cpml_t

contains

  subroutine cpml_init(cp, nx, ny, dx, dy, dt, npml, R0, m, kappa_max, alpha_max)
    type(cpml_t), intent(inout) :: cp
    integer, intent(in) :: nx, ny, npml
    real(dp), intent(in) :: dx, dy, dt, R0, m, kappa_max, alpha_max
    integer :: i, j, np
    real(dp) :: s, kap, alp, sx, sy, sigma_max_x, sigma_max_y, eta0, denom

    cp%nx = nx; cp%ny = ny
    ! profondeur PML effective, bornée au domaine
    cp%npml = min( max(npml,0), min(nx,ny) )
    np = cp%npml

    eta0 = sqrt(mu_0/epsilon_0)

    allocate(cp%kappa_x(0:nx), cp%sigma_x(0:nx), cp%alpha_x(0:nx), &
             cp%a_ex(0:nx), cp%b_ex(0:nx), cp%a_hx(0:nx), cp%b_hx(0:nx))
    allocate(cp%kappa_y(0:ny), cp%sigma_y(0:ny), cp%alpha_y(0:ny), &
             cp%a_ey(0:ny), cp%b_ey(0:ny), cp%a_hy(0:ny), cp%b_hy(0:ny))
    allocate(cp%psi_ez_x(0:nx,0:ny), cp%psi_ez_y(0:nx,0:ny), &
             cp%psi_hx_y(0:nx,0:ny), cp%psi_hy_x(0:nx,0:ny))

    cp%psi_ez_x = 0.0_dp; cp%psi_ez_y = 0.0_dp
    cp%psi_hx_y = 0.0_dp; cp%psi_hy_x = 0.0_dp

    sigma_max_x = 0.0_dp; sigma_max_y = 0.0_dp
    if (np > 0) then
      sigma_max_x = -(m+1.0_dp) * log(R0) / (2.0_dp * eta0 * dx * real(np,dp))
      sigma_max_y = -(m+1.0_dp) * log(R0) / (2.0_dp * eta0 * dy * real(np,dp))
    end if

    cp%kappa_x = 1.0_dp; cp%alpha_x = 0.0_dp; cp%sigma_x = 0.0_dp
    cp%kappa_y = 1.0_dp; cp%alpha_y = 0.0_dp; cp%sigma_y = 0.0_dp

    if (np > 0) then
      ! X : remplir symétriquement gauche (i=0..np-1) et droite (i=nx..nx-np+1)
      do i = 0, np-1
        ! distance normalisée depuis le bord (1 à 0 vers l'intérieur)
        s   = real(np - i, dp) / real(np, dp)
        kap = 1.0_dp + (kappa_max-1.0_dp) * s**m
        alp = alpha_max * (1.0_dp - s)
        sx  = sigma_max_x * s**m
        cp%kappa_x(i)      = kap; cp%alpha_x(i)      = alp; cp%sigma_x(i)      = sx
        cp%kappa_x(cp%nx-i)= kap; cp%alpha_x(cp%nx-i)= alp; cp%sigma_x(cp%nx-i)= sx
      end do

      ! Y : bas (j=0..np-1) et haut (j=ny..ny-np+1)
      do j = 0, np-1
        s   = real(np - j, dp) / real(np, dp)
        kap = 1.0_dp + (kappa_max-1.0_dp) * s**m
        alp = alpha_max * (1.0_dp - s)
        sy  = sigma_max_y * s**m
        cp%kappa_y(j)      = kap; cp%alpha_y(j)      = alp; cp%sigma_y(j)      = sy
        cp%kappa_y(cp%ny-j)= kap; cp%alpha_y(cp%ny-j)= alp; cp%sigma_y(cp%ny-j)= sy
      end do
    end if

    ! Coeffs en X (élec/mag)
    do i = 0, nx
      cp%b_ex(i) = exp( - (cp%sigma_x(i)/max(cp%kappa_x(i),1.0e-30_dp) + cp%alpha_x(i)) * dt / epsilon_0 )
      denom = cp%sigma_x(i) + cp%kappa_x(i)*cp%alpha_x(i)
      if (denom > 0.0_dp) then
        cp%a_ex(i) = cp%sigma_x(i) * (cp%b_ex(i) - 1.0_dp) / ( denom * max(cp%kappa_x(i),1.0e-30_dp) )
      else
        cp%a_ex(i) = 0.0_dp
      end if

      cp%b_hx(i) = exp( - (cp%sigma_x(i)/max(cp%kappa_x(i),1.0e-30_dp) + cp%alpha_x(i)) * dt / mu_0 )
      denom = cp%sigma_x(i) + cp%kappa_x(i)*cp%alpha_x(i)
      if (denom > 0.0_dp) then
        cp%a_hx(i) = cp%sigma_x(i) * (cp%b_hx(i) - 1.0_dp) / ( denom * max(cp%kappa_x(i),1.0e-30_dp) )
      else
        cp%a_hx(i) = 0.0_dp
      end if
    end do

    ! Coeffs en Y (élec/mag)
    do j = 0, ny
      cp%b_ey(j) = exp( - (cp%sigma_y(j)/max(cp%kappa_y(j),1.0e-30_dp) + cp%alpha_y(j)) * dt / epsilon_0 )
      denom = cp%sigma_y(j) + cp%kappa_y(j)*cp%alpha_y(j)
      if (denom > 0.0_dp) then
        cp%a_ey(j) = cp%sigma_y(j) * (cp%b_ey(j) - 1.0_dp) / ( denom * max(cp%kappa_y(j),1.0e-30_dp) )
      else
        cp%a_ey(j) = 0.0_dp
      end if

      cp%b_hy(j) = exp( - (cp%sigma_y(j)/max(cp%kappa_y(j),1.0e-30_dp) + cp%alpha_y(j)) * dt / mu_0 )
      denom = cp%sigma_y(j) + cp%kappa_y(j)*cp%alpha_y(j)
      if (denom > 0.0_dp) then
        cp%a_hy(j) = cp%sigma_y(j) * (cp%b_hy(j) - 1.0_dp) / ( denom * max(cp%kappa_y(j),1.0e-30_dp) )
      else
        cp%a_hy(j) = 0.0_dp
      end if
    end do
  end subroutine cpml_init


  subroutine cpml_update_ez(cp, Ez, Hx, Hy, dx, dy, dt)
    type(cpml_t), intent(inout) :: cp
    real(dp), intent(inout) :: Ez(:,:), Hx(:,:), Hy(:,:)
    real(dp), intent(in)    :: dx, dy, dt
    integer :: i, j, nx, ny, np
    real(dp) :: dHy_dx, dHx_dy, term_x, term_y, kx, ky

    nx = cp%nx; ny = cp%ny; np = min(cp%npml, min(nx,ny))
    if (np <= 0) return

    ! vérif des formes (attend 0..nx, 0..ny)
    if ( size(Ez,1) /= nx+1 .or. size(Ez,2) /= ny+1 ) then
      print *, 'cpml_update_ez: shape(Ez)=', size(Ez,1), size(Ez,2), ' attendu=', nx+1, ny+1
      stop 1
    end if
    if ( size(Hx,1) /= nx+1 .or. size(Hx,2) /= ny+1 ) then
      print *, 'cpml_update_ez: shape(Hx) incorrecte'
      stop 1
    end if
    if ( size(Hy,1) /= nx+1 .or. size(Hy,2) /= ny+1 ) then
      print *, 'cpml_update_ez: shape(Hy) incorrecte'
      stop 1
    end if

    ! Bandes X (gauche)
    do j = 1, max(1, ny-1)
      do i = 1, min(np, nx-1)
        dHy_dx = (Hy(i,j) - Hy(i-1,j)) / dx
        cp%psi_ez_x(i,j) = cp%b_ex(i) * cp%psi_ez_x(i,j) + cp%a_ex(i) * dHy_dx
        kx = max(cp%kappa_x(i), 1.0e-30_dp)
        term_x = (dHy_dx / kx) + cp%psi_ez_x(i,j)

        dHx_dy = (Hx(i,j) - Hx(i,j-1)) / dy
        cp%psi_ez_y(i,j) = cp%b_ey(j) * cp%psi_ez_y(i,j) + cp%a_ey(j) * dHx_dy
        ky = max(cp%kappa_y(j), 1.0e-30_dp)
        term_y = (dHx_dy / ky) + cp%psi_ez_y(i,j)

        Ez(i,j) = Ez(i,j) + (dt/epsilon_0) * ( term_x - term_y )
      end do
    end do

    ! Bandes X (droite)
    do j = 1, max(1, ny-1)
      do i = max(1, nx-np+1), max(1, nx-1)
        dHy_dx = (Hy(i,j) - Hy(i-1,j)) / dx
        cp%psi_ez_x(i,j) = cp%b_ex(i) * cp%psi_ez_x(i,j) + cp%a_ex(i) * dHy_dx
        kx = max(cp%kappa_x(i), 1.0e-30_dp)
        term_x = (dHy_dx / kx) + cp%psi_ez_x(i,j)

        dHx_dy = (Hx(i,j) - Hx(i,j-1)) / dy
        cp%psi_ez_y(i,j) = cp%b_ey(j) * cp%psi_ez_y(i,j) + cp%a_ey(j) * dHx_dy
        ky = max(cp%kappa_y(j), 1.0e-30_dp)
        term_y = (dHx_dy / ky) + cp%psi_ez_y(i,j)

        Ez(i,j) = Ez(i,j) + (dt/epsilon_0) * ( term_x - term_y )
      end do
    end do

    ! Bandes Y (bas)
    do i = 1, max(1, nx-1)
      do j = 1, min(np, ny-1)
        dHy_dx = (Hy(i,j) - Hy(i-1,j)) / dx
        cp%psi_ez_x(i,j) = cp%b_ex(i) * cp%psi_ez_x(i,j) + cp%a_ex(i) * dHy_dx
        kx = max(cp%kappa_x(i), 1.0e-30_dp)
        term_x = (dHy_dx / kx) + cp%psi_ez_x(i,j)

        dHx_dy = (Hx(i,j) - Hx(i,j-1)) / dy
        cp%psi_ez_y(i,j) = cp%b_ey(j) * cp%psi_ez_y(i,j) + cp%a_ey(j) * dHx_dy
        ky = max(cp%kappa_y(j), 1.0e-30_dp)
        term_y = (dHx_dy / ky) + cp%psi_ez_y(i,j)

        Ez(i,j) = Ez(i,j) + (dt/epsilon_0) * ( term_x - term_y )
      end do
    end do

    ! Bandes Y (haut)
    do i = 1, max(1, nx-1)
      do j = max(1, ny-np+1), max(1, ny-1)
        dHy_dx = (Hy(i,j) - Hy(i-1,j)) / dx
        cp%psi_ez_x(i,j) = cp%b_ex(i) * cp%psi_ez_x(i,j) + cp%a_ex(i) * dHy_dx
        kx = max(cp%kappa_x(i), 1.0e-30_dp)
        term_x = (dHy_dx / kx) + cp%psi_ez_x(i,j)

        dHx_dy = (Hx(i,j) - Hx(i,j-1)) / dy
        cp%psi_ez_y(i,j) = cp%b_ey(j) * cp%psi_ez_y(i,j) + cp%a_ey(j) * dHx_dy
        ky = max(cp%kappa_y(j), 1.0e-30_dp)
        term_y = (dHx_dy / ky) + cp%psi_ez_y(i,j)

        Ez(i,j) = Ez(i,j) + (dt/epsilon_0) * ( term_x - term_y )
      end do
    end do
  end subroutine cpml_update_ez


  subroutine cpml_update_h(cp, Ez, Hx, Hy, dx, dy, dt)
    type(cpml_t), intent(inout) :: cp
    real(dp), intent(in)    :: dx, dy, dt
    real(dp), intent(in)    :: Ez(:,:)
    real(dp), intent(inout) :: Hx(:,:), Hy(:,:)
    integer :: i, j, nx, ny, np
    real(dp) :: dEz_dx, dEz_dy, kx, ky

    nx = cp%nx; ny = cp%ny; np = min(cp%npml, min(nx,ny))
    if (np <= 0) return

    if ( size(Ez,1) /= nx+1 .or. size(Ez,2) /= ny+1 ) then
      print *, 'cpml_update_h: shape(Ez) incorrecte'
      stop 1
    end if
    if ( size(Hx,1) /= nx+1 .or. size(Hx,2) /= ny+1 ) then
      print *, 'cpml_update_h: shape(Hx) incorrecte'
      stop 1
    end if
    if ( size(Hy,1) /= nx+1 .or. size(Hy,2) /= ny+1 ) then
      print *, 'cpml_update_h: shape(Hy) incorrecte'
      stop 1
    end if

    ! X -> Hy (gauche)
    do j = 0, max(0, ny-1)
      do i = 0, min(np-1, nx-1)
        dEz_dx = (Ez(i+1,j) - Ez(i,j)) / dx
        cp%psi_hy_x(i,j) = cp%b_hx(i) * cp%psi_hy_x(i,j) + cp%a_hx(i) * dEz_dx
        kx = max(cp%kappa_x(i), 1.0e-30_dp)
        Hy(i,j) = Hy(i,j) + (dt/mu_0) * ( (dEz_dx / kx) + cp%psi_hy_x(i,j) )
      end do
    end do

    ! X -> Hy (droite)
    do j = 0, max(0, ny-1)
      do i = max(0, nx-np), max(0, nx-1)
        dEz_dx = (Ez(i+1,j) - Ez(i,j)) / dx
        cp%psi_hy_x(i,j) = cp%b_hx(i) * cp%psi_hy_x(i,j) + cp%a_hx(i) * dEz_dx
        kx = max(cp%kappa_x(i), 1.0e-30_dp)
        Hy(i,j) = Hy(i,j) + (dt/mu_0) * ( (dEz_dx / kx) + cp%psi_hy_x(i,j) )
      end do
    end do

    ! Y -> Hx (bas)
    do i = 0, max(0, nx-1)
      do j = 0, min(np-1, ny-1)
        dEz_dy = (Ez(i,j+1) - Ez(i,j)) / dy
        cp%psi_hx_y(i,j) = cp%b_hy(j) * cp%psi_hx_y(i,j) + cp%a_hy(j) * dEz_dy
        ky = max(cp%kappa_y(j), 1.0e-30_dp)
        Hx(i,j) = Hx(i,j) - (dt/mu_0) * ( (dEz_dy / ky) + cp%psi_hx_y(i,j) )
      end do
    end do

    ! Y -> Hx (haut)
    do i = 0, max(0, nx-1)
      do j = max(0, ny-np), max(0, ny-1)
        dEz_dy = (Ez(i,j+1) - Ez(i,j)) / dy
        cp%psi_hx_y(i,j) = cp%b_hy(j) * cp%psi_hx_y(i,j) + cp%a_hy(j) * dEz_dy
        ky = max(cp%kappa_y(j), 1.0e-30_dp)
        Hx(i,j) = Hx(i,j) - (dt/mu_0) * ( (dEz_dy / ky) + cp%psi_hx_y(i,j) )
      end do
    end do
  end subroutine cpml_update_h

end module cpml_2d
