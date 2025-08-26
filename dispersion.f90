module dispersion_tmz2d_cg_fg
  use lesconstantes_numeriques
  implicit none
contains

subroutine dispersion_numerique(dx, dy, r, dx_sm, dy_sm, dt, dt_prime)
    integer, intent(in) :: r
    real(dp), intent(in) :: dx, dy, dx_sm, dy_sm, dt, dt_prime
  integer, parameter :: Nk = 200
  integer, parameter :: Ntheta = 13
  real(dp) :: k, kmax, dk, theta, dtheta, kx, ky
  real(dp) :: omega_ex, omega_GM, omega_PM
  real(dp) :: vp_GM, vp_PM, err_GM, err_PM
  integer :: i, j, u
 ! character(len=*), parameter :: fname = "dispersion_tmz2d_CG_FG.dat"
  !real(dp) :: dt_prime

 ! dt_prime = dt / r
  kmax = 0.5_dp * (2.0_dp*acos(-1.0_dp)) / max(dx, dy)
  dk   = kmax / real(Nk, dp)
  dtheta = (0.5_dp*acos(-1.0_dp)) / real(Ntheta-1, dp)

  open(unit=20, file="dispersion.dat", status= "replace", action="write")
  write(20,'(A)') "# k  theta(rad)  omega_exact  omega_CG  omega_FG  vp_CG  vp_FG  err_CG  err_FG"

  do j = 0, Ntheta-1
     theta = real(j,dp) * dtheta
     do i = 1, Nk
        k = real(i,dp) * dk
        kx = k * cos(theta)
        ky = k * sin(theta)

        omega_ex = c * k

        omega_GM = vp_over_c(dx, dy, dt, k*dx, theta)!omega_fdtd_tmz(kx, ky, dx, dy, dt, c)
        omega_PM =  vp_over_c(dx_sm, dy_sm, dt_prime, k*dx_sm, theta) !omega_fdtd_tmz(kx, ky, dx_sm, dy_sm, dt_prime, c)

        if (k > 0.0_dp) then
           vp_GM = omega_GM / k
           vp_PM = omega_PM / k
        else
           vp_GM = c
           vp_PM = c
        end if
        err_GM = (vp_GM - c) / c
        err_PM = (vp_PM - c) / c

     !   write(20,'(9e16.8)') k, theta, omega_ex, omega_GM, omega_PM, vp_GM, vp_PM, err_GM, err_PM
     end do
  end do
  close(20)

  !print *, "OK: dispersion 2D TMz écrite dans "
 ! print *, "Colonnes: k, theta, omega_exact, omega_CG, omega_FG, vp_CG, vp_FG, err_CG, err_FG"
end subroutine dispersion_numerique

pure function vp_over_c(dx, dy, dt, kDelta, theta) result(vp_c)
  use lesconstantes_numeriques
  implicit none
  real(dp), intent(in) :: dx, dy, dt, kDelta, theta
  real(dp) :: vp_c, kx, ky, omega, sx, sy, s, pi

  kx = (kDelta/dx) * cos(theta)
  ky = (kDelta/dy) * sin(theta)

  sx = sin(0.5_dp * kx * dx)
  sy = sin(0.5_dp * ky * dy)
  s  = (c*dt/dx)**2 * sx*sx + (c*dt/dy)**2 * sy*sy
  if (s > 1.0_dp) s = 1.0_dp
  if (s < 0.0_dp) s = 0.0_dp

  omega = (2.0_dp/dt) * asin(sqrt(s))
  vp_c  = (omega / ( (kDelta/dx) * c ))   ! car k = kDelta/dx si dx=dy
end function vp_over_c

end module dispersion_tmz2d_cg_fg