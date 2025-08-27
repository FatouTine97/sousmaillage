module lesreference
  use lesconstantes_numeriques   ! doit fournir dp, c, epsilon_0, mu_0
  use lesfonction                ! doit fournir le type derive "tableau" avec Ez_r1, Hx_r1, Hy_r1, c_E, c_H, etc.
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

  integer  :: i, j, n, isc_x, isc_y
  real(dp) :: alpha_x, alpha_y, S
  real(dp), allocatable, save :: Ez_prev(:,:)
  real(dp) :: t_start, t_end
  integer, parameter :: nobs = 3
  integer :: xobs(nobs), yobs(nobs)
  integer :: ix_fix, iy_fix
  real(dp) :: Ez_fix

  ! Allocation robuste d'Ez_prev
  if (.not. allocated(Ez_prev)) then
     allocate(Ez_prev(0:Nx,0:Ny), source=0.0_dp)
  else if ( size(Ez_prev,1) /= Nx+1 .or. size(Ez_prev,2) /= Ny+1 ) then
     deallocate(Ez_prev)
     allocate(Ez_prev(0:Nx,0:Ny), source=0.0_dp)
  end if

  ! Source au centre
  isc_x = Nx/2
  isc_y = Ny/2

  ! Coefficient de Courant (info only)
  S = c*dt*sqrt( (1.0_dp/dx**2) + (1.0_dp/dy**2) )
  if (S >= 1.0_dp) then
     print *, 'Avertissement: CFL (r1) >= 1 (S = ', S, ') — instable en 2D'
  end if

  ! ABC Mur 1er ordre : coefficients distincts
  alpha_x = (c*dt - dx) / (c*dt + dx)
  alpha_y = (c*dt - dy) / (c*dt + dy)

  ! Points d'observation
  xobs = [ Nx/2, 3*Nx/4, Nx/2 ]
  yobs = [ Ny/2, 3*Ny/4, Ny/2 ]

  ! Un 4e point "fixe" optionnel (protégé)
  ix_fix = 400
  iy_fix = 400

  call execute_command_line('mkdir -p observation', wait=.true.)
  open(unit=16, file="observation/Ez_r1.dat", status='replace', action='write')

  call cpu_time(t_start)

  do n = 0, Nt - 1
     ! --- sauvegarde Ez(n) pour l'ABC (temps n) ---
     Ez_prev = fd%Ez_r1

     ! --- 1) Update Ez intérieur -> temps n+1 ---
     do i = 1, Nx-1
        do j = 1, Ny-1
           fd%Ez_r1(i,j) = fd%Ez_r1(i,j) + fd%c_E(i,j) * ( (fd%Hy_r1(i,j) - fd%Hy_r1(i-1,j))/dx &
                                                        - (fd%Hx_r1(i,j) - fd%Hx_r1(i,j-1))/dy )
        end do
     end do

     ! --- 2) ABC de Mur (utilise Ez_prev = temps n) ---
     ! Bords gauche/droit : j = 1..Ny-1
     fd%Ez_r1(0, 1:Ny-1)   = Ez_prev(1, 1:Ny-1)    + alpha_x * ( fd%Ez_r1(1, 1:Ny-1)    - Ez_prev(0, 1:Ny-1) )
     fd%Ez_r1(Nx,1:Ny-1)   = Ez_prev(Nx-1,1:Ny-1)  + alpha_x * ( fd%Ez_r1(Nx-1,1:Ny-1)  - Ez_prev(Nx,1:Ny-1) )

     ! Bords bas/haut : i = 1..Nx-1
     fd%Ez_r1(1:Nx-1,0)    = Ez_prev(1:Nx-1,1)     + alpha_y * ( fd%Ez_r1(1:Nx-1,1)     - Ez_prev(1:Nx-1,0) )
     fd%Ez_r1(1:Nx-1,Ny)   = Ez_prev(1:Nx-1,Ny-1)  + alpha_y * ( fd%Ez_r1(1:Nx-1,Ny-1)  - Ez_prev(1:Nx-1,Ny) )

     ! Coins (moyenne simple)
     fd%Ez_r1(0,0)         = 0.5_dp*( fd%Ez_r1(1,0)     + fd%Ez_r1(0,1) )
     fd%Ez_r1(Nx,0)        = 0.5_dp*( fd%Ez_r1(Nx-1,0)  + fd%Ez_r1(Nx,1) )
     fd%Ez_r1(0,Ny)        = 0.5_dp*( fd%Ez_r1(1,Ny)    + fd%Ez_r1(0,Ny-1) )
     fd%Ez_r1(Nx,Ny)       = 0.5_dp*( fd%Ez_r1(Nx-1,Ny) + fd%Ez_r1(Nx,Ny-1) )

     ! --- 3) Source (au centre) -> temps n+1 ---
     fd%Ez_r1(isc_x, isc_y) = fd%Ez_r1(isc_x, isc_y) + (Esrc(n))
     ! --- 4) Update Hx, Hy en utilisant Ez(n+1) ---
     do i = 0, Nx-1
        do j = 0, Ny-1
           fd%Hx_r1(i,j) = fd%Hx_r1(i,j) - (fd%c_H(i,j)/dy) * ( fd%Ez_r1(i,j+1) - fd%Ez_r1(i,j) )
           fd%Hy_r1(i,j) = fd%Hy_r1(i,j) + (fd%c_H(i,j)/dx) * ( fd%Ez_r1(i+1,j) - fd%Ez_r1(i,j) )
        end do
     end do

     ! Ecriture observation (4e point sûr)
     if (ix_fix>=0 .and. ix_fix<=Nx .and. iy_fix>=0 .and. iy_fix<=Ny) then
        Ez_fix = fd%Ez_r1(ix_fix,iy_fix)
     else
        Ez_fix = 0.0_dp
     end if
     write(16,*) n*dt, fd%Ez_r1(xobs(1),yobs(1)), fd%Ez_r1(xobs(2),yobs(2)), Ez_fix
  end do

  call cpu_time(t_end)
  print *, 'Temps CPU (reference1) :', t_end - t_start, ' secondes'
  close(16)

end subroutine reference1


!======================================================================
! Subroutine mise à jour des champs - Domaine raffiné (r2)
!======================================================================
subroutine reference2(fd, Nx_r, Ny_r, Nt, dx_r, dy_r, dt_r, Esrc)
  use lesconstantes_numeriques
    use cpml_2d
  implicit none
  class(tableau), intent(inout) :: fd
  integer, intent(in)           :: Nx_r, Ny_r, Nt
  real(dp), intent(in)          :: dx_r, dy_r, dt_r
  real(dp), intent(in)          :: Esrc(0:Nt-1)

  integer  :: i, j, n, ics_x, ics_y
  real(dp) :: alpha_x, alpha_y, S
  real(dp), allocatable, save :: Ez_prev(:,:)
  real(dp) :: t_start, t_end
  real(dp) :: cE_r, cHx_r, cHy_r
  integer, parameter :: nobs = 3
  integer :: xobs(nobs), yobs(nobs)
  integer :: ix_fix, iy_fix
  real(dp) :: Ez_fix

  ! (Re)allocation robuste d'Ez_prev
  if (.not. allocated(Ez_prev)) then
     allocate(Ez_prev(0:Nx_r,0:Ny_r), source=0.0_dp)
  else if ( size(Ez_prev,1) /= Nx_r+1 .or. size(Ez_prev,2) /= Ny_r+1 ) then
     deallocate(Ez_prev)
     allocate(Ez_prev(0:Nx_r,0:Ny_r), source=0.0_dp)
  end if

  ! Points d'observation
  xobs = [ Nx_r/2, 3*Nx_r/4, Nx_r/2 ]
  yobs = [ Ny_r/2, 3*Ny_r/4, Ny_r/2 ]

  ! Un 4e point "fixe" optionnel (protégé)
  ix_fix = 1200
  iy_fix = 1200

  ics_x = Nx_r/2
  ics_y = Ny_r/2

  ! ABC Mur 1er ordre : coefficients distincts
  alpha_x = (c*dt_r - dx_r) / (c*dt_r + dx_r)
  alpha_y = (c*dt_r - dy_r) / (c*dt_r + dy_r)

  ! Coeffs (vide). Si milieux, remplace comme en r1 avec fd%c_E / fd%c_H spécifiques.
  cE_r  = dt_r / epsilon_0
  cHx_r = dt_r / (mu_0 * dy_r)
  cHy_r = dt_r / (mu_0 * dx_r)

  ! Courant
  S = c*dt_r*sqrt( (1.0_dp/dx_r**2) + (1.0_dp/dy_r**2) )
  if (S >= 1.0_dp) then
     print *, 'Avertissement: CFL (r2) >= 1 (S = ', S, ') — instable en 2D'
  end if

  call execute_command_line('mkdir -p observation', wait=.true.)
  open(unit=19, file='observation/Ez_r2.dat', status='replace', action='write')

  call cpu_time(t_start)

  do n = 0, Nt - 1
     ! --- sauvegarde Ez(n) pour l'ABC ---
     Ez_prev = fd%Ez_r2

     ! --- 1) Update Ez intérieur -> temps n+1 ---
     do i = 1, Nx_r-1
        do j = 1, Ny_r-1
           fd%Ez_r2(i,j) = fd%Ez_r2(i,j) + cE_r * ( (fd%Hy_r2(i,j) - fd%Hy_r2(i-1,j))/dx_r &
                                                  - (fd%Hx_r2(i,j) - fd%Hx_r2(i,j-1))/dy_r )
        end do
     end do

     ! --- 2) ABC de Mur (utilise Ez_prev = temps n) ---
     ! Bords gauche/droit : j = 1..Ny_r-1
     fd%Ez_r2(0, 1:Ny_r-1)    = Ez_prev(1, 1:Ny_r-1)     + alpha_x * ( fd%Ez_r2(1, 1:Ny_r-1)     - Ez_prev(0, 1:Ny_r-1) )
     fd%Ez_r2(Nx_r,1:Ny_r-1)  = Ez_prev(Nx_r-1,1:Ny_r-1) + alpha_x * ( fd%Ez_r2(Nx_r-1,1:Ny_r-1) - Ez_prev(Nx_r,1:Ny_r-1) )

     ! Bords bas/haut : i = 1..Nx_r-1
     fd%Ez_r2(1:Nx_r-1,0)     = Ez_prev(1:Nx_r-1,1)      + alpha_y * ( fd%Ez_r2(1:Nx_r-1,1)      - Ez_prev(1:Nx_r-1,0) )
     fd%Ez_r2(1:Nx_r-1,Ny_r)  = Ez_prev(1:Nx_r-1,Ny_r-1) + alpha_y * ( fd%Ez_r2(1:Nx_r-1,Ny_r-1) - Ez_prev(1:Nx_r-1,Ny_r) )

     ! Coins (moyenne simple)
     fd%Ez_r2(0,0)            = 0.5_dp*( fd%Ez_r2(1,0)        + fd%Ez_r2(0,1) )
     fd%Ez_r2(Nx_r,0)         = 0.5_dp*( fd%Ez_r2(Nx_r-1,0)   + fd%Ez_r2(Nx_r,1) )
     fd%Ez_r2(0,Ny_r)         = 0.5_dp*( fd%Ez_r2(1,Ny_r)     + fd%Ez_r2(0,Ny_r-1) )
     fd%Ez_r2(Nx_r,Ny_r)      = 0.5_dp*( fd%Ez_r2(Nx_r-1,Ny_r)+ fd%Ez_r2(Nx_r,Ny_r-1) )

     ! --- 3) Source (au centre) -> temps n+1 ---
     ! (gain 3x conservé comme dans ton code)
     fd%Ez_r2(ics_x, ics_y) = fd%Ez_r2(ics_x, ics_y) + (dt_r/epsilon_0) * (Esrc(n)*3.0_dp) / (dx_r*dy_r)

     ! --- 4) Mise a jour de Hx en utilisant Ez(n+1) ---
     do i = 0, Nx_r-1
        do j = 0, Ny_r-1
           fd%Hx_r2(i,j) = fd%Hx_r2(i,j) - cHx_r * ( fd%Ez_r2(i,j+1) - fd%Ez_r2(i,j) )
           end do
           end do

         !Mise ajour de Hy en utilisant Ez(n+1) 
        do i = 0, Nx_r-1    
        do j = 0, Ny_r-1
           fd%Hy_r2(i,j) = fd%Hy_r2(i,j) + cHy_r * ( fd%Ez_r2(i+1,j) - fd%Ez_r2(i,j) )
        end do
     end do

     ! Ecriture observation (4e point sûr)
     if (ix_fix>=0 .and. ix_fix<=Nx_r .and. iy_fix>=0 .and. iy_fix<=Ny_r) then
        Ez_fix = fd%Ez_r2(ix_fix,iy_fix)
     else
        Ez_fix = 0.0_dp
     end if
     write(19,*) n*dt_r, fd%Ez_r2(xobs(1),yobs(1)), fd%Ez_r2(xobs(2),yobs(2)), Ez_fix

  end do

  call cpu_time(t_end)
  print *, 'Temps CPU (reference2) :', t_end - t_start, ' secondes'
  close(19)

end subroutine reference2

end module lesreference
