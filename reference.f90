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
  integer :: iPobs(2,Npt)
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

  !
print '("Point src gros maillage en cellules n°",2I5)', isc_x,isc_y

do i=1,Npt
  iPobs(:,i) = INT(Pobs(:,i)/(/dx,dy/));
  print '("Point obs gros maillage en cellules n°",2I5)', iPobs(:,i)
end do

print *, dt
  

 ! call execute_command_line('mkdir -p observation', wait=.true.)
  open(unit=16, file="observation/Ez_r1.dat", status='replace', action='write')

  call cpu_time(t_start)
  fd%c_E = dt / EPSILON_0
  fd%c_H = dt / MU_0
  print '("dt/eps0, dt/mu0 =",2E15.7)', fd%c_E(250,250),fd%c_H(250,250)
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
     fd%Ez_r1(isc_x, isc_y) = fd%Ez_r1(isc_x, isc_y) +  (dt/epsilon_0) *Esrc(n) /(dx*dy)

     ! --- 4) Update Hx, Hy en utilisant Ez(n+1) ---
     do i = 0, Nx-1
        do j = 0, Ny-1
           fd%Hx_r1(i,j) = fd%Hx_r1(i,j) - (fd%c_H(i,j)/dy) * ( fd%Ez_r1(i,j+1) - fd%Ez_r1(i,j) )
           end do
           end do 

    do i = 0, Nx-1
        do j = 0, Ny-1
           fd%Hy_r1(i,j) = fd%Hy_r1(i,j) + (fd%c_H(i,j)/dx) * ( fd%Ez_r1(i+1,j) - fd%Ez_r1(i,j) )
        end do
     end do

     
     do i=1,Npt
       EPtobs(i) = fd%Ez_r1(iPobs(1,i),iPobs(2,i))
     end do
     !write(16,*) n*dt, fd%Ez_r1(xobs(1),yobs(1)), fd%Ez_r1(xobs(2),yobs(2)), Ez_fix
     write(16,*) n*dt, EPtobs(1:Npt) ;
     
  end do

  call cpu_time(t_end)
  print *, 'Temps CPU (reference1) :', t_end - t_start, ' secondes'
  close(16)

end subroutine reference1


!======================================================================
! Subroutine mise à jour des champs - Domaine raffiné (r2)
!======================================================================
!subroutine reference2(fd, Nx_r, Ny_r, Nt, dx_r, dy_r, dt_r, Esrc)
 subroutine reference2(fd, Nx, Ny, Nt, dx, dt, dy, Esrc)
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
  integer :: iPobs(2,Npt)

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
     print *, 'Avertissement: CFL (r2) >= 1 (S = ', S, ') — instable en 2D'
  end if
print *, dt
  ! ABC Mur 1er ordre : coefficients distincts
  alpha_x = (c*dt - dx) / (c*dt + dx)
  alpha_y = (c*dt - dy) / (c*dt + dy)



print '("Point src petit maillage en cellules n°",2I5)', isc_x,isc_y

do i=1,Npt
  iPobs(:,i) = INT(Pobs(:,i)/(/dx,dy/));
  print '("Point obs petit maillage en cellules n°",2I5)', iPobs(:,i)
end do
  ! Un 4e point "fixe" optionnel (protégé)
  ix_fix = 750
  iy_fix = 750

 ! call execute_command_line('mkdir -p observation', wait=.true.)
  open(unit=19, file="observation/Ez_r2.dat", status='replace', action='write')

  call cpu_time(t_start)

  do n = 0, Nt - 1
     ! --- sauvegarde Ez(n) pour l'ABC (temps n) ---
     Ez_prev = fd%Ez_r2

     ! --- 1) Update Ez intérieur -> temps n+1 ---
     do i = 1, Nx-1
        do j = 1, Ny-1
           fd%Ez_r2(i,j) = fd%Ez_r2(i,j) + (dt/epsilon_0)* ( (fd%Hy_r2(i,j) - fd%Hy_r2(i-1,j))/dx &
                                                        - (fd%Hx_r2(i,j) - fd%Hx_r2(i,j-1))/dy )
        end do
     end do

     ! --- 2) ABC de Mur (utilise Ez_prev = temps n) ---
     ! Bords gauche/droit : j = 1..Ny-1
     fd%Ez_r2(0, 1:Ny-1)   = Ez_prev(1, 1:Ny-1)    + alpha_x * ( fd%Ez_r2(1, 1:Ny-1)    - Ez_prev(0, 1:Ny-1) )
     fd%Ez_r2(Nx,1:Ny-1)   = Ez_prev(Nx-1,1:Ny-1)  + alpha_x * ( fd%Ez_r2(Nx-1,1:Ny-1)  - Ez_prev(Nx,1:Ny-1) )

     ! Bords bas/haut : i = 1..Nx-1
     fd%Ez_r2(1:Nx-1,0)    = Ez_prev(1:Nx-1,1)     + alpha_y * ( fd%Ez_r2(1:Nx-1,1)     - Ez_prev(1:Nx-1,0) )
     fd%Ez_r2(1:Nx-1,Ny)   = Ez_prev(1:Nx-1,Ny-1)  + alpha_y * ( fd%Ez_r2(1:Nx-1,Ny-1)  - Ez_prev(1:Nx-1,Ny) )

     ! Coins (moyenne simple)
     fd%Ez_r2(0,0)         = 0.5_dp*( fd%Ez_r2(1,0)     + fd%Ez_r2(0,1) )
     fd%Ez_r2(Nx,0)        = 0.5_dp*( fd%Ez_r2(Nx-1,0)  + fd%Ez_r2(Nx,1) )
     fd%Ez_r2(0,Ny)        = 0.5_dp*( fd%Ez_r2(1,Ny)    + fd%Ez_r2(0,Ny-1) )
     fd%Ez_r2(Nx,Ny)       = 0.5_dp*( fd%Ez_r2(Nx-1,Ny) + fd%Ez_r2(Nx,Ny-1) )

     ! --- 3) Source (au centre) -> temps n+1 ---
     fd%Ez_r2(isc_x, isc_y) = fd%Ez_r2(isc_x, isc_y) +  (dt/epsilon_0) *Esrc(n) /(dx*dy)

     ! --- 4) Update Hx, Hy en utilisant Ez(n+1) ---
     do i = 1, Nx-1
        do j = 1, Ny-1
           fd%Hx_r2(i,j) = fd%Hx_r2(i,j) - ((dt/mu_0)/dy) * ( fd%Ez_r2(i,j+1) - fd%Ez_r2(i,j) )
           end do 
           end do

        do i = 1, Nx-1
        do j = 1, Ny-1
           fd%Hy_r2(i,j) = fd%Hy_r2(i,j) +  ((dt/mu_0)/dx) * ( fd%Ez_r2(i+1,j) - fd%Ez_r2(i,j) )
        end do
     end do

     
     do i=1,Npt
       EPtobs(i) = fd%Ez_r2(iPobs(1,i),iPobs(2,i))
     end do
     write(19,*) n*dt, EPtobs(1:Npt) ;
  end do

  call cpu_time(t_end)
  print *, 'Temps CPU (reference1) :', t_end - t_start, ' secondes'
  close(19)



end subroutine reference2


end module lesreference
