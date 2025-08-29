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
  real(dp):: dt_prime

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

  ! Coefficients de mise à jour (vide)
  allocate(this%c_H(0:Nx, 0:Ny), this%c_E(0:Nx, 0:Ny))
  this%c_E = dt / EPSILON_0
  this%c_H = dt / MU_0
end subroutine initialiser_champs

!=======================================================================
!   MISE À JOUR DES CHAMPS (CG + FG, ABC Mur corrigées)
!=======================================================================
subroutine mise_a_jour_champs(fd, Nx, Ny, Nt, dx, dt, dy, Nx_sm, Ny_sm, Esrc)
  use lesconstantes_numeriques
  implicit none
  class(tableau), intent(inout) :: fd
  integer, intent(in) :: Nx, Ny, Nt, Nx_sm, Ny_sm
  real(dp), intent(in) :: dx, dy, dt
  real(dp), intent(in) :: Esrc(0:Nt-1)

  !------------------ Réglages ------------------
  integer, parameter :: r = 3
  logical, parameter :: enable_fg = .false.     ! <== IMPORTANT: FG coupé pour Mur interne
  real(dp), parameter :: A_src = 3.0_dp         ! amplitude de la source "soft"
  integer, parameter :: snap_every_coarse = 100
  integer, parameter :: nobs = 3
  !----------------------------------------------

  integer :: i, j,jj, jr, ig, i0, j0, n
  real(dp) :: dx_sm, dy_sm, dt_prime, dx_prime, dy_prime
  real(dp) :: alpha_x, alpha_y
  integer :: Nx_eff
  logical :: do_snapshot
  real(dp) :: t1, t2
  integer :: xobs(nobs), yobs(nobs)
  real(dp), allocatable, save :: Ez_prev(:,:)
  real(dp), allocatable :: Haux_y(:)
  real  :: dist_interf
  ! vars FG (laissées mais inutilisées si enable_fg=.false.)
  real(dp) :: e1, e2, e3, e4, e5, eaux_z, Hz2, Hz3
  integer :: iPobs(2,Npt)

  ! ===== init divers =====
  if (.not. allocated(Ez_prev)) then
     allocate(Ez_prev(0:Nx,0:Ny), source=0.0_dp)
  else if ( size(Ez_prev,1) /= Nx+1 .or. size(Ez_prev,2) /= Ny+1 ) then
     deallocate(Ez_prev)
     allocate(Ez_prev(0:Nx,0:Ny), source=0.0_dp)
  end if

  i0 = Nx/2;  j0 = Ny/2

  dx_sm   = dx / r
  dy_sm   = dy / r
  dt_prime = dt / r
  dx_prime = 0.5_dp*dx + 0.5_dp*dx_sm
  dy_prime = 0.5_dp*dy + 0.5_dp*dy_sm

  ! ABC Mur 1er ordre
  alpha_x = (c*dt_prime  - dx) / (c*dt_prime + dx)
  alpha_y = (c*dt_prime - dy) / (c*dt_prime + dy)

  ! -------- Frontière interne --------
  if (i1 <= 0 .or. i1 >= Nx) error stop "i1 doit vérifier 0 < i1 < Nx."
  if (.not. allocated(fd%Ez)) error stop "Champs non initialisés."
  Nx_eff = i1   ! domaine calculé: 0..i1 (droite = frontière interne)

  

  do i=1,Npt
    iPobs(:,i) = INT(Pobs(:,i)/(/dx,dy/));
    if (iPobs(1,i) >= i1) then
      dist_interf = (Pobs(1,i)-i1*dx)
      iPobs(1,i)  = dist_interf /dx_sm  
      iPobs(2,i)  =  Pobs(2,i)/dy_sm
    end if
    print '("Point obs maillage mixte en cellules n°",2I5)', iPobs(:,i)
  end do

  allocate(Haux_y(0:Ny_sm-1))

  !call execute_command_line('mkdir -p observation', wait=.true.)
  open(10,file="observation/Ez_t.dat",    status='replace', action='write')
  open(11,file="observation/Hx_t.dat",    status='replace', action='write')
  open(12,file="observation/Hy_t.dat",    status='replace', action='write')
  open(13,file="observation/carto_t.dat", status='replace', action='write')
  open(14,file="observation/carto_t1.dat",status='replace', action='write')
  open(15,file="observation/ez_s.dat",    status='replace', action='write')

  call cpu_time(t1)

  do n = 0, Nt-1
     if (mod(n,100)==0) print *, "Itération n°", n

     ! ===== sauvegarde Ez^n pour ABC =====
     Ez_prev = fd%Ez

     !============================================================
     ! 1) Update Ez intérieur (sans i=0, i=Nx_eff, j=0, j=Ny)
     !============================================================
     do i = 1, Nx_eff-1
        do j = 1, Ny-1
           fd%Ez(i,j) = fd%Ez(i,j) + (dt_prime/epsilon_0)  * ( (fd%Hy(i,j) - fd%Hy(i-1,j))/dx &
                                                   - (fd%Hx(i,j) - fd%Hx(i,  j-1))/dy )
        end do
     end do

     !============================================================
     ! 2) ABC Mur (gauche, droite interne, bas, haut) au même temps
     !============================================================
     ! gauche (i=0)
     fd%Ez(0,      1:Ny-1) = Ez_prev(1,       1:Ny-1) + alpha_x * ( fd%Ez(1,       1:Ny-1) - Ez_prev(0,      1:Ny-1) )
     ! droite interne (i=i1 = Nx_eff)
     fd%Ez(Nx_eff, 1:Ny-1) = Ez_prev(Nx_eff-1,1:Ny-1) + alpha_x * ( fd%Ez(Nx_eff-1,1:Ny-1) - Ez_prev(Nx_eff, 1:Ny-1) )
     ! bas/haut
     fd%Ez(1:Nx_eff-1, 0)  = Ez_prev(1:Nx_eff-1, 1)    + alpha_y * ( fd%Ez(1:Nx_eff-1, 1)    - Ez_prev(1:Nx_eff-1, 0) )
     fd%Ez(1:Nx_eff-1, Ny) = Ez_prev(1:Nx_eff-1, Ny-1) + alpha_y * ( fd%Ez(1:Nx_eff-1, Ny-1) - Ez_prev(1:Nx_eff-1, Ny) )
     ! coins
     fd%Ez(0,      0 ) = 0.5_dp*( fd%Ez(1,      0) + fd%Ez(0,      1) )
     fd%Ez(Nx_eff, 0 ) = 0.5_dp*( fd%Ez(Nx_eff-1,0) + fd%Ez(Nx_eff, 1) )
     fd%Ez(0,      Ny) = 0.5_dp*( fd%Ez(1,      Ny)+ fd%Ez(0,      Ny-1) )
     fd%Ez(Nx_eff, Ny) = 0.5_dp*( fd%Ez(Nx_eff-1,Ny)+fd%Ez(Nx_eff, Ny-1) )

     !============================================================
     ! 3) Source "soft" (évite la sur-injection)
     !============================================================
     fd%Ez(i0,j0) = fd%Ez(i0,j0) +  (dt_prime/epsilon_0) *Esrc(n) /(dx*dy) ! Esrc = Courant

     !============================================================
     ! 4) Update Hx, Hy (bornes CORRECTES près de i1)
     !============================================================
     ! Hx : i = 0 .. Nx_eff (on lit Ez(i, j+1)-Ez(i,j))
     do i = 1, i1-1
        do j = 1, Ny-1
           fd%Hx(i,j) = fd%Hx(i,j) - ((dt_prime/mu_0) /dy) * ( fd%Ez(i,j+1) - fd%Ez(i,j) )
        end do
     end do
     ! Hy : i = 0 .. Nx_eff-1 (on lit Ez(i+1,j)-Ez(i,j), PAS AU-DELÀ de i1)
     do i = 1, i1-1
        do j = 1, Ny-1
           fd%Hy(i,j) = fd%Hy(i,j) + ((dt_prime/mu_0) /dx) * ( fd%Ez(i+1,j) - fd%Ez(i,j) )
        end do
     end do

     !============================================================
     ! 5 & 6) Couplage / Mise à jour FG  -- **désactivés ici**
     !============================================================
    ! if (enable_fg) then
        do j = 3, Ny_sm-3,3
           jr = int(j/r)
           e1 = fd%ez_s(0,j-2); e2 = fd%ez_s(0,j-1); e3 = fd%ez_s(0,j)
           e4 = fd%ez_s(0,j+1); e5 = fd%ez_s(0,j+2)
           eaux_z = ((1.0_dp/3.0_dp)*e1 + (2.0_dp/3.0_dp)*e2 + e3 + (2.0_dp/3.0_dp)*e4 + (1.0_dp/3.0_dp)*e5)/3.0_dp
           fd%Hy(i1-1,jr) = fd%Hy(i1-1,jr) + ((dt_prime/mu_0) /dx) * (eaux_z - fd%Ez(i1-1,jr))
        end do

        do jj = 1, Ny_sm-1,3
           ig = int(jj/r)
           Hz2 = fd%Hy(i1-1, ig)
           Hz3 = fd%Hy(i1-1, ig+1)
           j=jj
           Haux_y(j) = (2.0_dp/3.0_dp)*Hz2 + (1.0_dp/3.0_dp)*Hz3
           fd%ez_s(0,j) = fd%ez_s(0,j) + (dt_prime/epsilon_0) * ( (fd%hy_s(0,j) - Haux_y(j))/dx_prime &
                                                                - (fd%hx_s(0,j) - fd%hx_s(0,j-1))/dy_sm )
            j=jj+1                                                                
           Haux_y(j) = (1.0_dp/3.0_dp)*Hz2 + (2.0_dp/3.0_dp)*Hz3
           fd%ez_s(0,j) = fd%ez_s(0,j) + (dt_prime/epsilon_0) * ( (fd%hy_s(0,j) - Haux_y(j))/dx_prime &
                                                               - (fd%hx_s(0,j) - fd%hx_s(0,j-1))/dy_sm )
            j=jj+2
           Haux_y(j) = Hz3
           fd%ez_s(0,j) = fd%ez_s(0,j) + (dt_prime/epsilon_0) * ( (fd%hy_s(0,j) - Haux_y(j))/dx_prime &
                                                                - (fd%hx_s(0,j) - fd%hx_s(0,j-1))/dy_sm )
        end do

        do i = 1, Nx_sm-1
           do j = 1, Ny_sm-1
              fd%ez_s(i,j) = fd%ez_s(i,j) + (dt_prime/epsilon_0) * ( (fd%hy_s(i,j) - fd%hy_s(i-1,j))/dx_sm &
                                                                   - (fd%hx_s(i,j) - fd%hx_s(i,j-1))/dy_sm )
           end do
        end do

        do i = 0, Nx_sm-1
           do j = 0, Ny_sm-1
              fd%hx_s(i,j) = fd%hx_s(i,j) - (dt_prime/mu_0) * (fd%ez_s(i,j+1) - fd%ez_s(i,j)) / dy_sm
              end do
              end do

        do i = 0, Nx_sm-1
        do j = 0, Ny_sm-1
              fd%hy_s(i,j) = fd%hy_s(i,j) + (dt_prime/mu_0) * (fd%ez_s(i+1,j) - fd%ez_s(i,j)) / dx_sm
           end do
        end do
     !end if

     !============================================================
     ! 7) Logs / snapshots
     !============================================================
   
     do i=1,Npt-2
       EPtobs(i) = fd%Ez(iPobs(1,i),iPobs(2,i))
     end do
     EPtobs(Npt-1) = fd%ez_s(iPobs(1,Npt-1),iPobs(2,Npt-1))
     EPtobs(Npt)   = fd%ez_s(iPobs(1,Npt),  iPobs(2,Npt))

     write(10,*) n*dt_prime, EPtobs(1:Npt) ;     

     do_snapshot = (mod(n, r*snap_every_coarse) == 0)
     if (do_snapshot) then
        write(14,'(A,F12.6)') '# t = ', real(n,dp)*dt_prime
        do i = 0, Nx, 4
           write(14,*) (fd%Ez(i,j), j=0, Ny, 4)
        end do
        write(14,*); call flush(14)

        write(13,'(A,F12.6)') '# t = ', real(n,dp)*dt_prime
        do i = 0, Nx_sm, 4
           write(13,*) (fd%ez_s(i,j), j=0, Ny_sm, 4)
        end do
        write(13,*); call flush(13)
     end if

  end do

  close(10); close(11); close(12); close(13); close(14); close(15)
  deallocate(Haux_y)

  call cpu_time(t2)
  print *, "Temps (mise à jour champs) =", t2 - t1, " seconde"
end subroutine mise_a_jour_champs


!=======================================================================
!   Moyenne conservatrice 3x3 (utilitaire)
!=======================================================================
real(dp) function compute_ez_aux(fd, i_f, j_f) result(ez_aux)
  use lesconstantes_numeriques
  implicit none
  type(tableau), intent(in) :: fd
  integer, intent(in)       :: i_f, j_f
  integer :: ii, jj
  real(dp) :: sum
  real(dp), dimension(-1:1, -1:1) :: w

  w = reshape([ &
      1.0_dp/16.0_dp, 2.0_dp/16.0_dp, 1.0_dp/16.0_dp, &
      2.0_dp/16.0_dp, 4.0_dp/16.0_dp, 2.0_dp/16.0_dp, &
      1.0_dp/16.0_dp, 2.0_dp/16.0_dp, 1.0_dp/16.0_dp], [3,3])

  sum = 0.0_dp
  do ii = -1, 1
     do jj = -1, 1
        if ( i_f+ii >= lbound(fd%ez_s,1) .and. i_f+ii <= ubound(fd%ez_s,1) .and. &
             j_f+jj >= lbound(fd%ez_s,2) .and. j_f+jj <= ubound(fd%ez_s,2) ) then
           sum = sum + w(ii+2, jj+2) * fd%ez_s(i_f+ii, j_f+jj)
        end if
     end do
  end do
  ez_aux = sum
end function compute_ez_aux

end module lesfonction

 