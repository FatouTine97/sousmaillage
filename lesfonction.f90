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
       ! procedure :: reference1
    end type tableau

contains

    ! Sous-routine pour initialiser les champs
    subroutine initialiser_champs(this, Nx, Ny, Nx_sm, Ny_sm, Nx_r, Ny_r, dt)
        use lesconstantes_numeriques
        implicit none
        class(tableau), intent(inout) :: this
        integer, intent(in) :: Nx, Ny, Nx_sm, Ny_sm, Nx_r, Ny_r
        real(dp), intent(in) :: dt

        allocate(this%hx_s(0:Nx_sm, 0:Ny_sm), this%ez_s(0:Nx_sm, 0:Ny_sm), this%hy_s(0:Nx_sm, 0:Ny_sm))
        this%hx_s = 0.0_dp
        this%ez_s = 0.0_dp
        this%hy_s = 0.0_dp
!allocation des champs pour les gros maillage
       allocate( this%Ez(0:Nx, 0:Ny), this%Hx(0:Nx, 0:Ny), this%Hy(0:Nx, 0:Ny))
        this%Ez = 0.0_dp
        this%Hy = 0.0_dp
        this%Hx = 0.0_dp

!allocation des champs pour le reference 1
        allocate(this%Hx_r1(0:Nx, 0:Ny), this%Ez_r1(0:Nx, 0:Ny), this%Hy_r1(0:Nx, 0:Ny)) 
        this%Ez_r1 = 0.0_dp  
        this%Hy_r1 = 0.0_dp
        this%Hx_r1 = 0.0_dp

! allocation des champs pour le reference 2
        allocate(this%Ez_r2(0:Nx_r, 0:Ny_r), this%Hx_r2(0:Nx_r, 0:Ny_r), this%Hy_r2(0:Nx_r, 0:Ny_r))
        this%Ez_r2 = 0.0_dp
        this%Hy_r2 = 0.0_dp
        this%Hx_r2 = 0.0_dp

        !coefficients pour les champs
allocate(this%c_H(0:Nx, 0:Ny),  this%c_E(0:Nx, 0:Ny))
        this%c_E = dt / EPSILON_0
        this%c_H = dt / MU_0
        
    end subroutine initialiser_champs


    ! Subroutine pour la mise à jour des champs
    subroutine mise_a_jour_champs(fd, Nx, Ny, Nt, n, dx, dt, dy, Nx_sm, Ny_sm, Esrc)
        use lesconstantes_numeriques
        implicit none
    class(tableau), intent(inout)                              :: fd
    integer, intent(in)                                        :: Nx, Ny, Nt, Nx_sm, Ny_sm, n
    real(dp), intent(in)                                       :: dx, dy, dt
    real(dp), intent(in)                                       :: Esrc(0:Nt-1)
    integer, parameter                                         :: r = 3
    integer                                                    :: i, j,  i_0, i_1
    real(dp)                                                   :: Hz2, Hz3, dx_sm, dy_sm, dt_prime
    real(dp)                                                   :: Haux_y, dx_prime, dy_prime, f_Hy, fez
   real(dp) :: eaux_z
real(dp) :: coefMur
        dx_sm        = dx/ r
        dy_sm        = dy/ r
        dx_prime     = 0.5*dx + 0.5*dx_sm
        dy_prime     = 0.5*dy + 0.5*dy_sm
        f_Hy         = dy/dy_sm
        fez          = 1.0
        dt_prime     = dt/r


i_0 = Nx/2
i_1 = Ny/2

         coefMur = (c * dt - dx) / (c * dt + dx)
     
do i = 1, Nx - 1
 if (i < i1-1) then
do j = 1, Ny - 1
fd%Ez(i,j) = fd%Ez(i,j) + fd%c_E(i,j) * ((fd%Hy(i,j) - fd%Hy(i-1,j))/ dx  &
                        - (fd%Hx(i,j) - fd%Hx(i,j-1)) / dy)
!end if
end do
end if
end do

 
  ! Bord gauche i=0
  !do j = 1, Ny-1
    fd%Ez(0,:) = fd%Ez(1,j) + coefMur * (fd%Ez(1,j) - fd%Ez(0,j))
  !end do

  ! Bord droit i=Nx
 ! do j = 1, Ny-1
    fd%Ez(Nx,:) = fd%Ez(Nx-1,j) + coefMur * (fd%Ez(Nx-1,j) - fd%Ez(Nx,j))
 ! end do

  ! Bord bas j=0
  !do i = 1, Nx-1
    fd%Ez(:,0) = fd%Ez(i,1) + coefMur * (fd%Ez(i,1) - fd%Ez(i,0))
  !end do

  ! Bord haut j=Ny
  !do i = 1, Nx-1
    fd%Ez(:,Ny) = fd%Ez(i,Ny-1) + coefMur * (fd%Ez(i,Ny-1) - fd%Ez(i,Ny))
  !end do


            

             ! Source dans domaine grossier
fd%Ez(i_0,i_1) = fd%Ez(i_0, i_1) + Esrc(n)




            ! Mise à jour Hx dans le domaine grossier       
do i = 0, Nx - 1
 if (i < i1-1) then
do j = 0, Ny - 1
     fd%Hx(i,j) = fd%Hx(i,j) - fd%c_H(i,j) / dy * (fd%Ez(i,j+1) - fd%Ez(i,j))

end do
end if
end do

            ! Mise à jour  Hy dans le domaine grossier
do i = 0, Nx - 1
 if (i < i1-1) then
do j = 0, Ny - 1
fd%Hy(i,j) = fd%Hy(i,j) + fd%c_H(i,j) / dx * (fd%Ez(i+1,j) - fd%Ez(i,j))
 end do
    end if
end do


do j = 0, Ny_sm - 1
   Hz2 = fd%Hy(i1-1, j-1)
   Hz3 = fd%Hy(i1-1, j+1)
   Haux_y = (2.0_dp/3.0_dp)*Hz2 + (1.0_dp/3.0_dp)*Hz3

   fd%ez_s(0,j) = fd%ez_s(0,j) + (dt_prime/ epsilon_0) * ( &
        (fd%hy_s(0,j) - f_Hy * Haux_y)/dx_prime - &
        (fd%hx_s(0,j) - fd%hx_s(0,j-1))/dy_sm )
end do


! Ez raffiné
do i = 1, Nx_sm - 1
do j = 1, Ny_sm- 1
    fd%ez_s(i,j) = fd%ez_s(i,j) + (dt_prime/ epsilon_0 ) * ((fd%hy_s(i,j) - &
    fd%hy_s(i-1,j))/dx_sm - (fd%hx_s(i,j) - fd%hx_s(i,j-1))/dy_sm)
   end do
  end do
 
 

 !le champs magnetique(Hx) à l'intérieur du sous-maillage
do i = 0, Nx_sm - 1
  do j = 0, Ny_sm - 1
     fd%hx_s(i,j) = fd%hx_s(i,j) - (dt_prime/mu_0) * &
                    (fd%ez_s(i,j+1) - fd%ez_s(i,j)) / dy_sm
  end do
end do

do i = 0, Nx_sm - 1
  do j = 0, Ny_sm - 1
     fd%hy_s(i,j) = fd%hy_s(i,j) + (dt_prime/mu_0) * &
                    (fd%ez_s(i+1,j) - fd%ez_s(i,j)) / dx_sm
  end do
end do

do j = 1, Ny_sm - 1
  eaux_z = compute_ez_aux(fd, 0, j)
end do


do j = 0, Ny- 1
!(1.0_dp/3.0_dp) *(1.0* fd%ez_s(0, j-2)/3.0_dp + 2.0_dp * 2*fd%ez_s(0, j-1)/3.0_dp + &
      !       fd%ez_s(0, j) + 2.0_dp * fd%ez_s(0, j+1)/3.0_dp + 1.0_dp *fd%ez_s(0, j+2)/3.0_dp) 
 
 fd%Hy(i1-1, j) = fd%Hy(i1-1, j) + (fd%c_H(i1-1,j)/dx) * (fez * eaux_z - fd%Ez(i1-1, j))

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