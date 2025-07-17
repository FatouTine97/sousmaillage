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
    subroutine mise_a_jour_champs(fd, Nx, Ny, Nt, n,dx, dt, dy, Nx_sm, Ny_sm, Esrc)
        use lesconstantes_numeriques
        implicit none
    class(tableau), intent(inout)                              :: fd
    integer, intent(in)                                        :: Nx, Ny, Nt, Nx_sm, Ny_sm
    real(dp), intent(in)                                       :: dx, dy, dt
    real(dp), intent(in)                                       :: Esrc(0:Nt-1)
    integer, parameter                                         :: r = 3
    integer                                                    :: i, j, n
    real(dp)                                                   :: Hz2, Hz3, dx_sm, dy_sm, dt_prime
    real(dp)                                                   :: Haux_y, dx_prime, dy_prime, f_Hy, fez
    real(dp), dimension(0 : 1       , 0 : Ny)                  :: Ezx0_n1 , Ezx0_n2                    
    real(dp), dimension(Nx - 1 : Nx , 0 : Ny)                  :: Ezx_n1 , Ezx_n2                 
    real(dp), dimension(0 : 1       , 0 : Nx)                  :: Ezy0_n1 , Ezy0_n2                     
    real(dp), dimension(Ny - 1 : Ny , 0 : Nx)                  :: Ezy_n1 , Ezy_n2
   ! real(dp), dimension(0 : 1       , 0 : Ny)                  :: Hx0_n1, Hx0_n2                    
   ! real(dp), dimension(Nx - 1 : Nx , 0 : Ny)                  :: Hx_n1, Hx_n2
  !  real(dp), dimension(0 : 1       , 0 : Nx)                  :: Hy0_n1, Hy0_n2                    
  !  real(dp), dimension(Ny - 1 : Ny , 0 : Nx)                  :: Hy_n1, Hy_n2
    real(dp)                                                   :: coef_mur1, coef_mur2, coef_mur3
   real(dp), allocatable :: eaux_z(:)
   allocate(eaux_z(Ny_sm))

        dx_sm        = dx/ r
        dy_sm        = dy/ r
        dx_prime     = 0.5*dx + 0.5*dx_sm
        dy_prime     = 0.5*dy + 0.5*dy_sm
        f_Hy         = dy/dy_sm
        fez          = 1.0
        dt_prime     = dt/r



Ezx0_n1 = 0.d0 ;    Ezx0_n2 = 0.0d0 
Ezx_n1  = 0.d0 ;    Ezx_n2  = 0.0d0
Ezy0_n1 = 0.d0 ;    Ezy0_n2 = 0.0d0
Ezy_n1  = 0.d0 ;    Ezy_n2  = 0.0d0
!Hx0_n1  = 0.d0 ;    Hx0_n2  = 0.0d0
!Hx_n1   = 0.d0 ;    Hx_n2   = 0.0d0
!Hy0_n1  = 0.d0 ;    Hy0_n1   = 0.0d0
!Hy_n1   = 0.d0 ;    Hy_n2   = 0.0d0

coef_mur1 = (c * dt - dx) / (c * dt + dx)
coef_mur2 = (2.d0 * dx) / (c * dt + dx)
coef_mur3 = (c * dt)**2 / ( 2 * dx * ( c * dt + dx ) )


! Condition de bord x = 0
Ezx0_n2        = Ezx0_n1                        ! temps n - 2 pas i = 0 et i = 1
Ezx0_n1(0 , :) = fd%Ez(0,:)                     ! temps n - 1 pas i = 0
Ezx0_n1(1 , :) = fd%Ez(1,:)                     ! temps n - 1 pas i = 1

            ! Condition de bord x = Nx
Ezx_n2             = Ezx_n1                     ! temps n - 2 pas i = Nx - 1 et i = Nx
Ezx_n1(Nx - 1 , :) = fd%Ez(Nx - 1 , :)          ! temps n - 1 pas i = Nx - 1
Ezx_n1(Nx , :)     = fd%Ez(Nx , :)              ! temps n - 1 pas i = Nx

            ! Condition de bord y = 0
Ezy0_n2       = Ezy0_n1                         ! temps n - 2 pas j = 0 et j = 1
Ezy0_n1(0, :) = fd%Ez(:, 0)                     ! temps n - 1 pas j = 0
Ezy0_n1(1, :) = fd%Ez(:, 1)                     ! temps n - 1 pas j = 1

            ! Condition de bord y = Ny
Ezy_n2            = Ezy_n1                      ! temps n - 2 pas j = Ny - 1 et j = Ny
Ezy_n1(Ny - 1, :) = fd%Ez(:, Ny - 1)            ! temps n - 1 pas j = Ny - 1
Ezy_n1(Ny, :)     = fd%Ez(:, Ny)                ! temps n - 1 pas j = Ny

            
do i = 1, Nx - 1
 if (i < i1-1 ) then
do j = 1, Ny - 1
fd%Ez(i,j) = fd%Ez(i,j) + fd%c_E(i,j) * ((fd%Hy(i,j) - fd%Hy(i-1,j))/ dx  &
                        - (fd%Hx(i,j) - fd%Hx(i,j-1)) / dy)
!end if
end do
end if
end do

fd%Ez(0,:) = 0.0_dp
fd%Ez(Nx,:) = 0.0_dp  
fd%Ez(:,0) = 0.0_dp
fd%Ez(:,Ny) = 0.0_dp

             ! Source dans domaine grossier
fd%Ez(250,250) = fd%Ez(250, 250) + Esrc(n)

            ! Mise à jour Hx dans le domaine grossier       
do i = 0, Nx - 1
 if (i < i1-1 ) then
do j = 0, Ny - 1
     fd%Hx(i,j) = fd%Hx(i,j) - fd%c_H(i,j) / dy * (fd%Ez(i,j+1) - fd%Ez(i,j))

end do
end if
end do

            ! Mise à jour  Hy dans le domaine grossier
do i = 0, Nx - 1
 if (i < i1-1 ) then
do j = 0, Ny - 1
fd%Hy(i,j) = fd%Hy(i,j) + fd%c_H(i,j) / dx * (fd%Ez(i+1,j) - fd%Ez(i,j))
 end do
    end if
end do



do j = 0, Ny_sm - 1
   Hz2 = fd%Hy(i1-1, j)
   Hz3 = fd%Hy(i1-1, j-1)
   Haux_y = (1.0_dp/3.0_dp)*Hz2 + (2.0_dp/3.0_dp)*Hz3

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
   eaux_z(j) = compute_ez_aux(fd, 0, j)
end do


do j = 1, Ny_sm - 1
!do i = 1, i1 - 2
  fd%Hy(i1-1, j) = fd%Hy(i1-1, j) + (fd%c_H(i1-1,j)/dx) * (fez * eaux_z(j) - fd%Ez(i1-1, j))
!end do
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
            sum = sum + w(ii+2, jj+2) * fd%ez_s(i_f , j_f + jj)
        end do
    end do

    ez_aux = sum
end function compute_ez_aux



end module lesfonction  