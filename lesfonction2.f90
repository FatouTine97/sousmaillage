module fdtd_d
  use constante_numerique
  implicit none

  type :: Fdtd2d  ! Correction du nom
  
    real(dp), dimension(:,:), allocatable :: Ez, Hx, Hy
    real(dp), dimension(:,:), allocatable :: c_E, c_H
  contains
    procedure, pass :: initial
    procedure, pass :: calcul
  end type Fdtd2d

contains

  ! Initialisation de la structure
  subroutine initial(this, Nx, Ny, dt)
    class(Fdtd2d), intent(inout) :: this
    integer, intent(in) :: Nx, Ny
    real(dp), intent(in) :: dt
    real(dp), parameter :: EPSILON_0 = 8.854187818d-12, MU_0 = 12.5663706144d-7
  allocate(this%c_E(0:Nx, 0:Ny))
  allocate(this%c_H(0:Nx, 0:Ny),this%Ez(0:Nx, 0:Ny), this%Hx(0:Nx, 0:Ny), this%Hy(0:Nx, 0:Ny))

   this%Ez = 0.0d0
    this%Hx = 0.0d0
    this%c_H = dt / MU_0 do i = 1, Nx - 1
      this%Ez(i,j) = this%Ey(i,j) - (this%c_E(i,j) 
   end do 
     print *, "c_E(250,250) = ", this%c_E(250,250)
    print *, "c_E(100,250) = ", this%c_E(100,250)
    print *, "c_H(250,250) = ", this%c_H(250,250)
    print *, "c_H(100,250) = ", this%c_H(100,250)
  end subroutine initial

  ! Calcul des champs
  subroutine calcul(this, Nx, Ny, Nt, Esrc, dx, dy, dt)
    class(Fdtd2d), intent(inout) :: this
    integer, intent(in) :: Nx, Ny, Nt
    real(dp), intent(in) :: dx, dy, dt
    real(dp), intent(in) :: Esrc(0:Nt-1)
    integer :: i, n, j, m
    real(dp), dimension(0 : 1       , 0 : Ny)                  :: Ezx0_n1 , Ezx0_n2                    
    real(dp), dimension(Nx - 1 : Nx , 0 : Ny)                  :: Ezx_n1 , Ezx_n2                 
    real(dp), dimension(0 : 1       , 0 : Nx)                  :: Ezy0_n1 , Ezy0_n2                     
    real(dp), dimension(Ny - 1 : Ny , 0 : Nx)                  :: Ezy_n1 , Ezy_n2
    real(dp), dimension(0 : 1       , 0 : Ny)                  :: Hx0_n1, Hx0_n2                    
    real(dp), dimension(Nx - 1 : Nx , 0 : Ny)                  :: Hx_n1, Hx_n2
    real(dp), dimension(0 : 1       , 0 : Nx)                  :: Hy0_n1, Hy0_n2                    
    real(dp), dimension(Ny - 1 : Ny , 0 : Nx)                  :: Hy_n1, Hy_n2
    real(dp)                                                   :: coef_mur1, coef_mur2, coef_mur3

    ! Ouverture des fichiers
    !open(unit=10, file="R1/Ez_t.dat", form="formatted")
    !open(unit=11, file="R1Hx_t.dat", form="formatted")
   ! open(unit=12, file="Hy_t.dat", form="formatted")
    !open(unit=13, file="carto_t.dat", form="formatted")

    m = 0
    Ezx0_n1 = 0.d0 ;    Ezx0_n2 = 0.0d0 
Ezx_n1  = 0.d0 ;    Ezx_n2  = 0.0d0
Ezy0_n1 = 0.d0 ;    Ezy0_n2 = 0.0d0
Ezy_n1  = 0.d0 ;    Ezy_n2  = 0.0d0
Hx0_n1  = 0.d0 ;    Hx0_n2  = 0.0d0
Hx_n1   = 0.d0 ;    Hx_n2   = 0.0d0
Hy0_n1  = 0.d0 ;    Hy0_n1   = 0.0d0
Hy_n1   = 0.d0 ;    Hy_n2   = 0.0d0

coef_mur1 = (c * dt - dx) / (c * dt + dx)
coef_mur2 = (2.d0 * dx) / (c * dt + dx)
coef_mur3 = (c * dt)**2 / ( 2 * dx * ( c * dt + dx ) )
    do n = 0, Nt-1
      if (mod(n, 10) == 0) print *, "Itération temporelle n°", n
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
      ! Mise à jour du champ électrique
     do i = 1, Nx - 1
      do j = 1, Ny - 1
        this%Ez(i, j) = this%Ez(i, j) + this%c_E(i, j) * ((this%Hy(i, j) - this%Hy(i-1, j)) / dx - &
           (this%Hx(i, j) - this%Hx(i, j-1)) / dy)
      end do
    end do
      !source 
       this%Ez(250, 250) = this%Ez(250, 250) + Esrc(n)

            !! Condition aux bord (Mur absorbant) x = 0
fd%Ez(0 , :)  =   - Ezx0_n2(1, :)                                             &
                  + coef_mur1 * (fd%Ez(1, :) + Ezx0_n2(0, :))                 &         ! n2 -> le temps n - 1    n1 -> n
                  + coef_mur2 * (Ezx0_n1(0, :) + Ezx0_n1(1, :))

            ! Condition aux bord (Mur absorbant) x = Nx
fd%Ez(Nx , :) =   - Ezx_n2(Nx - 1, :)                                          &
                  + coef_mur1 * (fd%Ez(Nx - 1, :) + Ezx_n2(Nx, :))             &
                  + coef_mur2 * (Ezx_n1(Nx, :) + Ezx_n1(Nx - 1, :))

            ! Condition aux bord (Mur absorbant) y = 0
fd%Ez(: , 0)  =   - Ezy0_n2(1, :)                                             &
                + coef_mur1 * (fd%Ez(:, 1) + Ezy0_n2(0, :))                   &
                + coef_mur2 * (Ezy0_n1(0, :) + Ezy0_n1(1, :))

            ! Condition aux bord (Mur absorbant) y = Ny
fd%Ez(: , Ny) =   - Ezy_n2(Ny - 1, :)                                         &
                  + coef_mur1 * (fd%Ez(:, Ny - 1) + Ezy_n2(Ny, :))            &
                  + coef_mur2 * (Ezy_n1(Ny, :) + Ezy_n1(Ny - 1,:))

             ! Source dans domaine grossier
fd%Ez(250, 250) = fd%Ez(250, 250) + Esrc(n)
       
      ! Mise à jour des champs magnétiques      
do i = 0, Nx - 1
   do j = 0, Ny - 1
          this%Hx(i, j) = this%Hx(i, j) - (this%c_H(i, j) / dy) *  (this%Ez(i, j+1) - this%Ez(i, j))
     end do
 end do
 
      do i = 0, Nx - 1
        do j = 0, Ny - 1
          this%Hy(i, j) = this%Hy(i, j) + (this%c_H(i, j) / dx) *  (this%Ez(i+1, j) - this%Ez(i, j))
        end do
      end do

     write(10,*) n*dt, this%Ez(250,250), this%Ez(150,250), this%Ez(50,250), this%Ez(1,250)
      write(11,*) n*dt, this%Hx(250,250), this%Hx(150,250), this%Hx(50,250), this%Hx(1,250)
      write(12,*) n*dt, this%Hy(250,250), this%Hy(150,250), this%Hy(50,250), this%Hy(1,250)


      if (mod(n, 100) == 0 .and. n > 0) then
        m = m + 1
        do i = 0, Nx, 2
          write(13,*) (this%Ez(i,j), j=0,Ny,2)
        end do
        write(13,*)
      end if
    end do
    print *, "Nombre de carto en temps = ", m
    close(10); close(11); close(12); close(13)

      
  end subroutine calcul

end module fdtd_d