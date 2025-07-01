module lesfonction
    use lesconstantes_numeriques
    implicit none

    type :: tableau
        real(dp), dimension(:,:), allocatable :: Hx, Ez, Hy
        real(dp), dimension(:,:), allocatable :: hx_s, ez_s, hy_s
        real(dp), dimension(:,:), allocatable :: c_E, c_H
        real(dp), dimension(:,:), allocatable :: sigma_x, sigma_y
    contains
        procedure :: initialiser_champs
        procedure :: mise_a_jour_champs
        procedure :: initialiser_PML
    end type tableau

contains

    subroutine initialiser_champs(this, Nx, Ny, Nx_sm, Ny_sm, dt)
        class(tableau), intent(inout) :: this
        integer, intent(in) :: Nx, Ny, Nx_sm, Ny_sm
        real(dp), intent(in) :: dt

        allocate(this%Hx(0:Nx,0:Ny), this%Hy(0:Nx,0:Ny), this%Ez(0:Nx,0:Ny))
        allocate(this%hx_s(0:Nx_sm,0:Ny_sm), this%hy_s(0:Nx_sm,0:Ny_sm), this%ez_s(0:Nx_sm,0:Ny_sm))
        allocate(this%c_E(0:Nx,0:Ny), this%c_H(0:Nx,0:Ny))
        allocate(this%sigma_x(0:Nx,0:Ny), this%sigma_y(0:Nx,0:Ny))

        this%Ez = 0.0_dp; this%Hx = 0.0_dp; this%Hy = 0.0_dp
        this%ez_s = 0.0_dp; this%hx_s = 0.0_dp; this%hy_s = 0.0_dp
        this%c_E = dt / epsilon_0
        this%c_H = dt / mu_0
        this%sigma_x = 0.0_dp
        this%sigma_y = 0.0_dp
    end subroutine initialiser_champs

    subroutine initialiser_PML(this, Nx, Ny, npml, sigma_max)
        class(tableau), intent(inout) :: this
        integer, intent(in) :: Nx, Ny, npml
        real(dp), intent(in) :: sigma_max
        integer :: i, j
        real(dp) :: x_rel

        ! Profil sigma_x
        do i = 0, Nx
            x_rel = 0.0_dp
            if (i < npml) then
                x_rel = real(npml - i, dp) / npml
            else if (i > Nx - npml) then
                x_rel = real(i - (Nx - npml), dp) / npml
            end if
            do j = 0, Ny
                this%sigma_x(i,j) = sigma_max * x_rel**3
            end do
        end do

        ! Profil sigma_y
        do j = 0, Ny
            x_rel = 0.0_dp
            if (j < npml) then
                x_rel = real(npml - j, dp) / npml
            else if (j > Ny - npml) then
                x_rel = real(j - (Ny - npml), dp) / npml
            end if
            do i = 0, Nx
                this%sigma_y(i,j) = sigma_max * x_rel**3
            end do
        end do
    end subroutine initialiser_PML

    subroutine mise_a_jour_champs(this, Nx, Ny, Nt, dx, dt, dy, Nx_sm, Ny_sm, Esrc)
        class(tableau), intent(inout) :: this
        integer, intent(in) :: Nx, Ny, Nt, Nx_sm, Ny_sm
        real(dp), intent(in) :: dx, dy, dt
        real(dp), intent(in) :: Esrc(0:Nt-1)

        integer :: i, j, n
        real(dp) :: sigma_eff, cE, cH

        do n = 0, Nt - 1
            ! --- Mise à jour Ez avec PML ---
            do i = 1, Nx - 1
                do j = 1, Ny - 1
                    sigma_eff = this%sigma_x(i,j) + this%sigma_y(i,j)
                    cE = this%c_E(i,j) / (1.0_dp + dt * sigma_eff / (2.0_dp * epsilon_0))
                    this%Ez(i,j) = (1.0_dp - dt * sigma_eff / (2.0_dp * epsilon_0)) / &
                                   (1.0_dp + dt * sigma_eff / (2.0_dp * epsilon_0)) * this%Ez(i,j) + &
                                   cE * ((this%Hy(i,j) - this%Hy(i-1,j)) / dx - &
                                         (this%Hx(i,j) - this%Hx(i,j-1)) / dy)
                end do
            end do

            ! --- Source ---
       !d     this%Ez(Nx/2, Ny/2) = this%Ez(Nx/2, Ny/2) + Esrc(n)

            ! --- Conditions aux limites Dirichlet (optionnel) ---
            this%Ez(0,:) = 0.0_dp
            this%Ez(Nx,:) = 0.0_dp
            this%Ez(:,0) = 0.0_dp
            this%Ez(:,Ny) = 0.0_dp

            ! --- Mise à jour Hx, Hy avec PML ---
            do i = 0, Nx - 1
                do j = 0, Ny - 1
                    sigma_eff = this%sigma_x(i,j)
                    cH = this%c_H(i,j) / (1.0_dp + dt * sigma_eff / (2.0_dp * mu_0))
                    this%Hx(i,j) = (1.0_dp - dt * sigma_eff / (2.0_dp * mu_0)) / &
                                   (1.0_dp + dt * sigma_eff / (2.0_dp * mu_0)) * this%Hx(i,j) - &
                                   cH * (this%Ez(i,j+1) - this%Ez(i,j)) / dy

                    sigma_eff = this%sigma_y(i,j)
                    cH = this%c_H(i,j) / (1.0_dp + dt * sigma_eff / (2.0_dp * mu_0))
                    this%Hy(i,j) = (1.0_dp - dt * sigma_eff / (2.0_dp * mu_0)) / &
                                   (1.0_dp + dt * sigma_eff / (2.0_dp * mu_0)) * this%Hy(i,j) + &
                                   cH * (this%Ez(i+1,j) - this%Ez(i,j)) / dx
                end do
            end do


        end do
    end subroutine mise_a_jour_champs

end module lesfonction
