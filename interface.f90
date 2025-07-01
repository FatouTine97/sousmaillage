 subroutine interface_mag_coarse_fine()
    ! Couplage des champs Hz (coarse) -> ez (fine) via interpolation symétrique
    integer :: i_f, j_f
    integer :: i_c, j_c
    real :: Hz2, Hz3, Hz2c, Hz3c, Hauxz
    real :: LCC1, LCC2, LIC1, LIC2, fHz

    do j_f = 1, size(ez,2)-2
      do i_f = 1, size(ez,1)-2
        ! Identifie les indices grossiers correspondants
        i_c = sx + i_f / r
        j_c = sy + j_f / r

        ! Coefficients de longueur et d'interpolation pour r = 3
        LCC1 = 1.0 / 3.0
        LCC2 = 1.0 / 3.0
        LIC1 = 2.0 / 3.0
        LIC2 = 1.0 / 3.0

        ! Champs Hz à interpoler depuis le maillage grossier
        Hz2 = Hx(i_c, j_c)   ! Ce sont des identifications à affiner selon orientation
        Hz3 = Hx(i_c+1, j_c)

        Hz2c = LCC1 * Hz2
        Hz3c = LCC2 * Hz3

        Hauxz = Hz2c * LIC1 + Hz3c * LIC2

        ! Injection dans la mise à jour du champ ez en bordure
        ! ez(i_f,j_f) = ... + fHz * Hauxz  (intégré dans maj_e_fine en réalité)

      end do
    end do
  end subroutine

  subroutine interface_elec_fine_coarse()
    ! Couplage des champs ez (fine) -> Hz (coarse)
    integer :: i_c, j_c
    integer :: i_f, j_f
    real :: eAuxx
    real :: LIC(3), LCC(3)
    real :: sum_ez
    real :: fex

    do j_c = sx, ex
      do i_c = sy, ey
        ! Exemple : moyenne pondérée de 3x3 valeurs ez autour du centre
        sum_ez = 0.0
        LIC = (/1.0/3.0, 1.0/3.0, 1.0/3.0/)
        LCC = (/1.0/3.0, 1.0/3.0, 1.0/3.0/)

        do j_f = -1,1
          do i_f = -1,1
            sum_ez = sum_ez + ez(r*(i_c-sx)+i_f+1, r*(j_c-sy)+j_f+1) * LIC(i_f+2) * LCC(j_f+2)
          end do
        end do

        eAuxx = sum_ez
        fex = 1.0 / 3.0

        Hx(i_c,j_c) = Hx(i_c,j_c) + fex * eAuxx

      end do
    end do