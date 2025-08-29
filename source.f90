module lexcitation_donde
use lesconstantes_numeriques
implicit none 

contains 

subroutine excitation_donde(Esrc, Nt, dt)
    ! Subroutine to calculate the Gaussian excitation source
    real(dp), dimension(:), allocatable, intent(out) :: Esrc
    integer, intent(in) :: Nt
    real(dp), intent(in) :: dt
    real(dp) :: T, t0
    integer :: i
    integer :: unite
    integer,save :: icount=0
    character(LEN=20),dimension(2) :: fichier=(/'Esrc1.dat','Esrc2.dat'/)

   allocate(Esrc(0:Nt-1))
    ! Calculate the Gaussian parameters
    T = sqrt(log(10.0_dp)) / (pi * fmax)
    t0 = T * sqrt(log(1000.0_dp))

    ! Calculate the excitation source
    do i = 0, Nt-1
        Esrc(i) = exp(-((dt * i - t0) / T)**2)
    end do

 !stokage de l'excitation
    icount = icount+1;
    unite =10+icount;
    open(unit=unite, file=fichier(icount), status='replace')

    do i = 0, Nt-1
      write(unite, *) i*dt, Esrc(i)
    end do
    close(unite)
end subroutine excitation_donde
end module lexcitation_donde