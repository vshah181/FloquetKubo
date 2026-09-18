module statistical_distributions
use kinds, only: dp
implicit none
private
public :: fermi_dirac
contains
    pure function fermi_dirac(energy, temperature) result(fd)
    use constants, only: kB=>boltzmann_constant_ev
    implicit none
        real(dp), intent(in)           :: energy ! relative to fermi level
        real(dp), intent(in), optional :: temperature

        real(dp)                       :: fd

        if (present(temperature) .and. (temperature .gt. 0.0_dp)) then
            fd = (1.0_dp - tanh(energy / (kB * temperature * 2.0_dp))) / 2.0_dp
        else
            fd = 0.5_dp - sign(0.5_dp, energy)  ! branchless step-function
        endif
    end function fermi_dirac
end module statistical_distributions
