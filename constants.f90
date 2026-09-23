module constants
use kinds, only: dp
implicit none
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
    real(dp), parameter :: two_pi = 8.0_dp * atan(1.0_dp)
    real(dp), parameter :: deg2rad = atan(1.0_dp) / 45.0_dp
    real(dp), parameter :: rad2deg = 45.0_dp / atan(1.0_dp)
    real(dp), parameter :: root_two = sqrt(2.0_dp)
    real(dp), parameter :: root_half = sqrt(0.5_dp)
    real(dp), parameter :: planck_constant = 6.62607015E-34_dp
    real(dp), parameter :: elementary_charge = 1.602176634E-19_dp
    real(dp), parameter :: speed_of_light = 299792458.0_dp
    real(dp), parameter :: boltzmann_constant = 1.380649E-23_dp
    real(dp), parameter :: vacuum_permittivity = 8.8541878188E-12_dp
    real(dp), parameter :: planck_constant_ev = planck_constant / elementary_charge
    real(dp), parameter :: reduced_planck_constant = planck_constant / two_pi
    real(dp), parameter :: reduced_planck_constant_ev = planck_constant_ev / two_pi
    real(dp), parameter :: boltzmann_constant_ev = boltzmann_constant / elementary_charge

    complex(dp), parameter :: cmplx_i = cmplx(0.0_dp, 1.0_dp, kind=dp)
    complex(dp), parameter :: cmplx_0 = cmplx(0.0_dp, 0.0_dp, kind=dp)
    complex(dp), parameter :: cmplx_1 = cmplx(1.0_dp, 0.0_dp, kind=dp)
end module constants

