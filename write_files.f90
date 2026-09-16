module write_files
use kinds, only: dp
implicit none
private
public :: write_energies, write_spin, write_oam, write_weights, write_kdists
public :: write_surface_projections
contains
    subroutine write_spin(spins)
    use parameters, only: nf_bands, nkp, seedname
    implicit none
        real(dp), intent(in) :: spins(nkp, nf_bands, 3)
        character(len=99) :: filename

        write(filename, fmt="(2A)") trim(adjustl(seedname)), "_spin.bin"
        open(210, file=trim(adjustl(filename)), form="unformatted",            &
            access="stream", status="replace")
        write(210) spins
        close(210)
    end subroutine write_spin
!******************************************************************************
    subroutine write_oam(angular_momenta)
    use parameters, only: nf_bands, nkp, seedname
    implicit none
        real(dp), intent(in) :: angular_momenta(nkp, nf_bands, 3)
        character(len=99) :: filename

        write(filename, fmt="(2A)") trim(adjustl(seedname)), "_oam.bin"
        open(211, file=trim(adjustl(filename)), form="unformatted",            &
            access="stream", status="replace")
        write(211) angular_momenta
        close(211)
    end subroutine write_oam
!******************************************************************************
    subroutine write_energies(energies)
    use parameters, only: nkp, nf_bands, seedname
    implicit none
        real(dp), intent(in) :: energies(nkp, nf_bands)
        character(len=99) :: filename

        write(filename, fmt="(2A)") trim(adjustl(seedname)), "_eigval.bin"
        open(212, file=trim(adjustl(filename)), form="unformatted",            &
            access="stream", status="replace")
        write(212) energies
        close(212)
    end subroutine write_energies
!******************************************************************************
    subroutine write_weights(spectral_weights)
    use parameters, only: nkp, nf_bands, seedname
    implicit none
        real(dp), intent(in) :: spectral_weights(nkp, nf_bands)
        character(len=99) :: filename

        write(filename, fmt="(2A)") trim(adjustl(seedname)), "_weights.bin"
        open(212, file=trim(adjustl(filename)), form="unformatted",            &
            access="stream", status="replace")
        write(212) spectral_weights
        close(212)
    end subroutine write_weights
!******************************************************************************
    subroutine write_kdists(kdists)
    use parameters, only: nkp, seedname
    implicit none
        real(dp), intent(in) :: kdists(nkp)
        character(len=99) :: filename

        write(filename, fmt="(2A)") trim(adjustl(seedname)), "_kdists.bin"
        open(212, file=trim(adjustl(filename)), form="unformatted",            &
            access="stream", status="replace")
        write(212) kdists
        close(212)
    end subroutine write_kdists
!******************************************************************************
    subroutine write_surface_projections(greens_function)
    use parameters, only: nlayers, nene, nkp, seedname
    use constants, only: pi
    implicit none
        complex(dp), intent(in) :: greens_function(nlayers, nene, nkp)
        real(dp),   allocatable :: spectral_function(:, :, :)
        character(len=99)       :: filename

        allocate(spectral_function(nlayers, nene, nkp))
        spectral_function = (-1.0_dp / pi) * aimag(greens_function)

        write(filename, fmt="(2A)") trim(adjustl(seedname)), "_layer_proj.bin"
        open(212, file=trim(adjustl(filename)), form="unformatted",            &             
            access="stream", status="replace")
        write(212) spectral_function
        close(212)
    end subroutine write_surface_projections
end module write_files
