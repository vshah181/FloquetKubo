module parameters
use kinds, only: dp
implicit none
private
public :: seedname, basis, nkp, num_bands, num_r_pts, avec, bvec, r_list, nk
public :: initialise_parameters, electric_field_si, omega, soc, nlayers, k_frac
public :: weights, energy_list, r_ham_list, num_photon, phase_shift, nf_bands
public :: a_0, projection_centres, broadening_factor, direction, k_shift, nene
public :: rlist_cart, fermi_energy, intersite_diffs, slab_area
    character(len=99), protected        :: seedname
    character(len=4),  protected        :: basis
    integer,           protected        :: num_bands, nf_bands, nkp, direction 
    integer,           protected        :: num_r_pts, num_photon, nlayers, nene
    integer,           protected        :: nk(3)
    real(dp),          protected        :: avec(3, 3), bvec(3, 3), phase_shift
    real(dp),          protected        :: energy_range(2), omega, slab_area
    real(dp),          protected        :: a_0, k_shift(3), k_frac(3)
    real(dp),          protected        :: electric_field_si, energy_step
    real(dp),          protected        :: broadening_factor, fermi_energy
    logical,           protected        :: soc

    integer,     protected, allocatable :: r_list(:, :), weights(:)
    real(dp),    protected, allocatable :: energy_list(:), rlist_cart(:, :)
    real(dp),    protected, allocatable :: projection_centres(:, :)
    real(dp),    protected, allocatable :: intersite_positions(:, :)
    real(dp),    protected, allocatable :: intersite_diffs(:, :, :)
    complex(dp), protected, allocatable :: r_ham_list(:, :, :)
contains
    subroutine initialise_parameters
    use array_utilities, only: arange, cross_product
    use read_files, only: read_input, read_hr, read_nnkp, read_kpoints,        &
        read_vector_potential
    use constants, only: hbar=>reduced_planck_constant_ev
    implicit none
        integer :: i, j, orb_i, orb_j

        call read_input(seedname, basis, soc, nlayers, energy_range,           &
            energy_step, broadening_factor, direction, fermi_energy)
        call read_hr(seedname, num_bands, num_r_pts, r_list, r_ham_list,       &
            weights)
        allocate(projection_centres(num_bands, 3))
        call read_nnkp(seedname, avec, bvec, num_bands, basis,                 &
            projection_centres)
        rlist_cart = make_rlist_cart(num_r_pts, r_list, avec)

        call read_kpoints(nk, k_shift, k_frac)
        nkp = product(nk)

        num_photon = 0
        omega = 0.0_dp
        phase_shift = 0.0_dp
        a_0 = 0.0_dp
        call read_vector_potential(num_photon, omega, phase_shift, a_0)
        nf_bands = (1 + (2 * num_photon)) * num_bands

        electric_field_si = hbar * omega * a_0 * 1.0E10_dp
        energy_list = arange(energy_range(1), energy_range(2), energy_step)
        nene = size(energy_list)

        allocate(intersite_positions(num_bands, 3))
        allocate(intersite_diffs(nf_bands, nf_bands, 2))

        ! Compute Cartesian projection centres for each orbital.
        intersite_positions = 0.0_dp
        do i = 1, num_bands
            do j = 1, 3
                intersite_positions(i, :) = intersite_positions(i, :)          &
                                          + (projection_centres(i, j)          &
                                          * avec(j, :))
            enddo
        enddo

        ! Build intersite difference matrices for the full Floquet basis
        intersite_diffs = 0.0_dp
        do i = 1, nf_bands
            orb_i = mod(i - 1, num_bands) + 1
            do j = 1, nf_bands
                orb_j = mod(j - 1, num_bands) + 1
                intersite_diffs(i, j, 1) = intersite_positions(orb_j, 1)       &
                                         - intersite_positions(orb_i, 1)
                intersite_diffs(i, j, 2) = intersite_positions(orb_j, 2)       &
                                         - intersite_positions(orb_i, 2)
            enddo
        enddo

        select case(direction)
            case(1)
                slab_area = norm2(cross_product(avec(2, :), avec(3, :)))
            case(2)
                slab_area = norm2(cross_product(avec(1, :), avec(3, :)))
            case(3)
                slab_area = norm2(cross_product(avec(1, :), avec(2, :)))
            case default
                slab_area = 0.0_dp
        end select
    end subroutine initialise_parameters
!******************************************************************************
    pure function make_rlist_cart(num_r_pts, frac_rlist, avec) result(cart_rlist)
    implicit none
        integer,   intent(in) :: num_r_pts
        integer,   intent(in) :: frac_rlist(num_r_pts, 3)
        real(dp),  intent(in) :: avec(3, 3)

        integer               :: i, j

        real(dp), allocatable :: cart_rlist(:, :)

        allocate(cart_rlist(num_r_pts, 3))
        cart_rlist = 0.0_dp

        do i = 1, num_r_pts
            do j = 1, 3
                cart_rlist(i, :) = cart_rlist(i, :)                            &
                                 + (frac_rlist(i, j) * avec(j, :))
            enddo
        enddo
        end function make_rlist_cart
end module parameters
