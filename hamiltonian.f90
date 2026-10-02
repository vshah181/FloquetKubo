module hamiltonian
use kinds, only: dp
implicit none
private
public :: bulk_hamiltonian, bulk_velocities_xy
public :: slab_hamiltonian, slab_velocities_xy
contains
    pure subroutine bulk_hamiltonian(k, r_ham_list, n_bands, hk)
    use parameters, only: num_r_pts, r_list, weights
    use constants, only: cmplx_0, cmplx_i, two_pi
    implicit none
        integer,      intent(in) :: n_bands
        real(dp),     intent(in) :: k(3)
        complex(dp),  intent(in) :: r_ham_list(num_r_pts, n_bands, n_bands)

        real(dp)                 :: phase
        integer                  :: ir

        complex(dp), intent(out) :: hk(n_bands, n_bands)

        hk = cmplx_0
        do ir = 1, num_r_pts
            phase = dot_product(k * two_pi, r_list(ir, :))
            associate(hr=>r_ham_list(ir, :, :), w=>real(weights(ir), kind=dp))
                hk = hk + (hr * exp(cmplx_i * phase) / w)
            end associate
        enddo
    end subroutine bulk_hamiltonian
!******************************************************************************
    pure subroutine slab_hamiltonian(k, r_ham_list, n_bands, hk_slab)
    use parameters, only: num_r_pts, r_list, weights, nlayers, id=>direction
    use constants, only: cmplx_0, cmplx_i, two_pi
    implicit none
        integer,      intent(in) :: n_bands
        real(dp),     intent(in) :: k(3)
        complex(dp),  intent(in) :: r_ham_list(num_r_pts, n_bands, n_bands)

        real(dp)                 :: phase, k_scaled(3)
        integer                  :: ir, il, irow, icol

        complex(dp), intent(out) :: hk_slab(n_bands * nlayers,                 &
            n_bands * nlayers)

        hk_slab = cmplx_0
        k_scaled = k * two_pi

        do ir = 1, num_r_pts
            associate(rham=>r_ham_list(ir, :, :),                              &
                w=>real(weights(ir), kind=dp), r=>r_list(ir, :))
                phase=dot_product(k_scaled, r)
                ! irow = 1 if +ve z, else irow=|z| * n_bands
                irow = int((n_bands * (abs(r(id)) - r(id)) / 2) + 1)
                ! icol = 1 if -ve z, else icol=|z| * n_bands
                icol = int((n_bands * (abs(r(id)) + r(id)) / 2) + 1)
                do il=1, nlayers - int(abs(r(id)))
                    associate(block_ham=>hk_slab(irow:irow + n_bands - 1,      &
                        icol:icol + n_bands - 1))

                        block_ham = block_ham                                  &
                                  + (rham * exp(cmplx_i * phase) / w)
                        irow=irow+n_bands
                        icol=icol+n_bands

                    end associate
                enddo
            end associate
        enddo
    end subroutine slab_hamiltonian
!******************************************************************************
    pure subroutine bulk_velocities_xy(k, r_ham_list, n_bands, vxy_bulk) ! V_xy
    use parameters, only: num_r_pts, r_list, weights, rlist_cart,              &
        intersite_diffs
    use constants, only: cmplx_0, cmplx_i, two_pi,                             &
        hbar=>reduced_planck_constant_ev
    ! assume slab in z / a_3 parallel to cartesian z
    implicit none
        integer,      intent(in) :: n_bands
        real(dp),     intent(in) :: k(3)
        complex(dp),  intent(in) :: r_ham_list(num_r_pts, n_bands, n_bands)

        real(dp)                 :: phase, k_scaled(3)
        complex(dp)              :: prefac_matrix(n_bands, n_bands)
        integer                  :: ir

        complex(dp), intent(out) :: vxy_bulk(n_bands, n_bands, 2)

        vxy_bulk = cmplx_0

        k_scaled = two_pi * k
        do ir = 1, num_r_pts
            associate(rham=>r_ham_list(ir, :, :), r=>r_list(ir, :),            &
                rx=>rlist_cart(ir, 1), ry=>rlist_cart(ir, 2),                  &
                w=>real(weights(ir), kind=dp))
                
                phase = dot_product(k_scaled, r)

                prefac_matrix = cmplx_i * rham * exp(cmplx_i * phase) / (hbar * w)
                associate(vx=>vxy_bulk(:, :, 1), vy=>vxy_bulk(:, :, 2))
                    ! Remember, we cannot simply differentiate the 
                    ! Hamiltonain wrt k! We must add inter-site phases!
                    vx = vx + (rx * prefac_matrix)                             &
                             + (intersite_diffs(:, :, 1) * prefac_matrix)
                    vy = vy + (ry * prefac_matrix)                             &
                             + (intersite_diffs(:, :, 2) * prefac_matrix)
                end associate
            end associate
        enddo
    end subroutine bulk_velocities_xy
!******************************************************************************
    pure subroutine slab_velocities_xy(k, r_ham_list, n_bands, vxy_slab) ! V_xy
    use parameters, only: num_r_pts, r_list, weights, nlayers, rlist_cart,     &
        intersite_diffs
    use constants, only: cmplx_0, cmplx_i, two_pi,                             &
        hbar=>reduced_planck_constant_ev
    ! assume slab in z / a_3 parallel to cartesian z
    implicit none
        integer,      intent(in) :: n_bands
        real(dp),     intent(in) :: k(3)
        complex(dp),  intent(in) :: r_ham_list(num_r_pts, n_bands, n_bands)

        real(dp)                 :: phase, k_scaled(3)
        complex(dp)              :: prefac_matrix(n_bands, n_bands)
        integer                  :: ir, il, irow, icol

        complex(dp), intent(out) :: vxy_slab(n_bands * nlayers,                &
            n_bands * nlayers, 2)

        vxy_slab = cmplx_0

        k_scaled = two_pi * k
        do ir = 1, num_r_pts
            associate(rham=>r_ham_list(ir, :, :), r=>r_list(ir, :),            &
                rx=>rlist_cart(ir, 1), ry=>rlist_cart(ir, 2),                  &
                w=>real(weights(ir), kind=dp))
                
                phase = dot_product(k_scaled, r)
                irow = int((n_bands * (abs(r(3)) - r(3)) / 2) + 1)
                ! icol = 1 if -ve z, else icol=|z| * n_bands
                icol = int((n_bands * (abs(r(3)) + r(3)) / 2) + 1)

                prefac_matrix = cmplx_i * rham * exp(cmplx_i * phase) / (hbar * w)
                do il=1, nlayers - int(abs(r(3)))
                    associate(vx_block=>vxy_slab(irow:irow + n_bands - 1,      &
                            icol:icol + n_bands - 1, 1),                       &
                            vy_block=>vxy_slab(irow:irow + n_bands - 1,        &
                            icol:icol + n_bands - 1, 2))

                        ! Remember, we cannot simply differentiate the 
                        ! Hamiltonain wrt k! We must add inter-site phases!
                        vx_block = vx_block + (rx * prefac_matrix)             &
                                 + (intersite_diffs(:, :, 1) * prefac_matrix)
                        vy_block = vy_block + (ry * prefac_matrix)             &
                                 + (intersite_diffs(:, :, 2) * prefac_matrix)

                        irow = irow + n_bands
                        icol = icol + n_bands

                    end associate
                enddo
            end associate
        enddo
    end subroutine slab_velocities_xy
end module hamiltonian
