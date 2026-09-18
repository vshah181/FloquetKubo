module hamiltonian
use kinds, only: dp
implicit none
private
public :: ft_hr, slab_hamiltonian, slab_velocities_xy
contains
    pure function ft_hr(k, r_ham_list, n_bands) result(hk)
    use parameters, only: num_r_pts, r_list, weights
    use constants, only: cmplx_0, cmplx_i, two_pi
    implicit none
        integer,     intent(in) :: n_bands
        real(dp),    intent(in) :: k(3)
        complex(dp), intent(in) :: r_ham_list(num_r_pts, n_bands, n_bands)

        real(dp)                :: phase, k_scaled(3)
        integer                 :: ir

        complex(dp)             :: hk(n_bands, n_bands)

        hk = cmplx_0
        k_scaled = k * two_pi
        do ir = 1, num_r_pts
            phase = dot_product(r_list(ir, :), k_scaled)
            associate(hr=>r_ham_list(ir, :, :), w=>real(weights(ir), kind=dp))
                hk = hk + ((hr / w) * exp(cmplx_i * phase))
            end associate
        enddo
    end function ft_hr
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
    pure subroutine slab_velocities_xy(k, r_ham_list, n_bands, vxy_slab) ! V_xy
    use parameters, only: num_r_pts, r_list, weights, nlayers, rlist_cart 
    use constants, only: cmplx_0, cmplx_i, two_pi
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

                prefac_matrix = cmplx_i * rham * exp(cmplx_i * phase) / w
                do il=1, nlayers - int(abs(r(3)))
                    associate(vx_block=>vxy_slab(irow:irow + n_bands - 1,      &
                            icol:icol + n_bands - 1, 1),                       &
                            vy_block=>vxy_slab(irow:irow + n_bands - 1,        &
                            icol:icol + n_bands - 1, 2))

                        vx_block = vx_block + (rx * prefac_matrix)
                        vy_block = vy_block + (ry * prefac_matrix)

                        irow=irow+n_bands
                        icol=icol+n_bands

                    end associate
                enddo
            end associate
        enddo
    end subroutine slab_velocities_xy
end module hamiltonian
