module transport
use kinds, only: dp
implicit none
private
public :: floquet_kubo_greenwood
contains
!******************************************************************************
    pure function floquet_kubo_greenwood(floquet_r_ham_list, static_r_ham_list,&
            photon_0_start, photon_0_end, klist) result(conductivity_tensor)
    use hamiltonian, only: slab_hamiltonian, slab_velocities_xy
    use parameters, only: num_r_pts, nf_bands, nkp, num_bands, nkp, nene,      &
        nlayers
    implicit none
    external ZHEEVD
        complex(dp), intent(in) :: static_r_ham_list(num_r_pts, num_bands, num_bands)
        complex(dp), intent(in) :: floquet_r_ham_list(num_r_pts, nf_bands, nf_bands)
        real(dp),    intent(in) :: kmesh(3, nkp)
        integer,     intent(in) :: photon_0_start, photon_0_end

!-------------------------------ZHEEVD Variables-------------------------------
        integer                  :: lwork, lrwork, liwork, stat
        real(dp),    allocatable :: rwork(:), trwork(:)
        complex(dp), allocatable :: work(:), twork(:)
        integer,     allocatable :: iwork(:), tiwork(:)

        integer                  :: ik, ts_bands, tf_bands
        complex(dp), allocatable :: static_eigenmat(:, :), floquet_eigenmat(:, :)
        real(dp),    allocatable :: static_eigvals(:), floquet_eigvals(:)
        real(dp),   dimension(3) :: k
        logical                  :: workspace_allocated

    end function floquet_kubo_greenwood
!******************************************************************************
    pure function floquet_fermionic_occ(static_energies, static_states,        &
        floquet_states, n_bands_floq, n_bands_stat) result(occupancy)
    use parameters, only: e_fermi=>fermi_energy, num_photon
    use statistical_distributions only: fermi_dirac
    implicit none
        integer,     intent(in) :: n_bands_floq, n_bands_stat
        complex(dp), intent(in) :: static_states(n_bands_stat, n_bands_stat)
        complex(dp), intent(in) :: floquet_states(n_bands_floq, n_bands_floq)
        real(dp),    intent(in) :: static_energies(n_bands_stat)
        real(dp),    intent(in) :: floquet_energies(n_bands_floq)

        integer                 :: mu, nu, photon_0_start, photon_0_end
        real(dp)                :: projection
        complex(dp)             :: inner_prod

        real(dp)                :: occupancy(n_bands_floq)

        photon_0_start = 1 + (num_photon * n_bands_stat)
        photon_0_end = (1 + num_photon) * n_bands_stat

        do nu = 1, n_bands_floq
            occupancy(nu) = 0.0_dp
            do mu = 1, n_bands_stat
                associate(f_nu=>occupancy(nu), g_mu=>static_states(:, mu),     &
                    ene_mu=>(static_energies(mu)-e_fermi),                     &
                    phi_nu_0=>floquet_states(photon_0_start:photon_0_end, nu))

                    inner_prod = dot_product(g_mu, phi_nu_0)
                    projection = real(conjg(inner_prod) * inner_prod, kind=dp)
                    f_nu = f_nu + (fermi_dirac(ene_mu) * projection)
                end associate
            enddo
        enddo
    end function floquet_fermionic_occ
end module transport
