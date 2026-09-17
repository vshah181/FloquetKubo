module transport
use kinds, only: dp
implicit none
private
public :: compute_conductivities
contains
!******************************************************************************
    function compute_conductivities(floquet_r_ham_list, static_r_ham_list,     &
            klist, ibeg, iend) result(conductivity_tensor)
    use hamiltonian, only: slab_hamiltonian, slab_velocities_xy
    use parameters, only: num_r_pts, nf_bands, nkp, num_bands, nkp, nene,      &
        nlayers, energy_list
    use constants, only: cmplx_0
    implicit none
    external ZHEEVD
        complex(dp), intent(in) :: static_r_ham_list(num_r_pts, num_bands, num_bands)
        complex(dp), intent(in) :: floquet_r_ham_list(num_r_pts, nf_bands, nf_bands)
        real(dp),    intent(in) :: kmesh(3, nkp)
        integer,     intent(in) :: photon_0_start, photon_0_end, ibeg, iend

!--------------------------ZHEEVD Variables (floquet)--------------------------
        integer                  :: flwork, flrwork, fliwork, fstat
        real(dp),    allocatable :: frwork(:), ftrwork(:)
        complex(dp), allocatable :: fwork(:), ftwork(:)
        integer,     allocatable :: fiwork(:), ftiwork(:)
!--------------------------ZHEEVD Variables (static)---------------------------
        integer                  :: slwork, slrwork, sliwork, sstat
        real(dp),    allocatable :: srwork(:), strwork(:)
        complex(dp), allocatable :: swork(:), stwork(:)
        integer,     allocatable :: siwork(:), stiwork(:)

        integer                  :: ik, ts_bands, tf_bands, iw
        complex(dp), allocatable :: static_eigenmat(:, :), floquet_eigenmat(:, :)
        complex(dp), allocatable :: velocity_operators(:, :, :)
        real(dp),    allocatable :: static_eigvals(:), floquet_eigvals(:)
        real(dp)                 :: k(3), occupations(nlayers * nf_bands)
        logical                  :: floquet_workspace_allocated
        logical                  :: static_workspace_allocated

        complex(dp)              :: conductivity_tensor(2, 2, nene)
!--------------------------Allocate arrays for ZHEEVD--------------------------
        if ((ibeg .lt. 1) .or. (iend .gt. nkp)) then
            error stop "Invalid ibeg/iend."
        endif

        ts_bands = nlayers * num_bands
        tf_bands = nlayers * nf_bands
        
        flwork  = -1
        flrwork = -1
        fliwork = -1
        allocate(floquet_eigenmat(tf_bands, tf_bands), floquet_eigvals(tf_bands))
        allocate(ftwork(max(flwork, 1)), ftrwork(max(flrwork, 1))) 
        allocate(ftiwork(max(fliwork, 1)))

        slwork  = -1
        slrwork = -1
        sliwork = -1
        allocate(static_eigenmat(ts_bands, ts_bands), static_eigvals(ts_bands))
        allocate(stwork(max(slwork, 1)), strwork(max(slrwork, 1))) 
        allocate(stiwork(max(sliwork, 1)))

        static_workspace_allocated = .false.
        floquet_workspace_allocated = .false.

        allocate(velocity_operators(tf_bands, tf_bands, 2))

        conductivity_tensor = cmplx_0

        do ik = ibeg, iend
            k = klist(:, ik)
            ! First, we need to allocate ZHEEVD arrays (static)
            call slab_hamiltonian(k, static_r_ham_list, num_bands, static_eigenmat)
            if(.not.(static_workspace_allocated)) then
                call ZHEEVD('V', 'L', ts_bands, static_eigenmat, ts_bands,     &
                    static_eigvals, stwork, slwork, strwork, slrwork, stiwork, &
                    sliwork, sstat)
                if (sstat .ne. 0) then
                    error stop "Failed to allocate workspace arrays for ZHEEVD"
                endif
                slwork = int(real(stwork(1), kind=dp))
                slrwork = int(strwork(1))
                sliwork = stiwork(1)
                allocate(swork(slwork), srwork(slrwork), siwork(sliwork))
                deallocate(stwork, strwork, stiwork)
                static_workspace_allocated = .true.
                ! recompute eigenmat, just in case it has been filled with 
                ! garbage
                call slab_hamiltonian(k, static_r_ham_list, num_bands,         &
                    static_eigenmat)
            endif

            ! Next, we need to allocate ZHEEVD arrays (floquet)
            call slab_hamiltonian(k, floquet_r_ham_list, nf_bands,             &
                floquet_eigenmat)
            if(.not.(floquet_workspace_allocated)) then
                call ZHEEVD('V', 'L', tf_bands, floquet_eigenmat, tf_bands,    &
                    floquet_eigvals, ftwork, flwork, ftrwork, flrwork, ftiwork,&
                    fliwork, fstat)
                if (fstat .ne. 0) then
                    error stop "Failed to allocate workspace arrays for ZHEEVD"
                endif
                flwork = int(real(ftwork(1), kind=dp))
                flrwork = int(ftrwork(1))
                fliwork = ftiwork(1)
                allocate(fwork(flwork), frwork(flrwork), fiwork(fliwork))
                deallocate(ftwork, ftrwork, ftiwork)
                floquet_workspace_allocated = .true.
                ! recompute eigenmat, just in case it has been filled with 
                ! garbage
                call slab_hamiltonian(k, floquet_r_ham_list, nf_bands,         &
                    floquet_eigenmat)
            endif

            ! Now, diagonalise static and Floquet Hamiltonian
            call ZHEEVD("V", "L", ts_bands, static_eigenmat, ts_bands,         &
                static_eigvals, swork, slwork, srwork, slrwork, siwork,        &
                sliwork, sstat)
            if (sstat .ne. 0) then
                error stop "ZHEEVD failed!"
            endif

            call ZHEEVD("V", "L", tf_bands, floquet_eigenmat, tf_bands,        &
                floquet_eigvals, fwork, flwork, frwork, flrwork, fiwork,       &
                fliwork, fstat)
            if (fstat .ne. 0) then
                error stop "ZHEEVD failed!"
            endif

            ! Now we need velocity matrices
            call slab_velocities_xy(k, floquet_r_ham_list, nf_bands,           &
                velocity_operators)

            ! Now we need the occupations
            occupations = floquet_fermionic_occ(static_eigvals,                &
                static_eigenmat, floquet_eigenmat, tf_bands, ts_bands)

            ! Finally, sum over probe energies and apply Kubo-Greenwood
            do iw = 1, nene
                conductivity_tensor(:, :, iw) = conductivity_tensor(:, :, iw)
                                              + kubo_greenwood(occupations,    &
                                                floquet_eigvals,               &
                                                velocity_operators,            &
                                                energy_list(iw), tf_bands)
            enddo
        enddo
    end function compute_conductivities
!******************************************************************************
    pure function floquet_fermionic_occ(static_energies, static_states,        &
        floquet_states, n_bands_floq, n_bands_stat) result(occupancy)
    use parameters, only: e_fermi=>fermi_energy, num_photon
    use statistical_distributions, only: fermi_dirac
    implicit none
        integer,     intent(in) :: n_bands_floq, n_bands_stat
        complex(dp), intent(in) :: static_states(n_bands_stat, n_bands_stat)
        complex(dp), intent(in) :: floquet_states(n_bands_floq, n_bands_floq)
        real(dp),    intent(in) :: static_energies(n_bands_stat)

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
!******************************************************************************
    pure function kubo_greenwood(occupations, energies, velocities, probe,     &
            n_bands) result(summand)
    use parameters, only: broadening_factor
    use constants, only: cmplx_0, cmplx_i
    implicit none
        integer,     intent(in) :: n_bands
        real(dp),    intent(in) :: energies(n_bands), occupations(n_bands), probe
        complex(dp), intent(in) :: velocities(n_bands, n_bands, 2)
        
        real(dp), parameter     :: tol = 1.0E-16_dp
        integer                 :: a, b, x, y
        real(dp)                :: ediff 
        complex(dp)             :: numerator, denominator, ieta

        complex(dp)             :: summand(2, 2)

        ieta = cmplx_i * broadening_factor
        summand = cmplx_0

        ! compute sigma_{x y}
        do x = 1, 2
            do y = 1, 2
                do a = 1, n_bands
                    do b = 1, n_bands
                        ediff = energies(b) - energies(a)
                        numerator = (occupations(a) - occupations(b))          &
                                  * velocities(a, b, x) * velocities(b, a, y)
                        denominator = (probe + ieta - ediff) * ediff
                        if (abs(ediff) .gt. tol) then
                            summand(x, y) = summand(x, y)                      &
                                          + (numerator / denominator)
                        else
                            summand(x, y) = summand(x, y) + cmplx_0
                        endif
                    enddo
                enddo
            enddo
        enddo
    end function kubo_greenwood

end module transport
