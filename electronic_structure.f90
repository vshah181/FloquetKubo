module electronic_structure
use kinds, only: dp
implicit none
private
public :: greens_function
contains
!******************************************************************************
    subroutine weighted_greens_function(eigenvalues, eigenstates, green_func)
    use parameters, only: nene, nlayers, nf_bands, energy_list, broadening_factor
    use array_utilities, only: inner_product
    use constants, only: cmplx_i, cmplx_0
    implicit none
        real(dp),    intent(in)  :: eigenvalues(nlayers * nf_bands)
        complex(dp), intent(in)  :: eigenstates(nlayers * nf_bands, nlayers * nf_bands)
        complex(dp), intent(out) :: green_func(nlayers, nene)

        integer                  :: ien, ib, il, tf_bands, idx_beg, idx_end
        real(dp)                 :: eigenval, layer_norm_sq, min_eps, sp
        complex(dp)              :: eigenstate(nlayers * nf_bands), numerator
        complex(dp)              :: local_layer_state(nf_bands), denominator

        min_eps = tiny(1.0_dp)
        green_func = cmplx_0

        tf_bands = nlayers * nf_bands
        do ien = 1, nene
            do ib = 1, nlayers * nf_bands
                eigenstate = eigenstates(:, ib) 
                eigenval = eigenvalues(ib)
                associate(ieta=>(broadening_factor * cmplx_i))
                    denominator = energy_list(ien) - eigenval + ieta
                end associate
                do il = 1, nlayers
                    idx_beg = ((il - 1) * nf_bands) + 1
                    idx_end = il * nf_bands
                    local_layer_state = eigenstate(idx_beg:idx_end)
                    associate(lls=>local_layer_state)
                        layer_norm_sq = real(dot_product(lls, lls), kind=dp)
                        if (layer_norm_sq .le. min_eps) then
                            numerator = 0.0_dp
                        else
                            sp = spectral_weight(lls / sqrt(layer_norm_sq))
                            numerator = sp * layer_norm_sq
                        endif
                    end associate
                    associate(gfloc=>green_func(il, ien))
                        gfloc = gfloc + (numerator / denominator)
                    end associate
                enddo
            enddo  
        enddo
    end subroutine weighted_greens_function
!******************************************************************************
    subroutine unweighted_greens_function(eigenvalues, eigenstates, green_func)
    use parameters, only: nene, nlayers, nf_bands, energy_list, broadening_factor
    use array_utilities, only: inner_product
    use constants, only: cmplx_i, cmplx_0
    implicit none
        real(dp),    intent(in)  :: eigenvalues(nlayers * nf_bands)
        complex(dp), intent(in)  :: eigenstates(nlayers * nf_bands, nlayers * nf_bands)
        complex(dp), intent(out) :: green_func(nlayers, nene)

        integer                  :: ien, ib, il, tf_bands, idx_beg, idx_end
        real(dp)                 :: eigenval
        complex(dp)              :: eigenstate(nlayers * nf_bands), numerator
        complex(dp)              :: local_layer_state(nf_bands), denominator

        green_func = cmplx_0

        tf_bands = nlayers * nf_bands
        do ien = 1, nene
            do ib = 1, nlayers * nf_bands
                eigenstate = eigenstates(:, ib) 
                eigenval = eigenvalues(ib)
                associate(ieta=>(broadening_factor * cmplx_i))
                    denominator = energy_list(ien) - eigenval + ieta
                end associate
                do il = 1, nlayers
                    idx_beg = ((il - 1) * nf_bands) + 1
                    idx_end = il * nf_bands
                    local_layer_state = eigenstate(idx_beg:idx_end)
                    numerator = dot_product(local_layer_state, local_layer_state)
                    associate(gfloc=>green_func(il, ien))
                        gfloc = gfloc + (numerator / denominator)
                    end associate
                enddo
            enddo  
        enddo
    end subroutine unweighted_greens_function
!******************************************************************************
    subroutine greens_function(floquet_r_ham_list, ibeg, iend, klist, &
            green_func)
    use hamiltonian, only: slab_hamiltonian
    use parameters, only: num_r_pts, nf_bands, nkp, nene, nlayers,
        floquet_switch
    implicit none
    external ZHEEVD
        complex(dp), intent(in)  :: floquet_r_ham_list(num_r_pts, nf_bands, nf_bands)
        real(dp),    intent(in)  :: klist(nkp, 3)
        integer,     intent(in)  :: ibeg, iend

        complex(dp), intent(out) :: green_func(nlayers, nene, ibeg:iend)
!-------------------------------ZHEEVD Variables-------------------------------
        integer                  :: lwork, lrwork, liwork, stat
        real(dp),    allocatable :: rwork(:), trwork(:)
        complex(dp), allocatable :: work(:), twork(:)
        integer,     allocatable :: iwork(:), tiwork(:)

        integer                  :: ik, t_bands
        complex(dp), allocatable :: eigenmat(:, :)
        real(dp),    allocatable :: eigvals(:)
        real(dp),   dimension(3) :: k
        logical                  :: workspace_allocated
!--------------------------Allocate arrays for ZHEEVD--------------------------
        if ((ibeg .lt. 1) .or. (iend .gt. nkp)) then
            error stop "Invalid ibeg/iend."
        endif

        t_bands = nlayers * nf_bands
        
        lwork  = -1
        lrwork = -1
        liwork = -1
        allocate(eigenmat(t_bands, t_bands), eigvals(t_bands))
        allocate(twork(max(lwork, 1)), trwork(max(lrwork, 1))) 
        allocate(tiwork(max(liwork, 1)))

        workspace_allocated = .false.

        do ik = ibeg, iend
            k = klist(ik, :)
            ! First, we need to get the eigenvalues and eigenstates
            call slab_hamiltonian(k, floquet_r_ham_list, eigenmat)
            if(.not.(workspace_allocated)) then
                call ZHEEVD('V', 'L', t_bands, eigenmat, t_bands, eigvals,     &
                    twork, lwork, trwork, lrwork, tiwork, liwork, stat)
                if (stat .ne. 0) then
                    error stop "Failed to allocate workspace arrays for zheevd"
                endif
                lwork=int(real(twork(1), kind=dp))
                lrwork=int(trwork(1))
                liwork=tiwork(1)
                allocate(work(lwork), rwork(lrwork), iwork(liwork))
                deallocate(twork, trwork, tiwork)
                workspace_allocated = .true.
                ! recompute eigenmat, just in case it has been filled with 
                ! garbage
                call slab_hamiltonian(k, floquet_r_ham_list, eigenmat)
            endif
            call ZHEEVD('V', 'L', t_bands, eigenmat, t_bands, eigvals, work,   &
                lwork, rwork, lrwork, iwork, liwork, stat)
            if (stat .ne. 0) then
                error stop "ZHEEVD failed!"
            endif
            ! Now, we must calculate the greens function (layer-resolved)
            associate(gf_slab=>green_func(:, :, ik))
                if (floquet_switch) then
                    call weighted_greens_function(eigvals, eigenmat, gf_slab)
                else
                    call unweighted_greens_function(eigvals, eigenmat, gf_slab)
                endif
            end associate
        enddo
    end subroutine greens_function
!******************************************************************************
    pure function spectral_weight(eigenstate) result(weight)
    use parameters, only: nf_bands, num_photon, num_bands
    implicit none
        complex(dp), intent(in) :: eigenstate(nf_bands)

        integer                 :: istart, iend

        real(dp)                :: weight

        istart = 1 + (num_bands * num_photon)
        iend   = num_bands * (num_photon + 1)

        if (iend .gt. nf_bands) then
            error stop "out of bound error for eigenstates!"
        endif

        weight = real(dot_product(eigenstate(istart:iend),                     &
            eigenstate(istart:iend)), kind=dp)
    end function spectral_weight
end module electronic_structure
