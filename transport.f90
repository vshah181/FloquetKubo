module transport
use kinds, only: dp
implicit none
private
public :: floquet_kubo_greenwood
contains
!******************************************************************************
    pure function floquet_kubo_greenwood(floquet_r_ham_list, static_r_ham_list,&
            ibeg, iend, klist) result(conductivity_tensor)
    use hamiltonian, only: slab_hamiltonian, slab_velocities_xy
    use parameters, only: num_r_pts, nf_bands, nkp, num_bands, nkp, nene,      &
        nlayers
    implicit none
    external ZHEEVD
        complex(dp), intent(in) :: static_r_ham_list(num_r_pts, num_bands, num_bands)
        complex(dp), intent(in) :: floquet_r_ham_list(num_r_pts, nf_bands, nf_bands)
        real(dp),    intent(in) :: kmesh(3, nkp)
        integer,     intent(in) :: ibeg, iend

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
    pure function floquet_fermionic_occ()
    end function floquet_fermionic_occ
end module transport
