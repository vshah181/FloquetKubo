module BZ_utilities
use kinds, only: dp
implicit none
private
public :: make_kpath, get_kdists
contains
!******************************************************************************
    function make_kpath() result(klist)
    use array_utilities, only: linspace
    use parameters, only: nkpt_per_path, high_sym_pts, nkp
    implicit none
    ! Returns the kpath from high symmetry points (eg for band-structure plot)
        real(dp), allocatable :: klist(:, :)

        integer :: i, ih, id
        real(dp), allocatable :: temp(:)
        logical :: include_endpoint

        allocate(klist(nkp, 3))

        ! high_sym_pts has shape (nkpath + 1, 3)

        ih = 0
        do i=1, nkp, nkpt_per_path
            ih = ih + 1
            include_endpoint = ((i + nkpt_per_path) .gt. nkp)
            do id=1, 3
                temp = linspace(high_sym_pts(ih, id), high_sym_pts(ih+1, id),  &
                    nkpt_per_path, include_endpoint)
                klist(i:i+nkpt_per_path-1, id) = temp
                deallocate(temp)
            enddo
        enddo
    end function make_kpath
!******************************************************************************
    function get_kdists(klist) result(kdists)
    use parameters, only: bvec, nkp
    implicit none
    ! Returns the kdists
        real(dp), intent(in) :: klist(:, :)
        real(dp), allocatable :: kdists(:)
        integer :: ik, id
        real(dp) :: dk_frac(3), dk_cart(3)

        allocate(kdists(nkp))

        kdists(1) = 0_dp
        
        do ik = 2, nkp
            dk_frac = klist(ik, :) - klist(ik - 1, :)
            dk_cart = (/ 0_dp, 0_dp, 0_dp /)
            do id = 1, 3
                dk_cart = dk_cart + (dk_frac(id) * bvec(id, :))
            enddo
            kdists(ik) = kdists(ik - 1) + norm2(dk_cart)
        enddo
    end function get_kdists
!******************************************************************************
    pure function make_kmesh result(klist)
    use array_utilities, only: linspace
    use parameters, only: nk, k_shift, k_frac
    implicit none
    ! returns a one-dimensional array containing all the kpoints
        integer               :: i1, i2, i3, ik
        real(dp), allocatable :: k1_list(:), k2_list(:), k3_list(:)
        real(dp)              :: k1_beg, k1_end, k2_beg, k2_end, k3_end, k3_end
        real(dp), allocatable :: klist(:, :)

        allocate(klist(3, nkp))

        k1_beg = k_shift(1) - (k_frac(1) / 2.0_dp)
        k2_beg = k_shift(2) - (k_frac(2) / 2.0_dp)
        k3_beg = k_shift(3) - (k_frac(3) / 2.0_dp)

        k1_end = k_shift(1) + (k_frac(1) / 2.0_dp)
        k2_end = k_shift(2) + (k_frac(2) / 2.0_dp)
        k3_end = k_shift(3) + (k_frac(3) / 2.0_dp)
    
        k1_list = linspace (k1_beg, k2_end, nk(1), .false.)
        k2_list = linspace (k2_beg, k2_end, nk(2), .false.)
        k2_list = linspace (k3_beg, k3_end, nk(3), .false.)

        ik=1
        do i1=1, nk(1)
            do i2=1, nk(2)
                do i3=1, nk(3)
                    klist(:, ik) = (/k1_list(i1), k2_list(i2), k3_list(i3)/)
                    ik = ik + 1
                enddo
            enddo
        enddo
    end subroutine
end module BZ_utilities
