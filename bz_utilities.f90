module BZ_utilities
use kinds, only: dp
implicit none
private
public :: make_kmesh
contains
!******************************************************************************
    pure function make_kmesh() result(klist)
    use array_utilities, only: linspace
    use parameters, only: nk, nkp, k_shift, k_frac
    implicit none
    ! returns a one-dimensional array containing all the kpoints
        integer               :: i1, i2, i3, ik
        real(dp), allocatable :: k1_list(:), k2_list(:), k3_list(:)
        real(dp)              :: k1_beg, k1_end, k2_beg, k2_end, k3_beg, k3_end
        real(dp), allocatable :: klist(:, :)

        allocate(klist(3, nkp))

        k1_beg = k_shift(1) - (k_frac(1) / 2.0_dp)
        k2_beg = k_shift(2) - (k_frac(2) / 2.0_dp)
        k3_beg = k_shift(3) - (k_frac(3) / 2.0_dp)

        k1_end = k_shift(1) + (k_frac(1) / 2.0_dp)
        k2_end = k_shift(2) + (k_frac(2) / 2.0_dp)
        k3_end = k_shift(3) + (k_frac(3) / 2.0_dp)
    
        k1_list = linspace(k1_beg, k1_end, nk(1), .false.)
        k2_list = linspace(k2_beg, k2_end, nk(2), .false.)
        k3_list = linspace(k3_beg, k3_end, nk(3), .false.)

        ik=1
        do i1=1, nk(1)
            do i2=1, nk(2)
                do i3=1, nk(3)
                    klist(:, ik) = (/k1_list(i1), k2_list(i2), k3_list(i3)/)
                    ik = ik + 1
                enddo
            enddo
        enddo
    end function make_kmesh
!******************************************************************************
end module BZ_utilities
