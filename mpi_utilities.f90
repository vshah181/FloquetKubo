module mpi_utilities
implicit none
private
public :: split_klist, prepare_gatherv_info_3d
contains
    subroutine split_klist(nkp, nprocs, pid, ibeg, iend, nkpar)
        implicit none
        integer, intent(in) :: nkp, nprocs, pid
        integer, intent(out) :: iend, ibeg, nkpar
        integer :: extra
    
        extra=mod(nkp, nprocs)
        ibeg=1+(pid*(nkp/nprocs))
        if (pid<extra) then
            ibeg=ibeg+pid
            iend=ibeg+(nkp/nprocs)
        else
            ibeg=ibeg+extra
            iend=ibeg+(nkp/nprocs)-1
        endif
        nkpar=iend-(ibeg-1)
    end subroutine split_klist
!******************************************************************************
    subroutine prepare_gatherv_info_3d(local_array, ndim_glob, nprocs,         &
            sendcounts, recvcounts, displs)
    use kinds, only: dp
    implicit none
        complex(dp),           intent(in) :: local_array(:, :, :) 
        integer,               intent(in) :: nprocs, ndim_glob
        integer, allocatable, intent(out) :: recvcounts(:), displs(:)
        integer,              intent(out) :: sendcounts
        integer                           :: ip, extra, ndim_1, ndim_2

        sendcounts = size(local_array)
        allocate(recvcounts(nprocs), displs(nprocs))
        recvcounts(1) = 0
        displs(1) = 0
        ndim_1 = size(local_array, dim=1)
        ndim_2 = size(local_array, dim=2)
        extra = mod(ndim_glob, nprocs)

        do ip=1, nprocs
            if(ip > extra) then
                recvcounts(ip) = (ndim_glob / nprocs) * ndim_1 * ndim_2
            else
                recvcounts(ip) = ((ndim_glob / nprocs) + 1) * ndim_1 * ndim_2
            endif
            if(ip < nprocs) then
                displs(ip + 1) = displs(ip) + recvcounts(ip)
            endif
        enddo
    end subroutine prepare_gatherv_info_3d
end module mpi_utilities
