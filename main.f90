program FloquetSurfaceWeighted
use kinds, only: dp
use mpi_f08
use mpi_utilities, only: split_klist, prepare_gatherv_info_3d
use BZ_utilities, only: make_kmesh
use hamiltonian, only: ft_hr
use parameters, only: nkp, nf_bands, initialise_parameters, num_r_pts,         &
    r_ham_list, floquet_switch, electric_field_si, nlayers, nene
use floquet, only: make_floquet_hamiltonian_real_space
use write_files, only: write_kdists, write_surface_projections
use electronic_structure, only: greens_function
implicit none
!--------------------------------MPI Variables---------------------------------
    integer                  :: ierr, nprocs, pid, ibeg, iend, sendcounts
    integer, allocatable     :: recvcounts(:), displs(:)
!-------------------------------Global Variables-------------------------------
    complex(dp), allocatable :: conductivity_tensor_glob(:, :)
!------------------------------------------------------------------------------
    integer                  :: nkpar
    real(dp), allocatable    :: klist(:, :), kdists(:)
    complex(dp), allocatable :: floquet_r_ham_list(:, :, :)
    complex(dp), allocatable :: conductivity_tensor(:, :, :)

    
!-----------------------------Initialise MPI here------------------------------
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
!------------------------------------------------------------------------------
    call initialise_parameters
    klist = make_kmesh()
    call split_klist(nkp, nprocs, pid, ibeg, iend, nkpar) ! for MPI
!------------------------------------------------------------------------------
    allocate(floquet_r_ham_list(num_r_pts, nf_bands, nf_bands))
    if (pid .eq. 0) then
        write(*, fmt="(A)", advance="no") "Chosen parameters correspond to"
        write(*, fmt="(A)", advance="no") " an electric field strength of "
        write(*, fmt="(F8.4,A)") electric_field_si * 1.0E-9_dp, " V/nm"
    endif
    call make_floquet_hamiltonian_real_space(floquet_r_ham_list)
    floquet_r_ham_list = r_ham_list
    if (pid .eq. 0) then
        write(*, fmt="(A)") "Floquet hamiltonian built!"
    endif

    if (pid .eq. 0) then
        allocate(conductivity_tensor_glob(2, 2, n_ene))
    else
        allocate(conductivity_tensor_glob(1, 1, 1))
    endif


    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    allocate(conductivity_tensor(2, 2, n_ene))
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    deallocate(klist)
    
    call prepare_gatherv_info_3d(green_func_loc, nkp, nprocs, sendcounts,      &
        recvcounts, displs)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (nprocs /= 1) then
        call MPI_Gatherv(green_func_loc, sendcounts, MPI_DOUBLE_COMPLEX,       &
            greens_function_glob, recvcounts, displs, MPI_DOUBLE_COMPLEX,      &
            0, MPI_COMM_WORLD, ierr)
    else
        greens_function_glob = green_func_loc
    endif

    

    deallocate(green_func_loc)
    if(pid .eq. 0) then
        call write_kdists(kdists)
        call write_surface_projections(greens_function_glob)
    endif
    call MPI_FINALIZE(ierr)
end program FloquetSurfaceWeighted
