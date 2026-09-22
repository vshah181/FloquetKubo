program FloquetKubo
use kinds, only: dp
use mpi_f08
use mpi_utilities, only: split_klist
use BZ_utilities, only: make_kmesh
use parameters, only: nkp, nf_bands, initialise_parameters, num_r_pts,         &
    r_ham_list, electric_field_si, nene
use floquet, only: make_floquet_hamiltonian_real_space
use write_files, only: write_conductivity_real, write_conductivity_imag
use transport, only: compute_conductivities
use constants, only: elementary_charge
implicit none
!--------------------------------MPI Variables---------------------------------
    integer                  :: ierr, nprocs, pid, ibeg, iend
!-------------------------------Global Variables-------------------------------
    complex(dp), allocatable :: conductivity_tensor_glob(:, :, :)
!------------------------------------------------------------------------------
    integer                  :: nkpar
    real(dp), allocatable    :: klist(:, :)
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
    if (pid .eq. 0) then
        write(*, fmt="(A)") "Floquet hamiltonian built!"
    endif

    if (pid .eq. 0) then
        allocate(conductivity_tensor_glob(2, 2, nene))
    else
        allocate(conductivity_tensor_glob(1, 1, 1))
    endif


    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    allocate(conductivity_tensor(2, 2, nene))
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    conductivity_tensor = compute_conductivities(floquet_r_ham_list,           &
        r_ham_list, klist, ibeg, iend) * elementary_charge  ! convert to SI

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    deallocate(klist)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (nprocs /= 1) then
        call MPI_Reduce(conductivity_tensor, conductivity_tensor_glob,         &
            2 * 2 * nene, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    else
        conductivity_tensor_glob = conductivity_tensor
    endif

    

    if(pid .eq. 0) then
        call write_conductivity_real(conductivity_tensor_glob)
        call write_conductivity_imag(conductivity_tensor_glob)
    endif
    call MPI_FINALIZE(ierr)
end program FloquetKubo
