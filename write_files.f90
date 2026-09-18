module write_files
use kinds, only: dp
implicit none
private
public :: write_conductivity_real, write_conductivity_imag
contains
    subroutine write_conductivity_real(conductivity_tensor)
    use parameters, only: seedname, nene, energy_list
    implicit none
        complex(dp), intent(in) :: conductivity_tensor(2, 2, nene)
        character(len=99)       :: filename
        integer                 :: iu, ie
        logical                 :: file_exists

        write(filename, fmt="(2A)") trim(adjustl(seedname)),                   &
            "_real_conductivity.csv"

        inquire(file=trim(adjustl(filename)), exist=file_exists)
        if (file_exists) then
            open(newunit=iu, file=trim(adjustl(filename)), status="replace",   &
                action="write")
        else
            open(newunit=iu, file=trim(adjustl(filename)), status="new",       &
                action="write")
        endif

        write(iu, fmt="(2A)") "# probe energy,      sigma_xx,      sigma_xy,", &
            "      sigma_yx,      sigma_yy"
        do ie = 1, nene
            write(iu, fmt="(ES14.6,4(A,ES14.6))") energy_list(ie),             &
                ",", real(conductivity_tensor(1, 1, ie), kind=dp),             &
                ",", real(conductivity_tensor(1, 2, ie), kind=dp),             &
                ",", real(conductivity_tensor(2, 1, ie), kind=dp),             &
                ",", real(conductivity_tensor(2, 2, ie), kind=dp)
        enddo

        close(iu)
    end subroutine write_conductivity_real
!******************************************************************************
    subroutine write_conductivity_imag(conductivity_tensor)
    use parameters, only: seedname, nene, energy_list
    implicit none
        complex(dp), intent(in) :: conductivity_tensor(2, 2, nene)
        character(len=99)       :: filename
        integer                 :: iu, ie
        logical                 :: file_exists

        write(filename, fmt="(2A)") trim(adjustl(seedname)),                   &
            "_imag_conductivity.csv"

        inquire(file=trim(adjustl(filename)), exist=file_exists)
        if (file_exists) then
            open(newunit=iu, file=trim(adjustl(filename)), status="replace",   &
                action="write")
        else
            open(newunit=iu, file=trim(adjustl(filename)), status="new",       &
                action="write")
        endif

        write(iu, fmt="(2A)") "# probe energy,      sigma_xx,      sigma_xy,", &
            "      sigma_yx,      sigma_yy"
        do ie = 1, nene
            write(iu, fmt="(ES14.6,4(A,ES14.6))") energy_list(ie),             &
                ",", aimag(conductivity_tensor(1, 1, ie)),                     &
                ",", aimag(conductivity_tensor(1, 2, ie)),                     &
                ",", aimag(conductivity_tensor(2, 1, ie)),                     &
                ",", aimag(conductivity_tensor(2, 2, ie))
        enddo

        close(iu)
    end subroutine write_conductivity_imag
end module write_files
