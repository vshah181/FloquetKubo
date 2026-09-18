module read_files
use kinds, only: dp
use, intrinsic :: iso_fortran_env, only: iostat_end
implicit none
private
public :: read_input, read_hr, read_nnkp, read_kpoints, read_vector_potential
contains
!******************************************************************************
    subroutine read_input(seedname, basis, soc, nlayers, energy_range,         &
        energy_step, broadening_factor, direction, fermi_energy)
    implicit none
        character(len=99), intent(out) :: seedname
        character(len=4),  intent(out) :: basis
        integer,           intent(out) :: nlayers, direction
        real(dp),          intent(out) :: energy_range(2), energy_step
        real(dp),          intent(out) :: broadening_factor, fermi_energy
        logical,           intent(out) :: soc

        character(len=99) :: label, ival, line, temp_line
        integer :: i, eof, inp_unit

        open(newunit=inp_unit, file='INPUT', status='old', action='read', &
            iostat=eof)
        do while(eof .ne. iostat_end)
            read(inp_unit, '(a)', iostat=eof) line
            temp_line=adjustl(line)
            i=index(temp_line, ' ')
            label=temp_line(1:i)
            ival=temp_line(1+i:99)
            if(trim(adjustl(label)) .eq. 'seedname') then
                seedname=trim(adjustl(ival))
            else if(trim(adjustl(label)) .eq. 'basis') then
                read(ival, *) basis
            else if(trim(adjustl(label)) .eq. 'soc') then
                read(ival, *) soc
            else if(trim(adjustl(label)) .eq. 'nlayers') then
                read(ival, *) nlayers
            else if(trim(adjustl(label)) .eq. 'energy_range') then
                read(ival, *) energy_range(1), energy_range(2)
            else if(trim(adjustl(label)) .eq. 'energy_step') then
                read(ival, *) energy_step
            else if(trim(adjustl(label)) .eq. 'e_fermi') then
                read(ival, *) fermi_energy
            else if(trim(adjustl(label)) .eq. 'broadening_factor') then
                read(ival, *) broadening_factor
            else if(trim(adjustl(label)) .eq. 'direction') then
                read(ival, *) direction 
            endif
        enddo
        close(inp_unit)
    end subroutine read_input
!******************************************************************************
    subroutine read_hr(seedname, num_bands, num_r_pts, r_list, r_ham_list,     &
        weights)
    use constants, only: cmplx_0
    implicit none
        character(len=99), intent(in) :: seedname
        integer, intent(out) :: num_bands, num_r_pts
        integer, allocatable, intent(out) :: r_list(:,:), weights(:)
        complex(dp), allocatable, intent(out) :: r_ham_list(:,:,:)
        
        integer :: hi_row, hi_col, ir, o_i, o_j, inp_unit
        real(dp) :: rp, ip

        open(newunit=inp_unit, file=trim(adjustl(seedname))//"_hr.dat", &
            status='old', action='read')
        read(inp_unit, *)
        read(inp_unit, *) num_bands
        read(inp_unit, *) num_r_pts
        allocate(r_list(num_r_pts, 3), r_ham_list(num_r_pts, num_bands,        &
            num_bands), weights(num_r_pts))
        r_ham_list = cmplx_0
        r_list = 0
        weights = 1
        read(inp_unit, *) weights
        do ir=1, num_r_pts
            do o_i=1, num_bands
                do o_j=1, num_bands
                    read(inp_unit, *)r_list(ir, 1), r_list(ir, 2),             &
                        r_list(ir, 3), hi_row, hi_col, rp, ip
                    r_ham_list(ir, hi_row, hi_col)=cmplx(rp, ip, kind=dp)  
                enddo
            enddo
        enddo
        close(inp_unit)
    end subroutine read_hr
!******************************************************************************
    subroutine read_nnkp(seedname, avec, bvec, num_bands, basis,               &
        projection_centres)
    implicit none
        character(len=99), intent(in) :: seedname
        character(len=4), intent(in) :: basis
        integer, intent(in) :: num_bands
        real(dp), intent(out) :: avec(3,3), bvec(3,3),                         &
            projection_centres(num_bands, 3)
        integer :: eof, i, num_proj, inp_unit
        character(len=99) line, temp_line

        open(newunit=inp_unit, file=trim(adjustl(seedname))//'.nnkp',          &
            iostat=eof, status='old', action='read')

        do while(eof .ne. iostat_end)
            read(inp_unit, '(a)', iostat=eof) line
            if (trim(adjustl(line)) .eq. 'begin real_lattice') then
                do i=1, 3
                    read(inp_unit, *, iostat=eof) avec(i, :)
                enddo
            else if (trim(adjustl(line)) .eq. 'begin recip_lattice') then
                do i=1, 3
                    read(inp_unit, *, iostat=eof) bvec(i, :)
                enddo
            else if (trim(adjustl(line)) .eq. 'begin projections') then
                read(inp_unit, *) num_proj
                do i=1, num_proj
                    read(inp_unit, *) projection_centres(i, :)
                    read(inp_unit, *) temp_line
                enddo
                if ((num_proj*2 .eq. num_bands) .and. (basis.eq.'uudd')) then
                    do i=1, num_proj
                        projection_centres(i+num_proj, :)                      &
                            =projection_centres(i, :)
                    enddo
                endif
            else if (trim(adjustl(line)) .eq. "begin spinor_projections") then
                read(inp_unit, *) num_proj
                do i=1, num_proj
                    read(inp_unit, *) projection_centres(i,:)
                    read(inp_unit, *) temp_line
                    read(inp_unit, *) temp_line
                enddo
            endif
        enddo
        close(inp_unit)
    end subroutine read_nnkp
!******************************************************************************
    subroutine read_kpoints(nk, k_shift, k_frac)
        integer,  intent(out) :: nk(3)
        real(dp), intent(out) :: k_shift(3), k_frac(3)

        character(len=99) :: label, ival, line, temp_line
        integer :: eof, i, inp_unit

        k_shift = 0.0_dp
        k_frac  = 1.0_dp
        
        open(newunit=inp_unit, file='KPOINTS', status='old', action='read',    &
            iostat=eof)
        do while(eof .ne. iostat_end) 
            read(inp_unit, '(a)', iostat=eof) line
            temp_line=adjustl(line)
            i=index(temp_line, ' ')
            label=temp_line(1:i)
            ival=temp_line(1+i:99)
            if(adjustl(trim(label)) .eq. 'NK1') then
                read(ival, *) nk(1)
            else if(adjustl(trim(label)) .eq. 'NK2') then
                read(ival, *) nk(2)
            else if(adjustl(trim(label)) .eq. 'NK3') then
                read(ival, *) nk(3)
            else if(adjustl(trim(label)) .eq. 'K1_SHIFT') then
                read(ival, *) k_shift(1)
            else if(adjustl(trim(label)) .eq. 'K2_SHIFT') then
                read(ival, *) k_shift(2)
            else if(adjustl(trim(label)) .eq. 'K3_SHIFT') then
                read(ival, *) k_shift(3)
            else if(adjustl(trim(label)) .eq. 'K1_FRAC') then
                read(ival, *) k_frac(1)
            else if(adjustl(trim(label)) .eq. 'K2_FRAC') then
                read(ival, *) k_frac(2)
            else if(adjustl(trim(label)) .eq. 'K3_FRAC') then
                read(ival, *) k_frac(3)
            endif
        enddo
        close(inp_unit)
    end subroutine read_kpoints
!******************************************************************************
    subroutine read_vector_potential(num_photon, omega, phase_shift, a_0)
    use constants, only : hbar_ev=>reduced_planck_constant_ev, deg2rad
    implicit none
        integer, intent(out) :: num_photon
        real(dp), intent(out) :: omega, phase_shift, a_0
        integer :: eof, i, inp_unit
        real (dp) :: s
        character(len=99) :: label, ival, line, temp_line

        open(newunit=inp_unit, file='vector_potential.dat', iostat=eof)
        do while(eof .ne. iostat_end)
            read(inp_unit, '(a)', iostat=eof) line
            temp_line=adjustl(line)
            i=index(temp_line, ' ')
            label=temp_line(1:i)
            ival=temp_line(1+i:99)
            if(trim(adjustl(label)) .eq. 'phase_shift(deg)') then
                read(ival, *) phase_shift
            else if(trim(adjustl(label)) .eq. 's') then
                read(ival, *) s
            else if(trim(adjustl(label)) .eq. 'photon_energy') then
                read(ival, *) omega
            else if(trim(adjustl(label)) .eq. 'num_photons') then
                read(ival, *) num_photon
            endif
        enddo
        close(inp_unit)
        omega = omega / hbar_ev ! to get it in hertz
        phase_shift = phase_shift * deg2rad ! to get it in radians
        a_0 = s
        ! Ignore hbar, c and elementary charge.

    end subroutine read_vector_potential
end module read_files
