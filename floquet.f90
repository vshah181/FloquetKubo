module floquet
use kinds, only: dp
implicit none
private
public :: make_floquet_hamiltonian
contains
    subroutine make_floquet_hamiltonian_real_space(new_hr_list)
    use parameters, only: num_photon, nf_bands, omega, num_bands, num_r_pts,   &
            r_list, r_ham_list
    use constants, only: hbar=>reduced_planck_constant_ev, cmplx_0
    implicit none
        complex(dp), intent(out) :: new_hr_list(num_r_pts, nf_bands, nf_bands)
        complex(dp)              :: new_r_ham(nf_bands, nf_bands) 
        complex(dp)              :: ham_m(num_bands, num_bands)
        real(dp)                 :: frequency_array(nf_bands), hw, min_eps
        integer                  :: im, m, irow, icol, ir, i, ire, ice

        hw = omega * hbar
        frequency_array = 0.0_dp
        im = num_photon
        min_eps = tiny(1.0_dp)
        
        if (nf_bands /= (1 + (2 * num_photon)) * num_bands) then
            error stop "FATAL: Invalid number of floquet side bands!"
        endif

        do i = 1, nf_bands, num_bands
            frequency_array(i:i + num_bands - 1) = -real(im, kind=dp) * hw
            im = im - 1
        enddo

        do ir=1, num_r_pts
            if (hw .gt. min_eps) then
                new_r_ham = cmplx_0
                do m = -num_photon, num_photon
                    irow = (num_bands * (abs(m) - m) / 2) + 1  ! 1 if +ve m, else (|m|*nbands)+1
                    icol = (num_bands * (abs(m) + m) / 2) + 1  ! 1 if -ve m, else (|m|*nbands)+1
                    associate(hr=>r_ham_list(ir, :, :), r=>r_list(ir, :))
                        ham_m = fourier_coefficient(num_bands, hr, m, 0.0_dp, r)
                    end associate
                    do i=1, (1 + 2 * num_photon) - abs(m)
                        ire = irow + num_bands - 1
                        ice = icol + num_bands - 1
                        new_r_ham(irow:ire, icol:ice) = ham_m 
                        irow = irow + num_bands   ! this loop fills in the 
                        icol = icol + num_bands   ! block-diagonals of H_F(r)
                    enddo
                enddo
                if (all(r_list(ir, :) .eq. (/0, 0, 0/))) then
                    do i = 1, nf_bands
                        new_r_ham(i, i) = new_r_ham(i, i) + frequency_array(i)
                    enddo
                endif
                new_hr_list(ir, :, :) = new_r_ham
            else
                new_r_ham = cmplx_0
                do m = 0, 2 * num_photon
                    irow = 1 + (m * num_bands)
                    icol = 1 + (m * num_bands)

                    ire = irow + (num_bands - 1)
                    ice = icol + (num_bands - 1)

                    ham_m = r_ham_list(ir, :, :)
                    
                    new_r_ham(irow:ire, icol:ice) = ham_m
                enddo
                new_hr_list(ir, :, :) = new_r_ham
            endif
        enddo
    end subroutine make_floquet_hamiltonian_real_space

    function fourier_coefficient(nbands, r_ham, m, t0, r) result(ham_m)
        use constants, only: two_pi, cmplx_i
        use parameters, only: omega
        implicit none
            integer, parameter       :: n = 500

            integer, intent(in)      :: nbands, m, r(3)
            real (dp), intent(in)    :: t0
            complex (dp), intent(in) :: r_ham(nbands, nbands)

            integer                  :: i
            real (dp)                :: period, t_end, t_i, t_step

            complex (dp)             :: ham_m(nbands, nbands), imw

            period = two_pi / omega
            t_end = t0 + period
            t_step = period / real(n, kind=dp)
            imw = real(m, kind=dp) * omega * cmplx_i
            ! Integrate using the trapezium rule between t0 and t0+(2pi/omega)
            ham_m = 0.5_dp * (t_ham(r, t0, r_ham) / exp(t0 * imw)              &
                           +  t_ham(r, t_end, r_ham) / exp(t_end * imw))
            do i = 1, n - 1
                t_i = t0 + (i * t_step)
                ham_m = ham_m + (t_ham(r, t_i, r_ham) / exp(t_i * imw))
            enddo
            
            ham_m = ham_m / real(n, kind=dp)
    end function fourier_coefficient

    function t_ham(r, time, r_ham) result(ham_t)
    use parameters, only: projection_centres, num_bands, a_0, omega,           &
        phase_shift, avec
    use constants, only: cmplx_i
    implicit none
        integer, intent(in)     :: r(3)
        real(dp), intent(in)    :: time
        complex(dp), intent(in) :: r_ham(num_bands, num_bands)

        integer                 ::  i, j, k
        real (dp)               ::  phase, vector_potential(3), r_real(3)

        complex(dp)             :: ham_t(num_bands, num_bands)

!==============================================================================!
!                                  See HUUS19                                  !
!                phase_shift = 0.0 => anticlockwise polarisation               !
!                  phase_shift = pi => clockwise polarisation                  !
!==============================================================================!
        vector_potential(1) = a_0 * cos(time * omega)
        vector_potential(2) = a_0 * sin(phase_shift + (time * omega))
        vector_potential(3) = 0.0_dp

        ham_t=r_ham
        do i=1, num_bands
            do j=1, num_bands
                r_real=0.0_dp
                do k=1, 3
                    r_real = r_real + ((r(k) + projection_centres(j, k)        &
                           - projection_centres(i, k)) * avec(k, :))
                enddo
                phase = dot_product(vector_potential, r_real)
                ham_t(i, j) = ham_t(i, j) * exp(cmplx_i * phase)
            enddo
        enddo
    end function t_ham
end module floquet
