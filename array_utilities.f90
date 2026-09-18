module array_utilities
use kinds, only: dp
implicit none
private
public :: arange, linspace, cross_product
contains
!******************************************************************************
    pure function arange(first, last, step) result(answer)
        implicit none
        ! This is the equivalent of numpy.arange
        real(dp), intent(in)  :: first, last, step

        ! Local variables
        real(dp) :: test_val
        integer  :: num_steps, i

        real(dp), allocatable :: answer(:)
       
        test_val=(last-first)/step

        if(test_val .le. 0.0_dp) then
            error stop "FATAL ERROR: impossible to use arange with these values"
        else
            num_steps = int(test_val)
            allocate(answer(num_steps))

            do i = 0, num_steps - 1
                answer(i + 1) = first + (i * step)
            enddo
        endif
    end function arange
!******************************************************************************
    pure function linspace(first, last, num_steps, include_endpoint) result(answer)
        implicit none
        ! This is similar to numpy.linspace. Includes the endpoint by default.
        real(dp), intent(in)          :: first, last
        integer, intent(in)           :: num_steps
        logical, intent(in), optional :: include_endpoint

        ! Local variables
        logical  :: use_endpoint
        real(dp) :: step_size
        integer  :: i

        real(dp), allocatable         :: answer(:)

        if (present(include_endpoint)) then
            use_endpoint = include_endpoint
        else
            use_endpoint = .true.
        endif

        allocate(answer(num_steps))
        if (use_endpoint) then
            step_size = (last - first) / real(num_steps - 1, kind=dp)
        else
            step_size = (last - first) / real(num_steps, kind=dp)
        endif

        do i = 0, num_steps - 1
            answer(i + 1) = first + (i * step_size)
        enddo
    end function linspace
!******************************************************************************
    function cross_product(vec_1, vec_2) result(vec_cross)
        implicit none
        real(dp), intent(in) :: vec_1(3), vec_2(3)
        real(dp) :: vec_cross(3)

        vec_cross(1) = (vec_1(2) * vec_2(3)) - (vec_1(3) * vec_2(2))
        vec_cross(2) = (vec_1(3) * vec_2(1)) - (vec_1(1) * vec_2(3))
        vec_cross(3) = (vec_1(1) * vec_2(2)) - (vec_1(2) * vec_2(1))
    end function cross_product
end module array_utilities
