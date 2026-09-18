module gauge_transformation
use kinds, only: dp
implicit none
private
public :: transform_velocity
contains
    pure subroutine transform_velocity(n_bands, rotation_matrix, velocities)
    use matrix_utilities, only: similarity_transform
    implicit none
        integer,        intent(in) :: n_bands
        complex(dp),    intent(in) :: rotation_matrix(n_bands, n_bands)
        complex(dp), intent(inout) :: velocities(n_bands, n_bands, 2)

        complex(dp)                :: velocities_temp(n_bands, n_bands, 2)
        integer                    :: x

        velocities_temp = velocities


        do x = 1, 2
            associate(vx=>velocities_temp(:, :, x), u=>rotation_matrix)
                velocities(:, :, x) = similarity_transform(vx, u)
            end associate
        enddo

    end subroutine transform_velocity
end module gauge_transformation
