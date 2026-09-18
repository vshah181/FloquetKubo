module matrix_utilities
use kinds, only: dp
implicit none
private
public :: similarity_transform
contains
!******************************************************************************
    pure function similarity_transform(matrix, rotation) result(rotated)
    implicit none
        complex(dp), intent(in) :: matrix(:, :), rotation(:, :)

        complex(dp)             :: rotated(size(matrix, 1), size(matrix, 2))

        associate(hc=>conjg(transpose(rotation)))
            rotated = matmul(hc, matmul(matrix, rotation))
        end associate
    end function similarity_transform
end module matrix_utilities
