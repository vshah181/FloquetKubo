module lapack_interfaces
use, intrinsic :: iso_fortran_env, only: dp=>real64
implicit none
private
public :: ZHEEVD
    interface
        subroutine ZHEEVD(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork,&
                iwork, liwork, info)

            import :: dp

            character(len=1), intent(in) :: jobz, uplo
            integer,          intent(in) :: n, lda, lwork, lrwork, liwork
            complex(dp),   intent(inout) :: a(lda, *)
            real(dp),        intent(out) :: w(*), rwork(*)
            complex(dp),     intent(out) :: work(*)
            integer,         intent(out) :: iwork(*), info
        end subroutine ZHEEVD
    end interface
end module lapack_interfaces
