module kinds
use, intrinsic :: iso_fortran_env, only: real64, real32, int64
implicit none
private
public :: sp, dp, ilp
    integer, parameter :: sp = real32
    integer, parameter :: dp = real64
    integer, parameter :: ilp = int64
end module kinds
