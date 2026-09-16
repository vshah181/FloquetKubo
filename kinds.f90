module kinds
use, intrinsic :: iso_fortran_env, only: real64, real32
implicit none
private
public :: sp, dp
    integer, parameter :: sp = real32
    integer, parameter :: dp = real64
end module kinds
