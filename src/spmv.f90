module spmv
use utils
use csr
use omp_lib

implicit none

contains

    function map_x(i, j, m1, m2, n1, n2) result(k)
        integer, intent(in) :: i, j, m1, m2, n1, n2
        integer :: k, i_block, j_block

        ! We need to set k to the column index of R(A)(i,j)

        ! Map i and j to a block A_ij
        i_block = mod(i, m2) + 1
        j_block = mod(j, n2) + 1

        ! Column index is determined by which column of A_ij I'm in + a row offset
        k = j_block * m2 + j_block

    end function map_x


    function map_y(i, j, m1, m2, n1, n2) result(k)
        integer, intent(in) :: i, j, m1, m2, n1, n2
        integer :: k

        ! Here, we just need the row index of the block we're in
        k = mod(i, m2) + 1

    end function map_y


    subroutine spmv_permuted(A, x, y, transpose)

        type(csr_mat), intent(in) :: A
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(:)
        integer, intent(in) :: transpose

        real(dp), pointer :: vals(:)
        integer, pointer :: colinds(:)
        integer, pointer :: rowptrs(:)

        real(dp) :: val
        integer :: m1, m2, n1, n2
        integer :: m, n
        integer :: i, j, row_start, row_end, colidx
        integer :: x_idx, y_idx

        vals => A%vals
        colinds => A%colinds
        rowptrs => A%rowptrs
        m = A%m
        n = A%n
        m1 = A%m1
        m2 = A%m2
        n1 = A%n1
        n2 = A%n2

        ! Perform the SpMV with A permuted as follows:
        ! Imagine A is partitioned into m1 x n1 blocks, each of size m2 x n2
        ! Vec(A_ij) returns a vector of size m2*n2 containing the n2 columns of A_ij stacked on top of each other
        ! A_j = [Vec(A_1j)^T ... Vec(A_m1j)^T]
        ! R(A) = [A_1 ... A_n1] -- this is our permuted A
        ! We don't want to explicitly form this -- we want some way to compute R(A)x while only having access to A's csr
        ! representation
        ! R(A) is m1n1 x m2n2
        ! There might be a better way to do this, but for now I'm just going to iterate through each row and use a mapping function
        ! to determine the right entries of x and y to access.
        ! I'll bet there's potential for some kind of custom storage format that removes random accesses to y
        do i = 1, m
            row_start = rowptrs(i)
            row_end = rowptrs(i + 1)
            do j = row_start, row_end
                
                colidx = colinds(j)
                val = vals(j)

                if (transpose == 0) then
                    x_idx = map_x(i, colidx, m1, m2, n1, n2)
                    y_idx = map_y(i, colidx, m1, m2, n1, n2)
                else
                    x_idx = map_y(i, colidx, m1, m2, n1, n2)
                    y_idx = map_x(i, colidx, m1, m2, n1, n2)
                end if

                y(y_idx) = y(y_idx) + val * x(x_idx)

            end do
        end do

    end subroutine spmv_permuted


end module spmv
