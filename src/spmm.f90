module spmm_m
use csr
use utils
use spmv
implicit none


contains

    subroutine spmm_permuted(A, X, Y, transpose)
        type(csr_mat), intent(in) :: A
        real(dp), intent(in) :: X(:,:)
        real(dp), allocatable, intent(out) :: Y(:,:)
        integer, intent(in) :: transpose

        integer :: m, n, nnz
        integer :: j

        m = A%m
        n = size(X, 2)


        allocate(Y(m,n))


        ! Slow but easy to program
        do j = 1, n
            call spmv_permuted(A, X(:, j), Y(:, j), transpose)
        end do

    end subroutine spmm_permuted

end module spmm_m
