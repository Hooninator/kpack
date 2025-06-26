module lanczos
use utils
use csr
use spmv

implicit none

contains

    subroutine lanczos_bidiag(A, U, V, B, k)

        type(csr_mat), intent(in) :: A
        integer, intent(in) :: k
        real(dp), intent(inout), allocatable, dimension(:, :) :: U, V, B
        integer :: m, n, nnz, m1, m2, n1, n2, m1n1, m2n2
        integer :: j
        real(dp) :: beta, alpha


        m = A%m
        n = A%n
        nnz = A%nnz
        m1 = A%m1
        m2 = A%m2
        n1 = A%n1
        n2 = A%n2

        m1n1 = m1 * n1
        m2n2 = m2 * n2

        allocate(U(m1n1, k))
        allocate(V(m2n2, k + 1))
        allocate(B(2, k)) ! First row stores betas, second row stores alphas

        U = 0
        V = 0
        B = 0

        call rand_vec(V(:, 1))
        V(:, 1) = V(:, 1) / norm2(V(:, 1))

        beta = 1
        do j=1, k

            if (beta == 0) exit

            V(:, j) = V(:, j) / beta

            call spmv_permuted(A, V(:, j), U(:, j), 0)

            if (j > 1) then
                U(:, j) = U(:, j) - beta * U(:, j - 1)
            end if 

            alpha = norm2(U(:, j))
            U(:, j) = U(:, j) / alpha

            call spmv_permuted(A, U(:, j), V(:, j + 1), 1)
            V(:, j + 1) = V(:, j + 1) - alpha * V(:, j)
            beta = norm2(V(:, j + 1))

            B(1, j) = beta
            B(2, j) = alpha

        end do

    end subroutine 


    function measure_orthogonality(Q) result (err)
        real(dp), intent(in) :: Q(:,:)
        real(dp) :: err

        real(dp), allocatable :: I_tilde(:,:)

        integer :: m, n, i, j

        m = size(Q, 1)
        n = size(Q, 2)

        allocate(I_tilde(m,m))

        call DGEMM('N', 'T', m, m, n, 1.0, Q, m, Q, m, 0.0, I_tilde, m)

        err = 0
        do i=1, m
            do j = 1, m
                if (i == j) then
                    err = err + (I_tilde(i,i) - 1)**2
                else
                    err = err + (I_tilde(i,i))**2
                end if
            end do
        end do

        err = sqrt(err)
    end function

end module lanczos
