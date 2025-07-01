module lanczos
use utils
use csr
use spmv
use spmm_m

implicit none

contains

    subroutine lanczos_bidiag(A, U, V, B, k)

        type(csr_mat), intent(in) :: A
        integer, intent(in) :: k
        real(dp), intent(inout), allocatable, dimension(:, :) :: U, V, B
        integer :: m, n, nnz, m1, m2, n1, n2, m1n1, m2n2
        integer :: j
        real(dp) :: beta, alpha
        real(dp), allocatable :: p(:)


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
        allocate(V(m2n2, k))
        allocate(B(2, k)) ! First row stores betas, second row stores alphas
        allocate(p(m2n2))

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

            if (j == k) then
                call spmv_permuted(A, U(:, j), p, 1)
                p = p - alpha * V(:, j)
                beta = norm2(p)
            else
                call spmv_permuted(A, U(:, j), V(:, j + 1), 1)
                V(:, j + 1) = V(:, j + 1) - alpha * V(:, j)
                beta = norm2(V(:, j + 1))
            end if

            B(1, j) = beta
            B(2, j) = alpha

        end do

    end subroutine 


    function measure_orthogonality(Q) result (err)
        real(dp), intent(in) :: Q(:,:)
        real(dp) :: err

        real(dp) :: nrm, nrmg
        real(dp) :: alpha, beta

        real(dp), allocatable :: I_tilde(:,:)

        integer :: m, n, i, j

        m = size(Q, 1)
        n = size(Q, 2)
        alpha = 1.0
        beta = 0.0

        allocate(I_tilde(n,n))

        call DGEMM('T', 'N', n, n, m, alpha, Q, m, Q, m, beta, I_tilde, n)

        nrm = norm2(Q)
        nrm = nrm * nrm

        nrmg = norm2(I_tilde)
        nrmg = nrmg * nrmg

        err = n - 2 * nrm + nrmg 
    end function


    function measure_bidiag(A, U, V) result(err)
        type(csr_mat), intent(in) :: A
        real(dp), intent(in) :: U(:, :), V(:, :)
        real(dp) :: err
        real(dp), allocatable :: W(:,:), B(:,:)
        real(dp) :: alpha = 1.0
        real(dp) :: beta = 0.0

        integer :: i, j

        ! W = AV
        call spmm_permuted(A, V, W, 0)

        allocate(B(size(U, 2), size(W, 2)))

        ! B = U^TW
        call DGEMM('T', 'N', size(U, 2), size(W, 2), size(U, 1), alpha, U, size(U, 1), W, size(W, 1), beta, B, size(U, 2))


        err = 0.0
        do i=1, size(B, 1)
            do j=1, size(B, 2)
                if (((j - i) > 1).or.(i > j)) then
                    err = err + abs(B(i,j))
                end if
            end do
        end do


    end function

end module lanczos
