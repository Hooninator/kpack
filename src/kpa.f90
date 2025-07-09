module kpa
use utils
use csr
use lanczos
use power_method

implicit none


contains

    subroutine kpa_svd(A, Q, C, sing, maxiters)

        type(csr_mat), intent(in) :: A
        real(dp), intent(out), allocatable, dimension(:, :) :: Q, C
        real(dp), intent(out) :: sing
        integer, intent(in) :: maxiters

        real(dp), allocatable :: U(:,:), V(:,:), B(:,:)
        real(dp) :: work(5*(maxiters - 1))
        real(dp) :: N(maxiters, maxiters)
        real(dp) :: u_vec(maxiters), v_vec(maxiters)

        real(dp) :: ortho_err
        real(dp) :: bidiag_err

        integer :: info

        allocate(Q(A%m1, A%n1))
        allocate(C(A%m2, A%n2))

        print*, "Q : ",A%m1,"x",A%n1
        print*, "C : ",A%m2,"x",A%n2

        ! Bidiagonalize A
        print*, "Beginning lanczos..."
        call lanczos_bidiag(A, U, V, B, maxiters)
        print*, "Done with lanczos"

        ortho_err = 1 - dot_product(V(:, 1), V(:, 1))
        print*, "Orthogonality of V: ", ortho_err
        ortho_err = 1 - dot_product(U(:, 1), U(:, 1))
        print*, "Orthogonality of U: ", ortho_err

        ! Compute singular vectors and values of B
        print*, "Computing SVD..."

        call DBDSQR('U', maxiters, A%m2 * A%n2, A%m1 * A%n1, 0, B(2, 1:maxiters), B(1, 1:maxiters-1), transpose(V), maxiters, U,&
            A%m1 * A%n1, N, 1, work, info)
        call check_lapack(info)
        print*, "Done computing SVD -- Info: ", info

        sing = B(2, 1)

        ! Sanity check
        call singval_err(A, U(:, 1), V(:, 1), sing)

        ! Reshape the arrays into the approximation matrices
        Q = reshape(U(:, 1), (/A%m1, A%n1/))
        C = reshape(V(:, 1), (/A%m2, A%n2/))

    end subroutine


    subroutine kpa_als(A, Q, C, sing, maxiters)

        type(csr_mat), intent(in) :: A
        real(dp), intent(out), allocatable, dimension(:, :) :: Q, C
        real(dp), intent(out) :: sing
        integer, intent(in) :: maxiters

    end subroutine 


    subroutine singval_err(A, u, v, sing) 
        type(csr_mat), intent(in) :: A
        real(dp), intent(in) :: u(:), v(:)
        real(dp), intent(in) :: sing

        real(dp) :: neg_sing

        real(dp) :: err
        real(dp), allocatable :: y1(:), y2(:)
        real(dp), allocatable :: v2(:)
        allocate(v2(size(v)))

        allocate(y1(A%n))
        allocate(y2(A%m))
        y1 = 0
        y2 = 0

        ! U^T*A*V = S
        ! V^T * A^T * U = S
        ! V^T * A^T = SU^T
        ! AV = US
        ! Av = Us = u * sing
        ! ...
        ! A^Tu = v * sing

        neg_sing = sing * (-1)
        call spmv_permuted(A, u, y1, 1)
        call DAXPY(size(v), neg_sing, v, 1, y1, 1)
        err = norm2(y1)

        print*, "||A^Tu - sing * v||: ", err

        y1 = 0
        y2 = 0

        ! Av = sing * u
        call spmv_permuted(A, v, y2, 0)
        call DAXPY(size(u), neg_sing, u, 1, y2, 1)
        err = norm2(y2)

        print*, "||Av - sing * u||: ", err

    end subroutine


    function compute_err(A, Q, C, sing) result(err)
        type(csr_mat), intent(in) :: A
        real(dp), intent(in) :: Q(:,:), C(:,:)
        real(dp), intent(in) :: sing
        real(dp) :: err

        integer :: i, j
        integer :: colidx, Qi, Qj, Ci, Cj
        real(dp) :: kprod_val

        err = 0
        do i=1,A%m
            do j=A%rowptrs(i), A%rowptrs(i + 1)

                colidx = A%colinds(j)

                Qi = (i - 1) / A%m1
                Qj = (colidx - 1) / A%n1
                Ci = mod(i, A%m2) + 1
                Cj = mod(colidx, A%n2) + 1

                kprod_val = sing * Q(Qi+1, Qj+1) * C(Ci, Cj)

                err = err + ( A%vals(j) - kprod_val )**2
            end do
        end do

        err = sqrt(err) / norm2(A%vals)

    end function


    function count_nnz(X, sing) result(nnz)
        real(dp), intent(in) :: X(:,:)
        real(dp), intent(in) :: sing
        integer :: nnz
        integer :: i, j

        nnz = 0
        do i=1, size(X, 1)
            do j=1, size(X, 2)
                if (sing*X(i,j) > epsilon(X(i,j))) then
                    nnz = nnz + 1
                end if
            end do
        end do
    end function

end module kpa
