module kpa
use utils
use csr
use lanczos

implicit none


contains

    subroutine kpa_svd(A, Q, C, maxiters)

        type(csr_mat), intent(in) :: A
        real(dp), intent(out), allocatable, dimension(:, :) :: Q, C
        integer, intent(in) :: maxiters

        real(dp), allocatable :: U(:,:), V(:,:), B(:,:)
        real(dp)  :: D(maxiters)
        real(dp) :: work(4*(maxiters - 1))
        real(dp) :: N(0)

        integer :: info

        allocate(Q(A%m1, A%n1))
        allocate(C(A%m2, A%n2))

        ! Bidiagonalize A
        print*, "Beginning lanczos..."
        call lanczos_bidiag(A, U, V, B, maxiters)
        print*, "Done with lanczos"

        ! Compute singular vectors and values of B
        print*, "Computing SVD..."
        call DBDSQR('U', maxiters, A%m2 * A%n2, A%m1 * A%n1, 0, D, B(1, 1:maxiters-1), transpose(V), maxiters, U,&
            A%m1 * A%n1, N, 1, work, info)
        print*, "Done computing SVD"

        ! Reshape the arrays into the approximation matrices
        Q = reshape(U(:, 1), (/A%m1, A%n1/))
        C = reshape(V(:, 1), (/A%m2, A%n2/))

    end subroutine


    function compute_err(A, Q, C) result(err)
        type(csr_mat), intent(in) :: A
        real(dp), intent(in) :: Q(:,:), C(:,:)
        real(dp) :: err

        integer :: i, j
        integer :: colidx, Qi, Qj, Ci, Cj
        real(dp) :: kprod_val

        err = 0
        do i=1,A%m
            do j=A%rowptrs(i), A%rowptrs(i + 1)

                colidx = A%colinds(j)

                Qi = i / A%m1 
                Qj = colidx / A%n1
                Ci = mod(i, A%m2) + 1
                Cj = mod(colidx, A%n2) + 1

                kprod_val = Q(Qi, Qj) * C(Ci, Cj)

                err = err + abs( A%vals(j) - kprod_val )**2
            end do
        end do

        err = err / norm2(A%vals)

    end function


    function count_nnz(X) result(nnz)
        real(dp), intent(in) :: X(:,:)
        integer :: nnz
        integer :: i, j

        nnz = 0
        do i=1, size(X, 1)
            do j=1, size(X, 2)
                if (X(i,j) > epsilon(X(i,j))) then
                    nnz = nnz + 1
                end if
            end do
        end do
    end function

end module kpa
