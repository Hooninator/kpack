module utils
use, intrinsic :: iso_fortran_env, dp=>real64
implicit none

contains

    subroutine prefix_sum_integer(A, S)
        
        integer, intent(in) :: A(:)
        integer, pointer, intent(out) :: S(:)
        integer :: i

        allocate(S(size(A) + 1))

        S(1) = 0
        do i=2, size(A) + 1
            S(i) = A(i - 1) + S(i - 1)
        end do

    end subroutine prefix_sum_integer


    subroutine rand_vec(x)
        real(dp), intent(inout) :: x(:)

        integer :: i

        do i = 1,size(x)
            call random_number(x(i))
        end do

    end subroutine rand_vec


    subroutine rand_mat(X)
        real(dp), intent(inout) :: X(:, :)
        integer :: i, j

        do i = 1, size(X,1)
            do j = 1, size(X,2)
                call random_number(X(i,j))
            end do
        end do
    end subroutine


    subroutine check_lapack(info)
        integer :: info
        if (info /= 0) then
            print*, "LAPACK function at line ",__LINE__, " gave info ", info
            call ABORT()
        end if
    end subroutine



    subroutine print_mat(X)
        real(dp), intent(in) :: X(:,:)
        integer :: i, j


        do j=1, size(X, 2)*15
            write(*, '(A)', advance='no') "-"
        end do

        do i = 1, size(X, 1)
            print *
            do j = 1, size(X, 2)
                write(*, "(A)", advance='no') "|"
                write(*, '(F12.6)', advance='no') X(i,j)
                write(*, "(A)", advance='no') "|"
            end do
            print *
            do j=1, size(X, 2)*15
                write(*, '(A)', advance='no') "-"
            end do
        end do
        print *

    end subroutine

end module utils
