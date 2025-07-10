module csr
use utils
use stdlib_sparse
use stdlib_sorting
implicit none


    type :: csr_mat
        integer :: nnz
        integer :: m, m1, m2
        integer :: n, n1, n2
        real(dp), pointer :: vals(:)
        integer, pointer :: colinds(:)
        integer, pointer :: rowptrs(:)
    end type


contains

    subroutine csr_set_subdims(A)
        type(csr_mat), intent(inout) :: A

        A%m1 = floor(sqrt(real(A%m)))
        do while (mod(A%m, A%m1) /= 0)
            A%m1 = A%m1 + 1
        end do
        A%m2 = A%m / A%m1
        A%m2 = A%m2 + abs(A%m - A%m2 * A%m1)

        A%n1 = floor(sqrt(real(A%n)))
        do while (mod(A%n, A%n1) /= 0)
            A%n1 = A%n1 + 1
        end do
        A%n2 = A%n / A%n1
        A%n2 = A%n2 + abs(A%n - A%n2 * A%n1)

    end subroutine


    function csr_to_stdlib_coo(A_csr) result(A_coo)
        type(csr_mat), intent(in) :: A_csr
        type(COO_dp_type) :: A_coo

        integer :: nnz, i, j
        nnz = A_csr%nnz
        allocate(A_coo%data(nnz))
        allocate(A_coo%index(2, nnz))

        do i = 1, A_csr%m
            do j = A_csr%rowptrs(i), A_csr%rowptrs(i+1)
                A_coo%data(j) = A_csr%vals(j)
                A_coo%index(1, j) = i
                A_coo%index(2, j) = A_csr%colinds(j)
            end do
        end do

    end function


    function ijv_to_csr(vals, inds, m, n) result(A)
        real(dp), intent(inout) :: vals(:)
        integer, intent(inout) :: inds(:,:)
        type(csr_mat) :: A
        integer, intent(in) :: m, n

        integer :: nnz
        integer :: i, j
        integer, allocatable :: sorted_inds(:), rowcounts(:)

        nnz = size(vals,1)

        A%m = m
        A%n = n
        A%nnz = nnz

        call csr_set_subdims(A)

        allocate(A%vals(nnz))
        allocate(A%colinds(nnz))
        allocate(A%rowptrs(m + 1))

        allocate(sorted_inds(nnz))
        allocate(rowcounts(m))
        rowcounts=0

        call sort_index(inds(1, :), sorted_inds)

        inds(2, :) = inds(2, sorted_inds(:))
        vals(:) = vals(sorted_inds(:))

        do i = 1, nnz
            rowcounts(inds(1, i)) = rowcounts(inds(1, i)) + 1
            A%vals(i) = vals(i)
            A%colinds(i) = inds(2, i)
        end do

        call prefix_sum_integer(rowcounts, A%rowptrs)
        A%rowptrs = 1

    end function


    subroutine csr_deepcopy(src, dst)
        type(csr_mat), intent(in) :: src
        type(csr_mat), intent(out) :: dst

        dst%m = src%m
        dst%m1 = src%m1
        dst%m2 = src%m2
        dst%n = src%n
        dst%n1 = src%n1
        dst%n2 = src%n2
        dst%nnz = src%nnz

        allocate(dst%vals(src%nnz))
        dst%vals(:) = src%vals(:)

        allocate(dst%colinds(src%nnz))
        dst%colinds(:) = src%colinds(:)

        allocate(dst%rowptrs(src%m + 1))
        dst%rowptrs(:) = src%rowptrs(:)

    end subroutine

    subroutine free_csr_mat(A)
        type(csr_mat), intent(inout) :: A
        deallocate(A%vals)
        deallocate(A%colinds)
        deallocate(A%rowptrs)
    end subroutine 
end module csr
        
