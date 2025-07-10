module csp
use utils
use csr
implicit none

    type :: csp_mat
        real(dp), pointer :: vals(:)
        integer, pointer :: colinds(:), rowptrs(:)
        integer, pointer :: nzc(:)
        integer :: m, n, nnz, m1, m2, n1, n2
    contains
        procedure :: map_block
    end type

contains


    subroutine dist_csp_static(A_csp, A_csp_dist)
        class(csp_mat), intent(in) :: A_csp
        class(csp_mat), intent(inout) :: A_csp_dist
        integer :: i, j, iid, myim, my_nnz, offset, row_len, my_rows
        integer, allocatable :: rowcounts(:), colnz(:)

        myim = this_image()

        A_csp_dist%m1 = A_csp%m1
        A_csp_dist%n1 = A_csp%n1
        A_csp_dist%m2 = A_csp%m2
        A_csp_dist%n2 = A_csp%n2
        A_csp_dist%m = A_csp%m / num_images()
        A_csp_dist%n = A_csp%n

        ! Count nnz for my image
        my_nnz = 0 
        do i = 1, A_csp%m
            iid = (i - 1) / A_csp_dist%m + 1
            do j = A_csp%rowptrs(i), A_csp%rowptrs(i+1)
                if (iid == this_image()) then
                    my_nnz = my_nnz + 1
                end if
            end do
        end do

        A_csp_dist%nnz = my_nnz

        ! Allocate 
        allocate(A_csp_dist%vals(my_nnz))
        allocate(A_csp_dist%colinds(my_nnz))
        allocate(A_csp_dist%rowptrs(A_csp_dist%m + 1))
        allocate(rowcounts(A_csp_dist%m))
        allocate(colnz(A_csp_dist%n))
        rowcounts = 0
        colnz = 0

        ! Copy data
        offset = 1
        my_rows = 1
        do i = 1, A_csp%m

            iid = (i - 1) / A_csp_dist%m + 1

            row_len = A_csp%rowptrs(i + 1) - A_csp%rowptrs(i)

            if (iid == this_image()) then

                A_csp_dist%vals(offset:offset+row_len-1) = A_csp%vals(A_csp%rowptrs(i) + 1:A_csp%rowptrs(i+1))
                A_csp_dist%colinds(offset:offset+row_len-1) = A_csp%colinds(A_csp%rowptrs(i) + 1:A_csp%rowptrs(i+1))

                do j = A_csp%rowptrs(i)+1, A_csp%rowptrs(i + 1)
                    colnz(A_csp%colinds(j)) = 1
                end do

                rowcounts(my_rows) = row_len

                offset = offset + row_len
                my_rows = my_rows + 1

            end if

        end do

        call prefix_sum_integer(rowcounts, A_csp_dist%rowptrs)

        allocate(A_csp_dist%nzc(sum(colnz)))

        offset = 1
        do j = 1, A_csp_dist%n
            if (colnz(j) == 1) then
                A_csp_dist%nzc(offset) = j
                offset = offset + 1
            end if
        end do

    end subroutine


    function map_block(self, i, j) result(idx)
        class(csp_mat), intent(in) :: self
        integer, intent(in) :: i, j
        integer :: idx(2)
        
        idx(1) = ((i - 1) / self%m2) + 1
        idx(2) = ((j - 1) / self%n2) + 1

    end function


    function csr_to_csp(A_csr) result(A_csp)

        type(csr_mat), intent(in) :: A_csr
        type(csp_mat) :: A_csp

        real(dp), allocatable :: vals(:)
        integer, allocatable :: inds(:,:)
        integer :: i, j, k, l
        integer :: colidx, offset
        integer :: block_idx(2)
        integer :: idx

        ! Determine number of blocks and block size
        A_csp%m1 = A_csr%m1
        A_csp%n1 = A_csr%n1
        A_csp%m2 = A_csr%m2
        A_csp%n2 = A_csr%n2
        A_csp%m = A_csr%m1 * A_csr%n1
        A_csp%n = A_csr%m2 * A_csr%n2
        A_csp%nnz = A_csr%nnz

        allocate(A_csp%vals(A_csp%nnz))
        allocate(A_csp%colinds(A_csp%nnz))
        allocate(A_csp%rowptrs(A_csp%m + 1))

        allocate(vals(A_csp%nnz))
        allocate(inds(2, A_csp%nnz))

        do k = 1, A_csr%m
            do l = A_csr%rowptrs(k) + 1, A_csr%rowptrs(k + 1)
                
                colidx = A_csr%colinds(l)

                block_idx = A_csp%map_block(k, colidx)

                vals(l) = A_csr%vals(l)

                inds(1, l) = (block_idx(1) - 1) + (block_idx(2) - 1) * A_csp%m1 + 1
                inds(2, l) = mod((k - 1), A_csp%m2) + mod((l - 1), A_csp%n2) * A_csp%m2 + 1

            end do
        end do

        call ijv_to_csp(vals, inds, A_csp)
                
    end function


    subroutine ijv_to_csp(vals, inds, A)
        real(dp), intent(inout) :: vals(:)
        integer, intent(inout) :: inds(:,:)
        type(csp_mat), intent(inout) :: A

        integer, allocatable :: sorted_inds(:), rowcounts(:)
        integer :: i

        allocate(sorted_inds(A%nnz))
        allocate(rowcounts(A%m))
        rowcounts = 0

        call sort_index(inds(1, :), sorted_inds)

        inds(2, :) = inds(2, sorted_inds(:))
        vals(:) = vals(sorted_inds(:))

        do i = 1, A%nnz
            rowcounts(inds(1, i)) = rowcounts(inds(1, i)) + 1
            A%vals(i) = vals(i)
            A%colinds(i) = inds(2, i)
        end do

        call prefix_sum_integer(rowcounts, A%rowptrs)

    end subroutine




    subroutine free_csp_mat(A)
        type(csp_mat), intent(inout) :: A
        integer :: i

        deallocate(A%vals)
        deallocate(A%colinds)
        deallocate(A%rowptrs)

    end subroutine

end module csp
