module csp
use utils
use csr
implicit none

    type :: csp_mat
        integer :: nblocks
        integer :: bm, bn, m, n, nnz, bdimm, bdimn
        integer, pointer :: block_inds(:,:)
        type(csr_mat), pointer :: blocks(:)
    contains
        procedure :: map_block
    end type

contains


    subroutine dist_csp_static(A_csp, A_csp_dist)
        class(csp_mat), intent(in) :: A_csp
        class(csp_mat), intent(inout) :: A_csp_dist
        integer :: i, j, iid, myim, image_blockrows, my_nblocks, offset

        myim = this_image()

        A_csp_dist%bm = A_csp%bm
        A_csp_dist%bn = A_csp%bn
        A_csp_dist%bdimm = A_csp%m / A_csp%bm
        A_csp_dist%bdimn = A_csp%n / A_csp%bn

        A_csp_dist%m = A_csp%m
        A_csp_dist%n = A_csp%n
        A_csp_dist%nnz = A_csp%nnz

        image_blockrows = (A_csp%m / A_csp%bm) / num_images()
        my_nblocks = 0

        ! Count blocks for my image
        do i = 1, A_csp%nblocks
            iid = (A_csp%block_inds(1, i) - 1) / image_blockrows + 1
            if (iid == this_image()) then
                my_nblocks = my_nblocks + 1
            end if
        end do

        A_csp_dist%nblocks = my_nblocks

        ! Allocate 
        allocate(A_csp_dist%blocks(my_nblocks))
        allocate(A_csp_dist%block_inds(2, my_nblocks))

        ! Deep copy blocks into my image
        offset = 1
        do i = 1, A_csp%nblocks
            iid = (A_csp%block_inds(1, i) - 1) / image_blockrows + 1
            if (iid == myim) then
                A_csp_dist%block_inds(:, offset) = A_csp%block_inds(:, i)
                call csr_deepcopy(A_csp%blocks(i), A_csp_dist%blocks(offset))
                offset = offset + 1
            end if
        end do

    end subroutine


    function map_block(self, i, j) result(idx)
        class(csp_mat), intent(in) :: self
        integer, intent(in) :: i, j
        integer :: idx(2)
        
        idx(1) = ((i - 1) / self%bm) + 1
        idx(2) = ((j - 1) / self%bn) + 1

    end function


    function csr_to_csp(A_csr) result(A_csp)

        type(csr_mat), intent(in) :: A_csr
        type(csp_mat) :: A_csp

        integer, allocatable :: block_nnz(:, :), block_exists(:,:)

        real(dp), allocatable :: vals(:)
        integer, allocatable :: inds(:,:)
        integer :: i, j, k, l
        integer :: colidx, offset
        integer :: block_idx(2)
        integer :: idx

        logical :: marked


        ! Determine number of blocks and block size
        A_csp%bm = A_csr%m2
        A_csp%bn = A_csr%n2
        A_csp%m = A_csr%m
        A_csp%n = A_csr%n
        A_csp%nnz = A_csr%nnz

        ! Count nnz in each block
        allocate(block_nnz(A_csr%m1, A_csr%n1))
        allocate(block_exists(A_csr%m1, A_csr%n1))
        block_nnz = 0
        block_exists = 0

        do i = 1, A_csr%m
            do j = A_csr%rowptrs(i), A_csr%rowptrs(i + 1)
                block_idx = A_csp%map_block(i, A_csr%colinds(j))
                block_nnz(block_idx(1), block_idx(2)) = block_nnz(block_idx(1), block_idx(2)) + 1
                block_exists(block_idx(1), block_idx(2)) = 1
            end do
        end do


        A_csp%nblocks = sum(block_exists)
        allocate(A_csp%blocks(A_csp%nblocks))
        allocate(A_csp%block_inds(2, A_csp%nblocks))


        ! This is bad
        idx = 1
        do i = 1, A_csr%m1
            do j = 1, A_csr%n1

                if (block_nnz(i,j) == 0) then
                    cycle
                end if

                marked = .FALSE.


                allocate(vals(block_nnz(i,j)))
                allocate(inds(2, block_nnz(i,j)))
                offset = 1

                do k = 1, A_csr%m
                    do l = A_csr%rowptrs(k), A_csr%rowptrs(k + 1)
                        
                        colidx = A_csr%colinds(l)

                        if (all(A_csp%map_block(k, colidx)==[i, j])) then

                            if (.NOT.marked) then
                                A_csp%block_inds(1, idx) = k
                                A_csp%block_inds(2, idx) = colidx
                                marked = .TRUE.
                            end if

                            vals(offset) = A_csr%vals(l)

                            inds(1, offset) = mod((k - 1), A_csr%m1) + 1
                            inds(2, offset) = mod((colidx - 1), A_csr%n1) + 1

                            offset = offset + 1

                        end if

                    end do
                end do

                print*, "Done with block ", idx
                A_csp%blocks(idx) = ijv_to_csr(vals, inds, A_csp%bm, A_csp%bn)
                idx = idx + 1

                deallocate(vals)
                deallocate(inds)

            end do
        end do
                

    end function


    subroutine free_csp_mat(A)
        type(csp_mat), intent(inout) :: A
        integer :: i

        do i = 1, A%nblocks
            call free_csr_mat(A%blocks(i))
        end do

        deallocate(A%block_inds)
        deallocate(A%blocks)

    end subroutine

end module csp
