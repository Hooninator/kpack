module io

use utils
use csr
use csp

implicit none


contains

    function read_mm(fpath) result(A)
        implicit none

        character(256), intent(in) :: fpath

        character(512) :: line

        type(csr_mat) :: A

        real(dp), dimension(:), allocatable :: vals
        integer, dimension(:), allocatable :: colinds
        integer, dimension(:), allocatable :: rowptrs
        integer, dimension(:), allocatable :: rowcounts

        integer :: nnz, m, n

        integer :: fd, start_line

        integer :: i
        integer :: rowidx 

        integer :: pattern = 0, symmetric = 0


        fd = 1
        start_line = 1

        open(fd, file=fpath)


        !! Advance until end of header
        read(fd, '(A)') line
        if (index(line, "pattern") /= 0) then
            pattern = 1
            print*, "This is a pattern matrix"
        end if
        if (index(line, "symmetric") /= 0) then
            symmetric = 1
            print*, "This is a symmetric matrix, this can't handle those yet"
            call ABORT()
        end if

        do while (index(line, '%') /= 0)
            read(fd, '(A)') line
            start_line = start_line + 1
        end do
        
        !! Start over
        start_line = start_line - 1
        rewind(fd)
        do i=1, start_line
            read(fd,*) line
        end do

        read(fd, *) m, n, nnz
        print*, "M: ", m, " N: ", n, " NNZ: ", nnz

        A%m = m
        A%n = n
        A%nnz = nnz

        call csr_set_subdims(A)
        
        allocate(A%vals(nnz))
        allocate(A%colinds(nnz))
        allocate(rowcounts(m))
        rowcounts = 0

        !! Parse the actual values 
        do i=1, nnz
            if (pattern == 1) then
                read(fd,*) rowidx, A%colinds(i)
                A%vals(i) = 1
            else
                read(fd,*) rowidx, A%colinds(i), A%vals(i)
            end if
            rowcounts(rowidx) = rowcounts(rowidx) + 1
        end do


        call prefix_sum_integer(rowcounts, A%rowptrs)

        close(fd)

    end function read_mm


    subroutine dist_read_mm_csp(fpath, A, mapping)
        character(256), intent(in) :: fpath
        character, intent(in) :: mapping
        type(csp_mat), intent(inout) :: A

        type(csr_mat) :: A_csr
        type(csp_mat) :: A_csp

        A_csr = read_mm(fpath)
        print*, "Done reading mm"
        A_csp = csr_to_csp(A_csr)
        print*, "Done converting to csp"

        ! Distribute blocks to appropriate images
        if (mapping == "S") then
            call dist_csp_static(A_csp, A)
        else
            print*, "Invalid mapping ", mapping
            call abort()
        end if

        call free_csr_mat(A_csr)
        call free_csp_mat(A_csp)

    end subroutine

end module io
