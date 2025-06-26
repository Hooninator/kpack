module io

use utils
use csr

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


        fd = 1
        start_line = 1

        open(fd, file=fpath)


        !! Advance until end of header
        read(fd, *) line
        start_line = start_line + 1
        do while (index(line, '%') /= 0)
            read(fd, *) line
            start_line = start_line + 1
        end do
        
        !! Start over
        start_line = start_line - 2
        rewind(fd)
        do i=1, start_line
            read(fd,*) line
        end do

        read(fd, *) m, n, nnz
        print*, "M: ", m, " N: ", n, " NNZ: ", nnz

        A%m = m
        A%n = n
        A%nnz = nnz

        A%m1 = floor(sqrt(real(m)))
        A%m2 = m / A%m1
        A%m2 = A%m2 + abs(A%m - A%m2 * A%m1)

        A%n1 = floor(sqrt(real(n)))
        A%n2 = n / A%n1
        A%n2 = A%n2 + abs(A%n - A%n2 * A%n1)
        
        allocate(A%vals(nnz))
        allocate(A%colinds(nnz))
        allocate(rowcounts(m))
        rowcounts = 1

        !! Parse the actual values 
        do i=1, nnz
            read(fd,*) rowidx, A%colinds(i), A%vals(i)
            rowcounts(rowidx) = rowcounts(rowidx) + 1
        end do

        allocate(A%rowptrs(m+1))

        call prefix_sum_integer(rowcounts, A%rowptrs)

        close(fd)

    end function read_mm

end module io
