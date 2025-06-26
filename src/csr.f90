module csr
use utils
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
    subroutine free_csr_mat(A)
        type(csr_mat), intent(inout) :: A
        deallocate(A%vals)
        deallocate(A%colinds)
        deallocate(A%rowptrs)
    end subroutine 
end module csr
        
