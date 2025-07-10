module kpa_dist

use csp
use utils
implicit none

contains

    subroutine kpa_svd_dist(A, Q, C, sing, maxiters)
        type(csp_mat), intent(in) :: A
        real(dp), intent(out), allocatable :: Q(:,:), C(:,:)
        real(dp), intent(in) :: sing
        integer, intent(in) :: maxiters

        ! U^TAV = B
        real(dp), allocatable :: co_U(:,:)[:]
        real(dp), allocatable :: co_V(:,:)[:]
        real(dp), allocatable :: co_B(:,:)[:]
        
        !allocate(co_U(A%bdimm


    end subroutine


end module
