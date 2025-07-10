module kpa_dist

use csp
use utils
use lanczos
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
        real(dp), allocatable :: B(:,:)
        
        allocate(Q(A%m1, A%n1))
        allocate(C(A%m2, A%n2))

        if (this_image()==1) then
            print*, "Q : ",A%m1,"x",A%n1
            print*, "C : ",A%m2,"x",A%n2
        end if

        if (this_image()==1) then
            print*, "Running lanczos..."
        end if
        call lanczos_bidiag_dist(A, co_U, co_V, B, maxiters)
        if (this_image()==1) then
            print*, "Done!"
        end if


    end subroutine


end module
