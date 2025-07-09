module trace
use utils
use csr
implicit none

contains

    ! T_ij = trace(X^TA_ij) / gamma
    subroutine trace_permuted(A, X, T, hat)
        type(csr_mat) :: A
        real(dp), intent(in) :: X(:, :)
        real(dp), intent(inout) :: T(:,:)
        character, intent(in) :: hat

        real(dp) :: gamma, val
        integer :: i, j, colidx, block_i, block_j, loc_i, loc_j
        integer :: m2, n2

        m2 = A%m2
        n2 = A%n2

        T = 0

        gamma = norm2(X)
        gamma = gamma * gamma
        
        if (hat == 'N') then

            do i = 1, A%m
                do j = A%rowptrs(i), A%rowptrs(i + 1)

                    val = A%vals(j)
                    colidx = A%colinds(j)
                        
                    block_i = (i - 1) / m2 + 1
                    block_j = (colidx - 1) / n2 + 1

                    loc_i = mod(i - 1, m2) + 1
                    loc_j = mod(colidx - 1, n2) + 1

                    T(block_i, block_j) = T(block_i, block_j) + val * X(loc_i, loc_j)

                end do
            end do

        else if (hat == 'Y') then
            do i = 1, A%m
                do j = A%rowptrs(i), A%rowptrs(i + 1)

                    val = A%vals(j)
                    colidx = A%colinds(j)

                    block_i = mod(i - 1, m2) + 1
                    block_j = mod(colidx - 1, n2) + 1

                    loc_i = (i - 1) / m2 + 1
                    loc_j = (colidx - 1) / n2 + 1

                    T(block_i, block_j) = T(block_i, block_j) + val * X(loc_i, loc_j)

                end do
            end do
        else 
            print*, "Bad hat value: ", hat
            call abort()
        end if
        
        T = T / gamma

    end subroutine

end module trace
