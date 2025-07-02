module power_method

use utils
use csr
use spmv

implicit none


contains

    subroutine power_method_bidiag(B, u, sing, leftright)
        real(dp), intent(in) :: B(:,:)
        real(dp), intent(inout) :: u(:)
        real(dp), intent(in) :: sing
        character, intent(in) :: leftright 

        real(dp) :: delta, eig, eig_prev, nrm
        integer :: i = 1
        integer :: maxiters = size(B, 2)

        !! BB^Ty ---> u
        do while (delta < epsilon(delta).or.(i <= maxiters))
            if (leftright == "L") then
                call bmv(B, u, 'T')
                call bmv(B, u, 'N')
            else if (leftright == "R") then
                call bmv(B, u, 'N')
                call bmv(B, u, 'T')
            else
                print*, "Invalid leftright: ",leftright
            end if
            nrm = norm2(u)
            u = u / nrm
            eig = dot(u, u)
            if (i == 1) then
                delta = eig
            else 
                delta = abs( eig - eig_prev)
            end if
            eig_prev = eig
            i = i + 1
        end do

        sing = sqrt(eig)

    end subroutine


    subroutine bmv(B, x, trans)
        real(dp), intent(in) :: B(:,:)
        real(dp), intent(inout) :: x(:)
        character, intent(in) :: trans

        real(dp) :: tmp
        real(dp) :: acc
        
        integer :: n = size(B, 2)

        integer :: i, j

        do i = 1, n
            acc = 0.0_dp
            acc = B(2, i) * x(i)
            if (trans == 'N') then
                if (i < n) then
                    acc = acc + x(i + 1) * B(1, i)
                end if
            else if (trans =='T') then
                if (i > 1) then
                    acc = acc + tmp * B(1, i)
                end if
                tmp = x(i)
            else
                print*, "Invalid transpose code ", trans
                call abort()
            end if
            x(i) = acc
        end do

    end subroutine


end module power_method
