module power_method

use utils
use csr
use spmv

implicit none


contains

    subroutine power_method_bidiag(B, u, sing, leftright)
        real(dp), intent(in) :: B(:,:)
        real(dp), intent(inout) :: u(:)
        real(dp), intent(out) :: sing
        character, intent(in) :: leftright 

        real(dp) :: delta, eig, eig_prev, nrm
        integer :: i
        integer :: maxiters 

        real(dp), allocatable :: u_tmp(:)

        maxiters = size(B, 2)
        call rand_vec(u)
        allocate(u_tmp(maxiters))

        do while (delta < epsilon(delta).or.(i <= maxiters))

            if (leftright == "L") then
                ! BB^Ty ---> u
                call bmv(B, u, 'T')
                call bmv(B, u, 'N')
            else if (leftright == "R") then
                ! B^TBy ---> u
                call bmv(B, u, 'N')
                call bmv(B, u, 'T')
            else
                print*, "Invalid leftright: ",leftright
            end if

            nrm = norm2(u)
            u = u / nrm

            u_tmp = u

            if (leftright == "L") then
                ! eig = u^TBB^Tu
                call bmv(B, u_tmp, 'T')
            else if (leftright == "R") then
                ! eig = u^TB^TBu
                call bmv(B, u_tmp, 'N')
            end if
            eig = dot_product(u_tmp, u_tmp)


            if (i == 1) then
                delta = eig
            else 
                delta = abs( eig - eig_prev) / abs(eig)
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
        
        integer :: n

        integer :: i, j

        n = size(B, 2)

        do i = 1, n
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
