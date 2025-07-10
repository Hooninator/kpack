module co_utils
use utils
implicit none


contains

    function co_norm2_dp(arr) result(x)
        real(dp), intent(in) :: arr(:)
        real(dp), allocatable :: loc_x[:]
        real(dp) :: x
        allocate(loc_x[*])
        loc_x = norm2(arr)
        loc_x = loc_x * loc_x
        call co_sum(loc_x)
        x = sqrt(loc_x)
    end function

end module co_utils
