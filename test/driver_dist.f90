program driver_dist
    use kpack
    implicit none

    character(256) :: fpath
    character(5) :: maxiters_str
    character(5) :: mode_str
    integer :: maxiters = 20

    real(dp) :: err, sing
    integer :: nnz_Q, nnz_C

    type(csp_mat) :: A
    real(dp), allocatable :: Q(:, :), C(:, :)


    ! Command line args
    call get_command_argument(1, mode_str)
    call get_command_argument(2, fpath)
    call get_command_argument(3, maxiters_str)
    read(maxiters_str, *) maxiters


    ! Print preamble
    if (this_image()==1) then
        print *, "====Beginning run of kpa===="
        print *, "Mode: ", mode_str 
        print *, "Matrix path: ", fpath
        print *, "Maxiters: ", maxiters_str 
    end if


    ! Read input matrix
    if (this_image()==1) then
        print *, "Reading input..."
    end if
    call dist_read_mm_csp(fpath, A, "S")
    if (this_image()==1) then
        print *, "Done!"
    end if

    sync all

    ! Compute the KPA 
    call kpa_svd_dist(A, Q, C, sing, maxiters)


    ! Error


    ! NNZ


    ! Done

end program driver_dist
