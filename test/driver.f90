program driver
    use kpack
    implicit none

    character(256) :: fpath
    character(5) :: maxiters_str
    character(10) :: mode_str
    integer :: maxiters = 20

    real(dp) :: err, sing
    integer :: nnz_Q, nnz_C

    type(csr_mat) :: A
    real(dp), allocatable :: Q(:, :), C(:, :)


    ! Command line args
    call get_command_argument(1, mode_str)
    call get_command_argument(2, fpath)
    call get_command_argument(3, maxiters_str)
    read(maxiters_str, *) maxiters


    ! Print preamble
    print *, "====Beginning run of kpa===="
    print *, "Mode: ", mode_str 
    print *, "Matrix path: ", fpath
    print *, "Maxiters: ", maxiters_str 


    ! Read input matrix
    print *, "Reading input..."
    A = read_mm(fpath)
    print *, "Done!"


    ! Compute the KPA 
    print *, "Performing kpa..."
    if (mode_str == "als") then
        call kpa_als(A, Q, C, sing, maxiters)
    else if (mode_str == "svd") then
        call kpa_svd(A, Q, C, sing, maxiters)
    else
        print*, "Error: ", mode_str, " is not a valid mode"
        call abort()
    end if 
    print *, "Done!"


    ! Error
    err = compute_err(A, Q, C, sing)
    print*, "|| A - (Q x C)||_F / ||A||_F : ",err


    ! NNZ
    nnz_Q = count_nnz(Q, sing)
    nnz_C = count_nnz(C, sing)
    print*, "NNZ(Q): ", nnz_Q, " NNZ(C): ", nnz_C
    print*, "Sparsity(Q): ", real(nnz_Q) / (size(Q,1)*size(Q,2))
    print*, "Sparsity(C): ", real(nnz_C) / (size(C,1)*size(C,2))


    ! Done
    call free_csr_mat(A)

end program driver
