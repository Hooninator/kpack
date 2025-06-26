program driver
    use kpack
    implicit none

    character(256) :: fpath
    character(5) :: maxiters_str
    integer :: maxiters = 20

    real(dp) :: err
    integer :: nnz_Q, nnz_C

    type(csr_mat) :: A
    real(dp), allocatable :: Q(:, :), C(:, :)


    !! Command line args
    call get_command_argument(1, fpath)
    call get_command_argument(2, maxiters_str)
    read(maxiters_str, *) maxiters


    !! Print preamble
    print *, "Beginning run of kpa"
    print *, "Matrix path: ", fpath


    !! Read input matrix
    print *, "Reading input..."
    A = read_mm(fpath)
    print *, "Done!"


    !! Compute the KPA 
    print *, "Performing kpa..."
    call kpa_svd(A, Q, C, maxiters)
    print *, "Done!"


    !! Error
    err = compute_err(A, Q, C)
    print*, "|| A - (Q x C)||_F / ||A||_F : ",err


    !! NNZ
    nnz_Q = count_nnz(Q)
    nnz_C = count_nnz(C)
    print*, "NNZ(Q): ", nnz_Q, " NNZ(C): ", nnz_C
    print*, "Sparsity(Q): ", real(nnz_Q) / (size(Q,1)*size(Q,2))
    print*, "Sparsity(C): ", real(nnz_C) / (size(C,1)*size(C,2))


    !! Done
    call free_csr_mat(A)

end program driver
