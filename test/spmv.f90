program spmv_permuted_test
    use kpack
    implicit none

    character(256) :: fpath
    type(csr_mat) :: A
    real(dp), allocatable :: x(:), y(:)
    integer :: m, n

    call get_command_argument(1, fpath)

    A = read_mm(fpath)
    m = A%m
    n = A%n

    allocate(y(m))
    allocate(x(n))

    y = 0
    call rand_vec(x)

    call spmv_permuted(A, x, y, 0)
    call spmv_permuted(A, y, x, 1)

    call free_csr_mat(A)
    print*, "SpMV ran without crashing"

end program spmv_permuted_test
