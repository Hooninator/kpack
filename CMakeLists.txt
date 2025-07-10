cmake_minimum_required(VERSION 3.30)
project(kpack LANGUAGES Fortran)

set(CMAKE_Fortran_COMPILER="ftn")
set(CMAKE_Fortran_MODULE_DIRECTORY="./build")

find_package(OpenMP REQUIRED)
add_library(BLAS::BLAS INTERFACE IMPORTED)
add_library(LAPACK::LAPACK INTERFACE IMPORTED)
find_package(fortran_stdlib REQUIRED)

add_library(kpack SHARED
            src/matrix_market_io.f90
            src/csr.f90
            src/csp.f90
            src/spmv.f90
            src/utils.f90
            src/kpa.f90
            src/lanczos.f90
            src/spmm.f90
            src/power_method.f90
            src/trace.f90
            src/kpa_dist.f90
            src/kpack.f90
            src/co_utils.f90
        )


target_compile_options(kpack PUBLIC -O3 -cpp -fcoarray=lib -fcheck=all)
target_link_libraries(kpack PUBLIC OpenMP::OpenMP_Fortran fortran_stdlib::fortran_stdlib )


macro(add_test test_name)
    add_executable(${test_name} test/${test_name}.f90)
    target_compile_options(${test_name} PUBLIC -O3 -cpp -fcoarray=lib)
    target_link_libraries(${test_name} PUBLIC kpack)
endmacro()

add_test(driver)
add_test(driver_dist)
add_test(spmv)
