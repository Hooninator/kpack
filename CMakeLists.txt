cmake_minimum_required(VERSION 3.30)
project(kpack LANGUAGES Fortran)

set(CMAKE_Fortran_COMPILER="ftn")
set(CMAKE_Fortran_MODULE_DIRECTORY="./build")

add_library(kpack SHARED
            src/matrix_market_io.f90
            src/csr.f90
            src/spmv.f90
            src/utils.f90
            src/kpa.f90
            src/lanczos.f90
            src/kpack.f90
            src/lra.f90
            src/spmm.f90)


target_compile_options(kpack PUBLIC -O3 -std=gnu -cpp)


macro(add_test test_name)
    add_executable(${test_name} test/${test_name}.f90)
    target_compile_options(${test_name} PUBLIC -O3 -std=gnu -cpp)
    target_link_libraries(${test_name} PUBLIC kpack)
endmacro()

add_test(driver)
add_test(spmv)
