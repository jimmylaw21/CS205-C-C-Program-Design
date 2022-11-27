# gcc -o main main.c mat.c matmul.c -lm -fgnu89-inline -mavx2 -DWITH_AVX2 -fopenmp -L/opt/OpenBLAS/lib -lopenblas  && ./main
cmake_minimum_required(VERSION 3.16)

project(main)

add_executable(main.c 
mat.c 
mat.h 
matmul.c
)

target_link_libraries(main.c    /opt/OpenBLAS/lib/libopenblas.so)

target_compile_options(main.c PUBLIC -lm -fgnu89-inline -mavx2 -DWITH_AVX2 -fopenmp -L/opt/OpenBLAS/lib -lopenblas)

find_package(OpenMP REQUIRED)

if(OpenMP_FOUND)
    target_link_libraries(${main} OpenMP::OpenMP_CXX)
else()
    message(FATAL_ERROR "openmp not found!")
endif()

# target_compile_definitions(main.c PUBLIC -lm -fgnu89-inline -mavx2 -DWITH_AVX2 -fopenmp -L/opt/OpenBLAS/lib -lopenblas)

target_include_directories(main.c PUBLIC /opt/OpenBLAS/include) 

