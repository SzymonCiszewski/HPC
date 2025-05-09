cmake_minimum_required(VERSION 3.5)
project(projekt) # this does not matter for our simple case
add_executable(projekt
    /home/szymo/HPC/Projekt7v9/projekt_HPC/src/main.cpp
    
)
target_include_directories(projekt PUBLIC include)
target_compile_features(projekt PUBLIC cxx_std_23) # we need C++20
set(benchmark_DIR "/home/szymo/spack/opt/spack/linux-ubuntu24.04-skylake/gcc-13.3.0/benchmark-1.9.1-rgxmuidn2b5tpbmjrubguyyjfeb37ohj/lib/cmake/benchmark")
find_package(
    benchmark REQUIRED
)
find_package(Eigen3 REQUIRED)

target_link_libraries(projekt PUBLIC 
    benchmark::benchmark
    Eigen3::Eigen
    )