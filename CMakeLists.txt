cmake_minimum_required(VERSION 3.5)
project(projekt)

add_executable(projekt
    /home/szymo/HPC/Projekt7v9/projekt_HPC/src/main.cpp
    /home/szymo/HPC/Projekt7v9/projekt_HPC/src/functions.cpp
    /home/szymo/HPC/Projekt7v9/projekt_HPC/src/globals.cpp

)
target_include_directories(projekt PUBLIC include)
target_compile_features(projekt PUBLIC cxx_std_23) # potrzebujemy C++23

# Konfiguracja benchmark i Eigen3
set(benchmark_DIR "/home/szymo/spack/opt/spack/linux-ubuntu24.04-skylake/gcc-13.3.0/benchmark-1.9.1-rgxmuidn2b5tpbmjrubguyyjfeb37ohj/lib/cmake/benchmark")
add_library(openblas SHARED IMPORTED GLOBAL)
set_target_properties(openblas PROPERTIES
    IMPORTED_LOCATION_DEBUG  "/home/szymo/spack/opt/spack/linux-ubuntu24.04-skylate/gcc-13.3.0/openblas-0.3.29-qwjnrd5rdrrahbp572hatkoha4xf7gdx/lib/libopenblas.so"
    IMPORTED_LOCATION_RELEASE "/home/szymo/spack/opt/spack/linux-ubuntu24.04-skylake/gcc-13.3.0/openblas-0.3.29-qwjnrd5rdrrahbp572hatkoha4xf7gdx/lib/libopenblas.so"
)


find_package(benchmark REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(SFML REQUIRED COMPONENTS graphics window system)
find_package(OpenMP REQUIRED)

target_link_libraries(projekt PUBLIC 
    benchmark::benchmark
    Eigen3::Eigen
    openblas
    sfml-graphics sfml-window sfml-system
    OpenMP::OpenMP_CXX
)
