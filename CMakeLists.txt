cmake_minimum_required(VERSION 3.6)

project(netReg LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE OFF)


include(FeatureSummary)

find_package(OpenMP)
find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -O2 -pipe -ffast-math -fomit-frame-pointer -flto -fPIC -Wall -mtune=core2 -DARMA_DONT_USE_WRAPPER")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  message("setting intel cxx flags")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-report=5 -vec -simd -qopenmp")
endif()

if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(SOURCE_FILES
        src/config.h
        src/cv_fold.hpp
        src/cv_set.cpp
        src/edgenet_gaussian.cpp
        src/edgenet_gaussian.hpp
        src/edgenet_gaussian_loss_function.hpp
        src/edgenet_gaussian_model_selection.cpp
        src/edgenet_gaussian_model_selection.hpp
        src/family.hpp
        src/error_functions.hpp
        src/graph_functions.cpp
        src/graph_functions.hpp
        src/graph_penalized_linear_model_data.hpp
        src/graph_penalized_linear_model_cv_data.cpp
        src/graph_penalized_linear_model_cv_data.hpp
        src/illegal_argument_exception.hpp
        src/io_exception.hpp
        src/linear_model_data.hpp
        src/math_functions.cpp
        src/not_implemented_exception.hpp
        src/optim.hpp
        src/penalized_linear_model_data.hpp
        src/stat_functions.cpp
        src/stat_functions.hpp
        src/utility_functions.cpp
        src/utility_functions.hpp
        src/vector_functions.cpp
        src/vector_functions.hpp
        cmd_line/main.cpp)

include_directories(${Boost_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS})
link_libraries(${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(netReg ${SOURCE_FILES} ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

feature_summary(WHAT ALL)
