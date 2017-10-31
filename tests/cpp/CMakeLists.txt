cmake_minimum_required(VERSION 3.6)

project(netReg LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

set( CMAKE_VERBOSE_MAKEFILE on )
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    message("setting intel cxx flags")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -vec -simd -qopenmp")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE OFF)
enable_testing()

include(GNUInstallDirs)
include(FeatureSummary)

find_package(OpenMP)
find_package(Boost REQUIRED)
find_package(Armadillo REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(SOURCES
    ../../src/config.h
    ../../src/cv_fold.hpp
    ../../src/cv_set.cpp
    ../../src/edgenet_gaussian.cpp
    ../../src/edgenet_gaussian.hpp
    ../../src/edgenet_gaussian_loss_function.hpp
    ../../src/edgenet_gaussian_model_selection.cpp
    ../../src/edgenet_gaussian_model_selection.hpp
    ../../src/family.hpp
    ../../src/error_functions.hpp
    ../../src/graph_functions.cpp
    ../../src/graph_functions.hpp
    ../../src/graph_penalized_linear_model_data.hpp
    ../../src/graph_penalized_linear_model_cv_data.cpp
    ../../src/graph_penalized_linear_model_cv_data.hpp
    ../../src/illegal_argument_exception.hpp
    ../../src/io_exception.hpp
    ../../src/linear_model_data.hpp
    ../../src/math_functions.cpp
    ../../src/not_implemented_exception.hpp
    ../../src/optim.hpp
    ../../src/penalized_linear_model_data.hpp
    ../../src/stat_functions.cpp
    ../../src/stat_functions.hpp
    ../../src/vector_functions.cpp
    ../../src/vector_functions.hpp)

add_executable(
    test_netReg
    EXCLUDE_FROM_ALL
        ${SOURCES}
        test.cpp
        test_data.cpp
        test_functions.cpp
        test_model.cpp
        test_modelselection.cpp)

target_compile_definitions(
    test_netReg PRIVATE
        -DARMA_DONT_USE_WRAPPER)

target_compile_options(
    test_netReg PRIVATE
        -Wall)

target_include_directories(
    test_netReg PRIVATE
        ${ARMADILLO_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        ${BLAS_INCLUDE_DIRS}
        ${LAPACK_INCLUDE_DIRS})

target_link_libraries(
    test_netReg PRIVATE
        ${ARMADILLO_LIBRARIES}
        ${Boost_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES})

add_test(tests test_netReg)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS test_netReg)

feature_summary(WHAT ALL)
MESSAGE( STATUS "Using CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER} )
MESSAGE( STATUS "Using CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS} )
MESSAGE( STATUS "Using CMAKE_CXX_FLAGS_RELWITHDEBINFO: " ${CMAKE_CXX_FLAGS_RELWITHDEBINFO} )
