cmake_minimum_required(VERSION 3.6)

project(netReg LANGUAGES CXX C)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

set( CMAKE_VERBOSE_MAKEFILE on )
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    message("setting intel cxx flags")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -vec -simd -qopenmp")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE OFF)

include(GNUInstallDirs)
include(FeatureSummary)
include(CTest)

find_package(OpenMP)
find_package(Armadillo REQUIRED)
find_package(HDF5) # required because the FindArmadillo module doesn't look for HDF5
find_package(Boost COMPONENTS program_options system REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)


if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_library(
    netRegLib STATIC
        src/cv_set.cpp
        src/data_factory.cpp
        src/edgenet.cpp
        src/edgenet_model_selection.cpp
        src/graph_functions.cpp
        src/graph_model_cv_data.cpp
        src/stat_functions.cpp
        src/vector_functions.cpp)

target_compile_definitions(
    netRegLib PUBLIC
        -DARMA_DONT_USE_WRAPPER)

target_compile_options(
    netRegLib PUBLIC
        -Wall)

target_include_directories(
    netRegLib PUBLIC
        ${ARMADILLO_INCLUDE_DIRS}
        ${HDF5_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        ${BLAS_INCLUDE_DIRS}
        ${LAPACK_INCLUDE_DIRS})

target_link_libraries(
    netRegLib PUBLIC
        ${ARMADILLO_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${Boost_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES})

add_executable(
    netReg
        cpp/main.cpp)

target_link_libraries(
    netReg PRIVATE
        netRegLib)

install(
    TARGETS netReg
    DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(
    test_netReg
        tests/cpp/test.cpp
        tests/cpp/test_data.cpp
        tests/cpp/test_edgenet.cpp
        tests/cpp/test_edgenet_modelselection.cpp
        tests/cpp/test_error_functions.cpp
        tests/cpp/test_graph_functions.cpp
        tests/cpp/test_math_functions.cpp
        tests/cpp/test_model.cpp
        tests/cpp/test_stat_functions.cpp
        tests/cpp/test_vector_functions.cpp)

target_link_libraries(
    test_netReg PRIVATE
        netRegLib)

add_test(NAME main_test COMMAND test_netReg)

feature_summary(WHAT ALL)
MESSAGE( STATUS "Using CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER} )
MESSAGE( STATUS "Using CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS} )
MESSAGE( STATUS "Using CMAKE_CXX_FLAGS_RELWITHDEBINFO: " ${CMAKE_CXX_FLAGS_RELWITHDEBINFO} )
