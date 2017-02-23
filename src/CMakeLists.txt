cmake_minimum_required(VERSION 3.5)
project(netReg)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Armadillo REQUIRED)

include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fno-inline -g -pipe -ffast-math -fomit-frame-pointer -flto -fPIC -Wall -mtune=core2 -g -O2 -DARMA_DONT_USE_WRAPPER -framework Accelerate")

set(SOURCE_FILES
        config.h
        cv_fold.hpp
        cv_set.cpp
        edgenet_gaussian.cpp
        family.hpp
        error_functions.cpp
        error_functions.hpp
        graph_functions.cpp
        graph_functions.hpp
        graph_penalized_linear_model_data.hpp
        illegal_argument_exception.hpp
        io_exception.hpp
        linear_model_data.hpp
        math_functions.cpp
        not_implemented_exception.hpp
        penalized_linear_model_data.hpp
        stat_functions.cpp
        stat_functions.hpp
        types.hpp
        utility_functions.cpp
        utility_functions.hpp
        vector_functions.cpp
        vector_functions.hpp
        main.cpp)

include_directories("/opt/local/boost/boost_1_61_0/" ${ARMADILLO_INCLUDE_DIRS})

link_directories("/opt/local/armadillo/usr/local/lib")

add_executable(netReg ${SOURCE_FILES} ${ARMADILLO_LIBRARIES})
