cmake_minimum_required(VERSION 3.5)

project(netReg LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FeatureSummary)

find_package(OpenMP)
find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED)

include_directories(${ARMADILLO_INCLUDE_DIRS})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -pipe -ffast-math -fomit-frame-pointer -flto -fPIC -Wall -mtune=core2 -g -O2 -DARMA_DONT_USE_WRAPPER -framework Accelerate")

if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(SOURCE_FILES
        src/config.h
        src/cv_fold.hpp
        src/cv_set.cpp
        src/edgenet_gaussian.cpp
        src/family.hpp
        src/error_functions.cpp
        src/error_functions.hpp
        src/graph_functions.cpp
        src/graph_functions.hpp
        src/graph_penalized_linear_model_data.hpp
        src/illegal_argument_exception.hpp
        src/io_exception.hpp
        src/linear_model_data.hpp
        src/math_functions.cpp
        src/not_implemented_exception.hpp
        src/penalized_linear_model_data.hpp
        src/stat_functions.cpp
        src/stat_functions.hpp
        src/types.hpp
        src/utility_functions.cpp
        src/utility_functions.hpp
        src/vector_functions.cpp
        src/vector_functions.hpp
        src/main.cpp)

include_directories(${Boost_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS})
add_executable(netReg ${SOURCE_FILES} ${ARMADILLO_LIBRARIES})

feature_summary(WHAT ALL)
