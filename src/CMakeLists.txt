cmake_minimum_required(VERSION 3.5)
project(netReg)

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message("OpenMP could not be found. ${PROJECT_NAME} will still work, albeit in single-threaded mode")
endif()


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11  -g -pipe  -ffast-math -fomit-frame-pointer -flto -DARMA_DONT_USE_WRAPPER   -fPIC  -Wall -mtune=core2 -g -O2")

target_include_directories(ngshmmalign /Users/simondi/Library/R/3.3/library/Rcpp/include /Library/Frameworks/R.framework/Versions/3.3/Resources/library/RcppArmadillo/include)

set(SOURCE_FILES
        config.h
        cv_fold.hpp
        cv_set.cpp
        edgenet.cpp
        edgenet_binomial.cpp
        edgenet_binomial_likelihood_function.hpp
        edgenet_binomial_loss_function.hpp
        edgenet_gaussian.cpp
        edgenet_model_selection.cpp
        family.hpp
        error_functions.cpp
        graph_functions.cpp
        graph_penalized_linear_model_cv_data.cpp
        graph_penalized_linear_model_data.hpp
        illegal_argument_exception.hpp
        linear_model_data.hpp
        math_functions.cpp
        not_implemented_exception.hpp
        optim.hpp
        penalized_linear_model_data.hpp
        stat_functions.cpp
        stat_functions.hpp
        types.hpp
        utility_functions.cpp
        utility_functions.hpp
        vector_functions.cpp
        vector_functions.hpp
        main.cpp)
add_executable(netReg ${SOURCE_FILES})