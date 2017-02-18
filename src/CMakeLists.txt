cmake_minimum_required(VERSION 3.5)
project(netReg)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -pipe -ffast-math -fomit-frame-pointer -flto -fPIC -Wall -mtune=core2 -g -O2")

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

include_directories(netReg
        "/Users/simondi/Library/R/3.3/library/Rcpp/include"
        "/Library/Frameworks/R.framework/Versions/3.3/Resources/library/RcppArmadillo/include"
        "/Library/Frameworks/R.framework/Resources/include"
        "/usr/local/include" "/opt/X11/include" "/usr/local/include/freetype2" )

link_directories(
        "/Library/Frameworks/R.framework/Resources/lib"
        "/Library/Frameworks/R.framework/Versions/3.3/Resources/library/Rcpp/libs"
        "/Library/Frameworks/R.framework/Versions/3.3/Resources/library/RcppArmadillo/libs"
        "/usr/local/lib")

add_executable(netReg ${SOURCE_FILES})