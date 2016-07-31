AM_CXXFLAGS = @OPENMP_CXXFLAGS@

netreg_headers = \
	../config.h \
	../cv_set.hpp \
	../cv_fold.hpp \
	../edgenet.hpp \
	../edgenet_loss_function.hpp \
	../edgenet_model_selection.hpp \
	../error_functions.hpp \
	../graph_functions.hpp \
	../graph_penalized_linear_model_data.hpp \
	../illegal_argument_exception.hpp \
	../io_exception.hpp \
	../linear_model_data.hpp \
	../math_functions.hpp \
	../not_implemented_exception.hpp \
	../optim.hpp \
	../pareto_optimal_point.hpp \
	../penalized_linear_model_data.hpp \
	../stat_functions.hpp \
	../types.hpp \
	../utility_functions.hpp \
	../vector_functions.hpp

noinst_LIBRARIES = libnetreg.a
libnetreg_a_SOURCES = $(netreg_headers) \
	../cv_set.cpp \
	../edgenet.cpp \
	../error_functions.cpp \
	../graph_functions.cpp \
	../math_functions.cpp \
	../stat_functions.cpp \
	../utility_functions.cpp \
	../vector_functions.cpp

TESTS = netreg_tests
check_PROGRAMS = netreg_tests

netreg_tests_SOURCES = test.cpp
netreg_tests_LDLIBS = @netReg_LDLIBS@
netreg_tests_CPPFLAGS = @netReg_CPPFLAGS@
netreg_tests_LDADD = libnetreg.a
