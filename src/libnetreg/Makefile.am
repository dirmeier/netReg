AM_CXXFLAGS = @OPENMP_CXXFLAGS@

netreg_headers = \
	config.h \
	data/cvdata/cv_set.hpp \
	data/cvdata/cv_fold.hpp \
	data/modeldata/graph_penalized_linear_model_data.hpp \
	data/modeldata/linear_model_data.hpp \
	data/modeldata/penalized_linear_model_data.hpp \
	exceptions/illegal_argument_exception.hpp \
	exceptions/io_exception.hpp \
	exceptions/not_implemented_exception.hpp \
	func/error_functions.hpp \
	func/graph_functions.hpp \
	func/math_functions.hpp \
	func/stat_functions.hpp \
	func/utility_functions.hpp \
	func/vector_functions.hpp \
	models/edgenet.hpp \
	modelselection/loss_functions/edgenet_loss_functions.hpp \
	modelselection/optim/optim.hpp \
	modelselection/edgenet_model_selection.hpp \
	struc/types.hpp

noinst_LIBRARIES = libnetreg.a
libnetreg_a_SOURCES = $(netreg_headers) \
	data/cvdata/cv_set.cpp \
	func/error_functions.cpp \
	func/graph_functions.cpp \
	func/math_functions.cpp \
	func/stat_functions.cpp \
	func/utility_functions.cpp \
	func/vector_functions.cpp \
	models/edgenet.cpp \
	modelselection/edgenet_model_selection.cpp

libnetreg_a_CPPFLAGS = @netReg_CPPFLAGS@
