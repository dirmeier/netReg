dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT(
	[netReg],
	[1.7.1],
	[simon.dirmeier@web.de],
	[netReg],
	[https://github.com/dirmeier/netReg])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/edgenet_wrapper.cpp])

if test -z "${R_HOME}";
then
  AC_MSG_ERROR([no R_HOME found])
else
  AC_MSG_NOTICE([R home found in ${R_HOME}])

	CXX11=`"${R_HOME}/bin/R" CMD config CXX11`
	CXX11STD=`"${R_HOME}/bin/R" CMD config CXX11STD`
	CXX="${CXX11} ${CXX11STD}"
	CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXX11FLAGS`
	CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
  AC_MSG_NOTICE([setting CXX to ${CXX}])
fi

AC_LANG([C++])

AC_MSG_CHECKING([whether configure should try to set CXXFLAGS])
AS_IF(
	[test "x${CXXFLAGS+set}" = "xset"],
		[enable_cxxflags_setting=no],
		[enable_cxxflags_setting=yes]
)
AC_MSG_RESULT([$enable_cxxflags_setting])

AC_PROG_CXX
AC_PROG_SED
AC_PROG_RANLIB

AS_IF(
	[test "x${enable_cxxflags_setting}" = "xyes"],
	[
		AX_APPEND_COMPILE_FLAGS([-ffast-math], [CXXFLAGS])
		CXXFLAGS=$( echo $CXXFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
	]
)

AC_SEARCH_LIBS([cblas_dgemm], [openblas],
[
	AC_MSG_NOTICE([using openblas])
	netReg_LDLIBS="${netReg_LDLIBS} -lopenblas"
],
[
	AC_MSG_NOTICE([using rblas])
	RBLAS=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
	netReg_LDLIBS="${netReg_LDLIBS} $RBLAS"
])

RLAPACK=`"${R_HOME}/bin/R" CMD config LAPACK_LIBS`
RFORT=`"${R_HOME}/bin/R" CMD config FLIBS`
netReg_LDLIBS="${netReg_LDLIBS} $RLAPACK $RFORT"
netReg_CPPFLAGS="-DUSE_RCPPARMADILLO"

dnl Other dependencies
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])
AC_OPENMP

AC_SUBST([netReg_LDLIBS])
AC_SUBST([netReg_CPPFLAGS])
AC_SUBST([OPENMP_CXXFLAGS])

AC_OUTPUT([src/Makevars])

AC_MSG_RESULT([
	Have fun with: $PACKAGE_NAME

	CXX:        $CXX
	CXXFLAGS:   $CXXFLAGS
	CPPFLAGS:   $CPPFLAGS
])
