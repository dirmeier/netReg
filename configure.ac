dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT(
	[netReg],
	[0.3],
	[netreg@simon-dirmeier.net],
	[netReg],
	[https://github.com/rafstraumur/netreg])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/edgenet.cpp])

AC_ARG_ENABLE([console],
    AS_HELP_STRING([--enable-console], [Enables building the raw library.]))

AS_IF(
	[ test "x$enable_console" = "xyes"],
	[ AC_MSG_NOTICE([building from commandline])
   	],
	[
		enable_console=no	
		AC_MSG_NOTICE([building from R])
		if test -z "${R_HOME}"; 
		then
  			AC_MSG_NOTICE([no R_HOME found])
  			CXX=`"${R_HOME}/bin/R" CMD config CXX1X`
		fi	
	]
)

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
		AX_APPEND_COMPILE_FLAGS([-O3 -fPIC -pipe -Wall -ffast-math -fomit-frame-pointer -flto -DARMA_DONT_USE_WRAPPER], [CXXFLAGS])
		CXXFLAGS=$( echo $CXXFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
	]	
)

dnl Other dependencies
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])
AC_OPENMP

AC_CHECK_HEADER([armadillo],
[
	AC_DEFINE([HAVE_ARMADILLO], [1], [Define if Armadillo is installed])
	netReg_LDLIBS="-larmadillo ${LIBS}"
	netReg_CPPFLAGS="${CPPFLAGS}"
], 
[
	AC_MSG_NOTICE([no armadillo found. using bundled lib])
	netReg_CPPFLAGS="-I../inst/armadillo"
])

AC_SEARCH_LIBS([cblas_dgemm], [openblas],
[
	AC_MSG_NOTICE([using openblas])
	netReg_LDLIBS="${netReg_LDLIBS} -lopenblas -llapack"
],
[
 AC_SEARCH_LIBS([dgemm], [blas],
 [
   AC_MSG_NOTICE([using blas])
   netReg_LDLIBS="${netReg_LDLIBS} -lblas -llapack"
 ], 
 [ 
 	AS_IF(
		[ test "x${enable_console}" = "xyes" ],
		[AC_MSG_ERROR([no blas found])],
		[
			AC_MSG_NOTICE([using rblas]) 
			RBLAS=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
			RLAPACK=`"${R_HOME}/bin/R" CMD config LAPACK_LIBS`
			netReg_LDLIBS="${netReg_LDLIBS} $RBLAS $RLAPACK"
		]
	)
 ]) 
])

AC_SUBST([netReg_LDLIBS])
AC_SUBST([netReg_CPPFLAGS])
AC_SUBST([OPENMP_CXXFLAGS])

AM_INIT_AUTOMAKE([1.14.1 foreign dist-bzip2 no-dist-gzip subdir-objects])
AC_CONFIG_FILES([Makefile])
AC_SUBST([netReg_LDLIBS])
AC_SUBST([netReg_CPPFLAGS])
AC_SUBST([OPENMP_CXXFLAGS])

AC_OUTPUT([src/Makevars])

AC_MSG_RESULT([
	Have fun with: $PACKAGE_NAME

	CXX:        $CXX
	CXXFLAGS:   $CXXFLAGS
])
