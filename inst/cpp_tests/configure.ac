dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT(
	[netReg],
	[0.1],
	[netreg@simon-dirmeier.net],
	[netReg],
	[https://github.com/rafstraumur/netreg])
AC_CONFIG_HEADERS([../config.h])
AC_CONFIG_SRCDIR([../gaussian_edgenet.cpp])


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
		CXXFLAGS=""
		AX_APPEND_COMPILE_FLAGS([-O2 -g -pg -fPIC -pipe -Wall -ffast-math -fomit-frame-pointer -flto], [CXXFLAGS])
		CXXFLAGS=$( echo $CXXFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
	]
)

dnl Other dependencies
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])
AC_OPENMP

netReg_LDLIBS="${LIBS}"
netReg_CPPFLAGS="${CPPFLAGS} -DARMA_DONT_USE_WRAPPER"

AC_CHECK_HEADER([armadillo],
[
	AC_MSG_NOTICE([found armadillo])
	AC_DEFINE([HAVE_ARMADILLO], [1], [Define if Armadillo is installed])
	netReg_LDLIBS="-larmadillo ${netReg_LDLIBS}"
],
[
	AC_MSG_NOTICE([no armadillo found. using bundled lib])
	netReg_CPPFLAGS="${netReg_CPPFLAGS}  -I../../inst/armadillo"
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
   AC_MSG_NOTICE([cannot find blas or lapack])
 ])
])

AM_INIT_AUTOMAKE([1.14.1 foreign dist-bzip2 no-dist-gzip subdir-objects])
AC_SUBST([netReg_LDLIBS])
AC_SUBST([netReg_CPPFLAGS])
AC_SUBST([OPENMP_CXXFLAGS])

AC_OUTPUT([Makefile])

AC_MSG_RESULT([
	Have fun with: $PACKAGE_NAME

	CXX:        $CXX
	CXXFLAGS:   $CXXFLAGS
	CPPFLAGS:   $netReg_CPPFLAGS
	LDLIDS:     $netReg_LDLIBS
])