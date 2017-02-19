dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT(
	[netReg],
	[0.99.0],
	[simon.dirmeier@gmx.de],
	[netReg],
	[https://github.com/dirmeier/netReg])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/edgenet.cpp])

AC_MSG_NOTICE([${R_HOME}])
AC_MSG_NOTICE([building from R])
if test -z "${R_HOME}"; 
then
 	AC_MSG_NOTICE([no R_HOME found])
	CXX=`"${R_HOME}/bin/R" CMD config CXX1X`
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
		AX_APPEND_COMPILE_FLAGS([-pipe -ffast-math -fomit-frame-pointer -flto -DUSE_RCPPARMADILLO], [CXXFLAGS])
		CXXFLAGS=$( echo $CXXFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
	]	
)

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
])
