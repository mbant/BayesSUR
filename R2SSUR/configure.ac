## require at least autoconf 2.62 ( OMP MACRO )
AC_PREREQ(2.62)

## Process this file with autoconf to produce a configure script.
AC_INIT(R2SSUR, m4_esyscmd_s([awk -e '/^Version:/ {print $2}' DESCRIPTION]))

: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])
fi

CXX=$("${R_HOME}/bin/R" CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

AC_LANG(C++)
AC_PROG_CXX

# Check for c++11 support
CXXFLAGS="-Wall -std=c++11 -Wno-reorder"
AC_MSG_CHECKING([whether CXX supports -std=c++11])
AC_COMPILE_IFELSE(
	[AC_LANG_PROGRAM([])], 
	[AC_MSG_RESULT([yes]); CPPUTEST_CPP11FLAG="-std=c++11" ], 
	[AC_MSG_RESULT([no]); 
		AC_MSG_RESULT([]);
		AC_MSG_RESULT([rBSEM needs a c++11 compatible compiler!]);
		AC_MSG_RESULT([Please install a more recent gcc if on Linux, or see]);
		AC_MSG_RESULT([https://cran.r-project.org/bin/macosx/tools/,]);
		AC_MSG_RESULT([https://cran.r-project.org/bin/windows/Rtools/]);
		AC_MSG_RESULT([]);
		AC_MSG_ERROR([Compiler not supported])])
CXXFLAGS="$saved_cxxflags"

# Check for OMP support
AC_OPENMP

		AC_MSG_RESULT([]);
		AC_MSG_RESULT([If the previous line display OpenMP as unsupported, ]);
		AC_MSG_RESULT([in order to enable parall execution of the C++ library ]);
		AC_MSG_RESULT([please install a more recent gcc if on Linux, or see]);
		AC_MSG_RESULT([https://cran.r-project.org/bin/macosx/tools/,]);
		AC_MSG_RESULT([https://cran.r-project.org/bin/windows/Rtools/]);
		AC_MSG_RESULT([]);

AC_SUBST(OPENMP_FLAG,"$OPENMP_CXXFLAGS")

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
