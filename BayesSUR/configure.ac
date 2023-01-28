## require at least autoconf 2.62 ( OMP MACRO )
AC_PREREQ(2.62)

## Process this file with autoconf to produce a configure script.
AC_INIT(BayesSUR, m4_esyscmd_s([awk '/^Version:/ {print $2}' DESCRIPTION]))

: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])
fi
RBIN="${R_HOME}/bin/R"

## CC=c++ -fsanitize=address,undefined -fno-sanitize=float-divide-by-zero -fno-sanitize=alignment -fno-omit-frame-pointer
CXX=$("${R_HOME}/bin/R" CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

AC_LANG(C++)
AC_PROG_CXX

# Check for c++17 support
CXXFLAGS="-Werror -std=c++17 -Wno-reorder -fPIC"
AC_MSG_CHECKING([whether CXX supports -std=c++17])
AC_COMPILE_IFELSE(
	[AC_LANG_PROGRAM([])], 
	[AC_MSG_RESULT([yes]); CPPUTEST_CPP11FLAG="-std=c++17" ], 
	[AC_MSG_RESULT([no]); 
		AC_MSG_RESULT([]);
		AC_MSG_RESULT([BayesSUR needs a c++17 compatible compiler!]);
		AC_MSG_RESULT([Please install a more recent gcc if on Linux, or see]);
		AC_MSG_RESULT([https://cran.r-project.org/bin/macosx/tools/,]);
		AC_MSG_RESULT([https://cran.r-project.org/bin/windows/Rtools/]);
		AC_MSG_RESULT([]);
		AC_MSG_ERROR([Compiler not supported])])
CXXFLAGS="$saved_cxxflags"

# Check for OMP support
AC_OPENMP

# since some systems have broken OMP libraries
# we also check that the actual package will work
ac_pkg_openmp=no
if test -n "${OPENMP_CXXFLAGS}"; then
	AC_MSG_CHECKING([whether OpenMP will actually work in a package])
	AC_LANG_CONFTEST([AC_LANG_PROGRAM([[#include <omp.h>]], [[ return omp_get_num_threads (); ]])])
	PKG_CFLAGS="${OPENMP_CXXFLAGS}" PKG_LIBS="${OPENMP_CXXFLAGS}" "$RBIN" CMD SHLIB conftest.cpp 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp=yes
	AC_MSG_RESULT([${ac_pkg_openmp}])
fi

# if test "${ac_pkg_openmp}" = yes all is fine! 
if test "${ac_pkg_openmp}" = no; then 
  AC_MSG_WARN([No OpenMP support. If using GCC, upgrade to >= 4.2. If using clang, upgrade to >= 3.8.0]);
	AC_MSG_RESULT([]);
	AC_MSG_RESULT([It seems OpenMP is unsupported on your system,]);
	AC_MSG_RESULT([in order to enable parall execution of the C++ library]);
	AC_MSG_RESULT([please install a more recent gcc if on Linux, or see]);
	AC_MSG_RESULT([https://cran.r-project.org/bin/macosx/tools/ ,]);
	AC_MSG_RESULT([https://cran.r-project.org/bin/windows/Rtools/]);
	AC_MSG_RESULT([]);
fi 

AC_SUBST(OPENMP_FLAG,"$OPENMP_CXXFLAGS")

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
