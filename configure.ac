#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([mgcal], [0.1.1])
AC_CONFIG_SRCDIR([src/mgcal.c])
AC_CONFIG_HEADERS([config.h])

MGCAL_CURRENT=2
MGCAL_REVISION=0
MGCAL_AGE=2
MGCAL_LT_VERSION="${MGCAL_CURRENT}:${MGCAL_REVISION}:${MGCAL_AGE}"
AC_SUBST(MGCAL_LT_VERSION)

AM_INIT_AUTOMAKE([subdir-objects -Wall])

m4_ifdef([AM_PROG_AR], [AM_PROG_AR([ACT-IF-FAIL])])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_LIBTOOL

# Checks for libraries.
AC_CHECK_LIB([m], [floor], [])
AC_CHECK_LIB([gfortran], [main],[])

# Checks for header files.
AC_CHECK_HEADER([stdlib.h])
AC_CHECK_HEADER([string.h])
#AC_CHECK_HEADER([blas.h])
#AC_CHECK_HEADER([lapack.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt log strrchr])
#AC_CHECK_FUNCS([dpotrf_ dpotrs_])

AC_ARG_ENABLE(
	[openmp],
	[AC_HELP_STRING([--enable-openmp],
	[enable openMP support [default=yes]])],
	[],
	[enable_openmp=yes]
)

if test "x$enable_openmp" != xno; then
	case "$CC" in
		gcc )
			OMPFLAGS="-fopenmp"
			;;
		pgcc )
			OMPFLAGS="-mp"
			;;
		icc )
			OMPFLAGS="-openmp"
			;;
		* )
			OMPFLAGS=""
			;;
	esac
	CFLAGS="$CFLAGS $OMPFLAGS"
fi
AC_SUBST([OMPFLAGS])

# pkgconfig
MGCAL_DESCRIPTION='support library for magnetic field calculation'
AC_SUBST([MGCAL_DESCRIPTION])
AC_CONFIG_FILES([mgcal.pc:pkgconfig/mgcal.pc.in])

AC_CONFIG_FILES([Makefile
                 include/Makefile
                 src/Makefile
                 examples/Makefile])
AC_OUTPUT
