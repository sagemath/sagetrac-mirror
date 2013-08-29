#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([csage], [sage-builtin], [sage-support@googlegroups.com])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([-Wall -Werror foreign -Wno-portability])
AC_CONFIG_SRCDIR([include/ccobject.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
LT_INIT

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AM_PATH_PYTHON([2.7])
AM_CHECK_PYTHON_HEADERS()
PYTHON_INCLUDE_DIR=$PYTHON_PREFIX/include
PYTHON_LIB_DIR=$PYTHON_PREFIX/lib

# freeze environment
AC_ARG_VAR(CPPFLAGS, cppflags)
AC_ARG_VAR(LDFLAGS, ldflags)
AC_ARG_VAR(PATH, path)

AC_ARG_WITH(python-include-dir,
            [  --with-python-include-dir=DIR look in DIR for Python headers],
            [PYTHON_INCLUDE_DIR=$withval],)
AC_SUBST(PYTHON_INCLUDE_DIR)
AC_ARG_WITH(python-lib-dir,
            [  --with-python-lib-dir=DIR look in DIR for Python libraries],
            [PYTHON_LIB_DIR=$withval],)
AC_SUBST(PYTHON_LIB_DIR)

AC_ARG_WITH([gmp], [AS_HELP_STRING([--with-gmp], [Build with gmp library])],
            [ AS_IF([test "x$with_gmp" != "xno"], [ gmp_or_mpir=gmp ]) ])

AC_ARG_WITH([mpir], [AS_HELP_STRING([--with-mpir], [Build with mpir library])],
            [ AS_IF([ test x$gmp_or_mpir = xgmp ],
                    [ AC_MSG_ERROR([Cannot select both, gmp and mpir])]);
              AS_IF([ test "x$with_mpir" != "xno" ], [ gmp_or_mpir=mpir ])
            ])

# Checks for libraries.
AS_IF([ test x$gmp_or_mpir != xgmp ], [
       AC_CHECK_LIB([mpir], [__gmp_version],[gmp_or_mpir=mpir] ) ])

AS_IF([ test x$gmp_or_mpir != xmpir ], [
       AC_CHECK_LIB([gmp], [__gmp_version],[gmp_or_mpir=gmp],[AC_MSG_ERROR([there is no gmp lib])])])

AS_IF([test -z "$gmp_or_mpir" ], [AC_MSG_ERROR([there is no mpir or gmp lib])])

AS_IF([ test x$gmp_or_mpir = xgmp ],[
    AC_DEFINE(USE_GMP, [1], [use gmp (not mpir)]) ])
AS_IF([ test x$gmp_or_mpir = xmpir ],[
    AC_DEFINE(USE_MPIR, [1], [use mpir (not gmp)]) ])

GMP_OR_MPIR_LIBS=-l$gmp_or_mpir
AC_SUBST(GMP_OR_MPIR_LIBS)

# Checks for header files.
AC_CHECK_HEADERS([limits.h string.h sys/time.h])

AS_IF([test $gmp_or_mpir = gmp ],
       [AC_CHECK_HEADER([gmp.h],,[AC_MSG_ERROR([no gmp header found])])
       ])
AS_IF([test $gmp_or_mpir = mpir ],
       [AC_CHECK_HEADER([mpir.h],,[AC_MSG_ERROR([no mpir header found])])
       ])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gettimeofday memset])

AC_CONFIG_FILES([Makefile
                 include/Makefile
                 src/Makefile
                ])
AC_OUTPUT

AC_MSG_RESULT([
** Configuration summary for $PACKAGE $VERSION:

   prefix:       $prefix
   gmp_or_mpir:  $gmp_or_mpir
])

# vim:ts=8:et
