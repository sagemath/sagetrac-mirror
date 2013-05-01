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

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([limits.h string.h sys/time.h])

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

# vim:ts=8:et
