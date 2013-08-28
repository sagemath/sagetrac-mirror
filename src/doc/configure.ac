#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([sage-doc], [none], [sage-devel@googlegroups.com])
AM_INIT_AUTOMAKE([-Wall -Werror -Wno-portability foreign])

# Checks for programs.
AC_PROG_LN_S
AC_CHECK_PROG(SAGE, sage, yes)
if test ! x"$SAGE" = "xyes"; then
	AC_MSG_NOTICE([cannot find sage])
fi

python_version_required=2.7
## python stuff
AM_PATH_PYTHON([$python_version_required])
if ! test x$PYTHON_VERSION = x"$python_version_required"; then
   AC_MSG_ERROR([python $PYTHON_VERSION is not $python_version_required])
fi

SAGE_DOC=`pwd`
AC_SUBST([SAGE_DOC])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([common/build_options.py])
AC_OUTPUT([Makefile common/Makefile])
