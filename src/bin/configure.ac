#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([sage-scripts], [none], [sage-devel@googlegroups.com])
AM_INIT_AUTOMAKE([-Wall -Werror -Wno-portability foreign])

# Checks for programs.

AC_PATH_PROGS([SINGULAR],
              [Singular singular],
              [no])
AS_IF([test "$SINGULAR" = no], [AC_MSG_ERROR([cannot find singular.])])

# HACK: don't know how to retrieve share from Singular
SINGULAR_SHARE=`dirname $SINGULAR`/../share/singular
if test ! -e $SINGULAR_SHARE/algebra.lib; then
	AC_MSG_NOTICE([unexpected singular install, cannot find share])
fi
AC_SUBST([SINGULAR_SHARE])

AS_AC_EXPAND([SAGE_ETC], $sysconfdir)

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

## FIXME: implement
SAGE64=no
AC_SUBST([SAGE64])

AC_CONFIG_FILES([sage], [chmod +x sage])

AC_OUTPUT([Makefile env.sh])
