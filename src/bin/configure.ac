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

SAGE_PATH_REDUNDANCY=
AC_ARG_ENABLE([path_redundancy],
              AS_HELP_STRING([--enable-path_redundancy],
                             [double check whether commands belong to sage]),
 [case "${enableval}" in
  yes) SAGE_PATH_REDUNDANCY=true ;;
  no)  SAGE_PATH_REDUNDANCY=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-path_redundancy]) ;;
 esac],[PATH_REDUNDANCY=false] )

AC_SUBST([PATH_REDUNDANCY])

JUST_SAGE_LOCAL_BIN=
if test x$SAGE_PATH_REDUNDANCY = xtrue; then
	JUST_SAGE_LOCAL_BIN='$SAGE_LOCAL/bin'
fi
AC_SUBST([JUST_SAGE_LOCAL_BIN])

## FIXME: implement
SAGE64=no
AC_SUBST([SAGE64])

AC_CONFIG_FILES([sage], [chmod +x sage])

# needs sysconfdir...
AC_CONFIG_FILES([sage-spkg], [chmod +x sage-spkg])

AC_OUTPUT([Makefile env.sh])
