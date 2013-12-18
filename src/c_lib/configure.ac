AC_INIT([csage], [0], [sage-support@googlegroups.com])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([foreign -Wall -Werror])

AC_CONFIG_MACRO_DIR([m4])

# check for python
AX_PYTHON_DEVEL

# use "make V=1" if you want to see the long awful lines
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

# check for programs
AC_PROG_CC
AC_PROG_CXX
AM_PROG_AR

# initialize libtool
LT_INIT

AC_CONFIG_FILES([Makefile
				 src/Makefile
				 include/Makefile])

AC_OUTPUT
