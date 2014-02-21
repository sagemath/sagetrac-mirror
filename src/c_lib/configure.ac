AC_INIT([csage], [0], [sage-support@googlegroups.com])

dnl Share macro and aux dir with $SAGE_ROOT autotools project
AC_CONFIG_MACRO_DIR([../../m4])
AC_CONFIG_AUX_DIR([../../config])

dnl AX_PYTHON_DEVEL macro needs aclocal-1.13 to work
AM_INIT_AUTOMAKE([1.13 foreign -Wall])

# Set default prefix to $SAGE_LOCAL
AC_PREFIX_DEFAULT(["$SAGE_LOCAL"])

# Check for Python
AX_PYTHON_DEVEL

# Use "make V=1" if you want to see the long awful lines
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

# Check for programs
AC_PROG_CC
AC_PROG_CXX
AM_PROG_AR

# Initialize libtool
LT_INIT

AC_CONFIG_HEADERS([config.h])

AC_CONFIG_FILES([Makefile
    src/Makefile
    include/Makefile])

AC_OUTPUT
