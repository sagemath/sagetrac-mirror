# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#

# this macro is inspired by various implementations by
# Curtis C. Hovey
# Midokura KK
# Romain Lenglet
# Stephan Peijnik
#

# SAGE_PYTHON_MODULE([MODULE], [ACTION-IF-FOUND],
# [ACTION_IF_NOT_FOUND], [VERSION_VARIABLE])
# ----------------------------------------------------------
AC_DEFUN([SAGE_CHECK_PYTHON_MODULE],[dnl
  VERSION_VARIABLE=[__version__]
  AS_IF([ test -n "$4" ],VERSION_VARIABLE=$4)
  AC_CACHE_CHECK([for python module $1],
    [AS_TR_SH([sage_cv_python_$1])],[dnl
    AS_TR_SH([sage_cv_python_$1])=\
`$PYTHON -c "import [$1]; print $1.$VERSION_VARIABLE" 2>&AS_MESSAGE_LOG_FD`
    AS_IF([test $? -ne 0],
           [AS_VAR_SET([AS_TR_SH([sage_cv_python_$1])], [no])])])
  AS_VAR_IF([AS_TR_SH([sage_cv_python_$1])], [no], [$2], [$3])
  AS_TR_SH([PYTHON_$1_VERSION])=$AS_TR_SH([sage_cv_python_$1])
])
