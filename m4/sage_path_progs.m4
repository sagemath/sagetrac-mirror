# This file is part of Autoconf.                       -*- Autoconf -*-
# Checking for programs.

# Copyright (C) 1992-1996, 1998-2012 Free Software Foundation, Inc.

# This file is part of Autoconf.  This program is free
# software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Under Section 7 of GPL version 3, you are granted additional
# permissions described in the Autoconf Configure Script Exception,
# version 3.0, as published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# and a copy of the Autoconf Configure Script Exception along with
# this program; see the files COPYINGv3 and COPYING.EXCEPTION
# respectively.  If not, see <http://www.gnu.org/licenses/>.

# Written by David MacKenzie, with help from
# Franc,ois Pinard, Karl Berry, Richard Pixley, Ian Lance Taylor,
# Roland McGrath, Noah Friedman, david d zuhn, and many others.


## ----------------------------- ##
## Generic checks for programs.  ##
## ----------------------------- ##

# _SAGE_PATH_PROGS_FEATURE_CHECK(VARIABLE, PROGNAME-LIST, FEATURE-TEST,
#                              [ACTION-IF-NOT-FOUND], [PATH=$PATH])
# -------------------------------------------------------------------
# FEATURE-TEST is called repeatedly with $ac_path_VARIABLE set to the
# name of a program in PROGNAME-LIST found in PATH.  FEATURE-TEST must set
# $ac_cv_path_VARIABLE to the path of an acceptable program, or else
# ACTION-IF-NOT-FOUND is executed; the default action (for internal use
# only) issues a fatal error message.  If a suitable $ac_path_VARIABLE is
# found in the FEATURE-TEST macro, it can set $ac_path_VARIABLE_found=':'
# to accept that value without any further checks.
AC_DEFUN([_SAGE_PATH_PROGS_FEATURE_CHECK],
[if test -z "$$1"; then
  ac_path_$1_found=false
  for ac_prog in $2; do
  # Loop through the user's path and test for each of PROGNAME-LIST
    _AS_PATH_WALK([$5],
    [for ac_exec_ext in '' $ac_executable_extensions; do
      ac_path_$1="$as_dir/$ac_prog$ac_exec_ext"
      AS_EXECUTABLE_P(["$ac_path_$1"]) || continue
$3
      $ac_path_$1_found && break 3
    done dnl ac_exec loop
    ])dnl ac_dir loop
  done
dnl
  if test -z "$ac_cv_path_$1"; then
    m4_default([$4],
      [AC_MSG_ERROR([no acceptable m4_bpatsubst([$2], [ .*]) could be dnl
found in m4_default([$5], [\$PATH])])])
  fi
else
  ac_cv_path_$1=$$1
fi
])


# SAGE_PATH_PROGS_FEATURE_CHECK(VARIABLE, PROGNAME-LIST,
#                             FEATURE-TEST, [ACTION-IF-NOT-FOUND=:],
#                             [PATH=$PATH])
# ------------------------------------------------------------------
# Designed to be used inside AC_CACHE_VAL.  It is recommended,
# but not required, that the user also use AC_ARG_VAR([VARIABLE]).
# If VARIABLE is not empty, set the cache variable
# $ac_cv_path_VARIABLE to VARIABLE without any further tests.
# Otherwise, call FEATURE_TEST repeatedly with $ac_path_VARIABLE
# set to the name of a program in PROGNAME-LIST found in PATH.  If
# no invocation of FEATURE-TEST sets $ac_cv_path_VARIABLE to the
# path of an acceptable program, ACTION-IF-NOT-FOUND is executed.
# FEATURE-TEST is invoked even when $ac_cv_path_VARIABLE is set,
# in case a better candidate occurs later in PATH; to accept the
# current setting and bypass further checks, FEATURE-TEST can set
# $ac_path_VARIABLE_found=':'.  Note that, unlike AC_CHECK_PROGS,
# this macro does not have any side effect on the current value
# of VARIABLE.
AC_DEFUN([SAGE_PATH_PROGS_FEATURE_CHECK],
[_$0([$1], [$2], [$3], m4_default([$4], [:]), [$5])dnl
])
