# SYNOPSIS
#
#   AX_CHECK_ROOT
#
# DESCRIPTION
#
#   Check if the current UNIX user is root.
#
#   Example:
#
#     AX_CHECK_ROOT([action-if-yes], [action-if-no])
#     AX_CHECK_ROOT([], [AC_MSG_ERROR([We need root privileges])])
#     AX_CHECK_ROOT([AC_MSG_ERROR([Must not be root])], [])
#
# LICENSE
#
#   Copyright (c) 2014 Volker Braun <vbraun.name@gmail.com>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 5

AC_DEFUN([AX_CHECK_ROOT],[
AC_MSG_CHECKING([for root user])

uid=`id -u`
if test "x$uid" = "x0"; then
    AC_MSG_RESULT([yes])
    $1
else
    AC_MSG_RESULT([no])
    $2
fi;dnl

])
