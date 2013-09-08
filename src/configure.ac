#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
AC_PREREQ([2.69])
AC_INIT([sagelib], [5.12-beta-3], [sage-devel@googlegroups.com])

# do not mix up with sage-the-distribution parts
AC_CONFIG_AUX_DIR([.])

AM_INIT_AUTOMAKE([1.11 tar-ustar -Wall -Werror -Wno-portability foreign])
AC_CONFIG_SRCDIR([c_lib/configure.ac])

AC_CONFIG_SUBDIRS([bin c_lib sage ext doc])

distdir_subs='distdir: DIST_SUBDIRS+=$(CORE_MODULES)'
AC_SUBST([distdir_subs])
AM_SUBST_NOTMAKE([distdir_subs])

# freeze environment
AC_ARG_VAR(CPPFLAGS, cppflags)
AC_ARG_VAR(LDFLAGS, ldflags)
AC_ARG_VAR(LD_LIBRARY_PATH, ld_library_path)
AC_ARG_VAR(PATH, path)

AC_ARG_ENABLE([doc], AS_HELP_STRING([--enable-doc], [Build docs]),
 [case "${enableval}" in
  yes) doc=true ;;
  no)  doc=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-doc]) ;;
esac],[doc=no])
AM_CONDITIONAL([BUILD_DOC], [test x$doc = xtrue])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
