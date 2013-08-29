#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([sage], [5.12.beta3-autotools4], [sage-support@googlegroups.com])
# ustar: allow long filenames
AM_INIT_AUTOMAKE([tar-ustar -Wall -Werror -Wno-portability])

# Use the silent-rules feature when possible.
m4_ifdef([AM_SILENT_RULES],AM_SILENT_RULES([yes]))

AC_CONFIG_SRCDIR([README])
AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_BUILD

# `pwd` resolves $abs_top_builddir
# abs_top_builddir is not available here.
AC_PREFIX_DEFAULT(`pwd`/local)

# VPATH
AM_CONDITIONAL([VPATH_BUILD], [test ! x$srcdir = "x."])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_FC       # FIXME: not required, just check.
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_GREP     # required to build singular
AC_PROG_SED      # used here and there.

### stage patch when needed
AC_PATH_PROG([PATCH], [patch])
CHK_patch=negative
AS_IF([ test -n "$PATCH" ], [CHK_patch=positive])

### stage bzip2 when needed
AC_PATH_PROG([BZIP2], [bzip2])
CHK_bzip2=negative
AS_IF([ test -n "$BZIP2" ], [CHK_bzip2=positive])

### stage pkg
PKG_PROG_PKG_CONFIG
CHK_pkgconfig=negative
if test -n "$PKG_CONFIG"; then
   PKG_CONFIG_PATH=`$PKG_CONFIG --variable pc_path pkg-config`
   CHK_pkgconfig=positive
fi

### stage GCC when needed.
# FIXME: override this (former SAGE_INSTALL_GCC environment)
CHK_gcc=positive

# Check whether $CC is some version of GCC.  If it's a different
# compiler, install GCC locally.
# CCtype=`testcc.sh $CC`

AS_IF([test "x$CC" = x],
      [AC_MSG_ERROR([No C compiler found in your path.
         Install a C compiler or set the environment variable CC to point to one.])]
     )

AS_IF([test "x$GCC" != xyes],[
    AC_MSG_NOTICE([Rejecting '$CC', it's not GNU gcc.])
    CHK_gcc=negative
],[
    build_gcc_version=`$CC -dumpversion`

    dnl Install our own GCC if the system-provided one is older than gcc-4.4.
    dnl * gcc-4.2.4 compiles a slow IML:
    dnl https://groups.google.com/forum/?fromgroups#!topic/sage-devel/Ux3t0dW2FSI
    dnl * gcc-4.3 might have trouble building ATLAS:
    dnl https://groups.google.com/forum/?fromgroups#!topic/sage-devel/KCeFqQ_w2FE
    AS_VERSION_COMPARE(["$build_gcc_version"],[4.4],[
        AC_MSG_NOTICE([Rejecting GCC. Your $CC ($build_gcc_version) is too old.])
        CHK_gcc=negative
    ])

    AS_IF([test "x$build_cpu" = "xia64"],[
        dnl GCC 4.4.x and GCC 4.5.x fail to compile PARI/GP on ia64:
        dnl * http://gcc.gnu.org/bugzilla/show_bug.cgi?id=46044
        AS_VERSION_COMPARE(["$build_gcc_version"],[4.6],[
            AC_MSG_NOTICE([Rejecting GCC. Your $CC ($build_gcc_version)
                           will fail to compile PARI/GP on ia64.])
            CHK_gcc=negative
        ])
        dnl GCC 4.7.0 is very broken on ia64, see
        dnl http://gcc.gnu.org/bugzilla/show_bug.cgi?id=48496
        AS_VERSION_COMPARE(["$build_gcc_version"],[4.7.0],[
            AC_MSG_NOTICE([Rejecting GCC. Your $CC ($build_gcc_version)
                           is known not to work on ia64.])
            CHK_gcc=negative
        ])
    ])

    dnl Also install GCC if we have version 4.6.0 or 4.6.1, which is
    dnl known to give trouble within Sage:
    dnl * http://gcc.gnu.org/bugzilla/show_bug.cgi?id=48702
    dnl * http://gcc.gnu.org/bugzilla/show_bug.cgi?id=48774
    AS_CASE([$build_gcc_version],[4.6.@<:@01@:>@],[
        AC_MSG_NOTICE([Rejecting GCC. Your $CC ($build_gcc_version)
                       has known bugs affecting Sage.])
    ])

    dnl GCC 4.6.x doesn't compile ECL on Cygwin:
    dnl * http://gcc.gnu.org/bugzilla/show_bug.cgi?id=52061
    AS_CASE([$build_os],[*cygwin*],[
        AS_CASE([$build_gcc_version],[4.6*],[
            AC_MSG_NOTICE([Rejecting GCC. Your $CC ($build_gcc_version)
                           will fail to compile ECL on Cygwin.])
           CHK_gcc=negative
        ])
    ])
])

# package checks
CHK_cliquer=
AC_CHECK_HEADER([cliquer/cliquer.h],,[CHK_cliquer=negative])

python_version_required=2.7
AM_PATH_PYTHON([$python_version_required],,[:])
AS_VAR_IF([PYTHON_VERSION],[$python_version_required],[dnl
   # whitelist python if desired:
   # CHK_python=positive
   SAGE_CHECK_PYTHON_MODULE([Cython],[CHK_cython=negative])
   SAGE_CHECK_PYTHON_MODULE([something],[CHK_something=negative])
   # AS_VERSION_COMPARE(["$SAGE_Cython_Module"],[0.6],
   #                    [CHK_Cython=negative], dnl too old
   #                    [CHK_Cython=positive], dnl exact version match
   #                    [CHK_Cython=negative]) dnl version too big
],[dnl
   AC_MSG_NOTICE([python $PYTHON_VERSION is not $python_version_required])
   # blacklist everything
   CHK_python=negative
   CHK_cython=negative
])

# try to guess SAGE_* stuff
sage_local=`pwd`/local
if test ! "x$prefix" = xNONE; then
    AC_MSG_NOTICE([custom prefix is pretty much experimental.])
fi

AC_SUBST(SAGE_LOCAL, $sage_local)

SAGE_SRC=`cd $srcdir/src; pwd`
AC_SUBST(SAGE_SRC)

# Check C++ and Fortran compilers.
if test -z "$CXX"; then
  echo >&2 "Installing GCC because a C++ compiler is missing."
  need_to_install_gcc=yes
fi

AS_IF([test "x$FC" = x],[
    AC_MSG_NOTICE([Installing GCC because a Fortran compiler is missing.])
    CHK_gcc=negative
])

# FIXME: we have to set FC,GCC,... if need_to_install_gcc=yes - why should we only set FC?
if test "$need_to_install_gcc" = yes -o -x "$sage_local/bin/gfortran"; then
    FC=gfortran
    unset SAGE_FORTRAN_LIB
fi

# FIXME: what is this trying to implement?
# (it's deprecated anyway)
# Make the fortran library symlink if requested
#if [ -f "$SAGE_FORTRAN_LIB" ]; then
#    ( cd "$SAGE_LOCAL/lib" && ln -sf "$SAGE_FORTRAN_LIB" . )
#fi

AC_MSG_CHECKING(GNU tar)
AC_PATH_PROGS_FEATURE_CHECK(
  [TAR], [tar gtar],
  [[tarout=`$ac_path_TAR --version 2>&1 | grep GNU`
    test -n "$tarout" && ac_cv_path_TAR=$ac_path_TAR ac_path_TAR_found=:]],
  [ AC_MSG_RESULT([no]);
    AC_MSG_ERROR([could not find GNU tar])
  ]
)
AC_MSG_RESULT([$ac_cv_path_TAR])
AC_SUBST([TAR], [$ac_cv_path_TAR])

AC_MSG_CHECKING(GNU make)
AC_PATH_PROGS_FEATURE_CHECK(
  [GMAKE], [make gmake],
  [[makeout=`$ac_path_GMAKE --version 2>&1 | grep GNU`
    test -n "$makeout" && ac_cv_path_GMAKE=$ac_path_GMAKE ac_path_GMAKE_found=:]],
  [ AC_MSG_RESULT([no]);
    AC_MSG_ERROR([could not find GNU make])
  ]
)
AC_MSG_RESULT([$ac_cv_path_GMAKE])
AC_SUBST([GMAKE], [$ac_cv_path_GMAKE])

# download tool (to make dist*)

# (FIXME, later)
AC_SUBST(URL_GRABBER, "wget -c")
AC_SUBST(DOWNLOAD_URL, "http://www.sagemath.org/packages/upstream")
AC_SUBST(UPSTREAM, upstream)

# SAGE_ROOT contains .BUILDSTART, so it must translate to $abs_top_builddir
# and $abs_top_builddir is just `pwd`
AC_SUBST(SAGE_ROOT, `pwd`)

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_UID_T
AC_C_INLINE
AC_TYPE_MODE_T
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_CHECK_MEMBERS([struct stat.st_blksize])
AC_STRUCT_ST_BLOCKS
AC_CHECK_MEMBERS([struct stat.st_rdev])
AC_TYPE_UINT32_T

# Checks for library functions.
AC_FUNC_ALLOCA
AC_FUNC_CHOWN
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FORK
AC_FUNC_LSTAT_FOLLOWS_SLASHED_SYMLINK
AC_FUNC_MALLOC
AC_FUNC_MKTIME
AC_FUNC_REALLOC
AC_CHECK_FUNCS([endgrent endpwent floor getcwd inet_ntoa memset mkdir pathconf setenv sqrt strchr strdup strstr uname utime])

AC_CHECK_PROG(GIT,git,yes)
have_git_repo=no
if test x"$GIT" = "xyes" && test -d $srcdir/.git; then
   have_git_repo=yes
fi

## options
AC_ARG_ENABLE([keepbuilt], AS_HELP_STRING([--enable-keepbuilt], [Delete foreign package source after build]),
 [case "${enableval}" in
  yes) keepbuilt=yes ;;
  no)  keepbuilt=no ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-keepbuilt]) ;;
 esac],[keepbuilt=$have_git_repo])
AM_CONDITIONAL([KEEPBUILT],[test $keepbuilt = yes])

AC_ARG_ENABLE([distupstream], AS_HELP_STRING([--enable-distupstream], [include upstream tarball in dist, default=yes]),
 [case "${enableval}" in
  yes) distupstream=yes ;;
  no)  distupstream=no ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-distupstream]) ;;
 esac],[distupstream=yes])
AM_CONDITIONAL([DISTUPSTREAM],[test $distupstream = yes])

## foreign package options (use with care)
AC_ARG_WITH([gmp], [AS_HELP_STRING([--with-gmp], [Build with gmp library])],
            [ AS_IF([test "x$with_gmp" != "xno"], [ gmp_or_mpir=gmp ]) ])

AC_ARG_WITH([mpir], [AS_HELP_STRING([--with-mpir], [Build with mpir library])],
            [ AS_IF([ test x$gmp_or_mpir = xgmp ],
                    [ AC_MSG_ERROR([Cannot select both, gmp and mpir])]);
              AS_IF([ test "x$with_mpir" != "xno" ], [ gmp_or_mpir=mpir ])
            ])

FOREIGN_OPTS=
AS_IF([test -n "$gmp_or_mpir" ],[WITH_GMP_OR_MPIR="--with-$gmp_or_mpir"])

FOREIGN_OPTS+=" $WITH_GMP_OR_MPIR"
AC_SUBST([FOREIGN_OPTS])

## sage foreign package in/exclusion
# base packages
# checks for packages

AC_CHECK_LIB([atlas], [ATL_scpsc], [], [CHK_atlas=negative])

AC_MSG_CHECKING([whether iconv is required])
case "$UNAME" in
    CYGWIN|HP-UX|SunOS)
        AC_MSG_RESULT([yes])
        CHK_iconv=negative
        ;;
    *)
        AC_MSG_RESULT([no])
        CHK_iconv=positive
        ;;
esac

AC_ARG_ENABLE([all], AS_HELP_STRING([--disable-all], [disable all foreign packages]),
 [case "${enableval}" in
  yes) FOREIGN_DEFAULT=true ;;
  no)  FOREIGN_DEFAULT=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-all]) ;;
 esac],[FOREIGN_DEFAULT=auto])

AC_ARG_VAR([FOREIGN_DEFAULT])

dnl these packages are available on the host system
FOREIGN_AVAILABLE=
dnl these packages have been explicitely disabled
FOREIGN_DISABLED=
dnl these packages have been explicitely enabled
FOREIGN_ENABLED=
dnl all package dependencies
FOREIGN_DEP=
dnl all package dependencies
VERSION_ALL=
dnl generated version assignments
FP_VERSIONS=

SAGE_FOREIGN_PACKAGE([atlas])
SAGE_FOREIGN_PACKAGE([autotools])
SAGE_FOREIGN_PACKAGE([cliquer])
SAGE_FOREIGN_PACKAGE([blas])
SAGE_FOREIGN_PACKAGE([boehm_gc])
SAGE_FOREIGN_PACKAGE([boost_cropped])
SAGE_FOREIGN_PACKAGE([bzip2])
SAGE_FOREIGN_PACKAGE([cddlib])
SAGE_FOREIGN_PACKAGE([cephes])
SAGE_FOREIGN_PACKAGE([conway_polynomials])
SAGE_FOREIGN_PACKAGE([cvxopt])
SAGE_FOREIGN_PACKAGE([cython])
SAGE_FOREIGN_PACKAGE([docutils])
SAGE_FOREIGN_PACKAGE([ecl])
SAGE_FOREIGN_PACKAGE([eclib])
SAGE_FOREIGN_PACKAGE([ecm])
SAGE_FOREIGN_PACKAGE([elliptic_curves])
SAGE_FOREIGN_PACKAGE([fflas_ffpack])
SAGE_FOREIGN_PACKAGE([flint])
SAGE_FOREIGN_PACKAGE([flintqs])
SAGE_FOREIGN_PACKAGE([libfplll])
SAGE_FOREIGN_PACKAGE([freetype])
SAGE_FOREIGN_PACKAGE([genus2reduction])
SAGE_FOREIGN_PACKAGE([gap])
SAGE_FOREIGN_PACKAGE([gcc],[gcc])
SAGE_FOREIGN_PACKAGE([gd])
SAGE_FOREIGN_PACKAGE([gdmodule])
SAGE_FOREIGN_PACKAGE([gf2x])
SAGE_FOREIGN_PACKAGE([gfan])
SAGE_FOREIGN_PACKAGE([givaro])
# SAGE_FOREIGN_PACKAGE([git])
SAGE_FOREIGN_PACKAGE([glpk])
SAGE_FOREIGN_PACKAGE([graphs])
SAGE_FOREIGN_PACKAGE([gsl])
SAGE_FOREIGN_PACKAGE([iconv])
SAGE_FOREIGN_PACKAGE([iml])
SAGE_FOREIGN_PACKAGE([ipython])
SAGE_FOREIGN_PACKAGE([jinja2])
SAGE_FOREIGN_PACKAGE([jmol])
SAGE_FOREIGN_PACKAGE([lapack])
SAGE_FOREIGN_PACKAGE([lcalc])
SAGE_FOREIGN_PACKAGE([lrcalc])
SAGE_FOREIGN_PACKAGE([libgap])
SAGE_FOREIGN_PACKAGE([libpng])
SAGE_FOREIGN_PACKAGE([linbox])
SAGE_FOREIGN_PACKAGE([libm4ri])
SAGE_FOREIGN_PACKAGE([libm4rie])
SAGE_FOREIGN_PACKAGE([matplotlib])
SAGE_FOREIGN_PACKAGE([maxima],[maxima])
SAGE_FOREIGN_PACKAGE([mpc])
SAGE_FOREIGN_PACKAGE([mpfi])
SAGE_FOREIGN_PACKAGE([mpfr])
SAGE_FOREIGN_PACKAGE([mpir])
SAGE_FOREIGN_PACKAGE([mpmath])
SAGE_FOREIGN_PACKAGE([networkx])
SAGE_FOREIGN_PACKAGE([ntl])
SAGE_FOREIGN_PACKAGE([numpy])
SAGE_FOREIGN_PACKAGE([palp])
SAGE_FOREIGN_PACKAGE([pari])
SAGE_FOREIGN_PACKAGE([patch],[patch])
SAGE_FOREIGN_PACKAGE([pexpect])
SAGE_FOREIGN_PACKAGE([pil])
SAGE_FOREIGN_PACKAGE([pkgconfig])
SAGE_FOREIGN_PACKAGE([polybori])
SAGE_FOREIGN_PACKAGE([polytopes_db])
SAGE_FOREIGN_PACKAGE([ppl])
SAGE_FOREIGN_PACKAGE([pycrypto])
SAGE_FOREIGN_PACKAGE([pygments])
SAGE_FOREIGN_PACKAGE([pynac])
SAGE_FOREIGN_PACKAGE([python],[python python2.7])
SAGE_FOREIGN_PACKAGE([ratpoints])
SAGE_FOREIGN_PACKAGE([r])
SAGE_FOREIGN_PACKAGE([rpy2])
SAGE_FOREIGN_PACKAGE([readline])
SAGE_FOREIGN_PACKAGE([rubiks])
# sagelib here, remove CORE_MODULES in the end
# (not yet, requires proper staging)
SAGE_FOREIGN_PACKAGE([sagelib])
SAGE_FOREIGN_PACKAGE([sagenb])
SAGE_FOREIGN_PACKAGE([sagetex])
SAGE_FOREIGN_PACKAGE([scipy])
SAGE_FOREIGN_PACKAGE([scons],[scons])
SAGE_FOREIGN_PACKAGE([setuptools])
SAGE_FOREIGN_PACKAGE([singular])
SAGE_FOREIGN_PACKAGE([sphinx])
SAGE_FOREIGN_PACKAGE([sqlalchemy])
SAGE_FOREIGN_PACKAGE([sqlite])
SAGE_FOREIGN_PACKAGE([symmetrica])
SAGE_FOREIGN_PACKAGE([sympow])
SAGE_FOREIGN_PACKAGE([sympy])
SAGE_FOREIGN_PACKAGE([tachyon])
SAGE_FOREIGN_PACKAGE([termcap])
SAGE_FOREIGN_PACKAGE([zlib])
SAGE_FOREIGN_PACKAGE([zn_poly])

AC_SUBST([FP_VERSIONS])
AC_SUBST([VERSION_ALL])
# AM_SUBST_NOTMAKE([FP_VERSIONS])

# include submodules in dist tarball
distdir_subs='distdir: DIST_SUBDIRS+=$(CORE_MODULES)
distdir: $(CORE_MODULES:%=%-configure-stamp)'
AC_SUBST([distdir_subs])
AM_SUBST_NOTMAKE([distdir_subs])

# disable the installation target.
# it simply doesnt make any sense right now.
noinstall='install-am install-recursive: install-doesntwork'
AC_SUBST([noinstall])
AM_SUBST_NOTMAKE([noinstall])

# misc
dotPRECIOUS='.PRECIOUS'
AC_SUBST([dotPRECIOUS])
AM_SUBST_NOTMAKE([dotPRECIOUS])
AC_SUBST(UNAME, `uname`)

AC_CONFIG_FILES([GNUmakefile
                 deps
                 sagedist/spkg.mk
                 sagedist/Makefile
                 sagedist/bin/Makefile
                 sagedist/etc/Makefile
                 sagedist/etc/env.sh
                 sagedist/etc/dist.sh])

AC_CONFIG_FILES([sagedist/bin/sage-dist-install],
       [chmod +x sagedist/bin/sage-dist-install])
AC_CONFIG_FILES([sagedist/bin/sage-dist-make],
       [chmod +x sagedist/bin/sage-dist-make])
AC_CONFIG_FILES([sagedist/bin/sage_fortran],
       [chmod +x sagedist/bin/sage_fortran])

# this is currently not needed.
AC_CONFIG_FILES([sage],
       [chmod +x sage])

# not a good idea. packages are not installed yet.
# AC_CONFIG_SUBDIRS([src/bin])

AC_OUTPUT
AC_SUBST([FOREIGN_ALL])
AC_SUBST([FOREIGN_DEP])
AC_SUBST([FOREIGN_AVAILABLE])
AC_SUBST([FOREIGN_DISABLED])
AC_SUBST([FOREIGN_ENABLED])

if test $FOREIGN_DEFAULT = true; then
   packageinfo="installing all packages"
   if test -n "$FOREIGN_DISABLED"; then
      packageinfo+=" except$FOREIGN_DISABLED"
   fi
elif test $FOREIGN_DEFAULT = false; then
   packageinfo="installing no foreign packages"
   if test -n "$FOREIGN_ENABLED"; then
      packageinfo+=" except$FOREIGN_ENABLED"
   fi
else
   packageinfo="classic (automatic) mode"
   if test -n "$FOREIGN_AVAILABLE"; then
      packageinfo+=", unnecessary:$FOREIGN_AVAILABLE"
   fi
   if test -n "$FOREIGN_DISABLED"; then
      packageinfo+=", also skipping$FOREIGN_DISABLED"
   fi
   if test -n "$FOREIGN_ENABLED"; then
      packageinfo+=", installing anyway$FOREIGN_ENABLED"
   fi
fi

AC_MSG_RESULT([
** Configuration summary for $PACKAGE $VERSION:

   prefix:     $prefix
   keepbuilt:  $keepbuilt
   $packageinfo
])

# vim:et
