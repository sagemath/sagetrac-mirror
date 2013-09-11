
AC_CONFIG_SRCDIR([graphs/graph_bundle.py])
AC_CONFIG_MACRO_DIR([m4])
AC_INIT([python-sage], [none], [sage-support@googlegroups.com])
AM_SILENT_RULES
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([-Wall -Werror foreign -Wno-portability subdir-objects])

# Use the silent-rules feature when possible.
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

# this is not in stable automake yet
# AM_EXTRA_RECURSIVE_TARGETS([py])

AM_PROG_AR
LT_INIT([disable-static])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_LN_S

## environment,
# see autoconf/manual/autoconf-2.60/html_node/Setting-Output-Variables.html,
# particularly: "variable is kept during automatic reconfiguration"
AC_ARG_VAR(CPPFLAGS, cppflags)
AC_ARG_VAR(LDFLAGS, ldflags)
AC_ARG_VAR(PATH, path)
AC_ARG_VAR(PYTHONPATH, pythonpath)
AC_ARG_VAR(LD_LIBRARY_PATH, ld_library_path)
AC_ARG_VAR(LIBRARY_PATH, library_path)
# don't know how to get SINGULARPATH from executable, just take it from env for now
AC_ARG_VAR(SINGULARPATH, singularpath)

# use configure environment under all circumstances
# this way, the user has to care the least.
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH"
export LIBRARY_PATH="$LIBRARY_PATH"
export SINGULARPATH="$SINGULARPATH"
export PATH="$PATH"
# ...

# see common.mk, this is important when using sage-the-distribution toolchain
cfg_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
AC_SUBST([cfg_LD_LIBRARY_PATH])

AC_ARG_ENABLE([doc], AS_HELP_STRING([--enable-doc], [Build docs]),
 [case "${enableval}" in
  yes) doc=true ;;
  no)  doc=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-doc]) ;;
esac],[doc=no])
AM_CONDITIONAL([BUILD_DOC], [test x$doc = xtrue])

SAGE_LOCAL=NONE
AC_ARG_WITH(sage_local,
            AS_HELP_STRING([--with-sage_local=DIR].[set SAGE_LOCAL=DIR]),
            [SAGE_LOCAL=$withval])
AC_SUBST(SAGE_LOCAL)

## VPATH
AM_CONDITIONAL([VPATH_BUILD], [test ! x$srcdir = "x."])

# FIXME: use VPATH_BUILD_{TRUE/FALSE}
VPATH_TRUE=
VPATH_FALSE=#
if test x$srcdir = "x."; then
  VPATH_TRUE=#
  VPATH_FALSE=
fi
AC_SUBST(VPATH_TRUE)
AC_SUBST(VPATH_FALSE)
## defaults
AC_SUBST(AM_LDFLAGS, "-module -avoid-version")

# doesn't work... (how to set a default?!)
#AC_SUBST([AM_DEFAULT_SOURCE_EXT], [.pyx])

### arg parsing
AC_ARG_WITH([gmp], [AS_HELP_STRING([--with-gmp], [Build with gmp library])],
            [ AS_IF([test "x$with_gmp" != "xno"], [ gmp_or_mpir=gmp; GMP_or_MPIR=GMP ]) ])

AC_ARG_WITH([mpir], [AS_HELP_STRING([--with-mpir], [Build with mpir library])],
            [ AS_IF([ test x$gmp_or_mpir = xgmp ],
                    [ AC_MSG_ERROR([Cannot select both, gmp and mpir])]);
              AS_IF([ test "x$with_mpir" != "xno" ], [ gmp_or_mpir=mpir; GMP_or_MPIR=MPIR ])
            ])

## python stuff
python_version_required=2.7
AM_PATH_PYTHON([$python_version_required])
if ! test x$PYTHON_VERSION = x"$python_version_required"; then
   AC_MSG_ERROR([python $PYTHON_VERSION is not $python_version_required])
fi
AM_PATH_PYTHON([2.7])

AC_PATH_PROGS([CYTHON],
              [cython2.7 cython],
              [no])
AS_IF([test "$CYTHON" = no], [AC_MSG_ERROR([cannot find cython.])])

AM_CHECK_PYTHON_HEADERS()
PYTHON_INCLUDE_DIR=$PYTHON_PREFIX/include
PYTHON_LIB_DIR=$PYTHON_PREFIX/lib

AC_ARG_WITH(python-include-dir,
            [  --with-python-include-dir=DIR look in DIR for Python headers],
            [PYTHON_INCLUDE_DIR=$withval],)
AC_SUBST(PYTHON_INCLUDE_DIR)
my_CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"
AC_ARG_WITH(python-lib-dir,
            [  --with-python-lib-dir=DIR look in DIR for Python libraries],
            [PYTHON_LIB_DIR=$withval],)
AC_SUBST(PYTHON_LIB_DIR)

AC_MSG_CHECKING([whether jinja2 is available])
if ! $PYTHON -c import\ jinja2 &>/dev/null; then
   AC_MSG_RESULT([no])
   AC_MSG_ERROR([no jinja2])
fi
AC_MSG_RESULT([yes])

AC_MSG_CHECKING([whether $CYTHON supports dependencies])
if $CYTHON --help | grep -q 'dependency makefile'; then
   pydep=yes
else
   pydep=no
fi
AC_MSG_RESULT($pydep)
AM_CONDITIONAL([PYDEP], [test $pydep = yes])

AC_MSG_CHECKING([whether $CYTHON supports gdb])
if $CYTHON --help | grep -q 'gdb-outdir'; then
   cygdb=yes
else
   cygdb=no
fi
AC_MSG_RESULT($cygdb)
AM_CONDITIONAL([CYGDB], [test $cygdb = yes])

# find clib. might be obsolete
# AC_ARG_WITH(blas, --with-clib   choose c_lib, [clib=$with_clib], [])

AC_CHECK_HEADERS([math.h assert.h sys/sysmacros.h fcntl.h mach/mach.h sys/mman.h])
AC_CHECK_HEADERS([sgtty.h stdlib.h string.h sys/ioctl.h sys/param.h sys/time.h sys/resource.h])
AC_CHECK_HEADERS([termio.h termios.h unistd.h signal.h stdio.h libc.h sys/times.h sys/wait.h])
AC_CHECK_HEADERS([vfork.h dirent.h errno.h setjmp.h])

AC_FUNC_VFORK
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_CHECK_FUNCS([getpseudotty _getpty getpt ptsname_r])
AC_CHECK_FUNCS([setpgid setjmp sbrk])
AC_CHECK_FUNCS([rld_load])
AC_CHECK_FUNCS([times getrusage])
AC_CHECK_FUNCS([vm_allocate])
AC_CHECK_FUNCS([fork popen waitpid wait4])
AC_CHECK_FUNCS([dup2 floor memmove mkdir pow rint])
AC_CHECK_FUNCS([ttyname strerror select])
AC_CHECK_FUNCS([access lstat unlink mkdir mkdtemp mkstemp rmdir link rename])
AC_CHECK_FUNCS([symlink readlink chmod fchmod chown fchown lchown dup dup2 mknod mkfifo])
AC_CHECK_FUNCS([socket listen accept connect bind getsockname recv recvfrom recvmsg send])
AC_CHECK_FUNCS([sendto sendmsg getsockopt setsockopt getprotobyname gethostbyname])
AC_CHECK_FUNCS([opendir closedir dirfd readdir rewinddir seekdir telldir])

### singular
AC_PATH_PROGS([SINGULAR],
              [Singular singular],
              [no])
AS_IF([test "$SINGULAR" = no], [AC_MSG_ERROR([cannot find singular.])])
AS_IF([$SINGULAR -q],,[AC_MSG_ERROR([singular doesn't work.])])
AC_MSG_CHECKING([where Singular has its libs])
SINGULAR_LIB=`echo 'system("SingularLib");' | $SINGULAR -q --no-rc --no-stdlib -b`
AC_MSG_RESULT([$SINGULAR_LIB])
if test ! -e $SINGULAR_LIB/algebra.lib; then
    AC_MSG_ERROR([unexpected singular install, cannot find lib])
fi
AC_SUBST([SINGULAR_LIB])
SINGULAR_BIN=`dirname $SINGULAR`
AC_SUBST([SINGULAR_BIN])

# libsingular.so might not match singular executable (does it have to?)
# anything else (retrieving from `ldd $SINGULAR`?) would be a hack.
# (lets hope this is not necessary at all, one day)
SINGULAR_SO=`readlink -f \`$CC $LDFLAGS -print-file-name=libsingular.so\``

AC_SUBST([SINGULAR_SO])

# singular needs singular and factory
OLDCPPFLAGS=$CPPFLAGS
SINGULAR_INCLUDES="-I$SINGULAR_BIN/../include/singular -I$SINGULAR_BIN/../include/factory"
CPPFLAGS="$SINGULAR_INCLUDES $CPPFLAGS"
AC_CHECK_HEADER([si_gmp.h],,[AC_MSG_ERROR(cannot find singular header)])
AC_SUBST([SINGULAR_INCLUDES])
CPPFLAGS=$OLDCPPFLAGS

FACTORY_INCLUDES="-I$SINGULAR_BIN/../include/factory"
AC_SUBST([FACTORY_INCLUDES])

### header/library check
# this is superficial.
# the focus is on not breaking sagelib on sage (the distribution)
# in some cases, passing prefixes/include directories (AC_ARG_WITH) might be
# preferrable, but this will not make sense before sagelib/spkg-install is done
# right...

# gmp and mpir
AS_IF([ test x$gmp_or_mpir != xgmp ], [
       AC_CHECK_LIB([mpir], [__gmp_version],[gmp_or_mpir=mpir; GMP_or_MPIR=MPIR] ) ])

AS_IF([ test x$gmp_or_mpir != xmpir ], [
       AC_CHECK_LIB([gmp], [__gmp_version],[gmp_or_mpir=gmp; GMP_or_MPIR=GMP],[AC_MSG_ERROR([there is no gmp lib])])])

AS_IF([test -z "$gmp_or_mpir" ], [AC_MSG_ERROR([there is no mpir or gmp lib])])

SAGE_CFLAGS="-w -O2 -DUSE_"$GMP_or_MPIR
AC_SUBST([SAGE_CFLAGS])

# lets assume nobody wants to mix up {gmp,mpir} with {mpir,gmp}xx
GMP_OR_MPIR_LIBS=-l$gmp_or_mpir
GMP_OR_MPIR_XXLIBS="-l$gmp_or_mpir"xx
AC_SUBST(GMP_OR_MPIR_LIBS)
AC_SUBST(GMP_OR_MPIR_XXLIBS)

# just in case.
AM_CONDITIONAL([USE_GMP], [ test x$gmp_or_mpir = xgmp ])
AM_CONDITIONAL([USE_MPIR], [ test x$gmp_or_mpir = xmpir ])

# Checks for header files.
AS_IF([test $gmp_or_mpir = gmp ],
       [AC_CHECK_HEADER([gmp.h],,[AC_MSG_ERROR([no gmp header found])])
        AC_LANG_PUSH([C++])
        AC_CHECK_HEADER([gmpxx.h],,[AC_MSG_ERROR([no gmpxx header found])])
        AC_LANG_POP()
       ])
AS_IF([test $gmp_or_mpir = mpir ],
       [AC_CHECK_HEADER([mpir.h],,[AC_MSG_ERROR([no mpir header found])])
        AC_LANG_PUSH([C++])
        AC_CHECK_HEADER([mpirxx.h],,[AC_MSG_ERROR([no mpirxx header found])])
        AC_LANG_POP()
       ])

AS_IF([ test x$gmp_or_mpir = xgmp ],[
    AC_DEFINE(USE_GMP, [1], [use gmp (not mpir)]) ])
AS_IF([ test x$gmp_or_mpir = xmpir ],[
    AC_DEFINE(USE_MPIR, [1], [use mpir (not gmp)]) ])

# later: --with-ecl?
# AS_IF([test ! x$SAGE_LOCAL = xNONE],
#       [LIBECL_INCLUDES=-I$SAGE_LOCAL/include/ecl],
#       [LIBECL_INCLUDES=-I/usr/include/ecl])
# CPPFLAGS+=" $LIBECL_INCLUDES"
AC_CHECK_HEADER([ecl/ecl.h],,[AC_MSG_ERROR(cannot find ecl header)])
AC_SUBST(LIBECL_INCLUDES)
CPPFLAGS=$OLDCPPFLAGS

# # later: --with-flint?
# AS_IF([test ! x$SAGE_LOCAL = xNONE],
#       [LIBFLINT_INCLUDES=-I$SAGE_LOCAL/include/flint],
#       [LIBFLINT_INCLUDES=-I/usr/include/flint])
# CPPFLAGS+=" $LIBFLINT_INCLUDES"
AC_CHECK_HEADER([flint/flint.h],,[AC_MSG_ERROR(cannot find flint header)])
AC_SUBST(LIBFLINT_INCLUDES)

CPPFLAGS="$OLDCPPFLAGS $PYTHON_INCLUDES"
dnl FIXME: this check fails for some reason, should be AC_MSG_ERROR
AC_CHECK_HEADER([pynac/order.h],,[AC_MSG_NOTICE(cannot find pynac header)])

CPPFLAGS=$OLDCPPFLAGS

AC_LANG_PUSH([C++])
# note: if you have lfunction headers in a different location, you must take
# care of CPPFLAGS manually
AS_IF([test ! x$SAGE_LOCAL = xNONE -a -f $SAGE_LOCAL/include/libLfunction/L.h],
      [LIBLFUNCTION_INCLUDES=-I$SAGE_LOCAL/include/libLfunction],
      [test -f /usr/include/Lfunction/L.h],
      [LIBLFUNCTION_INCLUDES=-I/usr/include/Lfunction],
      [test -f /usr/local/include/Lfunction/L.h],
      [LIBLFUNCTION_INCLUDES=-I/usr/local/include/Lfunction])
CPPFLAGS="$OLDCPPFLAGS $LIBLFUNCTION_INCLUDES"
AC_CHECK_HEADER([L.h],,[AC_MSG_ERROR(cannot find Lfunction header)])
CPPFLAGS=$OLDCPPFLAGS
AC_LANG_POP
AC_SUBST([LIBLFUNCTION_INCLUDES])

# these will be used without *_INCLUDES, just check
AC_CHECK_HEADER([cliquer/cliquer.h],,[AC_MSG_ERROR(cannot find cliquer header)])
AC_CHECK_HEADER([zn_poly/zn_poly.h],,[AC_MSG_ERROR(cannot find zn_poly header)])
AC_CHECK_HEADER([symmetrica/macro.h],,[AC_MSG_ERROR(cannot find symmetrica header)])
AC_CHECK_HEADER([ecm.h],,[AC_MSG_ERROR(cannot find ecm header)])
AC_CHECK_HEADER([gd.h],,[AC_MSG_ERROR(cannot find gd header)])

AC_LANG_PUSH([C++])
AC_CHECK_HEADER([polybori/groebner/add_up.h],,[AC_MSG_ERROR(cannot find polybori groebner header)])

# linbox-sage.h doesn't work right fix later, maybe
# AC_CHECK_HEADER([linbox/linbox-sage.h],,[AC_MSG_ERROR(cannot find linbox-sage header)])

AC_CHECK_HEADER([fplll.h],,[AC_MSG_ERROR(cannot find fplll header)])
AC_CHECK_HEADER([fplll/fplll.h],,[AC_MSG_ERROR(cannot find fplll header)])
AC_LANG_POP
AC_SUBST([FPLLL_INCLUDES]) dnl currently unneeded

AC_CHECK_HEADER([coxeter3.h], [HAVE_LIBCOXETER3=yes], [HAVE_LIBCOXETER3=no])
AC_SUBST([HAVE_LIBCOXETER3])
AM_CONDITIONAL([USE_LIBCOXETER3], [ test $HAVE_LIBCOXETER3 = yes ])

AC_CHECK_HEADER([fes.h], [HAVE_LIBFES=yes], [HAVE_LIBFES=no])
AC_SUBST([HAVE_LIBFES])
AM_CONDITIONAL([USE_LIBFES], [ test $HAVE_LIBFES = yes ])

old_CPPFLAGS=$CPPFLAGS
HAVE_LIBLRCALC=no
# sage (the distribution) looks in $SAGE_LOCAL first.
# use --with-lrcalc to specify something different
AC_ARG_WITH(lrcalc, AS_HELP_STRING([--with-lrcalc=prefix].[use lrcalc from prefix]),
    [  case $with_lrcalc in
          yes)
            ;;
          *)
            LIBLRCALC_INCLUDES="-I$with_lrcalc/include/lrcalc"
            ;;
        esac
        ],[ ]
)
# use --with-lrcalc-include to override
AC_ARG_WITH(lrcalc-include, AS_HELP_STRING([--with-lrcalc-include=includedir].[use lrcalc from includedir]),
    [  case $with_lrcalc in
          yes)
            ;;
          *)
            LIBLRCALC_INCLUDES="-I$with_lrcalc/include"
            ;;
        esac
        ],[ ]
)

AS_IF([test -n "$LIBLRCALC_INCLUDES"],[],
      [test ! x$SAGE_LOCAL = xNONE -a -f $SAGE_LOCAL/include/lrcalc/symfcn.h],
      [LIBLRCALC_INCLUDES="-I$SAGE_LOCAL/include/lrcalc"],
      [LIBLRCALC_INCLUDES="-I/usr/include/lrcalc"])
CPPFLAGS="$CPPFLAGS $LIBLRCALC_INCLUDES"
AC_CHECK_HEADER([symfcn.h],
                [HAVE_LIBLRCALC=yes],
                [AC_MSG_ERROR([could not find lrcalc])])
CPPFLAGS=$old_CPPFLAGS
AC_SUBST([LIBLRCALC_INCLUDES])
AC_SUBST([HAVE_LIBLRCALC])
dnl AC_CHECK_LIB([lrcalc], [lrcalc], lrcalc=1))
AM_CONDITIONAL([USE_LIBLRCALC], [ test $HAVE_LIBLRCALC = yes ])

# FIXME: also check for gap
HAVE_LIBGAP=no
AC_CHECK_HEADER([libgap/libgap.h], [HAVE_LIBGAP=yes], [AC_MSG_NOTICE(cannot find libgap header);
   SAGE_CHECK_PYTHON_MODULE([ply.lex],[AC_MSG_ERROR([need ply to hack gap])]) ])
AM_CONDITIONAL([HAVE_LIBGAP], [ test $HAVE_LIBGAP = yes ])

AC_CHECK_HEADER([gurobi.h], [HAVE_LIBGUROBI=yes], [HAVE_LIBGUROBI=no])
AC_SUBST([HAVE_LIBGUROBI])
AM_CONDITIONAL([USE_LIBGUROBI], [ test $HAVE_LIBGUROBI = yes ])

AC_CHECK_HEADER([cbc.h], [HAVE_LIBCBC=yes], [HAVE_LIBCBC=no])
AC_SUBST([HAVE_LIBCBC])
AM_CONDITIONAL([USE_LIBCBC], [ test $HAVE_LIBCBC = yes ])

AC_CHECK_HEADER([cplex.h], [HAVE_LIBCPLEX=yes], [HAVE_LIBCPLEX=no])
AC_SUBST([HAVE_LIBCPLEX])
AM_CONDITIONAL([USE_LIBCPLEX], [ test $HAVE_LIBCPLEX = yes ])

AC_CHECK_HEADER([coin.h], [HAVE_LIBCOIN=yes], [HAVE_LIBCOIN=no])
AC_SUBST([HAVE_LIBCOIN])
AM_CONDITIONAL([USE_LIBCOIN], [ test $HAVE_LIBCOIN = yes ])

AC_CHECK_HEADER([Solver.h], [HAVE_LIBCRYPTOMINISAT=yes], [HAVE_LIBCRYPTOMINISAT=no])
AC_SUBST([HAVE_LIBCRYPTOMINISAT])
AM_CONDITIONAL([USE_LIBCRYPTOMINISAT], [ test $HAVE_LIBCRYPTOMINISAT = yes ])

AC_CHECK_LIB([readline], [main])
AC_SUBST(HAVE_LIBREADLINE)

# not fetching from environment. use --with-blas
blas=
blas2=

## Choose cblas library -- note -- make sure to update sage/misc/cython.py
## if you change this!!
AC_ARG_WITH(blas,   --with-blas   choose cblas, [blas=$with_blas], [])

if test -n "${blas}"; then
  AC_CHECK_LIB([$blas], [cblas_dgemv], [blas2=$blas],
               [ AC_MSG_ERROR(user specified blas doesnt exist) ]
              )
fi

if test -z "${blas}"; then
  # we cannot know whether this is sage-distributed atlas.
  # and it doesnt matter.
  AC_CHECK_LIB([atlas], [ATL_scpsc], [blas=cblas; blas2=atlas], [])
  # some random symbol   ^
fi

if test -z "${blas}"; then
  AC_CHECK_LIB([cblas], [cblas_dgemv], [blas=cblas; blas2=atlas], [])
fi

if test -z "${blas}"; then
  AC_CHECK_LIB([blas], [cblas_dgemv], [blas=gslcblas; blas2=gslcblas], [])
fi

if test -z "${blas}"; then
  AC_CHECK_LIB([gslcblas], [cblas_dgemv], [blas=gslcblas; blas2=gslcblas], [])
fi

if test -z "${blas}"; then
  AC_MSG_ERROR([couldn't find any blas])
fi

#if -z os.environ.has_key('SAGE_BLAS'):
#    BLAS=os.environ['SAGE_BLAS']
#    BLAS2=os.environ['SAGE_BLAS']
#   AC_MSG_RESULT(from environment: ${BLAS})
#elif os.path.exists('%s/lib/libatlas.so'%os.environ['SAGE_LOCAL']):
#    BLAS='cblas'
#    BLAS2='atlas'
#elif os.path.exists('/usr/lib/libcblas.dylib') or \
#     os.path.exists('/usr/lib/libcblas.so'):
#    BLAS='cblas'
#    BLAS2='atlas'
#elif os.path.exists('/usr/lib/libblas.dll.a'):
#    BLAS='gslcblas'
#    BLAS2='gslcblas'
#else:
#    # This is very slow  (?), but *guaranteed* to be available.
#    BLAS='gslcblas'
#    BLAS2='gslcblas'
#   AC_MSG_RESULT(fallback to gslcblas)

# random symbol (fixme?)
AC_CHECK_LIB([ec], [_Z9show_timeRSo], [],
             [ AC_MSG_ERROR([Could not find ec library]) ])

#########################################################
### M4RI flags
#########################################################

PKG_CHECK_MODULES([M4RI],[m4ri],,AC_MSG_ERROR([no m4ri config found]))
M4RI_CFLAGS+=" -std=c99"
AC_SUBST([M4RI_CFLAGS])
AC_SUBST([LIBM4RI_INCLUDES])

PKG_CHECK_MODULES([M4RIE],[m4rie],,
                  [AC_MSG_NOTICE([no m4rie config found, lets try without...])]
                  [M4RIE_CFLAGS=-I/usr/include/m4rie]
                 )
AC_SUBST([M4RIE_CFLAGS])
AC_SUBST([LIBM4RIE_INCLUDES])

old_CPPFLAGS=$CPPFLAGS
dnl hack: pass -c99 via CPPFLAGS...
CPPFLAGS+=" $M4RIE_CFLAGS -std=c99"
AC_CHECK_HEADER([mzed.h],,[AC_MSG_ERROR(cannot find m4rie header)])
CPPFLAGS=$old_CPPFLAGS

AC_SUBST(BLAS, ${blas})
AC_SUBST(BLAS2, ${blas2})

# ouch?
my_CPPFLAGS+=' -I$(top_srcdir)/..'
AC_SUBST([my_CPPFLAGS])

# more paths (compile time)
LIBCSAGE_INCLUDES='-I$(top_builddir)/../c_lib/include'
AS_IF([ test ! x$srcdir = "x."],
      [LIBCSAGE_INCLUDES+=' -I$(top_srcdir)/../c_lib/include'])
AC_SUBST([LIBCSAGE_INCLUDES])

# run time paths
SAGE_LIBCSAGE_INCLUDES='-I${prefix}/include/csage'
AC_SUBST([SAGE_LIBCSAGE_INCLUDES])
SAGE_SRC="${prefix}/share/sage/src"
AC_SUBST([SAGE_SRC])

# not necessary with $( ) support in env.py
# AS_AC_EXPAND(LIBCSAGE_INCLUDES, "-I$pkgincludedir")

# fixme: use AC_...
NUMPY_INCLUDEDIR=$SAGE_LOCAL/lib/python/site-packages/numpy/core/include
NUMPY_INCLUDES=-I$NUMPY_INCLUDEDIR
AC_SUBST(NUMPY_INCLUDEDIR)
AC_SUBST(NUMPY_INCLUDES)

AS_IF([test ! x$SAGE_LOCAL = xNONE],
      [SAGE_INC=$SAGE_LOCAL/include],
      [SAGE_INC=NONE])
AC_SUBST(SAGE_INC)

SAGE_DIR=$srcdir
AC_SUBST(SAGE_DIR)

AS_AC_EXPAND(SAGE_ETC, $sysconfdir)

AC_CHECK_SIZEOF(void *)
# AC_CHECK_ALIGNOF([long])
AC_HEADER_STDBOOL
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_HEADER_DIRENT
AC_HEADER_MAJOR
AC_HEADER_STAT
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_PID_T
AC_TYPE_SIGNAL
AC_CHECK_FUNCS([signal])
AC_TYPE_SIGNAL
# weird check for signal handling
if test "$ac_cv_type_signal" = "void"; then
    AC_DEFINE([HAVE_SIGNAL_VOID], [1] ,[void as signal type??])
fi

AC_DEFINE(HAVE_STRERROR, 1, [Everybody has that])
AC_DEFINE(HAVE_SYS_ERRLIST, 1, [Or GAP will define its own])

#### data packages ####
# default: sort of pkgdatadir, maybe just sage_share
AC_ARG_WITH(dbdir,
            AS_HELP_STRING([--with-dbdir=DIR].[use databases from DIR]),
            [dbdir=$withval],
            [dbdir=$SAGE_LOCAL/share])
SAGE_SHARE="$dbdir"

AS_IF([test -d $dbdir/conway_polynomials],
      [CONWAY_DATA_DIR=$dbdir/conway_polynomials],
      [test -d /usr/share/sage/conway_polynomials], #gentoo
      [CONWAY_DATA_DIR=/usr/share/sage/conway_polynomials],
      [test -d /usr/share/sagemath/conway_polynomials], #debian
      [CONWAY_DATA_DIR=/usr/share/sagemath/conway_polynomials],
      [CONWAY_DATA_DIR=$SAGE_LOCAL/share/conway_polynomials]) dnl hack!
AC_SUBST([CONWAY_DATA_DIR])
AM_SUBST_NOTMAKE([CONWAY_DATA_DIR])

AS_IF([test -d $dbdir/graphs],
      [GRAPHS_DATA_DIR=$SAGE_LOCAL/share/graphs],
      [test -d /usr/share/sage/graphs],
      [GRAPHS_DATA_DIR=/usr/share/sage/graphs],
      [test -d /usr/share/sagemath/graphs],
      [GRAPHS_DATA_DIR=/usr/share/sagemath/graphs],
      [AC_MSG_ERROR([cannot find graphs database])])
AC_SUBST([GRAPHS_DATA_DIR])
AM_SUBST_NOTMAKE([GRAPHS_DATA_DIR])

AS_IF([test -d $dbdir/reflexive_polytopes],
      [POLYTOPE_DATA_DIR=$SAGE_LOCAL/share/reflexive_polytopes],
      [test -d /usr/share/sage/reflexive_polytopes],
      [POLYTOPE_DATA_DIR=/usr/share/sage/reflexive_polytopes],
      [test -d /usr/share/sagemath/reflexive_polytopes],
      [POLYTOPE_DATA_DIR=/usr/share/sagemath/reflexive_polytopes],
      [AC_MSG_ERROR([cannot find polytopes database])])
AC_SUBST([POLYTOPE_DATA_DIR])
AM_SUBST_NOTMAKE([POLYTOPE_DATA_DIR])

AS_IF([test -d $dbdir/ellcurves],
      [ELLCURVE_DATA_DIR=$dbdir],
      [test -d /usr/share/sage/ellcurves],
      [ELLCURVE_DATA_DIR=/usr/share/sage],
      [test -d /usr/share/sagemath/ellcurves],
      [ELLCURVE_DATA_DIR=/usr/share/sagemath],
      [AC_MSG_ERROR([cannot find elliptic curves])])
AC_SUBST([ELLCURVE_DATA_DIR])
AM_SUBST_NOTMAKE([ELLCURVE_DATA_DIR])

# FIXME: gentoo?
# FIXME: require GAP_INCLUDEDIR unless HAVE_LIBGAP
GAP_INCLUDEDIR=/usr/include/gap
AS_IF([test -d $SAGE_LOCAL/gap],
      [GAP_DATA_DIR=$SAGE_LOCAL/gap/latest;
       GAP_INCLUDEDIR=$SAGE_LOCAL/gap/latest/src],
      [test -d /usr/share/gap],
      [GAP_DATA_DIR=/usr/share/gap],
      [AC_MSG_ERROR([cannot find gap])])
AC_SUBST([GAP_DATA_DIR])
#AC_SUBST([GAP_INCLUDES])
AC_SUBST([GAP_INCLUDEDIR])

AC_MSG_CHECKING([for gap libdir])
AS_IF([test -x $SAGE_LOCAL/gap/latest/bin/*/gap],
      [GAP_LIB_DIR=`dirname $SAGE_LOCAL/gap/latest/bin/*/gap`],
      [test -d /usr/lib/gap/bin],
      [GAP_LIB_DIR=/usr/lib/gap/bin],
      [AC_MSG_ERROR([cannot find gap])])
AC_SUBST([GAP_LIB_DIR])
AC_MSG_RESULT([$GAP_LIB_DIR])

#### output files

AC_CONFIG_FILES([
    setup.py
    module_list.py
    Makefile
    algebras/letterplace/Makefile
    algebras/Makefile
    algebras/quatalg/Makefile
    calculus/Makefile
    categories/examples/Makefile
    categories/Makefile
    coding/Makefile
    combinat/Makefile
    combinat/crystals/Makefile
    combinat/matrices/Makefile
    combinat/words/Makefile
    crypto/Makefile
    databases/Makefile
    doctest/Makefile
    dynamics/Makefile
    ext/interpreters/Makefile
    ext/Makefile
    finance/Makefile
    functions/Makefile
    games/Makefile
    geometry/Makefile
    geometry/triangulation/Makefile
    graphs/base/Makefile
    graphs/graph_decompositions/Makefile
    graphs/Makefile
    graphs/modular_decomposition/Makefile
    groups/abelian_gps/Makefile
    groups/Makefile
    groups/matrix_gps/Makefile
    groups/perm_gps/Makefile
    groups/perm_gps/partn_ref/Makefile
    gsl/Makefile
    homology/Makefile
    interacts/Makefile
    interfaces/Makefile
    lfunctions/Makefile
    libs/coxeter3/Makefile
    libs/cremona/Makefile
    libs/flint/Makefile
    libs/fplll/Makefile
    libs/gap/Makefile
    libs/lcalc/Makefile
    libs/lrcalc/Makefile
    libs/linbox/Makefile
    libs/Makefile
    libs/mpmath/Makefile
    libs/mwrank/Makefile
    libs/ntl/Makefile
    libs/pari/Makefile
    libs/singular/Makefile
    libs/symmetrica/Makefile
    logic/Makefile
    matrix/Makefile
    matroids/Makefile
    media/Makefile
    misc/Makefile
    modular/Makefile
    modular/arithgroup/Makefile
    modular/modform/Makefile
    modular/modsym/Makefile
    modules/fg_pid/Makefile
    modules/Makefile
    monoids/Makefile
    numerical/backends/Makefile
    numerical/Makefile
    parallel/Makefile
    plot/Makefile
    plot/plot3d/Makefile
    probability/Makefile
    quadratic_forms/Makefile
    rings/finite_rings/Makefile
    rings/function_field/Makefile
    rings/Makefile
    rings/number_field/Makefile
    rings/padics/Makefile
    rings/polynomial/Makefile
    rings/polynomial/padics/Makefile
    rings/semirings/Makefile
    rings/universal_cyclotomic_field/Makefile
    sandpiles/Makefile
    sat/Makefile
    sat/solvers/cryptominisat/Makefile
    sat/solvers/Makefile
    schemes/elliptic_curves/Makefile
    schemes/hyperelliptic_curves/Makefile
    schemes/Makefile
    schemes/toric/Makefile
    server/Makefile
    server/notebook/Makefile
    sets/Makefile
    stats/hmm/Makefile
    stats/Makefile
    structure/Makefile
    structure/proof/Makefile
    symbolic/integration/Makefile
    symbolic/Makefile
    tensor/Makefile
    tests/Makefile
    env.py
])

dnl create dummy files for cython file dependencies
dnl this is currently not natively supported by automake
AC_CONFIG_COMMANDS([cythondep],
                   [$ac_top_srcdir/cythondep AMDEP_TRUE=$AMDEP_TRUE])

AC_SUBST(PYTHONPATH)

AC_OUTPUT

AC_MSG_RESULT([
** Configuration summary for $PACKAGE $VERSION:

   prefix:     $prefix
   exec_prefix:$exec_prefix
   libdir:     $libdir
   CPPFLAGS:   $CPPFLAGS
   CFLAGS:     $CFLAGS
   CXXFLAGS:   $CXXFLAGS
   LDFLAGS:    $LDFLAGS
   LIBS:       $LIBS
   c_lib:      $clib
   BLAS:       $blas
   BLAS2:      $blas2
   PYTHON:     $PYTHON
   gmp or mpir:$gmp_or_mpir

   flags:
   M4RI_CFLAGS:  $M4RI_CFLAGS
   M4RIE_CFLAGS: $M4RIE_CFLAGS

   directories/stuff
   gap data      $GAP_DATA_DIR
   gap lib       $GAP_LIB_DIR
   singular lib  $SINGULAR_LIB
   singular so   $SINGULAR_SO

   optional packages
   lrcalc:       $HAVE_LIBLRCALC
   cbc:          $HAVE_LIBCBC
   coxeter:      $HAVE_LIBCOXETER3
   gurobi:       $HAVE_LIBGUROBI
   cryptominisat $HAVE_LIBCRYPTOMINISAT
   libgap        $HAVE_LIBGAP
])

# vim:et:sw=2:
