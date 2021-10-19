SAGE_SPKG_CONFIGURE([singular], [
  SAGE_SPKG_DEPCHECK([gmp mpir ntl flint readline mpfr cddlib], [

    AC_PATH_PROG([SINGULAR_BIN], [Singular])
    AS_IF([test -z "${SINGULAR_BIN}"], [sage_spkg_install_singular=yes])

    dnl Use pkg-config to ensure that Singular is new enough.
    PKG_CHECK_MODULES([SINGULAR],
                      [Singular >= 4.1.1],
                      [],
                      [sage_spkg_install_singular=yes])

    dnl The acl_shlibext variable is set in the top-level configure.ac,
    dnl and is ultimately substituted into sage_conf as SHLIBEXT.
    SINGULAR_LIB_BASE="libSingular"
    SINGULAR_LIB_FILE="${SINGULAR_LIB_BASE}.${acl_shlibext}"

    AC_MSG_CHECKING([if we can dlopen($SINGULAR_LIB_FILE)])
    ORIG_LIBS="${LIBS}"
    LIBS="${LIBS} -ldl"
    AC_LANG_PUSH(C)

    dnl if we can dlopen() it, substitute the name for sage_conf;
    dnl otherwise, fall back to using the SPKG.
    AC_RUN_IFELSE(
      [AC_LANG_PROGRAM(
        [[#include <dlfcn.h>]],
        [[void* h = dlopen("${SINGULAR_LIB_FILE}", RTLD_LAZY | RTLD_GLOBAL);
          if (h == 0) { return 1; } else { return dlclose(h); }]]
      )],
      [AC_SUBST(SINGULAR_LIB_BASE, "${SINGULAR_LIB_BASE}")],
      [sage_spkg_install_singular=yes]
    )

    AC_LANG_POP()
    LIBS="${ORIG_LIBS}"
  ])
],[],[],[
  dnl Post-check phase
  dnl We make the sage_conf substitutions here, because the "default"
  dnl substitution needs to be made even if we skipped the system-Singular
  dnl checks themselves.

  dnl If we're using the SPKG, we might as well use the FULL path to the
  dnl library, because we know it.
  AS_IF([test "x${sage_spkg_install_singular}" = "xyes"],
    [AC_SUBST(SINGULAR_LIB_BASE, '$SAGE_LOCAL/lib/libSingular')],
    [AC_SUBST(SINGULAR_LIB_BASE, "${SINGULAR_LIB_BASE}")]
  )
])
