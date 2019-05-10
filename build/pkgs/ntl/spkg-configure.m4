SAGE_SPKG_CONFIGURE([ntl], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_GMP])
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_GF2X])
    AC_MSG_CHECKING([installing gmp/mpir? ])
    if test x$sage_spkg_install_mpir = xyes -o x$sage_spkg_install_gmp = xyes; then
        AC_MSG_RESULT([yes; install ntl as well])
        sage_spkg_install_ntl=yes
    else
        AC_MSG_RESULT([no])
    fi

    AC_MSG_CHECKING([installing gf2x? ])
    if test x$sage_spkg_install_gf2x = xyes; then
        AC_MSG_RESULT([yes; install ntl as well])
        sage_spkg_install_ntl=yes
    else
        AC_MSG_RESULT([no])
    fi

    if test x$sage_spkg_install_gf2x != xyes; then
        AC_CHECK_HEADER([NTL/ZZ.h], [], [sage_spkg_install_ntl=yes])
        AC_LINK_IFELSE([
	        AC_LANG_PROGRAM([[#include <NTL/ZZ.h>]],
                            [[NTL::ZZ a;]]
            )], [], [sage_spkg_install_ntl=yes])
	    AC_RUN_IFELSE([
	        AC_LANG_PROGRAM(
            [[#include <NTL/version.h>]],
	        [[if (NTL_MAJOR_VERSION<$ntl_major_ver) return 1; \
              if (NTL_MINOR_VERSION<$ntl_minor_ver) return 1; \
			  return 0;]]
            )], [],  [sage_spkg_install_ntl=yes])
    fi
], [], [], [
    if test x$sage_spkg_install_ntl = xyes; then
        AC_SUBST(SAGE_NTL_PREFIX, ['$SAGE_LOCAL'])
        AC_MSG_RESULT([using Sage's ntl SPKG])
    else
        AC_SUBST(SAGE_NTL_PREFIX, [''])
        AC_MSG_RESULT([using ntl library from the system])
    fi
])

