SAGE_SPKG_CONFIGURE([ntl], [
    AC_ARG_WITH([ntl],
    [AS_HELP_STRING([--with-ntl=system],
        [use the system NTL, if possible (default)])]
    [AS_HELP_STRING([--with-ntl=install],
        [use the Sage SPKG for NTL])])
    
    ntl_major_ver=10 dnl   At least this major version
    ntl_minor_ver=3  dnl   At least this minor version

dnl Just part the options here
    case "$with_ntl" in
        system) ;;
        install) ;;
        "") with_ntl=system;;
        *)
            AC_MSG_ERROR([allowed values for --with-ntl are system and install]);;
    esac
    
    case "$with_ntl" in
        system)
dnl           LB_CHECK_NTL(10.3, [sage_spkg_install_ntl=no], [sage_spkg_install_ntl=yes])
            AC_CHECK_HEADER([NTL/ZZ.h], [], [sage_spkg_install_ntl=yes])
            AC_LINK_IFELSE([
	     AC_LANG_PROGRAM([[#include <NTL/ZZ.h>]],
                            [[NTL::ZZ a;]])],
			[], [sage_spkg_install_ntl=yes])
	    AC_RUN_IFELSE([
	     AC_LANG_PROGRAM([[#include <NTL/version.h>]],
	                  [[if (NTL_MAJOR_VERSION<$ntl_major_ver) return 1; \
			    if (NTL_MINOR_VERSION<$ntl_minor_ver) return 1; \
			    return 0;]])],
			[],  [sage_spkg_install_ntl=yes])

            if test x$sage_spkg_install_ntl = xyes; then
               AC_SUBST(SAGE_NTL_PREFIX, ['$SAGE_LOCAL'])
           dnl AC_SUBST(SAGE_NTL_INCLUDE, ['$SAGE_LOCAL/include']) --- not used
               AC_MSG_RESULT([using Sage's ntl SPKG])
            else
               AC_SUBST(SAGE_NTL_PREFIX, [''])
               AC_MSG_RESULT([using ntl library from the system])
            fi
            ;;
        install)
            sage_spkg_install_ntl=yes
            AC_SUBST(SAGE_NTL_PREFIX, ['$SAGE_LOCAL'])
        dnl AC_SUBST(SAGE_NTL_INCLUDE, ['$SAGE_LOCAL/include'])
            AC_MSG_RESULT([using Sage's ntl SPKG])
            ;;
    esac

])

