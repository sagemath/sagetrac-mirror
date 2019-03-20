SAGE_SPKG_CONFIGURE([mpfr], [
    AC_ARG_WITH([mpfr],
    [AS_HELP_STRING([--with-mpfr=system],
        [use the system MPFR, if possible (default)])]
    [AS_HELP_STRING([--with-mpfr=install],
        [use the Sage SPKG for MPFR])])

dnl Just part the options here
    case "$with_mpfr" in
        system) ;;
        install) ;;
        "") with_mpfr=system;;
        *)
            AC_MSG_ERROR([allowed values for --with-mpfr are system and install]);;
    esac
    
    case "$with_mpfr" in
        system)
            AC_CHECK_HEADER(mpfr.h, [], [sage_spkg_install_mpfr=yes])
        dnl mpfr_free_pool appeared in r11922 (Dec 2017) on MPFR svn 
            AC_SEARCH_LIBS([mpfr_free_pool], [mpfr], [break], [sage_spkg_install_mpfr=yes])

            if test x$sage_spkg_install_mpfr = xyes; then
               AC_SUBST(SAGE_MPFR_PREFIX, ['$SAGE_LOCAL'])
           dnl AC_SUBST(SAGE_MPFR_INCLUDE, ['$SAGE_LOCAL/include']) --- not used
               AC_MSG_RESULT([using Sage's mpfr SPKG])
            else
           dnl If found, we want to get the absolute path to where we
           dnl found it for use with some packages want
           dnl this information at configure time
               AX_ABSOLUTE_HEADER([mpfr.h])
               if test x$gl_cv_absolute_mpfr_h = x; then
                   AC_MSG_ERROR(m4_normalize([
                     failed to find absolute path to mpfr.h despite it being reported found
                   ]))
               fi
           dnl AC_SUBST(SAGE_MPFR_INCLUDE, [`AS_DIRNAME($gl_cv_absolute_mpfr_h)`])
               AC_SUBST(SAGE_MPFR_PREFIX, [''])
               AC_MSG_RESULT([using mpfr library from the system])
            fi
            ;;
        install)
            sage_spkg_install_mpfr=yes
            AC_SUBST(SAGE_MPFR_PREFIX, ['$SAGE_LOCAL'])
        dnl AC_SUBST(SAGE_MPFR_INCLUDE, ['$SAGE_LOCAL/include'])
            AC_MSG_RESULT([using Sage's mpfr SPKG])
            ;;
    esac

])
