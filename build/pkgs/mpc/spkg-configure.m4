SAGE_SPKG_CONFIGURE([mpc], [
    AC_ARG_WITH([mpc],
    [AS_HELP_STRING([--with-mpc=system],
        [use the system MPC, if possible (default)])]
    [AS_HELP_STRING([--with-mpc=install],
        [use the Sage SPKG for MPC])])

dnl Just part the options here
    case "$with_mpc" in
        system) ;;
        install) ;;
        "") with_mpc=system;;
        *)
            AC_MSG_ERROR([allowed values for --with-mpc are system and install]);;
    esac
    
    case "$with_mpc" in
        system)
            AC_CHECK_HEADER(mpc.h, [], [sage_spkg_install_mpc=yes])
        dnl mpc_cmp_abs appeared in MPC 1.1.0
            AC_SEARCH_LIBS([mpc_cmp_abs], [mpc], [break], [sage_spkg_install_mpc=yes])

            if test x$sage_spkg_install_mpc = xyes; then
               AC_SUBST(SAGE_MPC_PREFIX, ['$SAGE_LOCAL'])
           dnl AC_SUBST(SAGE_MPC_INCLUDE, ['$SAGE_LOCAL/include']) --- not used
               AC_MSG_RESULT([using Sage's mpc SPKG])
            else
           dnl If found, we want to get the absolute path to where we
           dnl found it for use with some packages want
           dnl this information at configure time
               AX_ABSOLUTE_HEADER([mpc.h])
               if test x$gl_cv_absolute_mpc_h = x; then
                   AC_MSG_ERROR(m4_normalize([
                     failed to find absolute path to mpc.h despite it being reported found
                   ]))
               fi
           dnl AC_SUBST(SAGE_MPC_INCLUDE, [`AS_DIRNAME($gl_cv_absolute_mpc_h)`])
               AC_SUBST(SAGE_MPC_PREFIX, [''])
               AC_MSG_RESULT([using mpc library from the system])
            fi
            ;;
        install)
            sage_spkg_install_mpc=yes
            AC_SUBST(SAGE_MPC_PREFIX, ['$SAGE_LOCAL'])
        dnl AC_SUBST(SAGE_MPC_INCLUDE, ['$SAGE_LOCAL/include'])
            AC_MSG_RESULT([using Sage's mpc SPKG])
            ;;
    esac

])
