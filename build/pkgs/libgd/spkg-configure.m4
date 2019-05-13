SAGE_SPKG_CONFIGURE([libgd], [
    dnl Check for libpng is done in the following macro, no need to repeat it here
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_FREETYPE])
    AC_MSG_CHECKING([Installing freetype? ])
    if test x$sage_spkg_install_freetype= xyes; then
      AC_MSG_RESULT([Yes. Install libgd as well.])
      sage_spkg_install_libgd=yes
    else
      AC_MSG_RESULT([No.])
      PKG_CHECK_MODULES([LIBGD], [gdlib >= 2.1], [], [sage_spkg_install_libgd=yes])
    fi
])


