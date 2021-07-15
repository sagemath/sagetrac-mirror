SAGE_SPKG_CONFIGURE([distlib], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_VIRTUALENV])
    sage_spkg_install_distlib=yes
  ], [dnl REQUIRED-CHECK
    dnl only needed as a dependency of virtualenv.
    AS_VAR_SET([SPKG_REQUIRE], [$sage_spkg_install_virtualenv])
  ])
