SAGE_SPKG_CONFIGURE([zstd], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_GCC])
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_GFORTRAN])
    sage_spkg_install_zstd=yes
  ], [dnl REQUIRED-CHECK
    dnl only needed as a dependency of gcc/gfortran.
    AS_VAR_SET([SPKG_REQUIRE], [no])
    AS_VAR_IF([sage_spkg_install_gfortran], [yes], [dnl
        AS_VAR_SET([SPKG_REQUIRE], [yes])])
    AS_VAR_IF([sage_spkg_install_gcc], [yes], [dnl
        AS_VAR_SET([SPKG_REQUIRE], [yes])])
  ])
