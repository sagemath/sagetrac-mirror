SAGE_SPKG_CONFIGURE([polymake], [
    sage_spkg_install_polymake=yes
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_GCC])
    if test $IS_REALLY_GCC = yes ; then
       AC_MSG_NOTICE([checking if gcc version is at least 5.1 for polymake])
       AS_CASE(["$GXX_VERSION.0"],
           [[[0-4]].*|5.0.*], [
               SAGE_MUST_INSTALL_GCC([you have $CXX version $GXX_VERSION, which is too old to build polymake])
                ])
    fi])
