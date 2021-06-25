SAGE_SPKG_CONFIGURE([libspatialindex], [
    dnl We check all executables that are tested by sage.features.libspatialindex
    AC_CHECK_PROGS([MVRTREE], [mvrtree])
    AS_IF([test x$MVRTREE = x], [sage_spkg_install_libspatialindex=yes])
    AC_CHECK_PROGS([TPRTREE], [tprtree])
    AS_IF([test x$TPRTREE = x], [sage_spkg_install_libspatialindex=yes])
])
