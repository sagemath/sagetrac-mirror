SAGE_SPKG_CONFIGURE([libpng], [
    m4_pushdef([SAGE_LIBPNG_MINVER],[1.2])
    SAGE_SPKG_DEPCHECK([zlib], [
      dnl First try checking for libpng with pkg-config
      PKG_CHECK_MODULES([LIBPNG], [libpng >= $SAGE_LIBPNG_MINVER], [], [
        dnl Fallback to manually grubbing around for headers and libs
        AC_CHECK_HEADERS([png.h], [
           AC_SEARCH_LIBS([png_get_io_ptr], [png], [
             dnl write out libpng.pc
           ], [sage_spkg_install_libpng=yes])
         ], [sage_spkg_install_libpng=yes])
      ])
    ])
    m4_popdef([SAGE_LIBPNG_MINVER])
])
