SAGE_SPKG_CONFIGURE([libpng], [
    m4_pushdef([SAGE_LIBPNG_MINVER],[1.2])
    SAGE_SPKG_DEPCHECK([zlib], [
      PKG_CHECK_MODULES([LIBPNG], [libpng >= $SAGE_LIBPNG_MINVER], [
      ], [sage_spkg_install_libpng=yes])
    ])
    m4_popdef([SAGE_LIBPNG_MINVER])
])
