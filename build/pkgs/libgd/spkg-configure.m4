SAGE_SPKG_CONFIGURE([libgd], [
    m4_pushdef([SAGE_LIBGD_MINVER],[2.1])
    SAGE_SPKG_DEPCHECK([libpng freetype], [
        PKG_CHECK_MODULES([LIBGD], [gdlib >= $SAGE_LIBGD_MINVER], [], [sage_spkg_install_libgd=yes])
    ])
    m4_popdef([SAGE_LIBGD_MINVER])
])


