SAGE_SPKG_CONFIGURE([freetype], [
    m4_pushdef([SAGE_FREETYPE_MINVER],[2.4])
    SAGE_SPKG_DEPCHECK([libpng], [
      PKG_CHECK_MODULES([FREETYPE], [freetype2 >= $SAGE_FREETYPE_MINVER], [], [sage_spkg_install_freetype=yes])
    ])
    AS_IF([test x$sage_spkg_install_freetype = xyes],
      [AC_SUBST(SAGE_FREETYPE_PREFIX, ['$SAGE_LOCAL'])],
      [AC_SUBST(SAGE_FREETYPE_PREFIX, [''])])
    m4_popdef([SAGE_FREETYPE_MINVER])
])


