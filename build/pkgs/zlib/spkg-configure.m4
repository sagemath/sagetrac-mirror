AC_DEFUN([SAGE_ZLIB_GEN_PC],[dnl generate zlib.pc
  AX_ABSOLUTE_HEADER([zlib.h])
  AS_IF([test x$gl_cv_absolute_zlib_h = x], [
     AC_MSG_ERROR(m4_normalize([
       failed to find absolute path to zlib.h despite it being reported found
     ]))
  ])
  AC_SUBST(SAGE_ZLIB_INCLUDE, [`AS_DIRNAME($gl_cv_absolute_zlib_h)`])
  dnl now the same for the full path to libz...
  dnl ...
])

SAGE_SPKG_CONFIGURE([zlib], [
    PKG_CHECK_MODULES([ZLIBANDLIBPNG], [zlib libpng >= 1.2], [
      AC_MSG_NOTICE([Use zlib and libpng from the system, pc files supplied.])
      ], [
      PKG_CHECK_MODULES([ZLIB], [zlib >= 1.2.9], [
        dnl inflateValidate is needed for Sage's libpng, newer than 1.2; this ensures
        dnl we have the minimum required for building zlib version
        AC_MSG_NOTICE([Use zlib the system, pc file supplied.])
        ], [dnl no zlib.pc file. try to see if this zlib is still good
        AC_CHECK_HEADER([zlib.h], [
           AC_CHECK_LIB([z], [inflateEnd], [
              PKG_CHECK_MODULES([LIBPNG], [libpng >= 1.2], [
                dnl hope we found the zlib of system's libpng
                SAGE_ZLIB_GEN_PC
                ], [
                dnl inflateValidate is needed for Sage's libpng (cf. above)
                AC_CHECK_LIB([z], [inflateValidate], [
                   dnl zlib is good enough to build Sage's libpng
                   SAGE_ZLIB_GEN_PC
                   ], [sage_spkg_install_zlib=yes])
                ])
           ], [sage_spkg_install_zlib=yes])
        ], [sage_spkg_install_zlib=yes])
      ])
    ])
])
