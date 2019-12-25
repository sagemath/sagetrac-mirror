AC_DEFUN([SAGE_ZLIB_GEN_PC],[dnl generate zlib.pc
  AX_ABSOLUTE_HEADER([zlib.h])
  AS_IF([test x$gl_cv_absolute_zlib_h = x], [
     AC_MSG_ERROR(m4_normalize([
       failed to find absolute path to zlib.h despite it being reported found
     ]))
  ])
  AC_SUBST(SAGE_ZLIB_INCLUDEDIR, [`AS_DIRNAME($gl_cv_absolute_zlib_h)`])
  AC_MSG_CHECKING([getting absolute path to libz...])
  AC_LANG_PUSH(C)
  ZLIB_SAVED_LIBS=$LIBS
  LIBS="-lz -ldl"
  AC_RUN_IFELSE([
            AC_LANG_PROGRAM(
            [[
              #include <stdio.h>
              #define __USE_GNU
              #include <dlfcn.h>
              #ifdef __APPLE__
                 #define EXT "dylib"
              #else
                 #define EXT "so"
              #endif
            ]],
            dnl run as $ cc -o conftest conftest.c -lz -ldl && ./conftest
            dnl output should be like  /usr/lib/x86_64-linux-gnu/libz.so
            [[ /* int main()... */
             Dl_info  info;
             int res;
             void *ha;
             void *z = dlopen("libz."EXT,  RTLD_LAZY); /* load zlib */
             if (z) {
               ha = dlsym(z, "inflateEnd");       /* get address of inflateEnd in libz */
               res = dladdr(ha, &info);           /* get info for the function */
               printf("%s\n", info.dli_fname);
               return 0; /* dladdr return value on success is platform-dependent */
             }
             printf("dlopen() call failed!\n");
             return 1;
            ]])], [
             computed_zlibdir=`./conftest$EXEEXT`
             AC_MSG_RESULT([ got it: "$computed_zlibdir"])
            ], [
             AC_MSG_RESULT([ failure.])
             computed_zlibdir="/usr/lib"dnl a pretty random choice
  ])
  LIBS=$ZLIB_SAVED_LIBS
  AC_LANG_POP(C)
  AC_SUBST(SAGE_ZLIB_LIBDIR, [`AS_DIRNAME($computed_zlibdir)`])
  cp build/pkgs/zlib/zlib.pc.in $SAGE_LOCAL/lib/pkgconfig/
])

SAGE_SPKG_CONFIGURE([zlib], [
    AC_CONFIG_FILES([$SAGE_LOCAL/lib/pkgconfig/zlib.pc])
    mkdir -p $SAGE_LOCAL/lib/pkgconfig/
    PKG_CHECK_MODULES([ZLIBANDLIBPNG], [zlib libpng >= 1.2], [
      AC_MSG_NOTICE([Use zlib and libpng from the system, pc files supplied.])
      PKG_CHECK_VAR([ZLIBPCDIR], [zlib], [pcfiledir], [
         dnl install zlib.pc as zlib.pc.in
         cp "$ZLIBPCDIR"/zlib.pc $SAGE_LOCAL/lib/pkgconfig/zlib.pc.in
         ], [
         AC_MSG_ERROR([Unable to locate the directory of zlib.pc. This should not happen!])
      ])
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
