SAGE_SPKG_CONFIGURE([gsl], [
    m4_pushdef([SAGE_GSL_MINVER],["2.4"])
    SAGE_SPKG_DEPCHECK([atlas openblas], [
      PKG_CHECK_MODULES([GSL], [gsl >= $SAGE_GSL_MINVER], [
        PKG_CHECK_VAR([GSLPCDIR], [gsl], [pcfiledir], [
          AC_CONFIG_FILES([$SAGE_LOCAL/lib/pkgconfig/gsl.pc:$GSLPCDIR/gsl.pc])
          AC_CONFIG_COMMANDS([GSLPCPROCESS], [
            sed -i'' -e 's/\${GSL_CBLAS_LIB}\ //' "$SAGE_LOCAL"/lib/pkgconfig/gsl.pc
            sed -i'' -e 's/GSL_CBLAS_LIB.*/Require: cblas/' "$SAGE_LOCAL"/lib/pkgconfig/gsl.pc
          ])
        ], [
        AC_MSG_WARN([Unable to locate the directory of gsl.pc. This should not happen!])
       sage_spkg_install_gsl=yes
       ])
      ], [sage_spkg_install_gsl=yes])
    ])
    m4_popdef([SAGE_GSL_MINVER])
])
