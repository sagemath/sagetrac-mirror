SAGE_SPKG_CONFIGURE([sage_system_blas_facade], [
  AC_SUBST([SAGE_SYSTEM_BLAS_FACADE_PC_FILES])
  AC_CONFIG_FILES([build/pkgs/sage_system_blas_facade/spkg-install], [chmod +x build/pkgs/sage_system_blas_facade/spkg-install])
])
