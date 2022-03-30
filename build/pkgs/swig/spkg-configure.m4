SAGE_SPKG_CONFIGURE([swig], [
   AC_PATH_PROG([SWIG], [swig])
   AS_IF([test -z "$ac_cv_path_SWIG"], [sage_spkg_install_swig=yes])
])
