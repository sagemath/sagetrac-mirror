SAGE_SPKG_CONFIGURE([swig], [
   AC_PATH_PROG([SWIG_CREATE_CMAKELISTS], [swig])
   AS_IF([test -z "$ac_cv_path_SWIG_CREATE_CMAKELISTS"], [sage_spkg_install_swig=yes])
])
