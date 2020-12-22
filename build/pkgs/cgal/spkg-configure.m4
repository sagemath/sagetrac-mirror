SAGE_SPKG_CONFIGURE([cgal], [
   AC_PATH_PROG([CGAL_CREATE_CMAKELISTS], [cgal_create_CMakeLists])
   AS_IF([test -z "$ac_cv_path_CGAL_CREATE_CMAKELISTS"], [sage_spkg_install_cgal=yes])
])
