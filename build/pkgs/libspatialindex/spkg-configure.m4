SAGE_SPKG_CONFIGURE([libspatialindex], [
    AC_LANG_PUSH(C++)
    AC_CHECK_HEADER([spatialindex/SpatialIndex.h], [], [sage_spkg_install_libspatialindex=yes])
    AC_MSG_CHECKING([whether we can link a program using libspatialindex])
    LIBSPATIALINDEX_SAVED_LIBS=$LIBS
    LIBS="$LIBS -llibspatialindex"
    AC_LINK_IFELSE([
        AC_LANG_PROGRAM([[#include <spatialindex/SpatialIndex.h>]],
                            [[SpatialIndex::Point p1 = SpatialIndex::Point(1.0,2.0);]]
            )], [AC_MSG_RESULT([yes])], 
	    [AC_MSG_RESULT([no]); sage_spkg_install_libspatialindex=yes])
            LIBS=$LIBSPATIALINDEX_SAVED_LIBS
    AC_LANG_POP(C++)
])