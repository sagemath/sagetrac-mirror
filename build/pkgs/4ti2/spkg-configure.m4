SAGE_SPKG_CONFIGURE([4ti2], [
    SAGE_SPKG_DEPCHECK([gmp mpir glpk zlib], [
        dnl Debian installs these programs with an executable prefix "4ti2-"
        dnl but polymake and our own code do not handle this yet (Singular does).
        m4_foreach([prog], [hilbert,markov,graver,zsolve,qsolve,rays,ppi,circuits,groebner], [
            AC_CHECK_PROGS([FOURTITWO_]m4_toupper(prog), prog)
            AS_VAR_IF([FOURTITWO_]m4_toupper(prog), [""], [sage_spkg_install_4ti2=yes])
        ])
        dnl Adapted from https://github.com/latte-int/latte/blob/master/m4/4ti2-check.m4
        AC_MSG_CHECKING(for library 4ti2gmp)
        BACKUP_CXXFLAGS=${CXXFLAGS}
        BACKUP_LIBS=${LIBS}
        FORTYTWO_CXXFLAGS="-D__STDC_LIMIT_MACROS -D_4ti2_GMP_"
        FORTYTWO_LIBS="-l4ti2gmp -lzsolve"
        CXXFLAGS="${BACKUP_CXXFLAGS} ${FORTYTWO_CXXFLAGS} ${GMP_CFLAGS}"
        LIBS="${BACKUP_LIBS} ${FORTYTWO_LIBS} ${GMP_LIBS}"
        AC_TRY_LINK([
#include "4ti2/4ti2.h"
],
[ _4ti2_rays_create_state(_4ti2_PREC_INT_ARB);
], [
        AC_MSG_RESULT([yes])
], [
        AC_MSG_RESULT([no])
        sage_spkg_install_4ti2=yes
])
        CXXFLAGS=${BACKUP_CXXFLAGS}
        LIBS=${BACKUP_LIBS}
    ])
])
