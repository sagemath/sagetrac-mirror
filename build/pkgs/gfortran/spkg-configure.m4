AC_DEFUN([SAGE_NEEDS_FORTRAN], [
        AC_MSG_NOTICE([$1])
        AC_MSG_ERROR([Sage needs working Fortran compiler.])
])

dnl This macro saves current FCFLAGS for later use.
AC_DEFUN([SAGE_SAVE_FCFLAGS], [
    sage_saved_fcflags="$FCFLAGS"
])

dnl This macro restores saved FCFLAGS.
AC_DEFUN([SAGE_RESTORE_FCFLAGS], [
    FCFLAGS="$sage_saved_fcflags"
])


SAGE_SPKG_CONFIGURE([gfortran], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_GCC])
    AC_REQUIRE([AC_PROG_FC])
    # Check that the Fortran compiler accepts free-format source code (as
    # opposed to the older fixed-format style from Fortran 77).
    # This helps verify the compiler works too, so if some idiot sets FC to
    # /usr/bin/ls, we will at least know it's not a working Fortran
    # compiler.
    AC_REQUIRE([SAGE_SAVE_FCFLAGS])
    AC_FC_FREEFORM([SAGE_HAVE_FC_FREEFORM=yes], [
        SAGE_HAVE_FC_FREEFORM=no
        AC_MSG_NOTICE([Your Fortran compiler does not accept free-format source code])
        AC_MSG_NOTICE([which means the compiler is either seriously broken, or])
        SAGE_NEEDS_FORTRAN([is too old to build Sage.])
	])

    # AC_FC_FREEFORM may have added flags.
    # However, it is up to the individual package how they invoke the
    # Fortran compiler.
    # We only check here, whether the compiler is suitable.
    AC_REQUIRE([SAGE_RESTORE_FCFLAGS])

    # AX_COMPILER_VENDOR does not work for Fortran. So we just match the name of the executable
    AS_CASE(["$FC"],
            [*gfortran*], [
                AC_MSG_CHECKING([the version of $FC])
                GFORTRAN_VERSION="`$FC -dumpversion`"
                AC_MSG_RESULT([$GFORTRAN_VERSION])
                # Add the .0 because Debian/Ubuntu gives version numbers like
                # 4.6 instead of 4.6.4 (Trac #18885)
                AS_CASE(["$GFORTRAN_VERSION.0"],
                    [[[0-3]].*|4.[[0-7]].*], [
                        SAGE_NEEDS_FORTRAN([$FC is version $GFORTRAN_VERSION, which is quite old])
                    ],
                    [1[[2-9]].*], [
                        # system-provided gfortran is newer than 11.x.
                        # See https://trac.sagemath.org/ticket/29456, https://trac.sagemath.org/ticket/31838
                        SAGE_NEEDS_FORTRAN([$FC is version $GFORTRAN_VERSION, which is too recent for this version of Sage])
                    ])
    ])
])
