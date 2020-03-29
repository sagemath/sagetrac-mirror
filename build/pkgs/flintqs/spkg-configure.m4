SAGE_SPKG_CONFIGURE([flintqs], [
  SAGE_SPKG_DEPCHECK([gmp], [
    # FlintQS only builds against GMP, not MPIR. If we're using the
    # system GMP, then we check for the QuadraticSieve program, which
    # is the only interface to FlintQS that sagelib uses.
    AC_CHECK_PROG(HAVE_QUADRATICSIEVE, QuadraticSieve, yes, no)
    AS_IF([test "x$HAVE_QUADRATICSIEVE" = "xno"],
          [sage_spkg_install_flintqs=yes])
  ],
  [ # If we're using sage's GMP, we have to use its FlintQS, too.
    sage_spkg_install_flintqs=yes
  ])
])
