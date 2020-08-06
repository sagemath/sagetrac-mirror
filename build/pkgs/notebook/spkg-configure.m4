dnl We check for all Jupyter components in one spkg-configure file
dnl and only accept system packages if ALL system packages are suitable.

SAGE_SPKG_CONFIGURE([notebook], [
   dnl SAGE_SPKG_DEPCHECK([], [
      dnl 4.4.0 is the version in Sage 9.1, https://trac.sagemath.org/ticket/24168
      m4_pushdef([JUPYTER_CORE_MIN_VERSION],   [4.4.0])
      m4_pushdef([JUPYTER_CORE_SUBCOMMAND],    [])
      dnl 5.7.6 is the version in Sage 9.1, updated because
      dnl 5.7.4 was broken, according to the https://trac.sagemath.org/ticket/27463
      m4_pushdef([NOTEBOOK_MIN_VERSION],       [5.7.6])
      m4_pushdef([NOTEBOOK_SUBCOMMAND],        [notebook])
      dnl 5.4.0 is the version in Sage 9.1, https://trac.sagemath.org/ticket/26969
      m4_pushdef([NBCONVERT_MIN_VERSION],      [5.4.0])
      m4_pushdef([NBCONVERT_SUBCOMMAND],       [nbconvert])
      dnl 5.2.4 is the version in Sage 9.1, https://trac.sagemath.org/ticket/26969
      m4_pushdef([JUPYTER_CLIENT_MIN_VERSION], [5.2.4])
      m4_pushdef([JUPYTER_CLIENT_SUBCOMMAND],  [kernelspec])

      AC_PATH_PROG([JUPYTER], [jupyter])



   dnl ])
])
