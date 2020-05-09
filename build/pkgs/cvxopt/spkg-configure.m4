SAGE_SPKG_CONFIGURE(
    [cvxopt], [dnl direct testing for import module
      AC_MSG_CHECKING(Import module cvxopt...)
      python3 -c "import cvxopt;" > /dev/null 2>&1
      if test $? -ne 0; then
        AC_MSG_RESULT(KO)
        sage_spkg_install_cvxopt=yes
      else
        AC_MSG_RESULT(OK)
      fi
])
