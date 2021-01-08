SAGE_SPKG_CONFIGURE([python3_tkinter], [
   SAGE_SPKG_DEPCHECK([python3], [
       dnl The _tkinter module is actually built and installed by python3.
       AS_IF([test -n "$PYTHON_FOR_VENV"], [
           AC_MSG_CHECKING([whether $PYTHON_FOR_VENV can import the _tkinter module])
           AS_IF(["$PYTHON_FOR_VENV" -m _tkinter >& AS_MESSAGE_LOG_FD 2>&1], [
               AC_MSG_RESULT([yes])
           ], [
               AC_MSG_RESULT([no])
               sage_spkg_install_python3_tkinter=yes
           ])
       ])
   ])
])
