SAGE_SPKG_CONFIGURE([lrslib], [
    # The lrs and lrsnash binaries are the only interfaces to lrslib that
    # sagelib uses. As a result, we don't need to call SAGE_SPKG_DEPCHECK
    # here because there's no possibility for a library conflict.
    AC_CHECK_PROGS([LRSNASH], [lrsnash])
    AS_IF([test -z "$LRSNASH"], [
        sage_spkg_install_lrslib=yes
    ], [
        AC_MSG_CHECKING([whether $LRSNASH can handle the new input format])
        cat > conftest.lrsnash <<EOF
1 1

0

0
EOF
        AS_IF([$LRSNASH conftest.lrsnash >& AS_MESSAGE_LOG_FD 2>&1], [
            AC_MSG_RESULT([yes])
        ], [
            AC_MSG_RESULT([no])
            sage_spkg_install_lrslib=yes
        ])
        rm -f conftest.lrsnash
    ])
])
