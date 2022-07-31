SAGE_SPKG_CONFIGURE([pari_nftables], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI])
    AC_MSG_CHECKING([installing pari? ])
    if test x$sage_spkg_install_pari = xyes; then
        AC_MSG_RESULT([yes; install pari_nftables as well])
        sage_spkg_install_pari_nftables=yes
    else
        AC_MSG_RESULT([no])
        AC_MSG_CHECKING([whether nftables is installed in pari/GP datadir ])
        gp_nft_check=`echo 'v=readvec(concat([default(datadir), "/nftables/T71.gp"])); -v[[1]][[1]]' | $GP -qf 2>> config.log`
        if test x$gp_nft_check = x184607; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; nftables are not there. ])
            sage_spkg_install_pari_nftables=yes
            AC_MSG_NOTICE([cf. http://pari.math.u-bordeaux1.fr/pub/pari/packages/nftables/README.txt])
            AC_MSG_NOTICE([For Sage, installing nftables should be done into pari/GP datadir,])
        fi
    fi
])
