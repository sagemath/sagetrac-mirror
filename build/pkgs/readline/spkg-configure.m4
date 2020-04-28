SAGE_SPKG_CONFIGURE([readline], [
  SAGE_SPKG_DEPCHECK([ncurses], [
  dnl First try checking for readline with pkg-config
    PKG_CHECK_MODULES([READLINE], [readline >= 6.0], [],
      [AC_CHECK_HEADERS([readline/readline.h],
  dnl rl_bind_keyseq is not present in macos's readline
  dnl and is not present in readline version 4 (the one in OpenBSD)
        [AC_SEARCH_LIBS([rl_bind_keyseq], [readline], [break],
                        [sage_spkg_install_readline=yes])],
        [sage_spkg_install_readline=yes])],
      [sage_spkg_install_readline=yes])
  ])
])
