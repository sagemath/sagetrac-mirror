SAGE_SPKG_CONFIGURE([fflas_ffpack], [
  dnl https://github.com/linbox-team/fflas-ffpack/blob/master/macros/instr_set.m4
  dnl discovers these flags from the processor but fails to check whether
  dnl compiler (and assembler) actually support these instruction sets.
  m4_foreach([ISFLAG], [avx512f, avx512vl, avx512dq, fma, fma4], [
    AX_CHECK_COMPILE_FLAG([-m]ISFLAG, [], [AS_VAR_APPEND]([SAGE_CONFIGURE_FFLAS_FFPACK], [" --disable-]ISFLAG[ "]))
  ])
  AC_SUBST([SAGE_CONFIGURE_FFLAS_FFPACK])
])
