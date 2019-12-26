AC_DEFUN([SAGE_ABSOLUTE_LIB], [
dnl
dnl ARG1 - the name of the dynamic library to find absolute name of
dnl ARG2 - a function in the library to test for
dnl
  AC_MSG_CHECKING([[absolute path to lib]$1[...]])
  AC_LANG_PUSH(C)
  ABS_LIB_SAVED_LIBS=$LIBS
  LIBS="[-l]$1[ -ldl]"
  AC_RUN_IFELSE([
            AC_LANG_PROGRAM(
            [[
              #include <stdio.h>
              #define __USE_GNU
              #include <dlfcn.h>
              #ifdef __APPLE__
                 #define EXT "dylib"
              #else
                 #define EXT "so"
              #endif
            ]],
            dnl run as $ cc -o conftest conftest.c -lARG1 -ldl && ./conftest
            dnl output should be like  /usr/lib/x86_64-linux-gnu/libARG1.so
            [[ /* int main()... */
             Dl_info  info;
             int res;
             void *ha;
             void *z = dlopen("lib$1."EXT,  RTLD_LAZY); /* load libARG1 */
             if (z) {
               ha = dlsym(z, "$2");       /* get address of inflateEnd in libz */
               res = dladdr(ha, &info);           /* get info for the function */
               printf("%s\n", info.dli_fname);
               return 0; /* dladdr return value on success is platform-dependent */
             }
             printf("dlopen() call failed!\n");
             return 1;
            ]])], [
             [computed_]$1[libdir]=`./conftest$EXEEXT`
             AC_MSG_RESULT([ got it: "$[computed_]$1[libdir]"])
            ], [
             AC_MSG_RESULT([ failure.])
             [computed_]$1[libdir]="/usr/lib"dnl a pretty random choice
  ])
  LIBS=$ABS_LIB_SAVED_LIBS
  AC_LANG_POP(C)
])
