SAGE_SPKG_CONFIGURE(
    [patch], [
	AC_CACHE_CHECK([for a POSIX patch], [ac_cv_path_PATCH], [
        AC_PATH_PROGS_FEATURE_CHECK([PATCH], [patch], [
            patch_not_posix=`$ac_path_PATCH --posix --version 2>&1 | grep unrecognized`
            AS_IF([test -n "$patch_not_posix"],
                  [sage_spkg_install_patch=yes],
                  [ac_cv_path_PATCH="$ac_path_PATCH"])
        ])
    ])
    AS_IF([test -z "$ac_cv_path_PATCH"], [sage_spkg_install_patch=yes])
])
