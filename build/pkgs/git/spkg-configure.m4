# Check whether git works by executing "git --version"
SAGE_SPKG_CONFIGURE([git], [
	AC_CACHE_CHECK([for git], [ac_cv_path_GIT], [
		AC_PATH_PROGS_FEATURE_CHECK([GIT], [git],
		[${ac_path_GIT} --version >/dev/null 2>/dev/null && ac_cv_path_GIT=${ac_path_GIT}],
		[sage_spkg_install_git=yes; ac_cv_path_GIT=no])])
])
