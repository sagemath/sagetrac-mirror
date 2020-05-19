SAGE_SPKG_CONFIGURE([palp], [
   # Sage prefers to execute this program using the name "class-4d.x",
   # because some distributions install a version of it that is
   # optimized for small dimensions under that name. If "class-4d.x"
   # isn't in the user's PATH, then the default upstream name of
   # "class.x" can be used, but only if we wind up using the system
   # copy of palp. We default the value of SAGE_PALP_CLASS_EXE to what
   # is appropriate if the palp SPKG is used, and override it if we
   # find a suitable system copy of palp.
   SAGE_PALP_CLASS_EXE="class-4d.x"

   AC_PATH_PROG([PALP], [poly.x])
   AS_IF([test -z "$ac_cv_path_PALP"],
         [sage_spkg_install_palp=yes],
	 [
           AC_PATH_PROG([PALP_CLASS4D_PATH], [class-4d.x])
           AS_IF([test -n "${PALP_CLASS4D_PATH}"],[
                   SAGE_PALP_CLASS_EXE="${PALP_CLASS4D_PATH}"
           ],[
                   AC_PATH_PROG([PALP_CLASS_PATH], [class.x])
                   AS_IF([test -z "${PALP_CLASS_PATH}"],
                         [sage_spkg_install_palp=yes],
                         [SAGE_PALP_CLASS_EXE="${PALP_CLASS_PATH}"])
           ])
         ])
   AC_SUBST(SAGE_PALP_CLASS_EXE)
])
