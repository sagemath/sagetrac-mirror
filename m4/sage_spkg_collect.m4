# SYNOPSIS
#
#   SAGE_SPKG_COLLECT
#
# DESCRIPTION
#
#   This macro gathers up information about SPKGs defined in the build/pkgs
#   directory of the Sage source tree, and generates variables to be
#   substitued into the build/make/Makefile.in template which list all the
#   SPKGs, their versions, their dependencies, and categorizes them based
#   on how they should be installed.
#
#   In particular, this generates the Makefile variables:
#
#      - SAGE_BUILT_PACKAGES - lists the names of SPKGs that should be built
#        and installed from source.
#
#      - SAGE_DUMMY_PACKAGES - lists the names of packages that are not built
#        as part of Sage--either they are not required for the current
#        platform, or the dependency on them is satisfied by an existing
#        system package.
#
#      - SAGE_STANDARD_PACKAGES - lists the names of all packages that have
#        the "standard" type.  All "standard" packages are installed by
#        default (if they are listed in SAGE_DUMMY_PACKAGES "installed" in
#        this case is a no-op).
#
#      - SAGE_OPTIONAL_PACKAGES - lists the names of packages with the
#        "optional" type that should be installed.
#
#      - SAGE_SDIST_PACKAGES - lists the names of all packages that should be
#        included in the source distribution.
#
#      - SAGE_PACKAGE_VERSIONS - this template variable defines multiple
#        Makefile variables in the format "vers_<packagename>" the value
#        of which is the current version of the SPKG <packagename>.
#
#      - SAGE_PACKAGE_DEPENDENCIES - this template variable defines multiple
#        Makefile variables in the format "deps_<packagename>" the value
#        of which is the names of the dependencies of <packagename> as read
#        from the build/<packagename>/dependencies file.
#
#      - SAGE_NORMAL_PACKAGES - lists the names of packages that are installed
#        by the "normal" method (using the sage-spkg program to download and
#        extract the source tarball, and run the relevant scripts from
#        build/<packagename>/spkg-*.
#
#      - SAGE_PIP_PACKAGES - lists the names of packages with the "pip" type
#        which are installed by directly invoking the pip command.
#
#      - SAGE_SCRIPT_PACKAGES - lists the names of packages with the "script"
#        type which are installed by running a custom script, which may
#        download additional source files.
#
AC_DEFUN_ONCE([SAGE_SPKG_COLLECT], [
# To deal with ABI incompatibilities when gcc is upgraded, every package
# (except gcc) should depend on gcc if gcc is already installed.
# See https://trac.sagemath.org/ticket/24703
if test x$SAGE_INSTALL_GCC = xexists; then
    SAGE_GCC_DEP='$(SAGE_LOCAL)/bin/gcc'
else
    SAGE_GCC_DEP=''
fi
AC_SUBST([SAGE_GCC_DEP])

AC_MSG_CHECKING([package versions])
AC_MSG_RESULT([])

# Usage: newest_version $pkg
# Print version number of latest package $pkg
newest_version() {
    PKG=$[1]
    if test -f "$SAGE_ROOT/build/pkgs/$PKG/package-version.txt" ; then
        cat "$SAGE_ROOT/build/pkgs/$PKG/package-version.txt"
    else
        echo "$PKG"
    fi
}

# Lists of packages that are actually built/installed and dummy packages
SAGE_BUILT_PACKAGES='\
'
SAGE_DUMMY_PACKAGES='\
'
# List of all standard packages
SAGE_STANDARD_PACKAGES='\
'
# List of all currently installed optional packages
SAGE_OPTIONAL_INSTALLED_PACKAGES='\
'
# List of all packages that should be downloaded
SAGE_SDIST_PACKAGES='\
'
# Generate package version and dependency lists
SAGE_PACKAGE_VERSIONS=""
SAGE_PACKAGE_DEPENDENCIES=""
# Lists of packages categorized according to their build rules
SAGE_NORMAL_PACKAGES='\
'
SAGE_PIP_PACKAGES='\
'
SAGE_SCRIPT_PACKAGES='\
'

# for each package in pkgs/, add them to the SAGE_PACKAGE_VERSIONS and
# SAGE_PACKAGE_DEPENDENCIES lists, and to one or more of the above variables
# depending on the package type and other criteria (such as whether or not it
# needs to be installed)

for DIR in $SAGE_ROOT/build/pkgs/*; do
    test -d "$DIR" || continue

    PKG_TYPE_FILE="$DIR/type"
    if test -f "$PKG_TYPE_FILE"; then
        PKG_TYPE=`cat $PKG_TYPE_FILE`
    else
        AC_MSG_ERROR(["$PKG_TYPE_FILE" is missing.])
    fi

    PKG_NAME=$(basename $DIR)
    PKG_VERSION=$(newest_version $PKG_NAME)

    in_sdist=no

    # Check consistency of 'DIR/type' file
    case "$PKG_TYPE" in
    base) ;;
    standard) 
        SAGE_STANDARD_PACKAGES+="    $PKG_NAME \\"$'\n'
        in_sdist=yes
        ;;
    optional)
        if test -f $SAGE_SPKG_INST/$PKG_NAME-*; then
            SAGE_OPTIONAL_INSTALLED_PACKAGES+="    $PKG_NAME \\"$'\n'
        fi;
        ;;
    experimental) ;;
    script) ;;
    pip) ;;
    *)
        AC_MSG_ERROR([The content of "$PKG_TYPE_FILE" must be 'base', 'standard', 'optional', 'experimental', 'script', or 'pip'])
        ;;
    esac

    SAGE_PACKAGE_VERSIONS+="vers_$PKG_NAME = $PKG_VERSION"$'\n'

    # If $need_to_install_{PKG_NAME} is set to no, then set inst_<pkgname> to
    # some dummy file to skip the installation. Note that an explicit
    # "./sage -i PKG_NAME" will still install the package.
    if test "$PKG_NAME" != "$PKG_VERSION"; then
        need_to_install="need_to_install_${PKG_NAME}"

        if test "${!need_to_install}" != no ; then
            SAGE_BUILT_PACKAGES+="    $PKG_NAME \\"$'\n'
            AC_MSG_RESULT([    $PKG_NAME-$PKG_VERSION])
        else
            SAGE_DUMMY_PACKAGES+="    $PKG_NAME \\"$'\n'
            AC_MSG_RESULT([    $PKG_NAME-$PKG_VERSION not installed (configure check)])
        fi
    fi

    # Packages that should be included in the source distribution
    # This includes all standard packages and two special cases
    case "$PKG_NAME" in
    mpir|python2)
        in_sdist=yes
        ;;
    esac

    if test "$in_sdist" = yes; then
        SAGE_SDIST_PACKAGES+="    $PKG_NAME \\"$'\n'
    fi

    # Determine package dependencies
    DEP_FILE="$SAGE_ROOT/build/pkgs/$PKG_NAME/dependencies"
    if test -f "$DEP_FILE"; then
        # - the # symbol is treated as comment which is removed
        DEPS=`sed 's/^ *//; s/ *#.*//; q' $DEP_FILE`
    else
        case "$PKG_TYPE" in
        optional)
            DEPS=' | $(STANDARD_PACKAGES)' # default for optional packages
            ;;
        script)
            DEPS=' | $(STANDARD_PACKAGES)' # default for script-only packages
            ;;
        pip)
            DEPS=' | pip'
            ;;
        *)
            DEPS=""
            ;;
        esac
    fi

    SAGE_PACKAGE_DEPENDENCIES+="deps_$PKG_NAME = $DEPS"$'\n'

    # Determine package build rules
    case "$PKG_TYPE" in
    pip)
        SAGE_PIP_PACKAGES+="    $PKG_NAME \\"$'\n'
        ;;
    script)
        SAGE_SCRIPT_PACKAGES+="    $PKG_NAME \\"$'\n'
        ;;
    *)
        SAGE_NORMAL_PACKAGES+="    $PKG_NAME \\"$'\n'
        ;;
    esac
done

AC_SUBST([SAGE_PACKAGE_VERSIONS])
AC_SUBST([SAGE_PACKAGE_DEPENDENCIES])
AC_SUBST([SAGE_NORMAL_PACKAGES])
AC_SUBST([SAGE_PIP_PACKAGES])
AC_SUBST([SAGE_SCRIPT_PACKAGES])
AC_SUBST([SAGE_BUILT_PACKAGES])
AC_SUBST([SAGE_DUMMY_PACKAGES])
AC_SUBST([SAGE_STANDARD_PACKAGES])
AC_SUBST([SAGE_OPTIONAL_INSTALLED_PACKAGES])
AC_SUBST([SAGE_SDIST_PACKAGES])
])
