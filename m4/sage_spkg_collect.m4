# SYNOPSIS
#
#   SAGE_SPKG_COLLECT
#
# DESCRIPTION
#
#   This macro gathers up information about SPKGs defined in the build/pkgs
#   directory of the Sage source tree, and generates variables to be
#   substituted into the build/make/Makefile.in template which list all the
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
#      - SAGE_SDIST_PACKAGES - lists the names of all packages whose sources
#        need to be downloaded to be included in the source distribution.
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

dnl ==========================================================================
dnl define PKG_CHECK_VAR for old pkg-config < 0.28; see Trac #29001
m4_ifndef([PKG_CHECK_VAR], [
AC_DEFUN([PKG_CHECK_VAR],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
AC_ARG_VAR([$1], [value of $3 for $2, overriding pkg-config])dnl

_PKG_CONFIG([$1], [variable="][$3]["], [$2])
AS_VAR_COPY([$1], [pkg_cv_][$1])

AS_VAR_IF([$1], [""], [$5], [$4])dnl
])dnl PKG_CHECK_VAR
])
dnl ==========================================================================

AC_DEFUN_ONCE([SAGE_SPKG_COLLECT], [
# Configure all spkgs with configure-time checks
m4_include([m4/sage_spkg_configures.m4])

# To deal with ABI incompatibilities when gcc is upgraded, every package
# (except gcc) should depend on gcc if gcc is already installed.
# See https://trac.sagemath.org/ticket/24703
if test x$SAGE_INSTALL_GCC = xexists; then
    SAGE_GCC_DEP='$(SAGE_LOCAL)/bin/gcc'
else
    SAGE_GCC_DEP=''
fi
AC_SUBST([SAGE_GCC_DEP])

AS_BOX([Build status for each package:                                         ]) >& AS_MESSAGE_FD
AS_BOX([Build status for each package:                                         ]) >& AS_MESSAGE_LOG_FD

# Usage: newest_version $pkg
# Print version number of latest package $pkg
newest_version() {
    SPKG=$[1]
    if test -f "$SAGE_ROOT/build/pkgs/$SPKG/package-version.txt" ; then
        cat "$SAGE_ROOT/build/pkgs/$SPKG/package-version.txt"
    else
        echo none
    fi
}

# Packages that are actually built/installed as opposed to packages that are
# not required on this platform or that can be taken from the underlying system
# installation. Note that this contains packages that are not actually going to
# be installed by most users because they are optional/experimental.
SAGE_BUILT_PACKAGES=''

# The complement of SAGE_BUILT_PACKAGES, i.e., packages that are not required
# on this platform or packages where we found a suitable package on the
# underlying system.
SAGE_DUMMY_PACKAGES=''

# Standard packages
SAGE_STANDARD_PACKAGES=''

# List of currently installed and to-be-installed optional packages - filled in SAGE_SPKG_ENABLE
#SAGE_OPTIONAL_INSTALLED_PACKAGES
# List of optional packages to be uninstalled - filled in SAGE_SPKG_ENABLE
#SAGE_OPTIONAL_CLEANED_PACKAGES

# List of all packages that should be downloaded
SAGE_SDIST_PACKAGES=''

# Generate package version and dependency lists
SAGE_PACKAGE_VERSIONS=""
SAGE_PACKAGE_DEPENDENCIES=""
# Lists of packages categorized according to their build rules
SAGE_NORMAL_PACKAGES=''
SAGE_PIP_PACKAGES=''
SAGE_SCRIPT_PACKAGES=''

SAGE_NEED_SYSTEM_PACKAGES=""

# for each package in pkgs/, add them to the SAGE_PACKAGE_VERSIONS and
# SAGE_PACKAGE_DEPENDENCIES lists, and to one or more of the above variables
# depending on the package type and other criteria (such as whether or not it
# needs to be installed)

for DIR in $SAGE_ROOT/build/pkgs/*; do
    test -d "$DIR" || continue

    SPKG_TYPE_FILE="$DIR/type"
    if test -f "$SPKG_TYPE_FILE"; then
        SPKG_TYPE=`cat $SPKG_TYPE_FILE`
    else
        AC_MSG_WARN(["$SPKG_TYPE_FILE" is missing.  Leftovers from another branch?])
        continue
    fi

    SPKG_NAME=$(basename $DIR)
    SPKG_VERSION=$(newest_version $SPKG_NAME)

    in_sdist=no

    uninstall_message=""
    # Check consistency of 'DIR/type' file
    case "$SPKG_TYPE" in
    base)
        message="came preinstalled with the SageMath tarball"
        ;;
    standard)
        SAGE_STANDARD_PACKAGES="${SAGE_STANDARD_PACKAGES} \\$(printf '\n    ')${SPKG_NAME}"
        in_sdist=yes
        message="will be installed as an SPKG"
        ;;
    optional|experimental)
        AS_VAR_IF([SAGE_ENABLE_]${SPKG_NAME}, [yes], [
            message="$SPKG_TYPE, will be installed as an SPKG"
        ], [
            message="$SPKG_TYPE, use \"$srcdir/configure --enable-$SPKG_NAME\" to install"
        ])
        uninstall_message=", use \"$srcdir/configure --disable-$SPKG_NAME\" to uninstall"
        ;;
    *)
        AC_MSG_ERROR([The content of "$SPKG_TYPE_FILE" must be 'base', 'standard', 'optional', or 'experimental'])
        ;;
    esac

    case "$SPKG_TYPE" in
    optional|experimental)
        stampfile=""
        for f in "$SAGE_SPKG_INST/$SPKG_NAME"-*; do
            AS_IF([test -r "$f"], [
                AS_IF([test -n "$stampfile"], [
                    AC_MSG_ERROR(m4_normalize([
                        multiple installation records for $SPKG_NAME:
                        m4_newline($(ls -l "$SAGE_SPKG_INST/$SPKG_NAME"-*))
                        m4_newline([only one should exist, so please delete some or all
                        of these files and re-run \"$srcdir/configure\"])
                    ]))
                ])
                stampfile=yes
            ])
        done
        ;;
    esac

    # Trac #29629: Temporary solution for Sage 9.1: Do not advertise installing pip packages
    # using ./configure --enable-SPKG
    if test -f "$DIR/requirements.txt"; then
        message="$SPKG_TYPE pip package; use \"./sage -i $SPKG_NAME\" to install"
        uninstall_message="$SPKG_TYPE pip package (installed)"
    fi

    SAGE_PACKAGE_VERSIONS="${SAGE_PACKAGE_VERSIONS}$(printf '\nvers_')${SPKG_NAME} = ${SPKG_VERSION}"

        AS_VAR_PUSHDEF([sage_spkg_install], [sage_spkg_install_${SPKG_NAME}])dnl
        AS_VAR_PUSHDEF([sage_require], [sage_require_${SPKG_NAME}])dnl
        AS_VAR_PUSHDEF([sage_use_system], [sage_use_system_${SPKG_NAME}])dnl

        # If $sage_spkg_install_{SPKG_NAME} is set to no, then set inst_<pkgname> to
        # some dummy file to skip the installation. Note that an explicit
        # "./sage -i SPKG_NAME" will still install the package.
        AS_VAR_IF([sage_spkg_install], [no], [
            dnl We will use the system package (or not required for this platform.)
            SAGE_DUMMY_PACKAGES="${SAGE_DUMMY_PACKAGES} \\$(printf '\n    ')${SPKG_NAME}"
            AS_VAR_IF([sage_require], [yes], [ message="using system package; SPKG will not be installed"
            ],                               [ message="not required on your platform; SPKG will not be installed"
            ])
        ], [
            dnl We won't use the system package.
            SAGE_BUILT_PACKAGES="${SAGE_BUILT_PACKAGES} \\$(printf '\n    ')${SPKG_NAME}"
            AS_VAR_SET_IF([sage_use_system], [
                AS_VAR_COPY([reason], [sage_use_system])
                AS_CASE([$reason],
                [yes],                       [ message="no suitable system package; $message"
                                               AS_VAR_APPEND([SAGE_NEED_SYSTEM_PACKAGES], [" $SPKG_NAME"])
                                             ],
                [installed],                 [ message="already installed as an SPKG$uninstall_message" ],
                                             [ message="$reason; $message" ])
            ], [
                # Package does not use spkg-configure.m4 yet
                message="does not support check for system package; $message"
            ])
        ])

    dnl Trac #29124: Do not talk about underscore club
    case "$SPKG_NAME" in
    _*)
        ;;
    *)
        formatted_message=$(printf '%-45s%s' "$SPKG_NAME-$SPKG_VERSION:" "$message")
        AC_MSG_RESULT([$formatted_message])
        ;;
    esac

        AS_VAR_POPDEF([sage_use_system])dnl
        AS_VAR_POPDEF([sage_require])dnl
        AS_VAR_POPDEF([sage_spkg_install])dnl

    # Packages that should be included in the source distribution
    # This includes all standard packages and two special cases
    case "$SPKG_NAME" in
    mpir)
        in_sdist=yes
        ;;
    esac

    # Determine package source
    #
    if test -f "$DIR/requirements.txt"; then
        SPKG_SOURCE=pip
        # Since pip packages are downloaded and installed by pip, we don't
        # include them in the source tarball. At the time of this writing,
        # all pip packages are optional.
        in_sdist=no
    elif test ! -f "$DIR/checksums.ini"; then
        SPKG_SOURCE=script
        # We assume that either (a) the sources for an optional script
        # package will be downloaded by the script, or (b) that a
        # standard script package's sources are already a part of the
        # sage repository (and thus the release tarball). As a result,
        # we don't need to download the sources, which is what
        # "in_sdist" really means. At the time of this writing, the
        # only standard script packages are sage_conf and sagelib.
        # The source of sage_conf is included under build/pkgs/sage_conf/src,
        # and the source of sagelib is provided by symlinks in
        # build/pkgs/sagelib/src.
        in_sdist=no
    else
        SPKG_SOURCE=normal
    fi

    if test "$in_sdist" = yes; then
        SAGE_SDIST_PACKAGES="${SAGE_SDIST_PACKAGES} \\$(printf '\n    ')${SPKG_NAME}"
    fi

    # Determine package dependencies
    #
    DEP_FILE="$DIR/dependencies"
    if test -f "$DEP_FILE"; then
        # - the # symbol is treated as comment which is removed
        DEPS=`sed 's/^ *//; s/ *#.*//; q' $DEP_FILE`
    else
        ORDER_ONLY_DEPS=""
        case "$SPKG_SOURCE" in
        pip)
            ORDER_ONLY_DEPS='pip'
            ;;
        esac
        if test -n "$ORDER_ONLY_DEPS"; then
            DEPS="| $ORDER_ONLY_DEPS"
        else
            DEPS=""
        fi
    fi

    SAGE_PACKAGE_DEPENDENCIES="${SAGE_PACKAGE_DEPENDENCIES}$(printf '\ndeps_')${SPKG_NAME} = ${DEPS}"

    # Determine package build rules
    case "$SPKG_SOURCE" in
    pip)
        SAGE_PIP_PACKAGES="${SAGE_PIP_PACKAGES} \\$(printf '\n    ')${SPKG_NAME}"
        ;;
    script)
        SAGE_SCRIPT_PACKAGES="${SAGE_SCRIPT_PACKAGES} \\$(printf '\n    ')${SPKG_NAME}"
        ;;
    normal)
        SAGE_NORMAL_PACKAGES="${SAGE_NORMAL_PACKAGES} \\$(printf '\n    ')${SPKG_NAME}"
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
AC_SUBST([SAGE_OPTIONAL_CLEANED_PACKAGES])
AC_SUBST([SAGE_SDIST_PACKAGES])
])

AC_DEFUN([SAGE_SYSTEM_PACKAGE_NOTICE], [
    AS_IF([test -n "$SAGE_NEED_SYSTEM_PACKAGES"], [
        AC_MSG_NOTICE([

    notice: the following SPKGs did not find equivalent system packages:

       $SAGE_NEED_SYSTEM_PACKAGES
        ])
        AC_MSG_CHECKING([for the package system in use])
        SYSTEM=$(build/bin/sage-guess-package-system 2>& AS_MESSAGE_FD)
        AC_MSG_RESULT([$SYSTEM])
        AS_IF([test $SYSTEM != unknown], [
            SYSTEM_PACKAGES=$(build/bin/sage-get-system-packages $SYSTEM $SAGE_NEED_SYSTEM_PACKAGES)
            AS_IF([test -n "$SYSTEM_PACKAGES"], [
                PRINT_SYS="build/bin/sage-print-system-package-command $SYSTEM --verbose=\"    \" --prompt=\"      \$ \" --sudo"
                COMMAND=$(eval "$PRINT_SYS" update && eval "$PRINT_SYS" install $SYSTEM_PACKAGES && SAGE_ROOT="$SAGE_ROOT" eval "$PRINT_SYS" setup-build-env )
                AC_MSG_NOTICE([

    hint: installing the following system packages, if not
    already present, is recommended and may avoid having to
    build them (though some may have to be built anyway):

$COMMAND

    After installation, re-run configure using:

      \$ ./config.status --recheck && ./config.status
                ])
            ], [
                AC_MSG_NOTICE([No equivalent system packages for $SYSTEM are known to Sage])
            ])
        ])
    ])
])
