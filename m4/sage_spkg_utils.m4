###############################################################################
# defines some utility functions for Sage spkgs for use in configure.ac
###############################################################################

AC_DEFUN_ONCE([SAGE_SPKG_UTILS], [
# Usage: newest_version $pkg
# Print version number of latest package $pkg
newest_version() {
    PKG=$[1]
    if test -f "$SAGE_ROOT/build/pkgs/$PKG/package-version.txt" ; then
        AS_ECHO_N(["$PKG-"])
        cat "$SAGE_ROOT/build/pkgs/$PKG/package-version.txt"
    else
        echo "$PKG"
    fi
}

# Outputs the list of packages, filtered by 'type', e.g.:
#
#     filtered_packages_list base
#     filtered_packages_list standard
#     filtered_packages_list optional
#     filtered_packages_list experimental
#
# Or, if you want all packages:
#
#     filtered_packages_list all
#
# Or, if you want all packages which should appear in the source tarball:
#
#     filtered_packages_list sdist
#
# The output consists of: 
#
# PKG_NAME PKG_VERSION PKG_VAR PKG_TYPE
#
# where outputting the package type can be useful when listing all packages

changequote(<,>)
filtered_packages_list() {
    FILTER=$<1>
    # for each package in pkgs/
    for DIR in $SAGE_ROOT/build/pkgs/*; do
        test -d "$DIR" || continue

        PKG_TYPE_FILE="$DIR/type"
        if [ -f "$PKG_TYPE_FILE" ]; then
            PKG_TYPE=`cat $PKG_TYPE_FILE`
        else
            # exit won't necessarily exit 'configure', so signal an error this way, too:
            echo INVALID
            changequote([,])
            AC_MSG_ERROR(["$PKG_TYPE_FILE" is missing.])
            changequote(<,>)
        fi

        # Check consistency of 'DIR/type' file
        case "$PKG_TYPE" in
            base) ;;
            standard) ;;
            optional) ;;
            experimental) ;;
            script) ;;
            pip) ;;
            abstract) ;;
            *)
                # exit won't necessarily exit 'configure', so signal an error this way, too:
                echo INVALID
                changequote([,])
                AC_MSG_ERROR([The content of "$PKG_TYPE_FILE" must be 'base', 'standard', 'optional', 'experimental', 'script' 'pip', or 'abstract'])
                changequote(<,>)
                ;;
        esac

        PKG_NAME=$(basename $DIR)

        # Filter
        selected=no
        if [ "$FILTER" = all ]; then
            selected=yes
        elif [ "$FILTER" = "$PKG_TYPE" ]; then
            selected=yes
        elif [ "$FILTER" = sdist ]; then
            # sdist packages are all standard packages, together with
            # mpir and python2.
            if [ "$PKG_TYPE" = standard ]; then
                selected=yes
            elif [ "$PKG_NAME" = mpir ]; then
                selected=yes
            elif [ "$PKG_NAME" = python2 ]; then
                selected=yes
            fi
        fi
        if [ $selected = yes ]; then
            PKG_VAR="$(echo $PKG_NAME | sed 's/^/inst_/')"
            PKG_VERSION=$(newest_version $PKG_NAME)
            echo "$PKG_NAME $PKG_VERSION $PKG_VAR $PKG_TYPE"
        fi
    done
}
changequote([,])
])
