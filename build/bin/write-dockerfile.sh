#! /usr/bin/env bash
## Write a Dockerfile to stdout that tests that the packages listed in the debian.txt/fedora.txt files of standard spkg exist
## and satisfy the requirements tested by spkg-configure.m4
## This is called by $SAGE_ROOT/tox.ini
set -e
SYSTEM="${1:-debian}"
shopt -s extglob
TYPE_PATTERN="${2:-standard}"
WITH_SYSTEM_SPKG="${3:-yes}"
IGNORE_MISSING_SYSTEM_PACKAGES="${4:-no}"
#
STRIP_COMMENTS="sed s/#.*//;"
SAGE_ROOT=.
SYSTEM_PACKAGES=$(echo $(${STRIP_COMMENTS} $SAGE_ROOT/build/pkgs/$SYSTEM{,-bootstrap}.txt))
CONFIGURE_ARGS="--enable-option-checking "
for PKG_SCRIPTS in build/pkgs/*; do
    if [ -d $PKG_SCRIPTS ]; then
        PKG_BASE=$(basename $PKG_SCRIPTS)
        SYSTEM_PACKAGES_FILE=$PKG_SCRIPTS/distros/$SYSTEM.txt
        PKG_TYPE=$(cat $PKG_SCRIPTS/type)
        if [ -f $SYSTEM_PACKAGES_FILE -a -f $PKG_SCRIPTS/spkg-configure.m4 ]; then
           PKG_SYSTEM_PACKAGES=$(echo $(${STRIP_COMMENTS} $SYSTEM_PACKAGES_FILE))
           if [ -n "PKG_SYSTEM_PACKAGES" ]; then
               case "$PKG_TYPE" in
                   $TYPE_PATTERN)
                       SYSTEM_PACKAGES+=" $PKG_SYSTEM_PACKAGES"
                       if [ -f $PKG_SCRIPTS/spkg-configure.m4 ]; then
                           CONFIGURE_ARGS+="--with-system-$PKG_BASE=${WITH_SYSTEM_SPKG} "
                       fi
                       ;;
               esac
           fi
        fi
    fi
done
echo "# Automatically generated by SAGE_ROOT/build/bin/write-dockerfile.sh"
echo "# the :comments: separate the generated file into sections"
echo "# to simplify writing scripts that customize this file"
case $SYSTEM in
    debian*|ubuntu*)
        cat <<EOF
ARG BASE_IMAGE=ubuntu:latest
FROM \${BASE_IMAGE}
EOF
        UPDATE="apt-get update &&"
        INSTALL="DEBIAN_FRONTEND=noninteractive apt-get install -qqq --no-install-recommends --yes"
        CLEAN="&& apt-get clean"
        ;;
    fedora*|redhat*|centos*)
        cat <<EOF
ARG BASE_IMAGE=fedora:latest
FROM \${BASE_IMAGE}
EOF
        INSTALL="yum install -y"
        ;;
    arch*)
        # https://hub.docker.com/_/archlinux/
        cat <<EOF
ARG BASE_IMAGE=archlinux:latest
FROM \${BASE_IMAGE}
EOF
        INSTALL="pacman -Syu --noconfirm"
        ;;
    conda*)
        cat <<EOF
ARG BASE_IMAGE=continuumio/miniconda3:latest
FROM \${BASE_IMAGE}
ARG USE_CONDARC=condarc.yml
ADD *condarc*.yml /tmp/
RUN echo \${CONDARC}; cd /tmp && conda config --stdin < \${USE_CONDARC}
RUN conda update -n base conda
EOF
        INSTALL="conda install --update-all --quiet --yes"
        EXISTS="2>/dev/null >/dev/null conda search -f"
        #EXISTS="conda search -f"
        ;;
    *)
        echo "Not implemented: package installation for SYSTEM=$SYSTEM" >&2
        exit 1
        ;;
esac
cat <<EOF
#:packages:
ENV PACKAGES="$SYSTEM_PACKAGES"
EOF
case "$IGNORE_MISSING_SYSTEM_PACKAGES" in
    no)
        cat <<EOF
RUN $UPDATE $INSTALL $SYSTEM_PACKAGES $CLEAN
EOF
        ;;
    yes)
        if [ -n "$EXISTS" ]; then
            # Filter by existing packages, try to install these in one shot; fall back to one by one.
            cat <<EOF
RUN $UPDATE EXISTING_PACKAGES=""; for pkg in \$PACKAGES; do echo -n .; if $EXISTS \$pkg; then EXISTING_PACKAGES="\$EXISTING_PACKAGES \$pkg"; echo -n "\$pkg"; fi; done; $INSTALL \$EXISTING_PACKAGES || (echo "Trying again one by one:"; for pkg in \$EXISTING_PACKAGES; do echo "Trying to install \$pkg"; $INSTALL \$pkg || echo "(ignoring error)"; done); : $CLEAN
EOF
        else
            # Try in one shot, fall back to one by one.  Separate "RUN" commands
            # for caching by docker.
            cat <<EOF
RUN $UPDATE $INSTALL \${PACKAGES} || echo "(ignoring error)"
EOF
            for pkg in $SYSTEM_PACKAGES; do
                cat <<EOF
RUN $INSTALL $pkg || echo "(ignoring error)"
EOF
            done
            if [ -n "$CLEAN" ]; then
                cat <<EOF
RUN : $CLEAN
EOF
            fi
        fi
        ;;
    *)
        echo "Argument IGNORE_MISSING_SYSTEM_PACKAGES must be yes or no"
        ;;
esac
cat <<EOF
#:bootstrapping:
RUN mkdir -p /sage
WORKDIR /sage
ADD Makefile VERSION.txt README.md bootstrap configure.ac sage ./
ADD m4 ./m4
ADD build ./build
ADD src/bin/sage-version.sh src/bin/sage-version.sh
RUN ./bootstrap
#:configuring:
ADD src/ext src/ext
ADD src/bin src/bin
ADD src/Makefile.in src/Makefile.in
ARG EXTRA_CONFIGURE_ARGS=""
EOF
if [ ${WITH_SYSTEM_SPKG} = "force" ]; then
    cat <<EOF
RUN echo "****** Configuring: ./configure --enable-build-as-root $CONFIGURE_ARGS \${EXTRA_CONFIGURE_ARGS} *******"; ./configure --enable-build-as-root $CONFIGURE_ARGS \${EXTRA_CONFIGURE_ARGS} || (echo "********** configuring without forcing ***********"; cat config.log; ./configure --enable-build-as-root; cat config.log; exit 1)
EOF
else
    cat <<EOF
RUN echo "****** Configuring: ./configure --enable-build-as-root $CONFIGURE_ARGS \${EXTRA_CONFIGURE_ARGS} *******"; ./configure --enable-build-as-root $CONFIGURE_ARGS \${EXTRA_CONFIGURE_ARGS} || (cat config.log; exit 1)
EOF
fi
cat <<EOF
# We first compile base-toolchain because otherwise lots of packages are missing their dependency on 'patch'
ARG NUMPROC=8
ENV MAKE="make -j\${NUMPROC}"
ARG USE_MAKEFLAGS="-k"
#:toolchain:
RUN make \${USE_MAKEFLAGS} base-toolchain
#:make:
# Avoid running the lengthy testsuite of the following.
RUN make \${USE_MAKEFLAGS} cython
# By default, compile something tricky but that does not take too long. scipy uses BLAS.
ARG TARGETS="scipy"
RUN SAGE_CHECK=yes SAGE_CHECK_PACKAGES="!r,!python3,!python2,!nose" make \${USE_MAKEFLAGS} \${TARGETS}
#:end:
EOF
