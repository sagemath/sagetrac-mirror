#! /usr/bin/env bash
## Write a Dockerfile to stdout that tests that the packages listed in the debian.txt/fedora.txt files of standard spkg exist
## and satisfy the requirements tested by spkg-configure.m4
## This is called by $SAGE_ROOT/tox.ini
set -e
SYSTEM="${1:-debian}"
shopt -s extglob
SAGE_PACKAGE_LIST_ARGS="${2:- --has-file=spkg-configure.m4 :standard:}"
WITH_SYSTEM_SPKG="${3:-yes}"
IGNORE_MISSING_SYSTEM_PACKAGES="${4:-no}"
EXTRA_SAGE_PACKAGES="${5:-_bootstrap}"
#
STRIP_COMMENTS="sed s/#.*//;"
SAGE_ROOT=.
export PATH="$SAGE_ROOT"/build/bin:$PATH
SYSTEM_PACKAGES=
CONFIGURE_ARGS="--enable-option-checking "
for PKG_BASE in $($SAGE_ROOT/sage -package list --has-file=distros/$SYSTEM.txt $SAGE_PACKAGE_LIST_ARGS) $EXTRA_SAGE_PACKAGES; do
    PKG_SCRIPTS="$SAGE_ROOT"/build/pkgs/$PKG_BASE
    if [ -d $PKG_SCRIPTS ]; then
       SYSTEM_PACKAGES_FILE=$PKG_SCRIPTS/distros/$SYSTEM.txt
       PKG_SYSTEM_PACKAGES=$(echo $(${STRIP_COMMENTS} $SYSTEM_PACKAGES_FILE))
       if [ -n "PKG_SYSTEM_PACKAGES" ]; then
           SYSTEM_PACKAGES+=" $PKG_SYSTEM_PACKAGES"
           if [ -f $PKG_SCRIPTS/spkg-configure.m4 ]; then
               CONFIGURE_ARGS+="--with-system-$PKG_BASE=${WITH_SYSTEM_SPKG} "
           fi
       fi
    fi
done
echo "# Automatically generated by SAGE_ROOT/build/bin/write-dockerfile.sh"
echo "# the :comments: separate the generated file into sections"
echo "# to simplify writing scripts that customize this file"
RUN=RUN
case $SYSTEM in
    debian*|ubuntu*)
        cat <<EOF
ARG BASE_IMAGE=ubuntu:latest
FROM \${BASE_IMAGE} as with-system-packages
EOF
        if [ -n "$DIST_UPGRADE" ]; then
            cat <<EOF
RUN sed -i.bak $DIST_UPGRADE /etc/apt/sources.list && apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y dist-upgrade
EOF
        fi
        EXISTS="2>/dev/null >/dev/null apt-cache show"
        UPDATE="apt-get update &&"
        INSTALL="DEBIAN_FRONTEND=noninteractive apt-get install -qqq --no-install-recommends --yes"
        CLEAN="&& apt-get clean"
        if [ -n "$EXTRA_REPOSITORY" ]; then
            cat <<EOF
RUN $UPDATE $INSTALL software-properties-common
RUN add-apt-repository $EXTRA_REPOSITORY
EOF
        fi
        ;;
    fedora*|redhat*|centos*)
        cat <<EOF
ARG BASE_IMAGE=fedora:latest
FROM \${BASE_IMAGE} as with-system-packages
EOF
        EXISTS="2>/dev/null >/dev/null yum install -y --downloadonly"
        INSTALL="yum install -y"
        ;;
    gentoo*)
        cat <<EOF
ARG BASE_IMAGE=sheerluck/sage-on-gentoo-stage4:latest
FROM \${BASE_IMAGE} as with-system-packages
EOF
        EXISTS="2>/dev/null >/dev/null emerge -f"
        UPDATE="" # not needed. "FROM gentoo/portage" used instead
        INSTALL="emerge -DNut --with-bdeps=y --complete-graph=y"
        ;;
    slackware*)
        # https://docs.slackware.com/slackbook:package_management
        cat <<EOF
ARG BASE_IMAGE=vbatts/slackware:latest
FROM \${BASE_IMAGE} as with-system-packages
EOF
        # slackpkg install ignores packages that it does not know, so we do not have to filter
        EXISTS="true"
        UPDATE="slackpkg update &&"
        INSTALL="slackpkg install"
        ;;
    arch*)
        # https://hub.docker.com/_/archlinux/
        cat <<EOF
ARG BASE_IMAGE=archlinux:latest
FROM \${BASE_IMAGE} as with-system-packages
EOF
        UPDATE="pacman -Sy &&"
        EXISTS="pacman -Si"
        INSTALL="pacman -Su --noconfirm"
        ;;
    nix*)
        # https://hub.docker.com/r/nixos/nix
        cat <<EOF
ARG BASE_IMAGE=nixos/nix:latest
FROM \${BASE_IMAGE} as with-system-packages
RUN nix-channel --add https://nixos.org/channels/nixpkgs-unstable nixpkgs
RUN nix-channel --update
EOF
        INSTALL="nix-env --install"
        RUN="RUN nix-shell --packages \$PACKAGES --run "\'
        ENDRUN=\'
        ;;
    void*)
	# https://hub.docker.com/r/voidlinux/masterdir-x86_64-musl
	cat <<EOF
ARG BASE_IMAGE=voidlinux:masterdir-x86_64-musl
FROM \${BASE_IMAGE} as with-system-packages
EOF
        UPDATE="xbps-install -Su &&"
        EXISTS="xbps-query"
        INSTALL="xbps-install --yes"
        ;;
    opensuse*)
	cat <<EOF
ARG BASE_IMAGE=opensuse/leap:latest
FROM \${BASE_IMAGE} as with-system-packages
EOF
        UPDATE="zypper refresh &&"
        INSTALL="zypper --ignore-unknown install --no-confirm --auto-agree-with-licenses --no-recommends --details"
        ;;
    conda*)
        cat <<EOF
ARG BASE_IMAGE=continuumio/miniconda3:latest
FROM \${BASE_IMAGE} as with-system-packages
ARG USE_CONDARC=condarc.yml
ADD *condarc*.yml /tmp/
RUN echo \${CONDARC}; cd /tmp && conda config --stdin < \${USE_CONDARC}
RUN conda update -n base conda
RUN ln -sf /bin/bash /bin/sh
EOF
        # On this image, /bin/sh -> /bin/dash;
        # but some of the scripts in /opt/conda/etc/conda/activate.d
        # from conda-forge (as of 2020-01-27) contain bash-isms:
        # /bin/sh: 5: /opt/conda/etc/conda/activate.d/activate-binutils_linux-64.sh: Syntax error: "(" unexpected
        # The command '/bin/sh -c . /opt/conda/etc/profile.d/conda.sh; conda activate base;  ./bootstrap' returned a non-zero code
        # We just change the link to /bin/bash.
        INSTALL="conda install --update-all --quiet --yes"
        EXISTS="2>/dev/null >/dev/null conda search -f"
        #EXISTS="conda search -f"
        CLEAN="&& conda info; conda list"
        RUN="RUN . /opt/conda/etc/profile.d/conda.sh; conda activate base; "  # to activate the conda env
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

FROM with-system-packages as bootstrapped
#:bootstrapping:
RUN mkdir -p /sage
WORKDIR /sage
ADD Makefile VERSION.txt COPYING.txt condarc.yml README.md bootstrap configure.ac sage .homebrew-build-env tox.ini Pipfile.m4 ./
ADD src/doc/bootstrap src/doc/bootstrap
ADD src/bin src/bin
ADD src/Pipfile.m4 src/pyproject.toml.m4 src/requirements.txt.m4 src/setup.cfg.m4 src/
ADD m4 ./m4
ADD pkgs pkgs
ADD build ./build
ARG BOOTSTRAP=./bootstrap
$RUN sh -x -c "\${BOOTSTRAP}" $ENDRUN

FROM bootstrapped as configured
#:configuring:
RUN mkdir -p logs/pkgs; ln -s logs/pkgs/config.log config.log
ARG EXTRA_CONFIGURE_ARGS=""
EOF
if [ ${WITH_SYSTEM_SPKG} = "force" ]; then
    cat <<EOF
$RUN echo "****** Configuring: ./configure --enable-build-as-root $CONFIGURE_ARGS \${EXTRA_CONFIGURE_ARGS} *******"; ./configure --enable-build-as-root $CONFIGURE_ARGS \${EXTRA_CONFIGURE_ARGS} || (echo "********** configuring without forcing ***********"; cat config.log; ./configure --enable-build-as-root; cat config.log; exit 1) $ENDRUN
EOF
else
    cat <<EOF
$RUN echo "****** Configuring: ./configure --enable-build-as-root $CONFIGURE_ARGS \${EXTRA_CONFIGURE_ARGS} *******"; ./configure --enable-build-as-root $CONFIGURE_ARGS \${EXTRA_CONFIGURE_ARGS} || (cat config.log; exit 1) $ENDRUN
EOF
fi
cat <<EOF

FROM configured as with-base-toolchain
# We first compile base-toolchain because otherwise lots of packages are missing their dependency on 'patch'
ARG NUMPROC=8
ENV MAKE="make -j\${NUMPROC}"
ARG USE_MAKEFLAGS="-k V=0"
ENV SAGE_CHECK=warn
ENV SAGE_CHECK_PACKAGES="!cython,!r,!python3,!gap,!cysignals,!linbox,!git,!ppl,!cmake,!rpy2,!sage_sws2rst"
#:toolchain:
$RUN make \${USE_MAKEFLAGS} base-toolchain $ENDRUN

FROM with-base-toolchain as with-targets-pre
ARG NUMPROC=8
ENV MAKE="make -j\${NUMPROC}"
ARG USE_MAKEFLAGS="-k V=0"
ENV SAGE_CHECK=warn
ENV SAGE_CHECK_PACKAGES="!gfan,!cython,!r,!python3,!gap,!cysignals,!linbox,!git,!ppl,!cmake,!rpy2,!sage_sws2rst"
#:make:
ARG TARGETS_PRE="all-sage-local"
$RUN make SAGE_SPKG="sage-spkg -y -o" \${USE_MAKEFLAGS} \${TARGETS_PRE} $ENDRUN

FROM with-targets-pre as with-targets
ARG NUMPROC=8
ENV MAKE="make -j\${NUMPROC}"
ARG USE_MAKEFLAGS="-k V=0"
ENV SAGE_CHECK=warn
ENV SAGE_CHECK_PACKAGES="!gfan,!cython,!r,!python3,!gap,!cysignals,!linbox,!git,!ppl,!cmake,!rpy2,!sage_sws2rst"
ADD src src
ARG TARGETS="build"
$RUN make SAGE_SPKG="sage-spkg -y -o" \${USE_MAKEFLAGS} \${TARGETS} $ENDRUN

FROM with-targets as with-targets-optional
ARG NUMPROC=8
ENV MAKE="make -j\${NUMPROC}"
ARG USE_MAKEFLAGS="-k V=0"
ENV SAGE_CHECK=warn
ENV SAGE_CHECK_PACKAGES="!gfan,!cython,!r,!python3,!gap,!cysignals,!linbox,!git,!ppl,!cmake,!rpy2,!sage_sws2rst"
ARG TARGETS_OPTIONAL="ptest"
$RUN make SAGE_SPKG="sage-spkg -y -o" \${USE_MAKEFLAGS} \${TARGETS_OPTIONAL} || echo "(error ignored)" $ENDRUN

#:end:
EOF
