#! /bin/bash
set -e # exit on errors

function finish {
    for a in `find logs -name "*.log"`; do
        echo "###### $a ######"
        cat $a
    done
}

#trap finish EXIT

export MAKE="make -j16"

export SAGE_INSTALL_CCACHE=yes

export V=0
export SAGE_PV="env PYTHONUNBUFFERED=1 sage-progress-meter"
make all-toolchain && make $TARGETS

