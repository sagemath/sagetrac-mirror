#!/bin/sh

# ****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

set -ex

# Determine the number of threads that can run simultaneously on this system
# (we might not have nproc available.)
# Note that this value is incorrect for some CI providers (notably CircleCI:
# https://circleci.com/docs/2.0/configuration-reference/#resource_class) which
# provision fewer vCPUs than shown in /proc/cpuinfo. Also, setting this value
# too high can lead to RAM being insufficient, so it's best to set the NTHREADS
# variable manually in your CI configuration.
if [ -z "$CPUTHREADS" ]; then
    CPUTHREADS=`grep -E '^processor' /proc/cpuinfo | wc -l`
fi
if [ -z "$RAMTHREADS" ]; then
    RAMTHREADS=$(( `grep MemTotal /proc/meminfo | awk '{ print $2 }'` / 1024 / 1024 / 2 ))
    if [ $RAMTHREADS = 0 ];then
        RAMTHREADS=1;
    fi
fi
if [ -z "$NTHREADS" ]; then
    NTHREADS=$([ $RAMTHREADS -le $CPUTHREADS ] && echo "$RAMTHREADS" || echo "$CPUTHREADS")
fi
export NTHREADS="$NTHREADS"

# Set -j and -l for make (though -l is probably stripped by Sage)
if [ -z "$MAKEOPTS" ]; then
    MAKEOPTS="-j $RAMTHREADS -l $((CPUTHREADS-1)).8"
fi
export MAKEOPTS="$MAKEOPTS"

# Not all parts of Sage seem to honor MAKEOPTS, so the current way of telling
# the system which concurrency to use, seems to be setting $MAKE.
if [ -z "$MAKE" ];then
    MAKE="make $MAKEOPTS"
fi
export MAKE="$MAKE"
