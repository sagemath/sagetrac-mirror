#!/bin/sh

# This script gets called from CI to run doctests in the sagemath build

# Usage: ./test-doctest.sh IMAGE-NAME --new|--short|--long

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

SETUP=":"

case "$2" in
    --new)
        SETUP='git reset `git log --author release@sagemath.org -1 --format=%H`'
        DOCTEST_PARAMETERS="--long --new"
        ;;
    --short)
        DOCTEST_PARAMETERS="--short --all"
        ;;
    --long)
        DOCTEST_PARAMETERS="--long --all"
        ;;
    *)
        exit 1
        ;;
esac

docker run "$1" "$SETUP; \
sage -tp $DOCTEST_PARAMETERS || \
sage -tp --failed $DOCTEST_PARAMETERS || \
sage -tp --failed $DOCTEST_PARAMETERS"
