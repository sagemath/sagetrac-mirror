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

# Some doctest flavours require additional setup in the docker image before
# they are ready to run.
export DOCTEST_SETUP=":"

case "$2" in
    --new)
        export DOCTEST_SETUP="git reset `git log --author release@sagemath.org -1 --format=%H`"
        export DOCTEST_PARAMETERS="--long --new"
        ;;
    --short)
        export DOCTEST_PARAMETERS="--short --all"
        ;;
    --long)
        export DOCTEST_PARAMETERS="--long --all"
        ;;
    *)
        exit 1
        ;;
esac

# Run tests once, and then try the failing files twice to work around flaky doctests.
docker run --entrypoint sh -e DOCTEST_PARAMETERS "$1" -c 'sh -c "$DOCTEST_SETUP"
                                                          sage -tp $DOCTEST_PARAMETERS ||
                                                          sage -tp --failed $DOCTEST_PARAMETERS ||
                                                          sage -tp --failed $DOCTEST_PARAMETERS'

