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
        # We need an image that contains a .git directory as this is not
        # contained in the sagemath image (and also not in sagemath-dev.)
        # Note that we can not mount our own .git as the docker daemon might
        # not run on the current host.
        docker create --name sagemath-git-build "$1"
        docker cp `pwd`/.git sagemath-git-build:/home/sage/sage/.git
        docker commit sagemath-git-build sagemath-git
        # Replace $1 so that the following code uses that image instead of the
        # original "$1"
        shift
        set -- sagemath-git "$@"

        SETUP='sudo apt-get update && sudo apt-get install -y git && \
               cd /home/sage/sage && \
               sudo chown -R sage:sage .git && \
               git reset --hard && \
               git reset `git describe --abbrev=0 --tags`'
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

docker run "$1" "$SETUP && \
                 (sage -tp $DOCTEST_PARAMETERS || \
                  sage -tp --failed $DOCTEST_PARAMETERS || \
                  sage -tp --failed $DOCTEST_PARAMETERS)"
