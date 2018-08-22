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

# The image's entrypoints might not accept arbitrary shell scripts, so we feed
# our commands into the container through a FIFO. (Note that we want the
# entrypoints to run and not override the entrypoint to be sh as we might need
# the setup steps performed by the entrypoint.)
mkfifo /tmp/stdin
# Pipe everything that comes into the FIFO into our image (without closing the FIFO)
(while true; do cat /tmp/stdin; done) | docker run -i "$1" bash &
# When this script exits the background job won't die, so the docker container
# is going to keep running. Make sure the entitre process tree goes away.
trap "exit" INT TERM
trap "kill 0" EXIT

case "$2" in
    --new)
        echo 'git reset `git log --author release@sagemath.org -1 --format=%H`' > /tmp/stdin
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
echo <<EOF > /tmp/stdin
sage -tp $DOCTEST_PARAMETERS ||
sage -tp --failed $DOCTEST_PARAMETERS ||
sage -tp --failed $DOCTEST_PARAMETERS
EOF
