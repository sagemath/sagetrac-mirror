#!/bin/sh

# This script gets called from CI to run doctests in the sagemath build

# Usage: ./test-doctest.sh IMAGE-NAME

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

# Run tests once, and then try the failing files twice to work around flaky doctests.
docker run --entrypoint sh "$1" -c 'sage -tp --all $DOCTEST_PARAMETERS ||
                                    sage -tp --all $DOCTEST_PARAMETERS --failed ||
                                    sage -tp --all $DOCTEST_PARAMETERS --failed'

