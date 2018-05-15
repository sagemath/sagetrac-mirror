#!/bin/sh

# This script gets called from CircleCI to setup some common environment
# variables and print some diagnostic messages.

# Source this script before the individual CI steps.

# ****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

# CircleCI has no mechanism to hide secret variables.
# Therefore we roll our own to protect $SECRET_* variables.
. .ci/protect-secrets.sh
# Collect debug infos about the system we are running on
.ci/describe-system.sh
# Set MAKE and NTHREADS according to the machine we are running on
. .ci/setup-make-parallelity.sh

# Set DOCKER_TAG according to the current branch/tag
export DOCKER_TAG=${CIRCLE_TAG:-$CIRCLE_BRANCH}
. .ci/update-env.sh

# Select ARTIFACT_BASE depending on the current branch/tag
TODO: port the changes from 24655 here to build tags from scratch.
