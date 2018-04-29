#!/bin/sh

# This script gets called from CI to pull the Sage docker images that were
# built during the "build" phase to pull to the connected docker daemon
# (likely a docker-in-docker.)
# This script expects a single parameter, the base name of the docker image
# such as sagemath or sagemath-dev.
# The variable $DOCKER_IMAGE is set to the full name of the pulled image;
# source this script to use it.

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

# Pull the built images from the dockerhub registry and give them the original
# names they had after built.
# We require $DOCKER_USER and $SECRET_DOCKER_PASS to be set. Otherwise we would
# be pulling some stale images here. (Sadly, CircleCI does not provide us with
# an integrated container registry like GitLab does.)
if [ -z "$DOCKER_USER" -o -z "$SECRET_DOCKER_PASS" ]; then
  echo "DOCKER_USER/SECRET_DOCKER_PASS variables have not been configured in your Continuous Integration setup. Not pulling as the images would not be the one that has just been built."
fi

export DOCKER_IMAGE="$DOCKER_USER/$1:$DOCKER_TAG"
docker pull $DOCKER_IMAGE
