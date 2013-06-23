#!/bin/sh
# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
# Run the various GNU autotools to bootstrap the build

# this is supposedly called by toplevel autogen.sh.

echo "Running aclocal..."
aclocal $ACLOCAL_FLAGS || exit 1

echo "Running automake..."
automake -a -c --gnu --add-missing  || exit 1

echo "Running autoconf..."
autoconf || exit 1
