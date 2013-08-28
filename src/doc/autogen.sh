#! /bin/sh
#
# $Id: autogen.sh,v 1.1 2009-10-23 12:01:23 felix Exp $
#
# Run the various GNU autotools to bootstrap the build
# system.  Should only need to be done once.

# for now avoid using bash as not everyone has that installed
CONFIG_SHELL=/bin/sh
export CONFIG_SHELL
# ACLOCAL_FLAGS=-Im4

# libtoolize -c -f -i

echo "Running aclocal..."
aclocal $ACLOCAL_FLAGS || exit 1

echo "Running automake..."
automake -a -c --gnu --add-missing  || exit 1

echo "Running autoconf..."
autoconf || exit 1

echo "not Running configure..."
##./configure $@ || exit 1


