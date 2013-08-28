#!/bin/sh
# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#

MODULES="src/sage src/c_lib src/bin src/doc src/ext"

# for now avoid using bash as not everyone has that installed
CONFIG_SHELL=/bin/sh
export CONFIG_SHELL
ACLOCAL_FLAGS=-Im4

echo "Running aclocal..."
aclocal $ACLOCAL_FLAGS || exit 1

echo "Running automake..."
automake -a -c --gnu --add-missing  || exit 1

echo "Running autoconf..."
autoconf || exit 1

if test "x$1" = xhere; then
	exit 0;
fi

for i in $MODULES; do
	( cd $i;
	  echo autogen in $( pwd )
	  ./autogen.sh || ( echo trouble in $i/autogen.sh; exit 1 )
	) || exit 1
done
