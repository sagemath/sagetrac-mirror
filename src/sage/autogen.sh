#!/bin/sh
ACLOCAL_FLAGS=-Im4

libtoolize
aclocal $ACLOCAL_FLAGS
autoheader
autoconf
automake --add-missing
