#!/bin/sh
# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#

MODULES=src/c_lib

for i in $MODULES; do
	( cd $i;
	  echo autogen in $( pwd )
	  ./autogen.sh || echo trouble in $i/autogen.sh
	)
done
