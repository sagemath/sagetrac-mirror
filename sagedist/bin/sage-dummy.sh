#!/bin/sh

a=`basename $1`
echo >&2 this is not really $a
echo >&2 $a is part of $2
echo >&2 $2 is neither staged nor disabled

exit 1
