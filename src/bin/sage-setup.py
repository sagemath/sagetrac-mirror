#!/usr/bin/env python

# this program is intended to be called instead of "python setup.py"
# see #14804

import distutils.file_util, distutils.dir_util, os
import distutils.command.install
#import distutils.util

destdir = ""

try:
	destdir = os.environ["DESTDIR"]
except:
	pass

prefix = os.environ["SAGE_PREFIX"]

import sys
sys.argv[0] = "setup.py"

newargv=[]

# distutils does not support destdir. hack it in.
# good idea? better put into spkg-install files later...
for i in sys.argv:
	newargv.append(i)
	if i == "install":
		newargv.append("--prefix=" + destdir + prefix)

sys.argv = newargv
c = os.getcwd()
s = os.path.join(c, "setup.py")
sys.path.reverse()
sys.path.append(os.getcwd())
sys.path.reverse()
__file__ = c + "/setup.py"

execfile("setup.py")
