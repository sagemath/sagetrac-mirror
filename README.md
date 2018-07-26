<<<<<<< HEAD
#  Sage: Open Source Mathematical Software

>   "Creating a Viable Open Source Alternative to
>    Magma, Maple, Mathematica, and MATLAB"

>   Copyright (C) 2005-2017 The Sage Development Team

   http://www.sagemath.org

The Sage Library is GPLv2+, and included packages have compatible OSS
licenses (see COPYING.txt). Over 400 people have contributed code to
Sage. Please see the following web page for a list:

   http://www.sagemath.org/development-map.html

In many cases, documentation for modules and functions list the
authors.


Getting Started
---------------

This README.txt contains build instructions for Sage. If you downloaded  
a binary, you do not need to do anything (as long as your system
provides a SSL library ; OpenSSL is recommended) ; just execute:  

    ./sage

from the command line. If you downloaded the sources, please read below  
on how to build Sage and work around common issues.  

If you have questions or encounter problems, please do not hesitate  
to email the sage-support mailing list:  

   http://groups.google.com/group/sage-support


Supported Platforms
-------------------

Sage fully supports several Linux distributions, recent versions of  
Mac OS X, Windows (using virtualization), as well as a number of  
Solaris and OpenSolaris releases.  

There should be no serious bugs in an officially released version of  
Sage on any of the fully supported platforms, but any major issues with  
a particular release will be documented on an errata page:  

   http://wiki.sagemath.org/errata

Ports are in progress to some other, less common platforms. The list of  
supported platforms and their current statuses are given at the  
following web page:  

   http://wiki.sagemath.org/SupportedPlatforms

If you are interested in helping port Sage to a new platform, please let  
us know at the sage-devel mailing list:  

   http://groups.google.com/group/sage-devel


Quick Instructions to Build from Source
---------------------------------------

The following steps briefly outline the process of building Sage from  
source. More detailed instructions, including how to build faster on  
multicore machines, are contained later in this README and in the  
Installation Guide:  

   http://doc.sagemath.org/html/en/installation

__1. Make sure your system has an SSL library and its development
files installed.__

   >* Like Python, on which it is based, Sage uses the OpenSSL library
   for added performance if made available by the operating system. It
   has been shown that Sage can be successfully built against other
   SSL libraries, with some of its features disabled.

__2. Make sure you have the dependencies and 5 GB of free disk space.__

   >* __All Linux versions:__ gcc, make, m4, perl, ranlib, and tar (a
   matching set of gcc, gfortran and g++ will avoid the compilation
   of Sage-specific compilers).__
   
   >* __Fedora or RedHat systems:__ the perl-ExtUtils-MakeMaker package.  
   (install these using your package manager)  

   >* __OS X:__ Xcode. Make sure you have installed the most recent version  
   of Xcode. With recent versions of OS X (OS X Lion or later), you  
   can install Xcode for free from the App Store. For pre-Lion  
   versions of OS X, you can download Xcode from  
   http://developer.apple.com/downloads/.  

   >* With OS X, you also need to install the "command line tools". When  
   using OS X Mavericks, after installing Xcode, run this command from  
   a terminal window:  

         xcode-select --install

   >* Then click "Install" in the pop-up window.

   >* When using OS X Mountain Lion or earlier, you need to install the  
   command line tools from Xcode: run Xcode; then from the File  
   menu, choose "Preferences", then the "Downloads" tab, and then  
   "Install" the Command Line Tools.  

   >* Other platforms: See detailed instructions below.

__3. Extract the tarball:__

         tar xvf sage-*.tar

__4. cd into the Sage directory and type make:__

         cd sage-*/
         make

   That's it! Everything is automatic and non-interactive. The build  
   should work fine on all fully supported platforms. If it does not, we  
   want to know!  

   If you'd like to contribute to Sage, be sure to read the
   Developer's Guide:

   http://doc.sagemath.org/html/en/developer/index.html


Environment Variables
---------------------

There are a lot of environment variables which control the install  
process of Sage, see:  

   http://doc.sagemath.org/html/en/installation/source.html#environment-variables


Implementation
--------------

Sage has significant components written in the following languages:  
   > C/C++, Python, Cython, Lisp, and Fortran. Lisp (ECL), Python, and Cython  
   are built as part of Sage and a GNU Fortran (gfortran) binary is  
   included (OS X only), so you do not need them in order to build Sage.  


More Detailed Instructions to Build from Source
-----------------------------------------------

__1.__ Make sure you have about __5 GB of free disk space__.

__2. Install build dependencies.__

   >* __Linux:__ See quick instructions above.

   >* __OS X:__ Make sure you have XCode version >= 2.4, i.e. "gcc -v" should  
   output build >= 5363. If you don't, go to:  
      http://developer.apple.com/

   > sign up, and download the free XCode package. Only OS X >= 10.4 is  
   supported.

   >* __Solaris and OpenSolaris:__ Building Sage on these platforms is more  
   tricky than on Linux or OS X. For details on how to build Sage on  
   these platforms, see:  
      http://wiki.sagemath.org/solaris

   >* __Windows:__ Download and install VirtualBox, and then download the  
   Sage virtual appliance. For details, see:  
      http://wiki.sagemath.org/SageAppliance

   >* __NOTE:__ On some operating systems, it might be necessary to install  
   gas/as, gld/ld, gnm/nm. On most platforms, these are automatically  
   installed when you install the programs listed above.  

__3. Extract the Sage source tarball__ and cd into a directory with no  
   spaces in it. If you have a machine with 4 processors, say, type  
   the following to configure the build script to perform a parallel  
   compilation of Sage using 4 jobs:

         export MAKE="make -j4"

   >* (With 4 processors, you might also consider "-j5" or "-j6" --  
   building with more jobs than CPU cores can speed things up.)  
   You might in addition pass a "-l" flag to "make": this  
   sets a load limit, so for example if you execute  

         export MAKE="make -j4 -l5.5"

   > then "make" won't start more than one job at a time if the system  
   load average is above 5.5. See:  
   http://www.gnu.org/software/make/manual/make.html#Options-Summary  
   http://www.gnu.org/software/make/manual/make.html#Parallel.  

   > *If you want to run the test suite for each individual spkg as it is  
   installed, type:  

         export SAGE_CHECK="yes"

   > before starting the Sage build. This will run each test suite and  
   will raise an error if any failures occur. Python's test suite has  
   been disabled by default, because it causes failures on most  
   systems. To renable the Python testsuite, set the environment  
   variable SAGE_CHECK_PACKAGES to "python".

   >* To start the build, type:  

         make

__4. Wait about 20 minutes to 14 days, depending on your computer__ (it took  
   about 2 weeks to build Sage on the T-Mobile G1 Android cell phone).  

__5. Type "./sage" to try it out.__  

__6. Optional__: Type "make ptest" to test all examples in the documentation  
   (over 200,000 lines of input!) -- this takes from 10 minutes to  
   several hours. Don't get too disturbed if there are 2 to 3 failures,  
   but always feel free to email the section of logs/ptest.log that  
   contains errors to the sage-support mailing list. If there are  
   numerous failures, there was a serious problem with your build.  

__7.__ The HTML version of the documentation is built during the compilation  
   process of Sage and resides in the directory:

         $SAGE_ROOT/local/share/doc/sage/html/

   >* __Optional__: If you want to build the PDF version (requires LaTeX)  
   of the documentation, run:  

         make doc-pdf

__8. Optional__: It is recommended that you install the optional GAP
   database by typing:

         ./sage -i database_gap

   This will download the package and install it. While you're at it,
   you might install other optional packages of interest to you: type
   "./sage --optional" to get a list.

__9. Optional__: It is recommended that you have both LaTeX and the  
   ImageMagick tools (e.g. the "convert" command) installed since some  
   plotting functionality benefits from it.  

__10. Optional__: Read this if you are intending to run a Sage notebook  
    server for multiple users. For security (i.e., to run  
    "notebook(secure=True)") you want to access the server using the  
    HTTPS protocol. First, install OpenSSL and the OpenSSL development  
    headers on your system if they are not already installed. Then  
    install pyOpenSSL by building Sage and then typing  

         ./sage -i pyopenssl

   >* Note that this command requires internet access. Alternatively,  
   "make ssl" builds Sage and installs pyOpenSSL.


Problems
--------

* If you have problems building Sage, check the Sage Installation Guide,  
and also note the following. Each separate component of Sage is  
contained in an spkg; these are stored in spkg/standard/. As each one  
is built, a build log is stored in logs/pkgs/, so you can browse these  
to find error messages. If an spkg fails to build, the whole build  
process will stop soon after, so check the most recent log files  
first, or run  

         grep -li "^Error" logs/pkgs/*

from the top-level Sage directory to find log files with error  
messages in them.  Send (a small part of) the relevant log file to the  
sage-devel mailing list, making sure to include at least some of the  
error messages; probably someone there will have some helpful  
suggestions.  


Supported Compilers
-------------------

* Sage includes a GCC (_GNU Compiler Collection_) package. In order to  
build Sage, you need a C compiler which can build GCC and its  
prerequisites. gcc version 4.0.1 or later should probably work. On  
Solaris or OpenSolaris, building with the Sun compiler should also work.  

* The GCC package in Sage is not always installed. It is determined  
automatically whether it needs to be installed. You can override this  
by setting the environment variable **SAGE_INSTALL_GCC=yes** (to force  
installation of GCC) or **SAGE_INSTALL_GCC=no** (to disable installation of  
GCC). If you don't want to install GCC, you need to have recent  
versions of gcc, g++ and gfortran; moreover, the versions must be equal.  

* There are some known problems with old assemblers, in particular when  
building the ECM package. You should ensure that your assembler  
understands all instructions for your processor. On Linux, this means  
you need a recent version of binutils; on OS X you need a recent version  
of XCode.  


Directory Layout
----------------

__Simplified directory layout (only essential files/directories):__  
```
SAGE_ROOT                 Root directory (sage-x.y.z in Sage tarball)
├── build
|   ├── deps              Dependency information of packages
|   └── pkgs              Every package is a subdirectory here
|       ├── atlas
|       ...
|       └── zn_poly
├── COPYING.txt           Copyright information
├── local                 Compiled packages are installed here
|   ├── bin               Executables
|   ├── include           C/C++ headers
|   ├── lib               Shared libraries
|   ├── share             Databases, architecture-independent data
|   └── var
|       ├── sage          List of installed packages
|       └── tmp           Temporary files when building Sage
├── logs
|   ├── dochtml.log       Log of the documentation build
|   ├── install.log       Full install log
|   └── pkgs              Build logs of individual packages
|       ├── atlas-3.10.1.p7.log
|       ...
|       └── zn_poly-0.9.p11.log
├── Makefile              Running "make" uses this file
├── README.md             This file
├── sage                  Script to start Sage
├── src                   All of Sage source (not third-party packages)
|   ├── bin               Scripts that Sage uses internally
|   ├── doc               Sage documentation
|   └── sage              The Sage library source code
├── upstream              Source tarballs of packages
|   ├── atlas-3.10.1.tar.bz2
|   ...
|   └── zn_poly-0.9.tar.bz2
└── VERSION.txt
```
For more details, see:

    http://doc.sagemath.org/html/en/developer/coding_basics.html#files-and-directory-structure


Build System
------------

This is a brief summary of the Sage software distribution's build system.
There are two components to the full Sage system--the Sage Python library
and its associated user interfaces, and the larger software distribution of
Sage's main dependencies (for those dependencies not supplied by the user's
system).

Sage's Python library is built and installed using a `setup.py` script as is
standard for Python packages (Sage's `setup.py` is non-trivial, but not
unusual).

Most of the rest of the build system is concerned with building all of Sage's
dependencies in the correct order in relation to each other.  The dependencies
included by Sage are referred to as SPKGs (i.e. "Sage Packages") and are listed
under `build/pkgs`.

The main entrypoint to Sage's build system is the top-level `Makefile` at the
root of the source tree.  Unlike most normal projects that use autoconf (Sage
does as well, as described below), this `Makefile` is not generated.  Instead,
it contains a few high-level targets and targets related to bootstrapping the
system.  Nonetheless, we still run `make <target>` from the root of the source
tree--targets not explicitly defined in the top-level `Makefile` are passed
through to another Makefile under `build/make/Makefile`.

The latter `build/make/Makefile` *is* generated by an autoconf-generated
`configure` script, using the template in `build/make/Makefile.in`.  This
includes rules for building the Sage library itself (`make sagelib`), and for
building and installing each of Sage's dependencies (e.g. `make python2`).

Although it's possible to manually run Sage's `configure` script if one wants
to provide some customizations (e.g. it is possible to select which BLAS
implementation to use), the top-level `Makefile` will run `configure` for you,
in order to build `build/make/Makefile` since it's a prerequisite for most of
Sage's make targets.

The `configure` script itself, if it is not already built, can be generated by
running the `bootstrap` script.  The top-level `Makefile` also takes care of
this automatically.

To summarize, running a command like `make python2` at the top-level of the
source tree goes something like this:

1. `make python2`
2. run `./bootstrap` if `configure` does not exist
3. run `./configure` if `build/make/Makefile` doe not exist
4. `cd` into `build/make` and run the `install` script--this is little more
   than a front-end to running `make -f build/make/Makefile python2`, which
   sets some necessary environment variables and logs some information
5. `build/make/Makefile` contains the actual rule for building `python2`; this
   includes building all of `python2`'s dependencies first (and their
   dependencies, recursively); the actual package installation is performed
   with the `sage-spkg` program


Relocation
----------

You *should* be able to move the sage-x.y.z/ directory anywhere you  
want. If you copy the sage script or make a symbolic link to it, you  
should modify the script to reflect this (as instructed at the top of  
the script). It is best if the path to Sage does not have any spaces in  
it.  

For a system-wide installation, as root you can move the sage-x.y.z/  
directory to a system-wide directory. Afterwards, you need to start up  
Sage as root at least once prior to using the system-wide Sage as a  
normal user. See the Installation Guide for further information on  
performing a system-wide installation:

   http://doc.sagemath.org/html/en/installation/source.html#installation-in-a-multiuser-environment

If you find anything that doesn't work correctly after you moved the  
directory, please email the sage-support mailing list.  


Redistribution
--------------

Your local Sage install is almost exactly the same as any 'developer'  
install. You can make changes to documentation, source, etc., and very  
easily package the complete results up for redistribution just like we  
do.  

1. To make your own source tarball of Sage, type:  

         sage --sdist

   The result is placed in the directory "$SAGE_ROOT/dist/".

2. To make a binary distribution with your currently installed packages,
   visit https://github.com/sagemath/binary-pkg


Changes to Included Software
----------------------------

All software included with Sage is copyrighted by the respective authors  
and released under an open source license that is __GPL version 3 or  
later__ compatible. See the file COPYING.txt for more details.  

Sources are in unmodified (as far as possible) tarballs in the  
__"$SAGE_ROOT/upstream"__ directory. The remaining description, version  
information, patches, and build scripts are in the accompanying  
__"$SAGE_ROOT/build/pkgs/<packagename>"__ directory. This directory is  
part of the Sage git repository.  
=======
Git-Trac Integration
====================

About
-----

This module implements a "git trac" subcommand of the git suite that
interfaces with trac over XMLRPC.

Included is a one-page git cheat sheet for Sage/Git/Trac: 
http://github.com/sagemath/git-trac-command/raw/master/doc/git-cheat-sheet.pdf

Installation
------------

The easiest way to just try out the code in this repo is to source the
``enable.sh`` script, which will prepend it to your PATH. This enables
the git trac subcommand until you close that shell, so no permanent
change is made:

    $ git clone https://github.com/sagemath/git-trac-command.git
    $ source git-trac-command/enable.sh
    Prepending the git-trac command to your search PATH

To permanently install the code from this repo, clone it and run
``setup.py``:

    $ git clone https://github.com/sagemath/git-trac-command.git
    $ cd git-trac-command
    $ python setup.py install

Alternatively you can just symlink ``git-trac`` to anywhere in your path:

    $ git clone https://github.com/sagemath/git-trac-command.git
    $ cd git-trac-command
    $ ln -s `pwd`/git-trac ~/bin/

On a Mac, which may not have a default in the home directory for commands,
it may be easier to add this to your ``.profile``: 

    $ git clone https://github.com/sagemath/git-trac-command.git
    $ cd git-trac-command
    $ pico/vim/emacs $HOME/.profile
    <add a line like "export PATH=$PATH:$HOME/Downloads/git-trac-command/bin">

Usage
-----

* Print the trac ticket information using ``git trac print
  <ticket_number>``. 

      $ git trac print 12345
      ==============================================================================
      Trac #12345: Title of ticket 12345
      ...
      ==============================================================================

  Alternatively, you can pass a remote branch name, in which case trac
  is searched for a ticket whose (remote) "Branch:" field equals the
  branch name.  If that fails, the ticket number will be deduced from
  the branch name by scanning for a number. If you neither specify a
  ticket number or branch name, the local git branch name is used:

      $ git branch
      /u/user/description
      $ git trac print
      ==============================================================================
      Trac #nnnnn: Title
      <BLANKLINE>
      Description
      Status: Status                          Component: Component                
      ...
      Branch: u/user/description
      ==============================================================================


* Checkout 
  a remote branch:

      $ git trac checkout 12345

  Will automatically pick a local branch name ``t/12345/description``
  based on the remote branch name. If you want a particular local
  branch name, you can specify it manually:

      $ git trac checkout -b my_branch 12345


* Create a new ticket on trac, and a new local branch 
  corresponding to it:

      $ git trac create "This is the summary"

  This will automatically create a local branch name
  ``t/12345/this_is_the_summary``. You can specify it manually if you
  prefer with:
  
      $ git trac create -b my_branch "This is the summary"


* Pull (= fetch + merge) from the branch
  on a ticket:

      $ git trac pull 12345

  You can omit the ticket number, in which case the script will try to
  search for the ticket having the local branch name attached. If that
  fails, an attempt is made to deduce the ticket number from the local
  branch name.


* Push (upload) to the branch
  on a ticket, and set the trac "Branch:" field accordingly:

      $ git trac push 12345

  You can omit the ticket number, in which case the script will try to
  search for the ticket having the local branch name attached. If that
  fails, an attempt is made to deduce the ticket number from the local
  branch name.


* Log of the commits for a
  ticket:

      $ git trac log 12345
    

* Find the trac ticket for a 
  commit, either identified by its SHA1 or branch/tag name.

      $ git log --oneline -1 ee5e39e
      ee5e39e Allow default arguments in closures
      $ git trac find ee5e39e
      Commit has been merged by the release manager into your current branch.
      commit 44efa774c5f991ea5f160646515cfe8d3f738479
      Merge: 5fd5442 679310b
      Author: Release Manager <release@sagemath.org>
      Date:   Sat Dec 21 01:16:56 2013 +0000

          Trac #15447: implement evaluation of PARI closures

* Review tickets with minimal recompiling. This assumes that you are
  currently on the "develop" branch, that is, the latest beta. Just
  checking out an older ticket would most likely reset the Sage tree
  to an older version, so you would have to compile older versions of
  packages to make it work. Instead, you can create an anonymous
  ("detached HEAD") merge of the ticket and the develop branch::

      $ git trac try 12345

  This will only touch files that are really modified by the
  ticket. In particular, if only Python files are changed by the
  ticket (which is true for most tickets) then you just have to run
  `sage -b` to rebuild the Sage library. When you are finished
  reviewing, just checkout a named branch. For example::

      $ git checkout develop
     
  If you want to edit the ticket branch (that is, add additional
  commits) you cannot use `git trac try`. You must use `git trac
  checkout` to get the actual ticket branch as a starting point.


Too Long, Didn't Read
---------------------

To fix a bug, start with

    $ git trac create "Fix foo"
    
This will open the ticket and create a new local branch
``t/<number>/fix_foo``. Then edit Sage, followed by 

    $ git add <filename>
    $ git commit

Repeat edit/commit as necessary. When you are finished, run

    $ git trac push

It will take the ticket number out of the branch name, so you don't
have to specify it.

    
Configuration
-------------

The scripts assume that the trac remote repository is set up as the
remote ``trac`` in the local repo. That is, you should have the
following for the Sage git server:

    $ git remote add trac https://trac.sagemath.org/sage.git      # read-only
    $ git remote add trac ssh://git@trac.sagemath.org/sage.git   # read-write
    $ git remote -v
    trac	ssh://git@trac.sagemath.org/sage.git (fetch)
    trac	ssh://git@trac.sagemath.org/sage.git (push)

Trac username and password are stored in the local repo (the
DOT_GIT/config file):
  
    $ git trac config --user=Myself --pass=s3kr1t
    Trac xmlrpc URL:
        https://trac.sagemath.org/xmlrpc (anonymous)
        https://trac.sagemath.org/login/xmlrpc (authenticated)
    Username: Myself
    Password: s3kr1t

If you do not want to store your trac username/password on disk you
can temporarily override it with the environment variables
``TRAC_USERNAME`` and ``TRAC_PASSWORD``. These take precedence over
any other configuration.


Sage-Trac Specifics
-------------------

Some of the functionality depends on the special trac plugins (see
https://github.com/sagemath/sage_trac), namely:

* Searching for a trac ticket by branch name requires the
  ``trac_plugin_search_branch.py`` installed in trac and a custom trac
  field named "Branch:":

      $ git trac search --branch=u/vbraun/toric_bundle
      15328

* SSH public key management requires the ``sshkeys.py`` trac 
  plugin:

      $ git trac ssh-keys
      $ git trac ssh-keys --add=~/.ssh/id_rsa.pub
      This is not implemented yet


Release Management
------------------

The Sage release management scripts are in the `git-trac.releasemgr`
subdirectory. They are probably only useful to the Sage release
manager.


Testing and Python Compatibility
--------------------------------

* The git-trac command supports Python 2.6, 2.7, 3.3, and 3.4.
* Most recent [Travis CI](https://travis-ci.org/sagemath/git-trac-command) test:
  [![Build Status](https://travis-ci.org/sagemath/git-trac-command.svg?branch=master)](https://travis-ci.org/sagemath/git-trac-command)
>>>>>>> 01cb18f9432a41ff32ca5922bff8fbed3156d284
