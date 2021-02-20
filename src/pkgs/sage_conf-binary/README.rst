sage_conf: Configuration module for the SageMath library (relocatable binary version)
=====================================================================================

Description
-----------

This package provides:

- a single Python module, ``sage_conf``, providing configuration information
  to the SageMath library at the time of its installation and at its runtime

- a console script ``sage-config``, for querying the variables of ``sage_conf``
  from the shell

- a sourcable shell script ``sage-env-config``, providing additional configuration
  information in the form of environment variables

This version of the package is suitable for shipping it with binary
distributions of Sage.

It implements a relocation mechanism that allows an unprivileged user
to install the distribution in any location.

The Sage distribution should be configured to use a
build and installation tree (``SAGE_ROOT``, ``SAGE_LOCAL``) in a
subdirectory of the directory ``/var/tmp/``, whose name is specific to
the version of the distribution and the version of Python in use,
such as ``/var/tmp/sage-9.3.beta7-cpython-39-darwin``.

Importing ``sage_conf`` (or using the installed
``sage-config`` script) makes sure that a symlink from the hardcoded
``/var/tmp`` location to the actual persistent installation location
is created.  As the relocated libraries and programs contain the
hardcoded path ``SAGE_LOCAL`` in various ways (including as rpaths),
this symlink is necessary for the prebuilt libraries and programs to
work.

``/var/tmp`` is a sticky directory on all Linux distributions
following the Filesystem Hierarchy Standard, as well as on macOS and
on Cygwin.  On multi-user systems, only one user can use a given
version of the distribution; other installation schemes are recommended
for systems with multiple Sage users.

License
-------

GNU General Public License (GPL) v3 or later

Upstream Contact
----------------

https://www.sagemath.org

This package is included in the source code of the Sage distribution,
in ``src/pkgs/sage_conf-binary/``.
