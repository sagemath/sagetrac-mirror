sage_conf: Configuration module for the SageMath library (relocatable wheel version)
====================================================================================

Description
-----------

This package provides:

- a single Python module, ``sage_conf``, providing configuration information
  to the SageMath library at the time of its installation and at its runtime

- a console script ``sage-config``, for querying the variables of ``sage_conf``
  from the shell

- a sourcable shell script ``sage-env-config``, providing additional configuration
  information in the form of environment variables

This version of the package is suitable to be put as wheels on PyPI.

On building a wheel, it invokes ``sage_bootstrap`` to establish a
build and installation tree (``SAGE_ROOT``, ``SAGE_LOCAL``) in a
subdirectory of the directory ``/var/tmp/``, whose name is specific to
the version of the distribution and the version of Python in use.

The wheel distributes a copy of the prebuilt ``SAGE_ROOT`` and
``SAGE_LOCAL``.  Importing ``sage_conf`` (or using the installed
``sage-config`` script), makes sure that a symlink from the
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
in ``src/pkgs/sage_conf-relocatable/``.
