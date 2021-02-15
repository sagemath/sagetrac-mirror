sage_conf: Configuration module for the SageMath library (sdist version)
========================================================================

Description
-----------

This package provides:

- a single Python module, ``sage_conf``, providing configuration information
  to the SageMath library at the time of its installation and at its runtime

- a console script ``sage-config``, for querying the variables of ``sage_conf``
  from the shell

- a sourcable shell script ``sage-env-config``, providing additional configuration
  information in the form of environment variables

This version of the package is suitable to be put as a source distribution on PyPI.

On installation (or building a wheel), it invokes ``sage_bootstrap`` to establish
a build tree (``SAGE_ROOT``) and installation tree (``SAGE_LOCAL``) for
the SageMath distribution.  By default, it uses a subdirectory of ``$HOME/.sage``
that is specific to the version of the distribution and the version of Python in
use.  If several virtual environments install ``sage_conf``, they will share these
trees.


License
-------

GNU General Public License (GPL) v3 or later

Upstream Contact
----------------

https://www.sagemath.org

This package is included in the source code of the Sage distribution,
in ``build/pkgs/sage_conf/src-pypi/``.
