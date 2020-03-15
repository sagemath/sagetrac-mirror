r"""
Listing Sage packages

This module can be used to see which Sage packages are installed
and which packages are available for installation.

For more information about creating Sage packages, see the "Packaging
Third-Party Code" section of the Sage Developer's Guide.

Actually installing the packages should be done via the command
line, using the following commands:

- ``sage -i PACKAGE_NAME`` -- install the given package

- ``sage -f PACKAGE_NAME`` -- re-install the given package, even if it
  was already installed

To list the packages available, either use in a terminal one of ``sage
-standard``, ``sage -optional`` or ``sage -experimental``. Or the following
command inside Sage::

    sage: from sage.misc.package import list_packages
    sage: pkgs = list_packages(local=True)  # optional - build
    sage: sorted(pkgs.keys())  # optional - build, random
    ['4ti2',
     'alabaster',
     'arb',
     ...
     'zlib',
     'zn_poly',
     'zope_interface']

Functions
---------
"""

# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import print_function

import sage.env

import json
import os
import subprocess
import sys
try:
    # Python 3.3+
    from urllib.request import urlopen
    from urllib.error import URLError
except ImportError:
    # Python 2.7
    from urllib2 import urlopen, URLError

DEFAULT_PYPI = 'https://pypi.python.org/pypi'

def pkgname_split(name):
    r"""
    Split a pkgname into a list of strings, 'name, version'.

    For some packages, the version string might be empty.

    EXAMPLES::

        sage: from sage.misc.package import pkgname_split
        sage: pkgname_split('hello_world-1.2')
        ['hello_world', '1.2']
    """
    return (name.split('-',1) + [''])[:2]

def pip_remote_version(pkg, pypi_url=DEFAULT_PYPI, ignore_URLError=False):
    r"""
    Return the version of this pip package available on PyPI.

    INPUT:

    - ``pkg`` -- the package

    - ``pypi_url`` -- (string, default: standard PyPI url) an optional Python
      package repository to use

    - ``ignore_URLError`` -- (default: ``False``) if set to ``True`` then no
      error is raised if the connection fails and the function returns ``None``

    EXAMPLES:

    The following test does fail if there is no TLS support (see e.g.
    :trac:`19213`)::

        sage: from sage.misc.package import pip_remote_version
        sage: pip_remote_version('beautifulsoup4') # optional - internet # not tested
        u'...'

    These tests are reliable since the tested package does not exist::

        sage: nap = 'hey_this_is_NOT_a_python_package'
        sage: pypi = 'http://this.is.not.pypi.com/'
        sage: pip_remote_version(nap, pypi_url=pypi, ignore_URLError=True) # optional - internet
        doctest:...: UserWarning: failed to fetch the version of
        pkg='hey_this_is_NOT_a_python_package' at http://this.is.not.pypi.com/
        sage: pip_remote_version(nap, pypi_url=pypi, ignore_URLError=False) # optional - internet
        Traceback (most recent call last):
        ...
        HTTPError: HTTP Error 404: Not Found
    """
    url = '{pypi_url}/{pkg}/json'.format(pypi_url=pypi_url, pkg=pkg)

    try:
        f = urlopen(url)
        text = f.read()
        f.close()
    except URLError:
        if ignore_URLError:
            import warnings
            warnings.warn("failed to fetch the version of pkg={!r} at {}".format(pkg, pypi_url))
            return
        else:
            raise

    info = json.loads(text)
    stable_releases = [v for v in info['releases'] if 'a' not in v and 'b' not in v]
    return max(stable_releases)

def pip_installed_packages():
    r"""
    Return a dictionary `name->version` of installed pip packages.

    This command returns *all* pip-installed packages. Not only Sage packages.

    EXAMPLES::

        sage: from sage.misc.package import pip_installed_packages
        sage: d = pip_installed_packages()  # optional - build
        sage: 'scipy' in d  # optional - build
        True
        sage: d['scipy']  # optional - build
        u'...'
        sage: d['beautifulsoup4']   # optional - build beautifulsoup4
        u'...'
    """
    with open(os.devnull, 'w')  as devnull:
        proc = subprocess.Popen(
            [sys.executable, "-m", "pip", "list", "--no-index", "--format", "json"],
            stdout=subprocess.PIPE,
            stderr=devnull,
        )
        stdout = proc.communicate()[0].decode()
        return {package['name'].lower():package['version']
                for package in json.loads(stdout)}

def list_packages(*pkg_types, **opts):
    r"""
    Return a dictionary of information about each package.

    The keys are package names and values are dictionaries with the following
    keys:

    - ``'type'``: either ``'base``, ``'standard'``, ``'optional'``, or ``'experimental'``
    - ``'source'``: either ``'normal', ``'pip'``, or ``'script'``
    - ``'installed'``: boolean
    - ``'installed_version'``: ``None`` or a string
    - ``'remote_version'``: string

    INPUT:

    - ``pkg_types`` -- (optional) a sublist of ``'base``, ``'standard'``, ``'optional'``,
      or ``'experimental'``.  If provided, list only the packages with the
      given type(s), otherwise list all packages.

    - ``pkg_sources`` -- (optional) a sublist of ``'normal', ``'pip'``, or ``'script'``.
      If provided, list only the packages with the given source(s), otherwise list all
      packages.

    - ``local`` -- (optional, default: ``False``) if set to ``True``, then do not
      consult remote (PyPI) repositories for package versions (only applicable for
      ``'pip'`` type)

    - ``exclude_pip`` -- (optional, default: ``False``) if set to ``True``, then
      pip packages are not considered.

    - ``ignore_URLError`` -- (default: ``False``) if set to ``True``, then
      connection errors will be ignored

    EXAMPLES::

        sage: from sage.misc.package import list_packages
        sage: L = list_packages('standard')  # optional - build
        sage: sorted(L.keys())  # optional - build, random
        ['alabaster',
         'arb',
         'babel',
         ...
         'zn_poly',
         'zope_interface']
        sage: L['ppl']  # optional - build
        {'installed': True,
         'installed_version': '...',
         'remote_version': '...',
         'type': 'standard'}

        sage: L = list_packages(pkg_sources=['pip'], local=True)  # optional - build
        sage: L['beautifulsoup4']                    # optional - build
        {'installed': ...,
         'installed_version': ...,
         'remote_version': None,
         'source': 'pip',
         'type': 'optional'}

        sage: L = list_packages(pkg_sources=['pip'])   # optional - build internet
        sage: L['beautifulsoup4']         # optional - build internet
        {'installed': ...,
         'installed_version': ...,
         'remote_version': u'...',
         'source': 'pip',
         'type': 'optional'}

    Check the option ``exclude_pip``::

        sage: [p for p, d in list_packages('optional', exclude_pip=True).items()  # optional - build
        ....:  if d['source'] == 'pip']
        []
    """
    if not pkg_types:
        pkg_types = ('base', 'standard', 'optional', 'experimental')
    elif any(pkg_type not in ('base', 'standard', 'optional', 'experimental') for pkg_type in pkg_types):
        raise ValueError("Each pkg_type must be one of 'base', 'standard', 'optional', 'experimental'")

    pkg_sources = opts.pop('pkg_sources',
                           ('normal', 'pip', 'script'))

    local = opts.pop('local', False)
    ignore_URLError = opts.pop('ignore_URLError', False)
    exclude_pip = opts.pop('exclude_pip', False)
    if opts:
        raise ValueError("{} are not valid options".format(sorted(opts)))

    installed = installed_packages(exclude_pip)

    pkgs = {}
    SAGE_PKGS = sage.env.SAGE_PKGS
    for p in os.listdir(SAGE_PKGS):
        try:
            f = open(os.path.join(SAGE_PKGS, p, "type"))
        except IOError:
            # Probably an empty directory => ignore
            continue

        with f:
            typ = f.read().strip()

        if typ not in pkg_types:
            continue

        if os.path.isfile(os.path.join(SAGE_PKGS, p, "requirements.txt")):
            src = 'pip'
        elif os.path.isfile(os.path.join(SAGE_PKGS, p, "checksums.ini")):
            src = 'normal'
        else:
            src = 'script'

        if src not in pkg_sources:
            continue

        pkg = {'name': p, 'type': typ, 'source': src, 'installed_version': installed.get(p)}
        pkg['installed'] = pkg['installed_version'] is not None

        if pkg['source'] == 'pip':
            if exclude_pip:
                continue
            if not local:
                pkg['remote_version'] = pip_remote_version(p, ignore_URLError=ignore_URLError)
            else:
                pkg['remote_version'] = None
        elif pkg['source'] == 'normal':
            # If package-version.txt does not exist, that is an error
            # in the build system => we just propagate the exception
            package_filename = os.path.join(SAGE_PKGS, p, "package-version.txt")
            with open(package_filename) as f:
                pkg['remote_version'] = f.read().strip()
            pkg['installed_version'] = installed.get(p)
        else:
            pkg['remote_version'] = 'none'

        pkgs[p] = pkg

    return pkgs


def installed_packages(exclude_pip=True):
    """
    Return a dictionary of all installed packages, with version numbers.

    INPUT:

    - ``exclude_pip`` -- (optional, default: ``True``) whether "pip" packages
      are excluded from the list

    EXAMPLES::

        sage: installed_packages()  # optional - build
        {...'brial': ...'pynac': ...}

    .. SEEALSO::

        :func:`sage.misc.package.list_packages`
    """
    installed = {}
    if not exclude_pip:
        installed.update(pip_installed_packages())
    # Sage packages should override pip packages (Trac #23997)
    installed.update(pkgname_split(pkgname)
                     for pkgname in os.listdir(sage.env.SAGE_SPKG_INST))
    return installed


def is_package_installed(package, exclude_pip=True):
    """
    Return whether (any version of) ``package`` is installed.

    INPUT:

    - ``package`` -- the name of the package

    - ``exclude_pip`` -- (optional, default: ``True``) whether to consider pip
      type packages


    EXAMPLES::

        sage: is_package_installed('gap')  # optional - build
        True

    Giving just the beginning of the package name is not good enough::

        sage: is_package_installed('matplotli')  # optional - build
        False

    Otherwise, installing "pillow" would cause this function to think
    that "pil" is installed, for example.

    Check that the option ``exclude_pip`` is turned on by default::

        sage: from sage.misc.package import list_packages
        sage: for pkg in list_packages(pkg_sources=('pip'), local=True):  # optional - build
        ....:     assert not is_package_installed(pkg), "pip package is installed: {}".format(pkg)

    .. NOTE::

        Do not use this function to check whether you can use a feature from an
        external library. This only checks whether something was installed with
        ``sage -i`` but it may have been installed by other means (for example
        if this copy of Sage has been installed as part of a distribution.)
        Use the framework provided by :mod:`sage.features` to check
        whether a library is installed and functional.
    """
    return any(p.split('-')[0] == package for p in installed_packages(exclude_pip))


def package_versions(package_type, local=False):
    r"""
    Return version information for each Sage package.

    INPUT:

    - ``package_type`` -- (string) one of ``"standard"``, ``"optional"`` or
      ``"experimental"``

    - ``local`` -- (boolean, default: ``False``) only query local data (no internet needed)

    For packages of the given type, return a dictionary whose entries
    are of the form ``'package': (installed, latest)``, where
    ``installed`` is the installed version (or ``None`` if not
    installed) and ``latest`` is the latest available version. If the
    package has a directory in ``SAGE_ROOT/build/pkgs/``, then
    ``latest`` is determined by the file ``package-version.txt`` in
    that directory.  If ``local`` is ``False``, then Sage's servers are
    queried for package information.

    .. SEEALSO:: :func:`sage.misc.package.list_packages`

    EXAMPLES::

        sage: std = package_versions('standard', local=True)  # optional - build
        sage: 'gap' in std  # optional - build
        True
        sage: std['zn_poly']  # optional - build, random
        ('0.9.p12', '0.9.p12')
    """
    return {pkg['name']: (pkg['installed_version'], pkg['remote_version']) for pkg in list_packages(package_type, local=local).values()}


def standard_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed standard packages that are available
    from the Sage repository.

    OUTPUT:

    - installed standard packages (as a list)

    - NOT installed standard packages (as a list)

    Run ``sage -i package_name`` from a shell to install a given
    package or ``sage -f package_name`` to re-install it.

    .. SEEALSO:: :func:`sage.misc.package.list_packages`

    EXAMPLES::

        sage: from sage.misc.package import standard_packages
        sage: installed, not_installed = standard_packages()  # optional - build
        sage: installed[0], installed[-1]  # optional - build
        ('alabaster', 'zope_interface')
    """
    pkgs = list_packages('standard', local=True).values()
    return (sorted(pkg['name'] for pkg in pkgs if pkg['installed']),
            sorted(pkg['name'] for pkg in pkgs if not pkg['installed']))


def optional_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed optional packages that are available
    from the Sage repository.

    OUTPUT:

    - installed optional packages (as a list)

    - NOT installed optional packages (as a list)

    Run ``sage -i package_name`` from a shell to install a given
    package or ``sage -f package_name`` to re-install it.

    .. SEEALSO:: :func:`sage.misc.package.list_packages`

    EXAMPLES::

        sage: from sage.misc.package import optional_packages
        sage: installed, not_installed = optional_packages()  # optional - build
        sage: 'beautifulsoup4' in installed+not_installed  # optional - build
        True

        sage: 'beautifulsoup4' in installed   # optional - build beautifulsoup4
        True
    """
    pkgs = list_packages('optional', local=True)
    pkgs = pkgs.values()
    return (sorted(pkg['name'] for pkg in pkgs if pkg['installed']),
            sorted(pkg['name'] for pkg in pkgs if not pkg['installed']))


def experimental_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed experimental packages that are available
    from the Sage repository.

    OUTPUT:

    - installed experimental packages (as a list)

    - NOT installed experimental packages (as a list)

    Run ``sage -i package_name`` from a shell to install a given
    package or ``sage -f package_name`` to re-install it.

    .. SEEALSO:: :func:`sage.misc.package.list_packages`

    EXAMPLES::

        sage: from sage.misc.package import experimental_packages
        sage: installed, not_installed = experimental_packages()  # optional - build
    """
    pkgs = list_packages('experimental', local=True).values()
    return (sorted(pkg['name'] for pkg in pkgs if pkg['installed']),
            sorted(pkg['name'] for pkg in pkgs if not pkg['installed']))

def package_manifest(package):
    """
    Return the manifest for ``package``.

    INPUT:

    - ``package`` -- package name

    The manifest is written in the file
    ``SAGE_SPKG_INST/package-VERSION``. It is a JSON file containing a
    dictionary with the package name, version, installation date, list
    of installed files, etc.

    EXAMPLES::

        sage: from sage.misc.package import package_manifest
        sage: sagetex_manifest = package_manifest('sagetex')  # optional - build
        sage: sagetex_manifest['package_name'] == 'sagetex'  # optional - build
        True
        sage: 'files' in sagetex_manifest  # optional - build
        True

    Test a nonexistent package::

        sage: package_manifest('dummy-package')  # optional - build
        Traceback (most recent call last):
        ...
        KeyError: 'dummy-package'
    """
    version = installed_packages()[package]
    stamp_file = os.path.join(sage.env.SAGE_SPKG_INST,
                              '{}-{}'.format(package, version))
    with open(stamp_file) as f:
        spkg_meta = json.load(f)
    return spkg_meta


class PackageNotFoundError(RuntimeError):
    """
    This class defines the exception that should be raised when a
    function, method, or class cannot detect a Sage package that it
    depends on.

    This exception should be raised with a single argument, namely
    the name of the package.

    When a ``PackageNotFoundError`` is raised, this means one of the
    following:

    - The required optional package is not installed.

    - The required optional package is installed, but the relevant
      interface to that package is unable to detect the package.

    EXAMPLES::

        sage: from sage.misc.package import PackageNotFoundError
        sage: raise PackageNotFoundError("my_package")
        Traceback (most recent call last):
        ...
        PackageNotFoundError: the package 'my_package' was not found. You can install it by running 'sage -i my_package' in a shell
    """
    def __str__(self):
        """
        Return the actual error message.

        EXAMPLES::

            sage: from sage.misc.package import PackageNotFoundError
            sage: str(PackageNotFoundError("my_package"))
            "the package 'my_package' was not found. You can install it by running 'sage -i my_package' in a shell"
        """
        return ("the package {0!r} was not found. "
            "You can install it by running 'sage -i {0}' in a shell"
            .format(self.args[0]))
