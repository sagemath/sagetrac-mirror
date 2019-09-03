r"""
Sage Runtime Environment

AUTHORS:

- \R. Andrew Ohana (2012): Initial version.

Verify that Sage can be started without any ``SAGE_`` environment
variables::

    sage: env = {k:v for (k,v) in os.environ.items() if not k.startswith("SAGE_")}
    sage: import subprocess
    sage: cmd = "from sage.all import SAGE_ROOT; print(SAGE_ROOT)"
    sage: res = subprocess.call([sys.executable, "-c", cmd], env=env)  # long time
    None
"""

# ****************************************************************************
#       Copyright (C) 2013 R. Andrew Ohana <andrew.ohana@gmail.com>
#       Copyright (C) 2019 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import absolute_import

import sage
import glob
import os
import socket
import sys
import sysconfig
from . import version


# All variables set by var() appear in this SAGE_ENV dict and also
# appear as module global (contained in __all__).
SAGE_ENV = dict()
__all__ = ['sage_include_directories', 'cython_aliases']


def join(*args):
    """
    Join paths like ``os.path.join`` except that the result is ``None``
    if any of the components is ``None``.

    EXAMPLES::

        sage: from sage.env import join
        sage: print(join("hello", "world"))
        hello/world
        sage: print(join("hello", None))
        None
    """
    if any(a is None for a in args):
        return None
    return os.path.join(*args)


def var(key, *fallbacks, **kwds):
    """
    Set ``SAGE_ENV[key]``.

    If ``key`` is an environment variable, this is the value.
    Otherwise, the ``fallbacks`` are tried until one is found which
    is not ``None``. If the environment variable is not set and all
    fallbacks are ``None``, then the final value is ``None``.

    INPUT:

    - ``key`` -- string.

    - ``fallbacks`` -- tuple containing ``str`` or ``None`` values.

    - ``force`` -- boolean (optional, default is ``False``). If
      ``True``, skip the environment variable and only use the
      fallbacks.

    EXAMPLES::

        sage: import os, sage.env
        sage: sage.env.SAGE_ENV = dict()
        sage: os.environ['SAGE_FOO'] = 'foo'
        sage: sage.env.var('SAGE_FOO', 'unused')
        sage: sage.env.SAGE_FOO
        'foo'
        sage: sage.env.SAGE_ENV['SAGE_FOO']
        'foo'

    If the environment variable does not exist, the fallbacks (if any)
    are used. In most typical uses, there is exactly one fallback::

        sage: _ = os.environ.pop('SAGE_BAR', None)  # ensure that SAGE_BAR does not exist
        sage: sage.env.var('SAGE_BAR', 'bar')
        sage: sage.env.SAGE_BAR
        'bar'
        sage: sage.env.SAGE_ENV['SAGE_BAR']
        'bar'

    Test multiple fallbacks::

        sage: sage.env.var('SAGE_BAR', None, 'yes', 'no')
        sage: sage.env.SAGE_BAR
        'yes'

    If all fallbacks are ``None``, the result is ``None``::

        sage: sage.env.var('SAGE_BAR')
        sage: print(sage.env.SAGE_BAR)
        None
        sage: sage.env.var('SAGE_BAR', None)
        sage: print(sage.env.SAGE_BAR)
        None

    Test the ``force`` keyword::

        sage: os.environ['SAGE_FOO'] = 'foo'
        sage: sage.env.var('SAGE_FOO', 'forced', force=True)
        sage: sage.env.SAGE_FOO
        'forced'
        sage: sage.env.var('SAGE_FOO', 'forced', force=False)
        sage: sage.env.SAGE_FOO
        'foo'
    """
    if kwds.get("force"):
        value = None
    else:
        value = os.environ.get(key)
    # Try all fallbacks in order as long as we don't have a value
    for f in fallbacks:
        if value is not None:
            break
        value = f
    SAGE_ENV[key] = value
    globals()[key] = value
    __all__.append(key)


# system info
var('UNAME',               os.uname()[0])
var('HOSTNAME',            socket.gethostname())
var('LOCAL_IDENTIFIER',    "{}.{}".format(HOSTNAME, os.getpid()))

# version info
var('SAGE_VERSION',        version.version)
var('SAGE_DATE',           version.date)
var('SAGE_VERSION_BANNER', version.banner)

# bunch of sage directories and files
var('SAGE_LOCAL',          os.path.abspath(sys.prefix))
var('SAGE_ETC',            join(SAGE_LOCAL, 'etc'))
var('SAGE_INC',            join(SAGE_LOCAL, 'include'))
var('SAGE_SHARE',          join(SAGE_LOCAL, 'share'))
var('SAGE_DOC',            join(SAGE_SHARE, 'doc', 'sage'))
var('SAGE_SPKG_INST',      join(SAGE_LOCAL, 'var', 'lib', 'sage', 'installed'))
var('SAGE_LIB',            os.path.dirname(os.path.dirname(sage.__file__)))

var('SAGE_ROOT')           # no fallback for SAGE_ROOT
var('SAGE_SRC',            join(SAGE_ROOT, 'src'), SAGE_LIB)
var('SAGE_DOC_SRC',        join(SAGE_ROOT, 'src', 'doc'), SAGE_DOC)
var('SAGE_PKGS',           join(SAGE_ROOT, 'build', 'pkgs'))
var('SAGE_ROOT_GIT',       join(SAGE_ROOT, '.git'))

var('DOT_SAGE',            join(os.environ.get('HOME'), '.sage'))
var('SAGE_STARTUP_FILE',   join(DOT_SAGE, 'init.sage'))

# installation directories for various packages
var('SAGE_EXTCODE',                  join(SAGE_SHARE, 'sage', 'ext'))
var('CONWAY_POLYNOMIALS_DATA_DIR',   join(SAGE_SHARE, 'conway_polynomials'))
var('GRAPHS_DATA_DIR',               join(SAGE_SHARE, 'graphs'))
var('ELLCURVE_DATA_DIR',             join(SAGE_SHARE, 'ellcurves'))
var('POLYTOPE_DATA_DIR',             join(SAGE_SHARE, 'reflexive_polytopes'))
var('GAP_ROOT_DIR',                  join(SAGE_SHARE, 'gap'))
var('THEBE_DIR',                     join(SAGE_SHARE, 'thebe'))
var('COMBINATORIAL_DESIGN_DATA_DIR', join(SAGE_SHARE, 'combinatorial_designs'))
var('CREMONA_MINI_DATA_DIR',         join(SAGE_SHARE, 'cremona'))
var('CREMONA_LARGE_DATA_DIR',        join(SAGE_SHARE, 'cremona'))
var('JMOL_DIR',                      join(SAGE_SHARE, 'jmol'))
var('JSMOL_DIR',                     join(SAGE_SHARE, 'jsmol'))
var('MATHJAX_DIR',                   join(SAGE_SHARE, 'mathjax'))
var('MTXLIB',                        join(SAGE_SHARE, 'meataxe'))
var('THREEJS_DIR',                   join(SAGE_SHARE, 'threejs'))
var('SINGULARPATH',                  join(SAGE_SHARE, 'singular'))
var('PPLPY_DOCS',                    join(SAGE_SHARE, 'doc', 'pplpy'))
var('MAXIMA_FAS')

# misc
var('SAGE_BANNER', '')
var('SAGE_IMPORTALL', 'yes')


def _get_shared_lib_filename(libname, *additional_libnames):
    """
    Return the full path to a shared library file installed in the standard
    location for the system within the ``LIBDIR`` prefix (or
    ``$SAGE_LOCAL/lib`` in the case of manual build of Sage).

    This can also be passed more than one library name (e.g. for cases where
    some library may have multiple names depending on the platform) in which
    case the first one found is returned.

    This supports most *NIX variants (in which ``lib<libname>.so`` is found
    under ``$SAGE_LOCAL/lib``), macOS (same, but with the ``.dylib``
    extension), and Cygwin (under ``$SAGE_LOCAL/bin/cyg<libname>.dll``,
    or ``$SAGE_LOCAL/bin/cyg<libname>-*.dll`` for versioned DLLs).

    For distributions like Debian that use a multiarch layout, we also try the
    multiarch lib paths (i.e. ``/usr/lib/<arch>/``).

    Returns ``None`` if the file does not exist.

    EXAMPLES::

        sage: import sys
        sage: from fnmatch import fnmatch
        sage: from sage.env import _get_shared_lib_filename
        sage: lib_filename = _get_shared_lib_filename("Singular",
        ....:                                         "singular-Singular")
        sage: if sys.platform == 'cygwin':
        ....:     pattern = "*/cygSingular-*.dll"
        ....: elif sys.platform == 'darwin':
        ....:     pattern = "*/libSingular.dylib"
        ....: else:
        ....:     pattern = "*/lib*Singular.so"
        sage: fnmatch(lib_filename, pattern)
        True
        sage: _get_shared_lib_filename("an_absurd_lib") is None
        True
    """

    for libname in (libname,) + additional_libnames:
        if sys.platform == 'cygwin':
            bindir = sysconfig.get_config_var('BINDIR')
            pats = ['cyg{}.dll'.format(libname), 'cyg{}-*.dll'.format(libname)]
            filenames = []
            for pat in pats:
                filenames += glob.glob(os.path.join(bindir, pat))

            # Note: This is not very robust, since if there are multi DLL
            # versions for the same library this just selects one more or less
            # at arbitrary.  However, practically speaking, on Cygwin, there
            # will only ever be one version
            if filenames:
                return filenames[-1]
        else:
            if sys.platform == 'darwin':
                ext = 'dylib'
            else:
                ext = 'so'

            libdirs = [sysconfig.get_config_var('LIBDIR')]
            multilib = sysconfig.get_config_var('MULTILIB')
            if multilib:
                libdirs.insert(0, os.path.join(libdirs[0], multilib))

            for libdir in libdirs:
                basename = 'lib{}.{}'.format(libname, ext)
                filename = os.path.join(libdir, basename)
                if os.path.exists(filename):
                    return filename

    # Just return None if no files were found
    return None


# locate singular shared object
# On Debian it's libsingular-Singular so try that as well
SINGULAR_SO = _get_shared_lib_filename('Singular', 'singular-Singular')
var('SINGULAR_SO', SINGULAR_SO)

# locate libgap shared object
GAP_SO= _get_shared_lib_filename('gap','')
var('GAP_SO', GAP_SO)

# post process
if ' ' in DOT_SAGE:
    if UNAME[:6] == 'CYGWIN':
        # on windows/cygwin it is typical for the home directory
        # to have a space in it.  Fortunately, users also have
        # write privileges to c:\cygwin\home, so we just put
        # .sage there.
        var('DOT_SAGE', "/home/.sage", force=True)
    else:
        print("Your home directory has a space in it.  This")
        print("will probably break some functionality of Sage.  E.g.,")
        print("the GAP interface will not work. A workaround")
        print("is to set the environment variable HOME to a")
        print("directory with no spaces that you have write")
        print("permissions to before you start sage.")


CYGWIN_VERSION = None
if UNAME[:6] == 'CYGWIN':
    import re
    _uname = os.uname()
    if len(_uname) >= 2:
        m = re.match(r'(\d+\.\d+\.\d+)\(.+\)', _uname[2])
        if m:
            CYGWIN_VERSION = tuple(map(int, m.group(1).split('.')))


def sage_include_directories(use_sources=False):
    """
    Return the list of include directories for compiling Sage extension modules.

    INPUT:

    -  ``use_sources`` -- (default: False) a boolean

    OUTPUT:

    a list of include directories to be used to compile sage code
    1. while building sage (use_sources='True')
    2. while using sage (use_sources='False')

    EXAMPLES:

    Expected output while using Sage::

        sage: import sage.env
        sage: sage.env.sage_include_directories()
        ['.../include',
        '.../python.../site-packages/sage/ext',
        '.../include/python...',
        '.../python.../numpy/core/include']

    To check that C/C++ files are correctly found, we verify that we can
    always find the include file ``sage/cpython/cython_metaclass.h``,
    with both values for ``use_sources``::

        sage: file = os.path.join("sage", "cpython", "cython_metaclass.h")
        sage: dirs = sage.env.sage_include_directories(use_sources=True)
        sage: any(os.path.isfile(os.path.join(d, file)) for d in dirs)
        True
        sage: dirs = sage.env.sage_include_directories(use_sources=False)
        sage: any(os.path.isfile(os.path.join(d, file)) for d in dirs)
        True
    """
    import numpy
    import distutils.sysconfig

    TOP = SAGE_SRC if use_sources else SAGE_LIB

    return [SAGE_INC,
            TOP,
            os.path.join(TOP, 'sage', 'ext'),
            distutils.sysconfig.get_python_inc(),
            numpy.get_include()]


def cython_aliases():
    """
    Return the aliases for compiling Cython code. These aliases are
    macros which can occur in ``# distutils`` headers.

    EXAMPLES::

        sage: from sage.env import cython_aliases
        sage: cython_aliases()
        {...}
        sage: sorted(cython_aliases().keys())
        ['ARB_LIBRARY',
         'FFLASFFPACK_CFLAGS',
         'FFLASFFPACK_INCDIR',
         'FFLASFFPACK_LIBDIR',
         'FFLASFFPACK_LIBRARIES',
         'GIVARO_CFLAGS',
         'GIVARO_INCDIR',
         'GIVARO_LIBDIR',
         'GIVARO_LIBRARIES',
         'GSL_CFLAGS',
         'GSL_INCDIR',
         'GSL_LIBDIR',
         'GSL_LIBRARIES',
         'LINBOX_CFLAGS',
         'LINBOX_INCDIR',
         'LINBOX_LIBDIR',
         'LINBOX_LIBRARIES',
         'SINGULAR_CFLAGS',
         'SINGULAR_INCDIR',
         'SINGULAR_LIBDIR',
         'SINGULAR_LIBRARIES']
    """
    import pkgconfig

    aliases = {}

    for lib in ['fflas-ffpack', 'givaro', 'gsl', 'linbox', 'Singular']:
        var = lib.upper().replace("-", "") + "_"
        aliases[var + "CFLAGS"] = pkgconfig.cflags(lib).split()
        pc = pkgconfig.parse(lib)
        # INCDIR should be redundant because the -I options are also
        # passed in CFLAGS
        aliases[var + "INCDIR"] = pc['include_dirs']
        aliases[var + "LIBDIR"] = pc['library_dirs']
        aliases[var + "LIBRARIES"] = pc['libraries']

    # LinBox needs special care because it actually requires C++11 with
    # GNU extensions: -std=c++11 does not work, you need -std=gnu++11
    # (this is true at least with GCC 7.2.0).
    #
    # Further, note that LinBox does not add any C++11 flag in its .pc
    # file (possibly because of confusion between CFLAGS and CXXFLAGS?).
    # This is not a problem in practice since LinBox depends on
    # fflas-ffpack and fflas-ffpack does add such a C++11 flag.
    aliases["LINBOX_CFLAGS"].append("-std=gnu++11")
    aliases["ARB_LIBRARY"] = os.environ.get('SAGE_ARB_LIBRARY', 'arb')
    return aliases
