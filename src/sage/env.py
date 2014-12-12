"""
Sage Runtime Environment

AUTHORS:

- \R. Andrew Ohana (2012): Initial version.

- Jeroen Demeyer (2014-12-12): add functions for Cython compiler flags,
  see :trac:`17484`.

"""

########################################################################
#       Copyright (C) 2013 R. Andrew Ohana <andrew.ohana@gmail.com>
#       Copyright (C) 2014 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
########################################################################

import os, socket, site
import version

opj = os.path.join

# set default values for sage environment variables
# every variable can be overwritten by os.environ
SAGE_ENV = dict()

# Helper to build the SAGE_ENV dictionary
def _add_variable_or_fallback(key, fallback, force=False):
    """
    Set ``SAGE_ENV[key]``.

    If ``key`` is an environment variable, this is the
    value. Otherwise, the ``fallback`` is used.

    INPUT:

    - ``key`` -- string.

    - ``fallback`` -- anything.

    - ``force`` -- boolean (optional, default is ``False``). Whether
      to always use the fallback, regardless of environment variables.

    EXAMPLES::

        sage: import os, sage.env
        sage: sage.env.SAGE_ENV = dict()
        sage: os.environ['SAGE_FOO'] = 'foo'
        sage: sage.env._add_variable_or_fallback('SAGE_FOO', '---$SAGE_URL---')
        sage: sage.env.SAGE_FOO
        'foo'
        sage: sage.env.SAGE_ENV['SAGE_FOO']
        'foo'

    If the environment variable does not exist, the fallback is
    used. Previously-declared variables are replaced if they are
    prefixed with a dollar sign::

        sage: _ = os.environ.pop('SAGE_BAR', None)  # ensure that SAGE_BAR does not exist
        sage: sage.env._add_variable_or_fallback('SAGE_BAR', '---$SAGE_FOO---')
        sage: sage.env.SAGE_BAR
        '---foo---'
        sage: sage.env.SAGE_ENV['SAGE_BAR']
        '---foo---'
    """
    global SAGE_ENV
    try:
        import os
        value = os.environ[key]
    except KeyError:
        value = fallback
    if force:
        value = fallback
    if isinstance(value, basestring):
        for k,v in SAGE_ENV.iteritems():
            if isinstance(v, basestring):
                value = value.replace('$'+k, v)
    SAGE_ENV[key] = value
    globals()[key] = value

# system info
_add_variable_or_fallback('UNAME',           os.uname()[0])
_add_variable_or_fallback('HOSTNAME',        socket.gethostname())
_add_variable_or_fallback('LOCAL_IDENTIFIER','$HOSTNAME.%s'%os.getpid())

# bunch of sage directories and files
_add_variable_or_fallback('SAGE_ROOT',       None)
_add_variable_or_fallback('SAGE_LOCAL',      opj('$SAGE_ROOT', 'local'))
_add_variable_or_fallback('SAGE_ETC',        opj('$SAGE_LOCAL', 'etc'))
_add_variable_or_fallback('SAGE_INC',        opj('$SAGE_LOCAL', 'include'))
_add_variable_or_fallback('SAGE_SHARE',      opj('$SAGE_LOCAL', 'share'))

_add_variable_or_fallback('SAGE_SRC',        opj('$SAGE_ROOT', 'src'))
_add_variable_or_fallback('SITE_PACKAGES',   site.getsitepackages())
_add_variable_or_fallback('SAGE_LIB',        SITE_PACKAGES[0])

_add_variable_or_fallback('SAGE_EXTCODE',    opj('$SAGE_SHARE', 'sage', 'ext'))
_add_variable_or_fallback('SAGE_LOGS',       opj('$SAGE_ROOT', 'logs', 'pkgs'))
_add_variable_or_fallback('SAGE_SPKG_INST',  opj('$SAGE_LOCAL', 'var', 'lib', 'sage', 'installed'))
_add_variable_or_fallback('SAGE_DOC',        opj('$SAGE_SRC', 'doc'))
_add_variable_or_fallback('DOT_SAGE',        opj(os.environ.get('HOME','$SAGE_ROOT'), '.sage'))
_add_variable_or_fallback('SAGE_DOT_GIT',    opj('$SAGE_ROOT', '.git'))

# misc
_add_variable_or_fallback('SAGE_URL',                'http://sage.math.washington.edu/sage/')
_add_variable_or_fallback('REALM',                   'sage.math.washington.edu')
_add_variable_or_fallback('TRAC_SERVER_URI',         'https://trac.sagemath.org')
_add_variable_or_fallback('SAGE_REPO_AUTHENTICATED', 'ssh://git@trac.sagemath.org:2222/sage.git')
_add_variable_or_fallback('SAGE_REPO_ANONYMOUS',     'git://trac.sagemath.org/sage.git')
_add_variable_or_fallback('SAGE_VERSION',            version.version)
_add_variable_or_fallback('SAGE_DATE',               version.date)
_add_variable_or_fallback('SAGE_IMPORTALL',          'yes')

# post process
if ' ' in DOT_SAGE:
    if UNAME[:6] == 'CYGWIN':
        # on windows/cygwin it is typical for the home directory
        # to have a space in it.  Fortunately, users also have
        # write privileges to c:\cygwin\home, so we just put
        # .sage there.
        _add_variable_or_fallback('DOT_SAGE', "/home/.sage", force=True)
    else:
        print("Your home directory has a space in it.  This")
        print("will probably break some functionality of Sage.  E.g.,")
        print("the GAP interface will not work. A workaround")
        print("is to set the environment variable HOME to a")
        print("directory with no spaces that you have write")
        print("permissions to before you start sage.")

# things that need DOT_SAGE
_add_variable_or_fallback('PYTHON_EGG_CACHE',   opj('$DOT_SAGE', '.python-eggs'))
_add_variable_or_fallback('SAGE_STARTUP_FILE',  opj('$DOT_SAGE', 'init.sage'))


#########################################
# Compiler/linker flags for Cython code
#########################################

def get_include_dirs():
    """
    Return a list of include directories, used to search for
    dependencies and add to gcc -I<path>.

    OUTPUT: a list of include directories

    EXAMPLES::

        sage: from sage.env import get_include_dirs
        sage: get_include_dirs()
        ['.../local/include',
         '.../src',
         '.../src/c_lib/include',
         '.../src/sage/ext']
    """
    from os.path import join as opj
    return [SAGE_INC,
            SAGE_SRC,
            opj(SAGE_SRC, 'c_lib', 'include'),
            opj(SAGE_SRC, 'sage', 'ext')]

def get_compile_args(debug=False, warn=False):
    """
    Return additional compile flags which should be used for Cython
    extensions.

    INPUT:

    - ``debug`` -- (default: ``False``) enable debugging with ``gdb``

    - ``warn`` -- (default: ``False``) enable compiler warnings

    OUTPUT: a list of compiler flags

    EXAMPLES::

        sage: from sage.env import get_compile_args
        sage: '-ggdb' in get_compile_args()
        False
        sage: '-ggdb' in get_compile_args(debug=True)
        True
    """
    import subprocess
    from distutils import sysconfig

    # Manually add -fno-strict-aliasing, which is needed to compile Cython
    # and disappears from the default flags if the user has set CFLAGS.
    extra_compile_args = [ "-fno-strict-aliasing" ]

    # comment these four lines out to turn on warnings from gcc
    if not warn:
        if sysconfig.get_config_var('CC').startswith("gcc"):
            extra_compile_args.append('-w')

    if debug:
        extra_compile_args.append('-ggdb')

    # Work around GCC-4.8.0 bug which miscompiles some sig_on() statements,
    # as witnessed by a doctest in sage/libs/gap/element.pyx if the
    # compiler flag -Og is used. See also
    # * http://trac.sagemath.org/sage_trac/ticket/14460
    # * http://gcc.gnu.org/bugzilla/show_bug.cgi?id=56982
    if subprocess.call("""$CC --version | grep -i 'gcc.* 4[.]8' >/dev/null """, shell=True) == 0:
        extra_compile_args.append('-fno-tree-dominator-opts')

    return extra_compile_args

def get_link_args():
    """
    Extra linker flags (currently nothing).

    OUTPUT: a list of linker flags

    EXAMPLES::

        sage: from sage.env import get_link_args
        sage: get_link_args()
        []
    """
    return []


# delete temporary variables used for setting up sage.env
del opj, os, socket, version, site
