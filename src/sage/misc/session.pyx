r"""
Loading and saving sessions and listing all variables

EXAMPLES:

We reset the current session, then define a rational number ``2/3``, and
verify that it is listed as a newly defined variable::

    sage: reset()
    sage: w = 2/3; w
    2/3
    sage: show_identifiers()
    ['w']

We next save this session. We are using a file in ``SAGE_TMP``. We do this
*for testing* only --- please do not do this, when you want to save your
session permanently, since ``SAGE_TMP`` will be removed when leaving Sage!

::

    sage: save_session(os.path.join(SAGE_TMP, 'session'))

This saves a dictionary with ``w`` as one of the keys::

    sage: z = load(os.path.join(SAGE_TMP, 'session'))
    sage: z.keys()
    ['w']
    sage: z['w']
    2/3

Next we reset the session, verify this, and load the session back.::

    sage: reset()
    sage: show_identifiers()
    []
    sage: load_session(os.path.join(SAGE_TMP, 'session'))

Indeed ``w`` is now defined again.::

    sage: show_identifiers()
    ['w']
    sage: w
    2/3

It is not needed to clean up the file created in the above code, since it
resides in the directory ``SAGE_TMP``.

AUTHOR:

- William Stein

"""

#############################################################################
#       Copyright (C) 2007,2010 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################
from __future__ import print_function, absolute_import

# Standard python imports
from six.moves import cPickle
import os
import types

# We want the caller's locals, but locals() is emulated in Cython
from six.moves import builtins
cdef caller_locals = builtins.locals

# Sage imports
from .misc import embedded
from sage.structure.sage_object import load, save

# This module-scope variables is used to save the
# global state of the sage environment at the moment
# before the user starts typing or running code.

state_at_init = None

def init(state=None):
    """
    Initialize some dictionaries needed by the :func:`show_identifiers`,
    :func:`save_session`, and :func:`load_session` functions.

    INPUT:

    - ``state`` -- a dictionary or ``None``; if ``None`` the :func:`locals()`
      of the caller is used.

    EXAMPLES::

        sage: reset()
        sage: w = 10
        sage: show_identifiers()
        ['w']

    When we call :func:`init()` below it reinitializes the internal table, so
    the ``w`` we just defined doesn't count as a new identifier::

        sage: sage.misc.session.init()
        sage: show_identifiers()
        []
    """
    if state is None: state = caller_locals()  # use locals() by default
    global state_at_init
    # Make a *copy* of the state dict, since it is mutable
    state_at_init = dict(state)

def _is_new_var(x, v, hidden):
    """
    Return whether or not the variable named ``x`` with value ``v`` is
    considered newly defined in the current session.

    INPUT:

    - ``x`` -- string

    - ``v`` -- object

    - ``hidden`` -- bool; if ``True``, always return ``False`` on variables
      that start with ``_``)

    OUTPUT:

    A bool

    EXAMPLES:

    We reset the session, then check whether the builtin factor function
    is newly defined (it isn't)::

        sage: reset()
        sage: sage.misc.session._is_new_var('factor', factor, True)
        False

    We then redefine factor, and find that it is newly defined::

        sage: factor = 10
        sage: sage.misc.session._is_new_var('factor', factor, True)
        True

    We define a new variable ``'blue'``, and test::

        sage: blue = 10
        sage: sage.misc.session._is_new_var('blue', blue, True)
        True
        sage: sage.misc.session._is_new_var('_blue', blue, True)
        True
        sage: sage.misc.session._is_new_var('_blue', blue, False)
        False
    """
    # We ignore all _ variable names unless hidden is True
    if not hidden and x.startswith('_'):
        return False
    # If a variable names was not there at init time then it is
    # definitely new.
    if x not in state_at_init:
        return True
    # A variable could also be new even if it was there at init, say if
    # its value changed.
    return x not in state_at_init or state_at_init[x] is not v

def show_identifiers(hidden=False):
    r"""
    Returns a list of all variable names that have been defined during
    this session.  By default, this returns only those identifiers
    that don't start with an underscore.

    INPUT:

    - ``hidden`` -- bool (Default: ``False``); If ``True``, also return
      identifiers that start with an underscore.

    OUTPUT:

    A list of variable names

    EXAMPLES:

    We reset the state of all variables, and see that none are defined::

        sage: reset()
        sage: show_identifiers()
        []

    We then define two variables, one which overwrites the default factor
    function; both are shown by :func:`show_identifiers()`::

        sage: a = 10
        sage: factor = 20
        sage: show_identifiers()
        ['a', 'factor']

    To get the actual value of a variable from the list, use the
    :func:`globals()` function.::

        sage: globals()['factor']
        20

    By default :func:`show_identifiers()` only returns variables that
    don't start with an underscore.  There is an option hidden that
    allows one to list those as well::

        sage: _hello = 10
        sage: show_identifiers()
        ['a', 'factor']
        sage: '_hello' in show_identifiers(hidden=True)
        True

    Many of the hidden variables are part of the IPython command history, at
    least in command line mode.::

        sage: show_identifiers(hidden=True)        # random output
        ['__', '_i', '_6', '_4', '_3', '_1', '_ii', '__doc__', '__builtins__', '___', '_9', '__name__', '_', 'a', '_i12', '_i14', 'factor', '__file__', '_hello', '_i13', '_i11', '_i10', '_i15', '_i5', '_13', '_10', '_iii', '_i9', '_i8', '_i7', '_i6', '_i4', '_i3', '_i2', '_i1', '_init_cmdline', '_14']
    """
    state = caller_locals()
    return [x for x, v in state.iteritems() if _is_new_var(x, v, hidden)]

def save_session(name='sage_session', verbose=False):
    r"""
    Save all variables that can be saved to the given filename.  The
    variables will be saved to a dictionary, which can be loaded using
    ``load(name)`` or :func:`load_session`.

    .. NOTE::

        1. Function and anything else that can't be pickled is not
           saved. This failure is silent unless you set
           ``verbose=True``.

        2. In the Sage notebook the session is saved both to the current
           working cell and to the ``DATA`` directory.

        3. One can still make sessions that can't be reloaded.  E.g., define
           a class with::

               class Foo: pass

          and make an instance with::

               f = Foo()

          Then :func:`save_session` followed by ``quit`` and
          :func:`load_session` fails. I doubt there is any good way to
          deal with this. Fortunately, one can simply re-evaluate the
          code to define ``Foo``, and suddenly :func:`load_session`
          works fine.

    INPUT:

        - ``name`` -- string (default: 'sage_session') name of ``sobj``
          to save the session to.

        - ``verbose`` -- bool (default: ``False``) if ``True``, print
          info about why certain variables can't be saved.

    OUTPUT:

        - Creates a file and returns silently.

    EXAMPLES:

    For testing, we use a temporary file that will be removed as soon
    as Sage is left. Of course, for permanently saving your session,
    you should choose a permanent file.

    ::

        sage: a = 5
        sage: tmp_f = tmp_filename()
        sage: save_session(tmp_f)
        sage: del a
        sage: load_session(tmp_f)
        sage: print(a)
        5

    We illustrate what happens when one of the variables is a function::

        sage: f = lambda x : x^2
        sage: save_session(tmp_f)
        sage: save_session(tmp_f, verbose=True)
        Saving...
        Not saving f: f is a function, method, class or type
        ...

    Something similar happens for cython-defined functions::

        sage: g = cython_lambda('double x', 'x*x + 1.5')
        sage: save_session(tmp_f, verbose=True)
        Not saving g: g is a function, method, class or type
        ...
    """
    state = caller_locals()
    # This dict D will contain the session -- as a dict -- that we will save to disk.
    D = {}
    # We iterate only over the new variables that were defined in this
    # session, since those are the only ones we will save.
    for k in show_identifiers(hidden = True):
        try:
            x = state[k]
            if isinstance(x, (types.FunctionType, types.BuiltinFunctionType, types.BuiltinMethodType, type)):
                raise TypeError('{} is a function, method, class or type'.format(k))

            # We attempt to pickle *and* unpickle every variable to
            # make *certain* that we can pickled D at the end below.
            # This seems wasteful, but it guarantees (I hope) that
            # D itself can be pickled and unpickled (assuming something
            # doesn't change in the Sage library itself).  Otherwise,
            # we could easily pickle whole sessions but get something
            # not at all useful.
            _ = cPickle.loads(cPickle.dumps(x, protocol=2))
            if verbose:
                print("Saving %s" % k)
            D[k] = x
        except Exception as msg:
            if verbose:
                print("Not saving {}: {}".format(k, msg))
            pass
    save(D, name)
    if embedded():
        # Also save D to the data directory if we're using the notebook.
        # This is broken for now. Simply print some information to the user
        # if the user does not save it in the DATA directory.
        # save(D, '../../data/' + name)
        if name.find('.sagenb/') <= 0 or name.find('/data/') <= 0:
            print("To store the session in a common directory that the "
                  "entire worksheet can access, save it using the command:\n"
                  "save_session(DATA + '{0}')\n"
                  "You can later load it by running in any cell:\n"
                  "load_session(DATA + '{0}')".format(name.rsplit('/', 1)[-1]))


def load_session(name='sage_session', verbose=False):
    r"""
    Load a saved session.

    This merges in all variables from a previously saved session.  It
    does not clear out the variables in the current sessions, unless
    they are overwritten.  You can thus merge multiple sessions, and
    don't necessarily loose all your current work when you use this
    command.

    .. NOTE::

        In the Sage notebook the session name is searched for both
        in the current working cell and the ``DATA`` directory.

    EXAMPLES::

        sage: a = 5
        sage: f = lambda x: x^2

    For testing, we use a temporary file, that will be removed as soon
    as Sage is left. Of course, for permanently saving your session,
    you should choose a permanent file.

    ::

        sage: tmp_f = tmp_filename()
        sage: save_session(tmp_f)
        sage: del a; del f
        sage: load_session(tmp_f)
        sage: print(a)
        5

    Note that ``f`` does not come back, since it is a function, hence
    couldn't be saved::

        sage: print(f)
        Traceback (most recent call last):
        ...
        NameError: name 'f' is not defined
    """
    state = caller_locals()

    if embedded():
        if not os.path.exists(name):
            nm = '../../data/' + name
            if not nm.endswith('.sobj'): nm += '.sobj'
            if os.path.exists(nm):
                name = nm
    D = load(name)
    for k, x in D.items():
        state[k] = x

