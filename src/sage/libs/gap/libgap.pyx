"""
Library Interface to GAP

This module implements a fast C library interface to GAP.  To use it, you
simply call ``libgap`` and use it to convert Sage objects into GAP objects.

EXAMPLES::

    sage: a = libgap(10)
    sage: a
    10
    sage: type(a)
    <class 'gappy.gapobj.GapInteger'>
    sage: a*a
    100
    sage: timeit('a*a')   # random output
    625 loops, best of 3: 898 ns per loop

Compared to the expect interface this is >1000 times faster::

    sage: b = gap('10')
    sage: timeit('b*b')   # random output; long time
    125 loops, best of 3: 2.05 ms per loop

If you want to evaluate GAP commands, use the :meth:`Gap.eval` method::

    sage: libgap.eval('List([1..10], i->i^2)')
    [ 1, 4, 9, 16, 25, 36, 49, 64, 81, 100 ]

not to be confused with the ``libgap`` call, which converts Sage
objects to GAP objects, for example strings to strings::

    sage: libgap('List([1..10], i->i^2)')
    "List([1..10], i->i^2)"
    sage: type(_)
    <class 'gappy.gapobj.GapString'>

You can usually use the ``.sage()`` method to convert the resulting GAP element
back to its Sage equivalent::

    sage: a.sage()
    10
    sage: type(_)
    <type 'sage.rings.integer.Integer'>

    sage: libgap.eval('5/3 + 7*E(3)').sage()
    7*zeta3 + 5/3

    sage: gens_of_group = libgap.AlternatingGroup(4).GeneratorsOfGroup()
    sage: generators = gens_of_group.sage()
    sage: generators   # a Sage list of Sage permutations!
    [[2, 3, 1], [1, 3, 4, 2]]
    sage: PermutationGroup(generators).cardinality()   # computed in Sage
    12
    sage: libgap.AlternatingGroup(4).Size()            # computed in GAP
    12

We can also specify which group in Sage the permutations should
consider themselves as elements of when converted to Sage::

    sage: A4 = groups.permutation.Alternating(4)
    sage: generators = gens_of_group.sage(parent=A4); generators
    [(1,2,3), (2,3,4)]
    sage: all(gen.parent() is A4 for gen in generators)
    True

So far, the following GAP data types can be directly converted to the
corresponding Sage datatype:

#. GAP booleans ``true`` / ``false`` to Sage booleans ``True`` /
   ``False``. The third GAP boolean value ``fail`` raises a
   ``ValueError``.

#. GAP integers to Sage integers.

#. GAP rational numbers to Sage rational numbers.

#. GAP cyclotomic numbers to Sage cyclotomic numbers.

#. GAP permutations to Sage permutations.

#. The GAP containers ``List`` and ``rec`` are converted to Sage containers
   ``list`` and ``dict``.  Furthermore, the ``.sage()`` method is applied
   recursively to the entries.

Special support is available for the GAP container classes. GAP lists
can be used as follows::

    sage: lst = libgap([1,5,7]);  lst
    [ 1, 5, 7 ]
    sage: type(lst)
    <class 'gappy.gapobj.GapList'>
    sage: len(lst)
    3
    sage: lst[0]
    1
    sage: [ x^2 for x in lst ]
    [1, 25, 49]
    sage: type(_[0])
    <class 'gappy.gapobj.GapInteger'>

Note that you can access the elements of GAP ``List`` objects as you would
expect from Python (with indexing starting at 0), but the elements are still of
type :class:`~gappy.gapobj.GapObj`. The other GAP container type are records,
which are similar to Python dictionaries. You can construct them directly from
Python dictionaries::

    sage: libgap({'a':123, 'b':456})
    rec( a := 123, b := 456 )

Or get them as results of computations::

    sage: rec = libgap.eval('rec(a:=123, b:=456, Sym3:=SymmetricGroup(3))')
    sage: rec['Sym3']
    Sym( [ 1 .. 3 ] )
    sage: dict(rec)
    {"Sym3": Sym( [ 1 .. 3 ] ), "a": 123, "b": 456}

The output is a Sage dictionary whose keys are GAP strings and whose Values are
instances of :meth:`~sage.gap.GapObj`. So, for example,
``rec['a']`` is not a Sage integer. To recursively convert the keys and entries
into Sage objects, you should use the
``GapObj.sage()`` method::

    sage: rec.sage()
    {'Sym3': NotImplementedError('cannot construct equivalent Sage object'...),
     'a': 123,
     'b': 456}

Now ``rec['a']`` is a Sage integer. We have not implemented the
conversion of the GAP symmetric group to the Sage symmetric group yet,
so you end up with a ``NotImplementedError`` exception object. The
exception is returned and not raised so that you can work with the
partial result.

While we don't directly support matrices yet, you can convert them to Gap List
of Lists. These lists are then easily converted into Sage using the recursive
expansion of the ``GapObj.sage`` method::

    sage: M = libgap.eval('BlockMatrix([[1,1,[[1, 2],[ 3, 4]]], [1,2,[[9,10],[11,12]]], [2,2,[[5, 6],[ 7, 8]]]],2,2)')
    sage: M
    <block matrix of dimensions (2*2)x(2*2)>
    sage: M.List()   # returns a GAP List of Lists
    [ [ 1, 2, 9, 10 ], [ 3, 4, 11, 12 ], [ 0, 0, 5, 6 ], [ 0, 0, 7, 8 ] ]
    sage: M.List().sage()   # returns a Sage list of lists
    [[1, 2, 9, 10], [3, 4, 11, 12], [0, 0, 5, 6], [0, 0, 7, 8]]
    sage: matrix(ZZ, _)
    [ 1  2  9 10]
    [ 3  4 11 12]
    [ 0  0  5  6]
    [ 0  0  7  8]


Using the GAP C library from Cython
===================================

.. TODO:: Expand the following text

   We are using the GAP API provided by the GAP project since
   GAP 4.10.

AUTHORS:

  - William Stein, Robert Miller (2009-06-23): first version
  - Volker Braun, Dmitrii Pasechnik, Ivan Andrus (2011-03-25, Sage Days 29):
    almost complete rewrite; first usable version.
  - Volker Braun (2012-08-28, GAP/Singular workshop): update to
    gap-4.5.5, make it ready for public consumption.
  - Dima Pasechnik (2018-09-18, GAP Days): started the port to native
    libgap API
"""

###############################################################################
#       Copyright (C) 2009, William Stein <wstein@gmail.com>
#       Copyright (C) 2012, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   https://www.gnu.org/licenses/
###############################################################################


##############################################################################
#
#  If you want to add support for converting some GAP datatype to its
#  Sage equivalent, you have two options. Either
#
#  1. add an if-clause to GapElement.sage(). This is the easiest
#     option and you should probably start with it first
#
#  2. Subclass GapElement to GapElement_Mydatatype. In that case you
#     need to write the derived class, a factory function to
#     instantiate it (GapElement cannot be instantiated by __init__),
#     and add an if-clause in make_GapElement. See GapElement_List,
#     for example. The advantage of this more complicated approach is
#     that you can then add extra methods to your data type. For
#     example, GapElement_List instances can access the individual
#     elements with the usual gapelement[i] syntax.
#
# TODO
#
# Not all GAP types have their own GapElement. We should wrap more
# stuff. Talk to me (Volker) if you want to work on that.
#
##############################################################################


import os

from gappy.core cimport Gap
from gappy.exceptions import GAPError
from gappy.gapobj cimport make_GapList
from gappy.gapobj import GapObj
from gappy.gap_includes cimport Obj
from gappy.gap_globals import common_gap_globals

from sage.interfaces.gap_workspace import prepare_workspace_dir
from sage.misc.cachefunc import cached_method
from sage.misc.randstate cimport current_randstate
from sage.misc.superseded import deprecation
from sage.misc.temporary_file import atomic_write
from sage.rings.all import ZZ
from sage.structure.parent cimport Parent
from sage.structure.sage_object import SageObject

from .saved_workspace import workspace as get_workspace
from .util import gap_root


############################################################################
### Gap  ###################################################################
############################################################################
# The libGap interpreter object Gap is the parent of the GapElements
# Provides a wrapper to gappy.core.Gap, which we can't subclass directly
# since it has an incompatible binary layout with Parent.
cdef class SageGap(Gap):
    """
    Subclasses `gappy.core.Gap` to provide some Sage-specific default
    behaviors.
    """

    def __init__(self):
        workspace, _ = get_workspace()
        Gap.__init__(self, gap_root=gap_root(), workspace=workspace,
                     autoload=True)

    cpdef initialize(self):
        initializing = Gap.initialize(self)

        if initializing:
            # These steps are only performed if SageGap has just been
            # initialized for the first time; if we don't check this then we'll
            # cause an infinite recursion

            # Import the converters module to register SageObject -> GapObj and
            # GapObj.sage() converters on the default libgap.
            from . import converters

            workspace, workspace_is_up_to_date = get_workspace()
            if self.workspace == os.path.normpath(workspace):
                # Save a new workspace if necessary
                if not workspace_is_up_to_date:
                    prepare_workspace_dir()
                    with atomic_write(workspace) as f:
                        self.SaveWorkspace(f.name)

        return initializing

    cpdef eval(self, gap_command):
        """
        Evaluate a gap command and wrap the result.

        INPUT:

        - ``gap_command`` -- a string containing a valid gap command
          without the trailing semicolon.

        OUTPUT:

        A :class:`~gappy.gapobj.GapObj`.

        EXAMPLES::

            sage: libgap.eval('0')
            0
            sage: libgap.eval('"string"')
            "string"
        """

        if hasattr(gap_command, '_libgap_init_'):
            gap_command = str(gap_command._libgap_init_())

        return Gap.eval(self, gap_command)

    def load_package(self, pkg):
        """
        If loading fails, raise a RuntimeError exception.

        TESTS::

            sage: libgap.load_package("chevie")
            Traceback (most recent call last):
            ...
            RuntimeError: Error loading GAP package chevie. You may want to
            install the gap_packages SPKG.
        """
        try:
            return super().load_package(pkg)
        except RuntimeError as exc:
            # Catch the exception from gappy and amend it with a Sage-specific
            # hint.
            if exc.args[0].startswith('Error loading GAP package'):
                raise RuntimeError(
                    exc.args[0] + ' You may want to install the gap_packages '
                    'SPKG.')
            else:
                raise exc

    @cached_method
    def function_factory(self, function):
        """
        Return a GAP function wrapper

        This is almost the same as calling ``libgap.eval(function)``, but
        faster and makes it obvious in your code that you are wrapping a
        function.

        INPUT:

        - ``function`` -- string. The name of a GAP function or a GAP
          function definition.

        OUTPUT:

        A function wrapper :class:`~gappy.gapobj.GapFunction` for the GAP
        function. Calling it from Sage is equivalent to calling the wrapped
        function from GAP.

        EXAMPLES::

            sage: libgap.function_factory('Print')
            doctest:warning
            ...
            DeprecationWarning: use libgap.Print instead...
            <GAP function "Print">
        """
        if ';' not in function:
            deprecation_msg = (
                f'use libgap.{function} instead; you can reduce overhead by '
                f'assigning {function} = libgap.{function} if the '
                f'function will be used repeatedly in your code'
            )
        else:
            # This method has also been used to create and cache anonymous
            # functions in which case there should be a ';' at least
            # somewhere in the string
            deprecation_msg = f'use the libgap.gap_function decorator instead'
        # TODO: Set correct issue number
        deprecation(31297, deprecation_msg)
        return self.eval(function)

    def unset_global(self, variable):
        """
        Remove a GAP global variable

        INPUT:

        - ``variable`` -- string. The variable name.

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: libgap.get_global('FooBar')
            1
            sage: libgap.unset_global('FooBar')
            sage: libgap.get_global('FooBar')
            Traceback (most recent call last):
            ...
            GAPError: no value bound to FooBar
        """
        if self.IsReadOnlyGlobal(variable):
            deprecation(31297,
                f'{variable} is a read-only global; unsetting of read-only '
                f'globals is deprecated and will be removed in a future '
                f'version; if you need to unset a read-only global manually '
                f'call libgap.MakeReadWriteGlobal({variable!r}) first')
            self.MakeReadWriteGlobal(variable)

        return super().unset_global(variable)

    cpdef get_global(self, variable):
        """
        Get a GAP global variable

        INPUT:

        - ``variable`` -- string. The variable name.

        OUTPUT:

        A :class:`~gappy.gapobj.GapObj` wrapping the GAP output. A
        ``ValueError`` is raised if there is no such variable in GAP.

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: libgap.get_global('FooBar')
            1
            sage: libgap.unset_global('FooBar')
            sage: libgap.get_global('FooBar')
            Traceback (most recent call last):
            ...
            GAPError: no value bound to FooBar
        """
        # TODO: Should this return None like gappy does or still raise a GAP
        # error?  Perhaps we could raise a DeprecationWarning on the GAPError
        # case?
        val = Gap.get_global(self, variable)
        if val is None:
            raise GAPError(f'no value bound to {variable}')
        return val

    def set_seed(self, seed=None):
        """
        Reseed the standard GAP pseudo-random sources with the given seed.

        Uses a random seed given by ``current_randstate().ZZ_seed()`` if
        ``seed=None``.  Otherwise the seed should be an integer.

        EXAMPLES::

            sage: libgap.set_seed(0)
            0
            sage: [libgap.Random(1, 10) for i in range(5)]
            [2, 3, 3, 4, 2]
        """
        if seed is None:
            seed = current_randstate().ZZ_seed()

        return super().set_seed(seed)

    def zero(self):
        """
        Return (integer) zero in GAP.

        OUTPUT:

        A :class:`~gappy.gapobj.GapInteger`.

        EXAMPLES::

            sage: libgap.zero()
            0
        """
        return self(0)

    def one(self):
        r"""
        Return (integer) one in GAP.

        OUTPUT:

        A :class:`~gappy.gapobj.GapInteger`.

        EXAMPLES::

            sage: libgap.one()
            1
        """
        return self(1)


libgap = SageGap()
