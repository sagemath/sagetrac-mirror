"""
Library Interface to GAP

This module implements a fast C library interface to GAP.
To use it, you simply call ``libgap`` (the parent of all
:class:`~sage.libs.gap.element.GapElement` instances) and use it to
convert Sage objects into GAP objects.

EXAMPLES::

    sage: a = libgap(10)
    sage: a
    10
    sage: type(a)
    <type 'sage.libs.gap.element.GapElement_Integer'>
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
    <type 'sage.libs.gap.element.GapElement_String'>

You can usually use the :meth:`~sage.libs.gap.element.GapElement.sage`
method to convert the resulting GAP element back to its Sage
equivalent::

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

#. The GAP containers ``List`` and ``rec`` are converted to Sage
   containers ``list`` and ``dict``.  Furthermore, the
   :meth:`~sage.libs.gap.element.GapElement.sage` method is applied
   recursively to the entries.

Special support is available for the GAP container classes. GAP lists
can be used as follows::

    sage: lst = libgap([1,5,7]);  lst
    [ 1, 5, 7 ]
    sage: type(lst)
    <type 'sage.libs.gap.element.GapElement_List'>
    sage: len(lst)
    3
    sage: lst[0]
    1
    sage: [ x^2 for x in lst ]
    [1, 25, 49]
    sage: type(_[0])
    <type 'sage.libs.gap.element.GapElement_Integer'>

Note that you can access the elements of GAP ``List`` objects as you
would expect from Python (with indexing starting at 0), but the
elements are still of type
:class:`~sage.libs.gap.element.GapElement`. The other GAP container
type are records, which are similar to Python dictionaries. You can
construct them directly from Python dictionaries::

    sage: libgap({'a':123, 'b':456})
    rec( a := 123, b := 456 )

Or get them as results of computations::

    sage: rec = libgap.eval('rec(a:=123, b:=456, Sym3:=SymmetricGroup(3))')
    sage: rec['Sym3']
    Sym( [ 1 .. 3 ] )
    sage: dict(rec)
    {'Sym3': Sym( [ 1 .. 3 ] ), 'a': 123, 'b': 456}

The output is a Sage dictionary whose keys are Sage strings and whose
Values are instances of :meth:`~sage.libs.gap.element.GapElement`. So,
for example, ``rec['a']`` is not a Sage integer. To recursively
convert the entries into Sage objects, you should use the
:meth:`~sage.libs.gap.element.GapElement.sage` method::

    sage: rec.sage()
    {'Sym3': NotImplementedError('cannot construct equivalent Sage object'...),
     'a': 123,
     'b': 456}

Now ``rec['a']`` is a Sage integer. We have not implemented the
conversion of the GAP symmetric group to the Sage symmetric group yet,
so you end up with a ``NotImplementedError`` exception object. The
exception is returned and not raised so that you can work with the
partial result.

While we don't directly support matrices yet, you can convert them to
Gap List of Lists. These lists are then easily converted into Sage
using the recursive expansion of the
:meth:`~sage.libs.gap.element.GapElement.sage` method::

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

from gappy.core cimport Gap as Gappy
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

from .element cimport GapElement, make_any_gap_element, make_GapElement_Function
from .saved_workspace import workspace as get_workspace
from .util import gap_root


############################################################################
### Gap  ###################################################################
############################################################################
# The libGap interpreter object Gap is the parent of the GapElements
# Provides a wrapper to gappy.core.Gap, which we can't subclass directly
# since it has an incompatible binary layout with Parent.
cdef class SageGappy(Gappy):
    """
    Subclasses `gappy.core.Gap` to provide some Sage-specific default
    behaviors.
    """

    def __init__(self):
        workspace, _ = get_workspace()
        Gappy.__init__(self, gap_root=gap_root(), workspace=workspace,
                       autoload=True)

    cpdef initialize(self):
        initializing = Gappy.initialize(self)

        if initializing:
            # These steps are only performed if Gappy has just been initialized
            # for the first time; if we don't check this then we'll cause an
            # infinite recursion
            workspace, workspace_is_up_to_date = get_workspace()
            if self.workspace == os.path.normpath(workspace):
                # Save a new workspace if necessary
                if not workspace_is_up_to_date:
                    prepare_workspace_dir()
                    with atomic_write(workspace) as f:
                        self.SaveWorkspace(f.name)

        return initializing


cdef class Gap(Parent):
    r"""
    The libgap interpreter object.

    .. NOTE::

        This object must be instantiated exactly once by the
        libgap. Always use the provided ``libgap`` instance, and never
        instantiate :class:`Gap` manually.

    EXAMPLES::

        sage: libgap.eval('SymmetricGroup(4)')
        Sym( [ 1 .. 4 ] )

    TESTS::

        sage: TestSuite(libgap).run(skip=['_test_category', '_test_elements', '_test_pickling'])
    """

    cdef readonly SageGappy gap

    Element = GapElement

    cpdef _coerce_map_from_(self, S):
        """
        Whether a coercion from `S` exists.

        INPUT / OUTPUT:

        See :mod:`sage.structure.parent`.

        EXAMPLES::

            sage: libgap.has_coerce_map_from(ZZ)
            True
            sage: libgap.has_coerce_map_from(CyclotomicField(5)['x','y'])
            True
        """
        # TODO: This seems wrong to me...
        return True

    def _element_constructor_(self, x):
        r"""
        Construct elements of this parent class.

        INPUT:

        - ``x`` -- anything that defines a GAP object.

        OUTPUT:

        A :class:`GapElement`.

        EXAMPLES::

            sage: libgap(0)   # indirect doctest
            0
            sage: libgap(ZZ(0))
            0
            sage: libgap(int(0))
            0
            sage: libgap(vector((0,1,2)))
            [ 0, 1, 2 ]
            sage: libgap(vector((1/3,2/3,4/5)))
            [ 1/3, 2/3, 4/5 ]
            sage: libgap(vector((1/3, 0.8, 3)))
            [ 0.333333, 0.8, 3. ]
            sage: v = _
            sage: libgap(v) is v
            True
        """
        # TODO: It might be good to implement a "fast lane" for some built-in
        # Sage types like Integer (e.g. currently Sage Integers are first
        # converted to Python ints, and then to GAP Integers, whereas it would
        # be much faster to convert Sage Integers directly to GAP Integers since
        # they are both basically mpz_t limbs under the hood).

        # If already a GapElement just return it directly
        if isinstance(x, GapElement):
            return x
        elif isinstance(x, SageObject):
            # Fast lane for all SageObjects
            x = x._libgap_()
            if isinstance(x, GapElement):
                return x
            # Otherwise, pass through self.gap() to convert

        return make_any_gap_element(self, self.gap(x))

    def eval(self, gap_command):
        """
        Evaluate a gap command and wrap the result.

        INPUT:

        - ``gap_command`` -- a string containing a valid gap command
          without the trailing semicolon.

        OUTPUT:

        A :class:`GapElement`.

        EXAMPLES::

            sage: libgap.eval('0')
            0
            sage: libgap.eval('"string"')
            "string"
        """

        if not isinstance(gap_command, basestring):
            gap_command = str(gap_command._libgap_init_())

        return make_any_gap_element(self, self.gap.eval(gap_command))

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
            return self.gap.load_package(pkg)
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

        A function wrapper
        :class:`~sage.libs.gap.element.GapElement_Function` for the
        GAP function. Calling it from Sage is equivalent to calling
        the wrapped function from GAP.

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

    def gap_function(self, func):
        """
        Create GAP functions from decorated Sage functions.

        EXAMPLES:

        The code for the GAP function is actually written in the Python
        function's docstring like so::

            sage: @libgap.gap_function
            ....: def one():
            ....:     '''
            ....:     Returns the multiplicative identity of the ring of integers.
            ....:
            ....:     function()
            ....:         return 1;
            ....:     end;
            ....:     '''
            sage: one
            <GAP function "one">
            sage: one()
            1

        Any text in the docstring before the first line beginning the text
        ``function()`` is used as the function's docstring.  Any following
        text is considered part of the function definition:

            sage: one.__doc__
            'Returns the multiplicative identity of the ring of integers.'

        Note that using this decorator does *not* cause the GAP interpreter
        to be initialized, so it can be used in module or class-level code.
        The GAP interpreter will only be initialized (if needed) the first time
        the function is called.

        Any Python code in the function's body will be disregarded, so this is
        in effect syntactic sugar for::

            sage: one = libgap.eval('function() return 1; end;')

        with the difference being that it can be used to pre-define GAP
        functions without invoking the GAP interpreter directly.

        This decorator may also be used on methods in classes.  In this case
        the ``self``--the instance of the class on which it is defined, is
        always passed as the first argument to the GAP function, *if* it has
        a conversion to a GAP type::

            sage: class MyInt(int):
            ....:     @libgap.gap_function
            ....:     def n_partitions(self):
            ....:         '''
            ....:         Compute the number of integer partitions.
            ....:
            ....:         function(n)
            ....:             local np;
            ....:             if n < 0 then
            ....:                 Error("must be a non-negative integer");
            ....:             fi;
            ....:             np:= function(n, m)
            ....:                local i, res;
            ....:                if n = 0 then
            ....:                   return 1;
            ....:                fi;
            ....:                res:= 0;
            ....:                for i in [1..Minimum(n,m)] do
            ....:                   res:= res + np(n-i, i);
            ....:                od;
            ....:                return res;
            ....:             end;
            ....:             return np(n,n);
            ....:         end;
            ....:         '''
            ....:
            sage: ten = MyInt(10)
            sage: ten.n_partitions()
            42
        """

        return make_GapElement_Function(self, self.gap.gap_function(func))

    def set_global(self, variable, value):
        """
        Set a GAP global variable

        INPUT:

        - ``variable`` -- string. The variable name.

        - ``value`` -- anything that defines a GAP object.

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
        return self.gap.set_global(variable, value)

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
        if self.gap.IsReadOnlyGlobal(variable):
            # TODO: Set correct ticket number
            deprecation(31297,
                f'{variable} is a read-only global; unsetting of read-only '
                f'globals is deprecated and will be removed in a future '
                f'version; if you need to unset a read-only global manually '
                f'call libgap.MakeReadWriteGlobal({variable!r}) first')
            self.gap.MakeReadWriteGlobal(variable)

        return self.gap.unset_global(variable)

    def get_global(self, variable):
        """
        Get a GAP global variable

        INPUT:

        - ``variable`` -- string. The variable name.

        OUTPUT:

        A :class:`~sage.libs.gap.element.GapElement` wrapping the GAP
        output. A ``ValueError`` is raised if there is no such
        variable in GAP.

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
        val = make_any_gap_element(self, self.gap.get_global(variable))
        if val is None:
            raise GAPError(f'no value bound to {variable}')
        return val

    def global_context(self, variable, value):
        """
        Temporarily change a global variable

        INPUT:

        - ``variable`` -- string. The variable name.

        - ``value`` -- anything that defines a GAP object.

        OUTPUT:

        A context manager that sets/reverts the given global variable.

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: with libgap.global_context('FooBar', 2):
            ....:     print(libgap.get_global('FooBar'))
            2
            sage: libgap.get_global('FooBar')
            1
        """
        return self.gap.global_context(variable, value)

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

        return self.gap.set_seed(seed)

    def _an_element_(self):
        r"""
        Return a :class:`GapElement`.

        OUTPUT:

        A :class:`GapElement`.

        EXAMPLES::

            sage: libgap.an_element()   # indirect doctest
            0
        """
        return self(0)

    def zero(self):
        """
        Return (integer) zero in GAP.

        OUTPUT:

        A :class:`GapElement`.

        EXAMPLES::

            sage: libgap.zero()
            0
        """
        return self(0)

    def one(self):
        r"""
        Return (integer) one in GAP.

        EXAMPLES::

            sage: libgap.one()
            1
            sage: parent(_)
            C library interface to GAP
        """
        return self(1)

    def __cinit__(self):
        self.gap = SageGappy()

    def __init__(self):
        r"""
        The Python constructor.

        EXAMPLES::

            sage: type(libgap)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: type(libgap._get_object())
            <class 'sage.libs.gap.libgap.Gap'>
        """
        Parent.__init__(self, base=ZZ)

    def __repr__(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: libgap
            C library interface to GAP
        """
        return 'C library interface to GAP'

    @cached_method
    def __dir__(self):
        """
        Customize tab completion

        EXAMPLES::

           sage: 'OctaveAlgebra' in dir(libgap)
           True
        """
        return dir(self.__class__) + sorted(common_gap_globals)

    def __getattr__(self, name):
        r"""
        The attributes of the Gap object are the GAP functions, and in some
        cases other global variables from GAP.

        INPUT:

        - ``name`` -- string. The name of the GAP function you want to
          call.

        OUTPUT:

        A :class:`GapElement`. A ``AttributeError`` is raised
        if there is no such function or global variable.

        EXAMPLES::

            sage: libgap.List
            <GAP function "List">
            sage: libgap.GlobalRandomSource
            <RandomSource in IsGlobalRandomSource>
        """

        # First try to get attributes from the category which is necessary
        # for coercion to work; for some reason when this is a pure Python
        # class we don't need this, but for cdef classes it doesn't go
        # through CategoryObject's __getattr__
        try:
            return self.getattr_from_category(name)
        except AttributeError:
            pass

        val = getattr(self.gap, name)
        if isinstance(val, GapObj):
            # Wrap GapObjs as GapElements
            val = make_any_gap_element(self, val)

        return val

    def show(self):
        """
        Return statistics about the GAP owned object list

        This includes the total memory allocated by GAP as returned by
        ``libgap.eval('TotalMemoryAllocated()'), as well as garbage collection
        / object count statistics as returned by
        ``libgap.eval('GasmanStatistics')``, and finally the total number of
        GAP objects held by Sage as :class:`~sage.libs.gap.element.GapElement`
        instances.

        The value ``livekb + deadkb`` will roughly equal the total memory
        allocated for GAP objects (see
        ``libgap.eval('TotalMemoryAllocated()')``).

        .. note::

            Slight complication is that we want to do it without accessing
            libgap objects, so we don't create new GapElements as a side
            effect.

        EXAMPLES::

            sage: a = libgap(123)
            sage: b = libgap(456)
            sage: c = libgap(789)
            sage: del b
            sage: libgap.collect()
            sage: libgap.show()  # random output
            {'gasman_stats': {'full': {'cumulative': 110,
               'deadbags': 321400,
               'deadkb': 12967,
               'freekb': 15492,
               'livebags': 396645,
               'livekb': 37730,
               'time': 110,
               'totalkb': 65536},
              'nfull': 1,
              'npartial': 1},
             'nelements': 23123,
             'total_alloc': 3234234}
        """
        return self.gap.show()

    def count_GAP_objects(self):
        """
        Return the number of GAP objects that are being tracked by
        GAP.

        OUTPUT:

        An integer

        EXAMPLES::

            sage: libgap.count_GAP_objects()   # random output
            5
        """
        return self.gap.count_GAP_objects()

    def collect(self):
        """
        Manually run the garbage collector

        EXAMPLES::

            sage: a = libgap(123)
            sage: del a
            sage: libgap.collect()
        """
        return self.gap.collect()


libgap = Gap()


@libgap.gap.convert_from(SageObject)
def _sageobject_to_gapobj(gap, obj):
    r"""
    gappy converter for converting generic `.SageObject`\s to their
    corresponding `.GapObj` if any.

    This implements the libgap conversion functions already documented for
    `.SageObject`\s: `.SageObject._libgap_` and `.SageObject._libgap_init_`.
    """

    # NOTE: In the default implementation of SageObject._libgap_ it defers
    # to _libgap_init_, so we just need to try calling _libgap_
    ret = obj._libgap_()
    if isinstance(ret, GapElement):
        return (<GapElement>ret).obj
    elif isinstance(ret, gap.supported_builtins):
        return gap(ret)
    elif isinstance(ret, GapObj):
        return ret
    else:
        raise RuntimeError(
            f'{type(obj).__name__}._libgap_ returned something that cannot '
            f'be converted to a GAP object: {ret!r}')


@libgap.gap.convert_from(GapElement)
def _gapelement_to_gapobj(gap, elem):
    r"""
    gappy converter function for converting `.GapElement`\s to their
    corresponding `.GapObj`.

    In this case it simply returns the `~gappy.gapobj.GapObj` wrapped by the
    `.GapElement`.  This allows seamlessly passing `.GapElement`\s to gappy.
    """

    return (<GapElement>elem).obj
