r"""
Enumerated set from iterator

EXAMPLES:

We build a set from the iterator ``graphs`` that returns a canonical
representative for each isomorphism class of graphs::

    sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
    sage: E = EnumeratedSetFromIterator(
    ...     graphs,
    ...     name = "Graphs",
    ...     category = InfiniteEnumeratedSets(),
    ...     cache = True)
    sage: E
    Graphs
    sage: E.unrank(0)
    Graph on 0 vertices
    sage: E.unrank(4)
    Graph on 3 vertices
    sage: E.cardinality()
    +Infinity
    sage: E.category()
    Category of facade infinite enumerated sets

The module also provides decorator for functions and methods::

    sage: from sage.sets.set_from_iterator import set_from_function
    sage: @set_from_function
    ... def f(n): return xsrange(n)
    sage: f(3)
    {0, 1, 2}
    sage: f(5)
    {0, 1, 2, 3, 4}
    sage: f(100)
    {0, 1, 2, 3, 4, ...}

    sage: from sage.sets.set_from_iterator import set_from_method
    sage: class A:
    ...    @set_from_method
    ...    def f(self,n):
    ...        return xsrange(n)
    sage: a = A()
    sage: a.f(3)
    {0, 1, 2}
    sage: f(100)
    {0, 1, 2, 3, 4, ...}
"""
#*****************************************************************************
#  Copyright (C) 2012 Vincent Delecroix <vincent.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.parent import Parent
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_cat import EmptySetError
from itertools import izip_longest
import os
from sage.misc.function_mangling import ArgumentFixer
from sage.misc.lazy_list import lazy_list

class EnumeratedSetFromIterator(Parent):
    """
    A class for enumerated set built from an iterator.

    INPUT:

    - ``f`` -- a function that returns an iterable from which the set is built from

    - ``args`` -- tuple -- arguments to be sent to the function ``f``

    - ``kwds`` -- dictionnary -- keywords to be sent to the function ``f``

    - ``name`` -- an optional name for the set

    - ``category`` -- (default: ``None``) an optional category for that
      enumerated set. If you know that your iterator will stop after a finite
      number of steps you should set it as :class:`FiniteEnumeratedSets`, conversly if
      you know that your iterator will run over and over you should set it as
      :class:`InfiniteEnumeratedSets`.

    - ``cache`` -- boolean (default: ``False``) -- Whether or not use a cache
      mechanism for the iterator. If ``True``, then the function ``f`` is called
      only once.


    EXAMPLES::

        sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
        sage: E = EnumeratedSetFromIterator(graphs, args = (7,))
        sage: E
        {Graph on 7 vertices, Graph on 7 vertices, Graph on 7 vertices, Graph on 7 vertices, Graph on 7 vertices, ...}
        sage: E.category()
        Category of facade enumerated sets

    The same example with a cache and a custom name::

        sage: E = EnumeratedSetFromIterator(
        ...      graphs,
        ...      args = (8,),
        ...      category = FiniteEnumeratedSets(),
        ...      name = "Graphs with 8 vertices",
        ...      cache = True)
        sage: E
        Graphs with 8 vertices
        sage: E.unrank(3)
        Graph on 8 vertices
        sage: E.category()
        Category of facade finite enumerated sets

    TESTS:

    The cache is compatible with multiple call to ``__iter__``::

        sage: from itertools import count
        sage: E = EnumeratedSetFromIterator(count, args=(0,), category=InfiniteEnumeratedSets(), cache=True)
        sage: e1 = iter(E)
        sage: e2 = iter(E)
        sage: next(e1), next(e1)
        (0, 1)
        sage: next(e2), next(e2), next(e2)
        (0, 1, 2)
        sage: next(e1), next(e1)
        (2, 3)
        sage: next(e2)
        3

    The following warning is due to ``E`` being a facade parent. For more,
    see the discussion on :trac:`16239`::

        sage: TestSuite(E).run()
        doctest:...: UserWarning: Testing equality of infinite sets which will not end in case of equality

        sage: E = EnumeratedSetFromIterator(xsrange, args=(10,), category=FiniteEnumeratedSets(), cache=True)
        sage: TestSuite(E).run()

    .. NOTE::

        In order to make the ``TestSuite`` works, the elements of the set
        should have parents.
    """
    def __init__(self, f, args=None, kwds=None, name=None, category=None, cache=False):
        """
        TESTS::

            sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
            sage: S = EnumeratedSetFromIterator(xsrange, (1,200,-1), category=FiniteEnumeratedSets())
            sage: TestSuite(S).run()
        """
        if category is not None:
            Parent.__init__(self, facade = True, category = category)
        else:
            Parent.__init__(self, facade = True, category = EnumeratedSets())


        if name is not None:
            self.rename(name)

        self._func = f

        if args is not None:
            self._args = args
        if kwds is not None:
            self._kwds = kwds

        if cache:
            self._cache = lazy_list(iter(self._func(
                                         *getattr(self, '_args', ()),
                                        **getattr(self, '_kwds', {}))))

    def __reduce__(self):
        r"""
        Support for pickle.

        TESTS::

            sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
            sage: from sage.graphs.graph_generators import graphs
            sage: E = EnumeratedSetFromIterator(graphs,
            ...      args=(3,),
            ...      category=FiniteEnumeratedSets(),
            ...      name="Graphs on 3 vertices")
            sage: E
            Graphs on 3 vertices
            sage: F = loads(dumps(E)); F
            Graphs on 3 vertices
            sage: E == F
            True
        """
        return (EnumeratedSetFromIterator,
                (self._func,                           # func
                 getattr(self, '_args', None),         # args
                 getattr(self, '_kwds', None),         # kwds
                 getattr(self, '__custom_name', None), # name
                 self.category(),                      # category
                 hasattr(self, '_cache'))              # cache
                )

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
            sage: E = EnumeratedSetFromIterator(Partitions(7,min_part=2).__iter__)
            sage: repr(E)    # indirect doctest
            '{[7], [5, 2], [4, 3], [3, 2, 2]}'
            sage: E = EnumeratedSetFromIterator(Partitions(9,min_part=2).__iter__)
            sage: repr(E)    # indirect doctest
            '{[9], [7, 2], [6, 3], [5, 4], [5, 2, 2], ...}'
            sage: E = EnumeratedSetFromIterator(Partitions(9,min_part=2).__iter__, name="Some partitions")
            sage: repr(E)    # indirect doctest
            'Some partitions'
        """
        l = []
        i = iter(self)
        for _ in xrange(6):
            try:
                l.append(next(i))
            except StopIteration:
                break
        if len(l) < 6:
            return '{' + ', '.join(repr(x) for x in l) + '}'
        l.pop(-1)
        return '{' + ', '.join(repr(x) for x in l) + ', ...}'

    def __contains__(self, x):
        r"""
        Test whether ``x`` is in ``self``.

        If the set is infinite, only the answer ``True`` should be expected in
        finite time.

        EXAMPLES::

            sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
            sage: P = Partitions(12,min_part=2,max_part=5)
            sage: E = EnumeratedSetFromIterator(P.__iter__)
            sage: P([5,5,2]) in E
            True
        """
        return any(x == y for y in self)

    is_parent_of = __contains__

    #TODO: what should we do for comparisons of infinite sets
    def __eq__(self, other):
        r"""
        Equality test.

        The function returns ``True`` if and only if other is an enumerated
        set and has the same element as ``self``.

        TESTS::

            sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
            sage: E4 = EnumeratedSetFromIterator(graphs, args=(4,), category=FiniteEnumeratedSets())
            sage: F4 = EnumeratedSetFromIterator(graphs, args=(4,), category=FiniteEnumeratedSets())
            sage: E5 = EnumeratedSetFromIterator(graphs, args=(5,), category=FiniteEnumeratedSets())
            sage: E4 == E4
            True
            sage: E4 == F4
            True
            sage: E4 == E5
            False
            sage: E5 == E4
            False
            sage: E5 == E5
            True
        """
        if isinstance(other, EnumeratedSetFromIterator):
            # trick to allow equality between infinite sets
            # this assume that the function does not return randomized data!
            if (self._func == other._func and
                getattr(self, '_args', None) == getattr(other, '_args', None) and
                getattr(self, '_kwds', None) == getattr(other, '_kwds', None)):
                return True

        if other in EnumeratedSets():
            #TODO: think about what should be done at that point
            if self not in FiniteEnumeratedSets() and other not in FiniteEnumeratedSets():
                import warnings
                warnings.warn("Testing equality of infinite sets which will not end in case of equality")

            i1 = iter(self)
            i2 = iter(other)
            while True:
                try:
                    x = next(i1)
                except StopIteration:
                    try:
                        next(i2)
                        return False
                    except StopIteration:
                        return True
                try:
                    y = next(i2)
                except StopIteration:
                    return False
                if x != y:
                    return False

    def __ne__(self,other):
        r"""
        Difference test.

        The function calls the ``__eq__`` test.

        TESTS::

            sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
            sage: E4 = EnumeratedSetFromIterator(graphs, args=(4,), category=FiniteEnumeratedSets())
            sage: F4 = EnumeratedSetFromIterator(graphs, args=(4,), category=FiniteEnumeratedSets())
            sage: E5 = EnumeratedSetFromIterator(graphs, args=(5,), category=FiniteEnumeratedSets())
            sage: E4 != E4
            False
            sage: E4 != F4
            False
            sage: E4 != E5
            True
            sage: E5 != E4
            True
            sage: E5 != E5
            False
        """
        return not self.__eq__(other)

    def __iter__(self):
        r"""
        Returns an iterator over the element of ``self``.

        EXAMPLES::

            sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
            sage: E = EnumeratedSetFromIterator(graphs, args=(8,))
            sage: g1 = next(iter(E)); g1
            Graph on 8 vertices
            sage: E = EnumeratedSetFromIterator(graphs, args=(8,), cache=True)
            sage: g2 = next(iter(E)); g2
            Graph on 8 vertices
            sage: g1 == g2
            True
        """
        if hasattr(self, '_cache'):
            return iter(self._cache)
        return iter(self._func(*getattr(self, '_args', ()), **getattr(self, '_kwds', {})))

    def unrank(self, i):
        r"""
        Returns the element at position ``i``.

        EXAMPLES::

            sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
            sage: E = EnumeratedSetFromIterator(graphs, args=(8,), cache=True)
            sage: F = EnumeratedSetFromIterator(graphs, args=(8,), cache=False)
            sage: E.unrank(2)
            Graph on 8 vertices
            sage: E.unrank(2) == F.unrank(2)
            True
        """
        if hasattr(self, '_cache'):
            return self._cache[i]
        return super(EnumeratedSetFromIterator,self).unrank(i)

    def _element_constructor_(self, el):
        """
        Construct an element from ``el``.

        TESTS::

            sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
            sage: S = EnumeratedSetFromIterator(xrange, args=(1,4))
            sage: S(1)  # indirect doctest
            1
            sage: S(0)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: 0 not in {1, 2, 3}
        """
        if el in self:
            return el
        else:
            raise ValueError("%s not in %s"%(el, self))

    def clear_cache(self):
        r"""
        Clear the cache.

        EXAMPLES::

            sage: from itertools import count
            sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
            sage: E = EnumeratedSetFromIterator(count, args=(1,), cache=True)
            sage: e1 = E._cache
            sage: e1
            lazy list [1, 2, 3, ...]
            sage: E.clear_cache()
            sage: E._cache
            lazy list [1, 2, 3, ...]
            sage: e1 is E._cache
            False
        """
        if hasattr(self, '_cache'):
            self._cache = lazy_list(iter(self._func(
                                         *getattr(self, '_args', ()),
                                        **getattr(self, '_kwds', {}))))

#
# Decorators
#

#TODO: move it in sage.misc ?
class Decorator:
    r"""
    Abstract class that manage documentation and sources of the wrapped object.

    The method needs to be stored in the attribute ``self.f``
    """
    def _sage_doc_(self):
        """
        Provide documentation for the wrapped function.

        TESTS::

            sage: from sage.misc.sageinspect import sage_getdoc
            sage: from sage.sets.set_from_iterator import Decorator
            sage: d = Decorator()
            sage: d.f = Integer.is_prime
            sage: print sage_getdoc(d)   # indirect doctest
               Test whether "self" is prime.
            ...
               IMPLEMENTATION: Calls the PARI "isprime" function.
            <BLANKLINE>
        """
        from sage.misc.sageinspect import sage_getsourcelines, sage_getfile, _extract_embedded_position
        f = self.f
        doc = f.__doc__ or ''
        if _extract_embedded_position(doc.splitlines()[0]) is None:
            try:
                sourcelines = sage_getsourcelines(f)
                from sage.env import SAGE_LIB, SAGE_SRC
                filename = sage_getfile(f)
                # The following is a heuristics to get
                # the file name of the cached function
                # or method
                if filename.startswith(SAGE_SRC):
                    filename = filename[len(SAGE_SRC):]
                elif filename.startswith(SAGE_LIB):
                    filename = filename[len(SAGE_LIB):]
                file_info = "File: %s (starting at line %d)\n"%(filename,sourcelines[1])
                doc = file_info+doc
            except IOError:
                pass
        return doc

    def _sage_src_(self):
        r"""
        Returns the source code for the wrapped function.

        TESTS::

            sage: from sage.misc.sageinspect import sage_getsource
            sage: from sage.sets.set_from_iterator import Decorator
            sage: d = Decorator()
            sage: d.f = Rational.is_square
            sage: print sage_getsource(d.f)   # indirect doctest
            def is_square(self):
            ...
                return mpq_sgn(self.value) >= 0 and mpz_perfect_square_p(mpq_numref(self.value)) and mpz_perfect_square_p(mpq_denref(self.value))
        """
        from sage.misc.sageinspect import sage_getsource
        return sage_getsource(self.f)

    def _sage_src_lines_(self):
        r"""
        Returns the list of source lines and the first line number
        of the wrapped function.

        TESTS::

            sage: from sage.misc.sageinspect import sage_getsourcelines
            sage: from sage.sets.set_from_iterator import Decorator
            sage: d = Decorator()
            sage: d.f = MathieuGroup.order
            sage: S = sage_getsourcelines(d)   # indirect doctest
            sage: S[0][2]
            '        Return the number of elements of this group.\n'
            sage: S[0][18]
            '            return Integer(1)\n'
        """
        from sage.misc.sageinspect import sage_getsourcelines
        return sage_getsourcelines(self.f)

    def _sage_argspec_(self):
        """
        Return the argument specification of the wrapped function or method.

        TESTS::

            sage: from sage.misc.sageinspect import sage_getargspec
            sage: from sage.sets.set_from_iterator import Decorator
            sage: d = Decorator()
            sage: d.f = find_local_minimum
            sage: sage_getargspec(d) # indirect doctest
            ArgSpec(args=['f', 'a', 'b', 'tol', 'maxfun'], varargs=None, keywords=None, defaults=(1.48e-08, 500))
        """
        from sage.misc.sageinspect import sage_getargspec
        return sage_getargspec(self.f)

    def __call__(self, *args, **kwds):
        r"""
        Call function.

        Needs to be implemented in derived subclass.

        TEST::

            sage: from sage.sets.set_from_iterator import Decorator
            sage: d = Decorator()
            sage: d()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

class EnumeratedSetFromIterator_function_decorator(Decorator):
    r"""
    Decorator for :class:`EnumeratedSetFromIterator`.

    Name could be string or a function ``(args,kwds) -> string``.

    .. WARNING::

        If you are going to use this with the decorator ``cached_function``,
        you must place the ``cached_function`` first. See the example below.

    EXAMPLES::

        sage: from sage.sets.set_from_iterator import set_from_function
        sage: @set_from_function
        ... def f(n):
        ...    for i in xrange(n):
        ...        yield i**2 + i + 1
        sage: f(3)
        {1, 3, 7}
        sage: f(100)
        {1, 3, 7, 13, 21, ...}

    To avoid ambiguity, it is always better to use it with a call which
    provides optional global initialization for the call to
    :class:`EnumeratedSetFromIterator`::

        sage: @set_from_function(category=InfiniteEnumeratedSets())
        ... def Fibonacci():
        ...    a = 1; b = 2
        ...    while True:
        ...       yield a
        ...       a,b = b,a+b
        sage: F = Fibonacci()
        sage: F
        {1, 2, 3, 5, 8, ...}
        sage: F.cardinality()
        +Infinity

    A simple example with many options::

        sage: @set_from_function(
        ...        name = "From %(m)d to %(n)d",
        ...        category = FiniteEnumeratedSets())
        ... def f(m,n): return xsrange(m,n+1)
        sage: E = f(3,10); E
        From 3 to 10
        sage: E.list()
        [3, 4, 5, 6, 7, 8, 9, 10]
        sage: E = f(1,100); E
        From 1 to 100
        sage: E.cardinality()
        100
        sage: f(n=100,m=1) == E
        True

    An example which mixes together ``set_from_function`` and
    ``cached_method``::

        sage: @cached_function
        ... @set_from_function(
        ...    name = "Graphs on %(n)d vertices",
        ...    category = FiniteEnumeratedSets(),
        ...    cache = True)
        ... def Graphs(n): return graphs(n)
        sage: Graphs(10)
        Graphs on 10 vertices
        sage: Graphs(10).unrank(0)
        Graph on 10 vertices
        sage: Graphs(10) is Graphs(10)
        True

    The ``cached_function`` must go first::

        sage: @set_from_function(
        ...    name = "Graphs on %(n)d vertices",
        ...    category = FiniteEnumeratedSets(),
        ...    cache = True)
        ... @cached_function
        ... def Graphs(n): return graphs(n)
        sage: Graphs(10)
        Graphs on 10 vertices
        sage: Graphs(10).unrank(0)
        Graph on 10 vertices
        sage: Graphs(10) is Graphs(10)
        False
    """
    def __init__(self, f=None, name=None, **options):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.sets.set_from_iterator import set_from_function
            sage: F = set_from_function(category=FiniteEnumeratedSets())(xsrange)
            sage: TestSuite(F(100)).run()
            sage: TestSuite(F(1,5,2)).run()
            sage: TestSuite(F(0)).run()
        """
        if f is not None:
            self.f = f
            if hasattr(f, "__name__"):
                self.__name__ = f.__name__
            else:
                self.__name__ = f.__name__
            self.__module__ = f.__module__
            self.af = ArgumentFixer(f)
        if name is not None:
            self.name = name
        self.options = options

    def __call__(self, *args, **kwds):
        r"""
        Build a new :class:`EnumeratedSet` by calling ``self.f`` with
        apropriate argument. If ``f`` is ``None``, then returns a new instance
        of :class:`EnumeratedSetFromIterator`.

        EXAMPLES::

            sage: from sage.sets.set_from_iterator import set_from_function
            sage: F = set_from_function(category=FiniteEnumeratedSets())(xsrange)
            sage: F(3)
            {0, 1, 2}
            sage: F(end=7,start=3)
            {3, 4, 5, 6}
            sage: F(10).cardinality()
            10
        """
        options = self.options

        if hasattr(self, 'f'): # yet initialized
            if hasattr(self,'name'):
                if isinstance(self.name,str):
                    if args or kwds:
                        _,kk = self.af.fix_to_named(*args,**kwds)
                        name = self.name%dict(kk)
                    else:
                        name = self.name
                else:
                    name = self.name(*args,**kwds)
                return EnumeratedSetFromIterator(self.f, args, kwds, name=name, **self.options)
            return EnumeratedSetFromIterator(self.f, args, kwds, **self.options)

        else: # potential global options
            if args == ():
                assert len(kwds.keys()) == 1
                f = kwds.values()[0]
            else:
                assert len(args) == 1
                f = args[0]
            return EnumeratedSetFromIterator_function_decorator(
                f,
                name=getattr(self,'name',None),
                **self.options)

set_from_function = EnumeratedSetFromIterator_function_decorator

class EnumeratedSetFromIterator_method_caller(Decorator):
    r"""
    Caller for decorated method in class.

    INPUT:

    - ``inst`` -- an instance of a class

    - ``f`` -- a method of a class of ``inst`` (and not of the instance itself)

    - ``name`` -- optional -- either a string (which may contains substitution
      rules from argument or a function args,kwds -> string.

    - ``options`` -- any option accepted by :class:`EnumeratedSetFromIterator`
    """
    def __init__(self, inst, f, name=None, **options):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.sets.set_from_iterator import DummyExampleForPicklingTest
            sage: d = DummyExampleForPicklingTest()
            sage: d.f()
            {10, 11, 12, 13, 14, ...}

        It is possible to pickle/unpickle the class and the instance::

            sage: loads(dumps(DummyExampleForPicklingTest))().f()
            {10, 11, 12, 13, 14, ...}
            sage: loads(dumps(d)).f()
            {10, 11, 12, 13, 14, ...}

        But not the enumerated set::

            sage: loads(dumps(d.f()))
            Traceback (most recent call last):
            ...
            PicklingError: Can't pickle <type 'function'>: attribute lookup __builtin__.function failed
        """
        self.inst = inst
        self.f = f
        self.af = ArgumentFixer(self.f)
        if hasattr(f, "__name__"):
            self.__name__ = f.__name__
        else:
            self.__name__ = f.__name__
        self.__module__ = f.__module__

        self.name = name
        self.options = options

    def __call__(self,*args,**kwds):
        r"""
        Returns an instance of :class:`EnumeratedSetFromIterator` with
        proper argument.

        TESTS::

            sage: from sage.sets.set_from_iterator import set_from_method
            sage: class A:
            ...    @set_from_method(name = lambda self,n: str(self)*n)
            ...    def f(self,n):
            ...        return xsrange(n)
            ...    def __repr__(self):
            ...        return "A"
            sage: a = A()
            sage: a.f(3)                         # indirect doctest
            AAA
            sage: A.f(a,3)                       # indirect doctest
            AAA
            sage: [x for x in a.f(6)]            # indirect doctest
            [0, 1, 2, 3, 4, 5]
        """
        if self.inst is not None:
            args = (self.inst,) + args
        if self.name:
            if isinstance(self.name,str):
                aa,kk = self.af.fix_to_named(*args,**kwds)
                name = self.name%dict(kk)
            else:
                name = self.name(*args, **kwds)
            return EnumeratedSetFromIterator(self.f, args, kwds, name, **self.options)
        return EnumeratedSetFromIterator(self.f, args, kwds, **self.options)

    def __get__(self, inst, cls):
        r"""
        Get a :class:`EnumeratedSetFromIterator_method_caller` bound to a
        specific instance of the class of the cached method.

        .. NOTE::

            :class:`EnumeratedSetFromIterator_method_caller` has a separate
            ``__get__`` because of the special behavior of category framework
            for element classes which are not of extension type (see
            :meth:`sage.structure.element.Element.__get__`).

        TESTS::

            sage: from sage.sets.set_from_iterator import set_from_method
            sage: from sage.structure.misc import getattr_from_other_class
            sage: class A:
            ...      stop = 10000
            ...      @set_from_method
            ...      def f(self,start):
            ...          return xsrange(start,self.stop)
            sage: a = A()
            sage: getattr_from_other_class(a, A, 'f')(4)
            {4, 5, 6, 7, 8, ...}

            sage: class B:
            ...      stop = 10000
            ...      @set_from_method(category=FiniteEnumeratedSets())
            ...      def f(self,start):
            ...          return xsrange(start,self.stop)
            sage: b = B()
            sage: getattr_from_other_class(b, B, 'f')(2)
            {2, 3, 4, 5, 6, ...}
        """
        return EnumeratedSetFromIterator_method_caller(
                inst, self.f,
                self.name,
                **self.options)

class EnumeratedSetFromIterator_method_decorator(object):
    r"""
    Decorator for enumerated set built from a method.

    INPUT:

    - ``f`` -- Optional function from which are built the enumerated sets at
      each call

    - ``name`` -- Optional string (which may contains substitution rules from
      argument) or a function ``(args,kwds) -> string``.

    - any option accepted by :class:`EnumeratedSetFromIterator`.

    EXAMPLES::

        sage: from sage.sets.set_from_iterator import set_from_method
        sage: class A():
        ...    def n(self): return 12
        ...    @set_from_method
        ...    def f(self): return xsrange(self.n())
        sage: a = A()
        sage: print a.f.__class__
        sage.sets.set_from_iterator.EnumeratedSetFromIterator_method_caller
        sage: a.f()
        {0, 1, 2, 3, 4, ...}
        sage: A.f(a)
        {0, 1, 2, 3, 4, ...}

    A more complicated example with a parametrized name::

        sage: class B():
        ...    @set_from_method(
        ...        name = "Graphs(%(n)d)",
        ...        category = FiniteEnumeratedSets())
        ...    def graphs(self, n): return graphs(n)
        sage: b = B()
        sage: G3 = b.graphs(3)
        sage: G3
        Graphs(3)
        sage: G3.cardinality()
        4
        sage: G3.category()
        Category of facade finite enumerated sets
        sage: B.graphs(b,3)
        Graphs(3)

    And a last example with a name parametrized by a function::

        sage: class D():
        ...    def __init__(self, name): self.name = str(name)
        ...    def __str__(self): return self.name
        ...    @set_from_method(
        ...        name = lambda self,n: str(self)*n,
        ...        category = FiniteEnumeratedSets())
        ...    def subset(self, n):
        ...        return xsrange(n)
        sage: d = D('a')
        sage: E = d.subset(3); E
        aaa
        sage: E.list()
        [0, 1, 2]
        sage: F = d.subset(n=10); F
        aaaaaaaaaa
        sage: F.list()
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    .. TODO::

        It is not yet possible to use ``set_from_method`` in conjunction with
        ``cached_method``.
    """
    def __init__(self, f=None, **options):
        r"""
        Initialize ``self``.

        TESTS:

        We test if pickling works correctly on the Permutation class (in
        :mod:`sage.combinat.permutation`) because its method ``bruhat_succ``
        and ``bruhat_pred`` are decorated with ``set_from_method``::

            sage: from sage.combinat.permutation import Permutation
            sage: loads(dumps(Permutation))
            <class 'sage.combinat.permutation.Permutation'>
            sage: p = Permutation([3,2,1])
            sage: loads(dumps(p)) == p
            True
        """
        if f is not None:
            import types
            self.f = f
            if hasattr(f,"__name__"):
                self.__name__ = f.__name__
                self.__module__ = f.__module__

            else:
                if hasattr(f, '__module__'):
                    self.__module__ = f.__module__
                elif hasattr(f, '__func__'):
                    self.__module__ = f.__func__.__module__

                if hasattr(f, '__name__'):
                    self.__name__ = f.__name__
                elif hasattr(f, '__func__'):
                    self.__name__ = f.__func__.__name__

        self.options = options

    def __call__(self, f):
        r"""
        Trick if :class:`EnumeratedSetFromIterator_method` was created with
        some options and is called with a function as argument.

        TESTS::

            sage: from sage.sets.set_from_iterator import set_from_method
            sage: class A:
            ...    @set_from_method()    # indirect doctest
            ...    def f(self):
            ...        return xsrange(3)
            sage: a = A()
            sage: a.f()
            {0, 1, 2}
        """
        return EnumeratedSetFromIterator_method_decorator(f,**self.options)

    def __get__(self, inst, cls):
        r"""
        TESTS::

            sage: from sage.sets.set_from_iterator import set_from_method
            sage: class A():
            ...    def n(self): return 12
            ...    @set_from_method
            ...    def f(self): return xsrange(self.n())
            sage: a = A()
            sage: print A.f.__class__
            sage.sets.set_from_iterator.EnumeratedSetFromIterator_method_caller
            sage: print a.f.__class__
            sage.sets.set_from_iterator.EnumeratedSetFromIterator_method_caller
        """
        # You would hardly ever see an instance of this class alive.
        return EnumeratedSetFromIterator_method_caller(inst, self.f, **self.options)

set_from_method = EnumeratedSetFromIterator_method_decorator

class DummyExampleForPicklingTest:
    r"""
    Class example to test pickling with the decorator :class:`set_from_method`.

    .. WARNING::

        This class is intended to be used in doctest only.

    EXAMPLES::

        sage: from sage.sets.set_from_iterator import DummyExampleForPicklingTest
        sage: DummyExampleForPicklingTest().f()
        {10, 11, 12, 13, 14, ...}
    """
    start = 10
    stop  = 100
    @set_from_method
    def f(self):
        r"""
        Returns the set between ``self.start`` and ``self.stop``.

        EXAMPLES::

            sage: from sage.sets.set_from_iterator import DummyExampleForPicklingTest
            sage: d = DummyExampleForPicklingTest()
            sage: d.f()
            {10, 11, 12, 13, 14, ...}
            sage: d.start = 4
            sage: d.stop = 200
            sage: d.f()
            {4, 5, 6, 7, 8, ...}
        """
        from sage.misc.misc import xsrange
        return xsrange(self.start, self.stop)
