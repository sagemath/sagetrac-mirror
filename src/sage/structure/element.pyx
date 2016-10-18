r"""
Elements

AUTHORS:

- David Harvey (2006-10-16): changed CommutativeAlgebraElement to
  derive from CommutativeRingElement instead of AlgebraElement

- David Harvey (2006-10-29): implementation and documentation of new
  arithmetic architecture

- William Stein (2006-11): arithmetic architecture -- pushing it
  through to completion.

- Gonzalo Tornaria (2007-06): recursive base extend for coercion --
  lots of tests

- Robert Bradshaw (2007-2010): arithmetic operators and coercion

- Maarten Derickx (2010-07): added architecture for is_square and sqrt

- Jeroen Demeyer (2016-08): moved all coercion to the base class
  :class:`Element`, see :trac:`20767`

The Abstract Element Class Hierarchy
====================================

This is the abstract class hierarchy, i.e., these are all
abstract base classes.

::

    SageObject
        Element
            ModuleElement
                RingElement
                    CommutativeRingElement
                        IntegralDomainElement
                            DedekindDomainElement
                                PrincipalIdealDomainElement
                                    EuclideanDomainElement
                        FieldElement
                        CommutativeAlgebraElement
                    AlgebraElement   (note -- can't derive from module, since no multiple inheritance)
                        Matrix
                    InfinityElement
                AdditiveGroupElement
                Vector

            MonoidElement
                MultiplicativeGroupElement
        ElementWithCachedMethod


How to Define a New Element Class
=================================

Elements typically define a method ``_new_c``, e.g.,

.. code-block:: cython

    cdef _new_c(self, defining data):
        cdef FreeModuleElement_generic_dense x
        x = FreeModuleElement_generic_dense.__new__(FreeModuleElement_generic_dense)
        x._parent = self._parent
        x._entries = v

that creates a new sibling very quickly from defining data
with assumed properties.

.. _element_arithmetic:

Arithmetic for Elements
-----------------------

Sage has a special system for handling arithmetic operations on Sage
elements (that is instances of :class:`Element`), in particular to
manage uniformly mixed arithmetic operations using the :mod:`coercion
model <sage.structure.coerce>`. We describe here the rules that must
be followed by both arithmetic implementers and callers.

A quick summary for the impatient
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To implement addition for any :class:`Element` subclass, override the
``def _add_(self, other)`` method instead of the usual Python
``__add__`` :python:`special method <reference/datamodel.html#special-method-names>`.
Within ``_add_(self, other)``, you may assume that ``self`` and
``other`` have the same parent.

If the implementation is generic across all elements in a given
category `C`, then this method can be put in ``C.ElementMethods``.

When writing *Cython* code, ``_add_`` should be a cpdef method:
``cpdef _add_(self, other)``.

When doing arithmetic with two elements having different parents,
the :mod:`coercion model <sage.structure.coerce>` is responsible for
"coercing" them to a common parent and performing arithmetic on the
coerced elements.

Arithmetic in more detail
^^^^^^^^^^^^^^^^^^^^^^^^^

The aims of this system are to provide (1) an efficient calling protocol
from both Python and Cython, (2) uniform coercion semantics across Sage,
(3) ease of use, (4) readability of code.

We will take addition as an example; all other operators are similar.
There are two relevant functions, with differing names
(single vs. double underscores).

-  **def Element.__add__(left, right)**

   This function is called by Python or Cython when the binary "+"
   operator is encountered. It assumes that at least one of its
   arguments is an :class:`Element`.

   It has a fast pathway to deal with the most common case where both
   arguments have the same parent. Otherwise, it uses the coercion
   model to work out how to make them have the same parent. The
   coercion model then adds the coerced elements (technically, it calls
   ``operator.add``). Note that the result of coercion is not required
   to be a Sage :class:`Element`, it could be a plain Python type.

   Note that, although this function is declared as ``def``, it doesn't
   have the usual overheads associated with Python functions (either
   for the caller or for ``__add__`` itself). This is because Python
   has optimised calling protocols for such special functions.

-  **def Element._add_(self, other)**

   This is the function that you should override to implement addition
   in a subclass of :class:`Element`.

   The two arguments to this function are guaranteed to have the **same
   parent**, but not necessarily the same Python type.

   When implementing ``_add_`` in a Cython extension type, use
   ``cpdef _add_`` instead of ``def _add_``.

   In Cython code, if you want to add two elements and you know that
   their parents are identical, you are encouraged to call this
   function directly, instead of using ``x + y``. This only works if
   Cython knows that the left argument is an ``Element``. You can
   always cast explicitly: ``(<Element>x)._add_(y)`` to force this.
   In plain Python, ``x + y`` is always the fastest way to add two
   elements because the special method ``__add__`` is optimized
   unlike the normal method ``_add_``.

The difference in the names of the arguments (``left, right``
versus ``self, other``) is intentional: ``self`` is guaranteed to be an
instance of the class in which the method is defined. In Cython, we know
that at least one of ``left`` or ``right`` is an instance of the class
but we do not know a priori which one.

For addition and multiplication (not for other operators), there is a
fast path for operations with a Python ``int`` (which corresponds
to a C ``long``). Implement ``cdef _add_long(self, long n)`` or
``cdef _mul_long(self, long n)`` with optimized code for ``self + n``
or ``self * n``. These are assumed to be commutative, so they are also
called for ``n + self`` or ``n * self``.
From Cython code, you can also call ``_add_long`` or ``_mul_long``
directly.

Examples
^^^^^^^^

We need some :class:`Parent` to work with::

    sage: from sage.structure.parent import Parent
    sage: class ExampleParent(Parent):
    ....:     def __init__(self, name, **kwds):
    ....:         Parent.__init__(self, **kwds)
    ....:         self.rename(name)

We start with a very basic example of a Python class implementing
``_add_``::

    sage: from sage.structure.element import Element
    sage: class MyElement(Element):
    ....:     def _add_(self, other):
    ....:         return 42
    sage: p = ExampleParent("Some parent")
    sage: x = MyElement(p)
    sage: x + x
    42

When two different parents are involved, this no longer works since
there is no coercion::

    sage: q = ExampleParent("Other parent")
    sage: y = MyElement(q)
    sage: x + y
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '+': 'Some parent' and 'Other parent'

If ``_add_`` is not defined, an error message is raised, referring to
the parents::

    sage: x = Element(p)
    sage: x._add_(x)
    Traceback (most recent call last):
    ...
    AttributeError: 'sage.structure.element.Element' object has no attribute '_add_'
    sage: x + x
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '+': 'Some parent' and 'Some parent'
    sage: y = Element(q)
    sage: x + y
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '+': 'Some parent' and 'Other parent'

We can also implement arithmetic generically in categories::

    sage: class MyCategory(Category):
    ....:     def super_categories(self):
    ....:         return [Sets()]
    ....:     class ElementMethods:
    ....:         def _add_(self, other):
    ....:             return 42
    sage: p = ExampleParent("Parent in my category", category=MyCategory())
    sage: x = Element(p)
    sage: x + x
    42

Implementation details
^^^^^^^^^^^^^^^^^^^^^^

Implementing the above features actually takes a bit of magic. Casual
callers and implementers can safely ignore it, but here are the
details for the curious.

To achieve fast arithmetic, it is critical to have a fast path in Cython
to call the ``_add_`` method of a Cython object. So we would like
to declare ``_add_`` as a ``cpdef`` method of class :class:`Element`.
Remember however that the abstract classes coming
from categories come after :class:`Element` in the method resolution
order (or fake method resolution order in case of a Cython
class). Hence any generic implementation of ``_add_`` in such an
abstract class would in principle be shadowed by ``Element._add_``.
This is worked around by defining ``Element._add_`` as a ``cdef``
instead of a ``cpdef`` method. Concrete implementations in subclasses
should be ``cpdef`` or ``def`` methods.

Let us now see what happens upon evaluating ``x + y`` when ``x`` and ``y``
are instances of a class that does not implement ``_add_`` but where
``_add_`` is implemented in the category.
First, ``x.__add__(y)`` is called, where ``__add__`` is implemented
in :class:`Element`.
Assuming that ``x`` and ``y`` have the same parent, a Cython call to
``x._add_(y)`` will be done.
The latter is implemented to trigger a Python level call to ``x._add_(y)``
which will succeed as desired.

In case that Python code calls ``x._add_(y)`` directly,
``Element._add_`` will be invisible, and the method lookup will
continue down the MRO and find the ``_add_`` method in the category.
"""

#*****************************************************************************
#       Copyright (C) 2006-2016 ...
#       Copyright (C) 2016 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, division, print_function

from libc.limits cimport LONG_MAX, LONG_MIN

from cpython cimport *
from sage.ext.stdsage cimport *

from cpython.ref cimport PyObject

import types
cdef add, sub, mul, div, truediv, floordiv, mod
cdef iadd, isub, imul, idiv, itruediv, ifloordiv
from operator import (add, sub, mul, div, truediv, floordiv, mod,
        iadd, isub, imul, idiv, itruediv, ifloordiv)
cdef dict _coerce_op_symbols = dict(
        add='+', sub='-', mul='*', div='/', truediv='/', floordiv='//', mod='%',
        iadd='+', isub='-', imul='*', idiv='/', itruediv='/', ifloordiv='//')

cdef MethodType
from types import MethodType

from sage.structure.sage_object cimport rich_to_bool
from sage.structure.coerce cimport py_scalar_to_element
from sage.structure.parent cimport Parent
from sage.structure.misc import is_extension_type
from sage.structure.misc cimport getattr_from_other_class
from sage.misc.lazy_format import LazyFormat
from sage.misc import sageinspect
from sage.misc.classcall_metaclass cimport ClasscallMetaclass
from sage.misc.superseded import deprecated_function_alias
from sage.arith.numerical_approx cimport digits_to_bits
from sage.misc.decorators import sage_wraps


def make_element(_class, _dict, parent):
    """
    This function is only here to support old pickles.

    Pickling functionality is moved to Element.{__getstate__,__setstate__}
    functions.
    """
    from sage.misc.pickle_old import make_element_old
    return make_element_old(_class, _dict, parent)


cdef unary_op_exception(op, x):
    try:
        op = op.__name__
        op = _coerce_op_symbols[op]
    except (AttributeError, KeyError):
        pass
    px = parent(x)
    return TypeError(f"unsupported operand parent for {op}: '{px}'")


cdef bin_op_exception(op, x, y):
    try:
        op = op.__name__
        op = _coerce_op_symbols[op]
    except (AttributeError, KeyError):
        pass
    px = parent(x)
    py = parent(y)
    return TypeError(f"unsupported operand parent(s) for {op!r}: '{px}' and '{py}'")


def is_Element(x):
    """
    Return ``True`` if x is of type Element.

    EXAMPLES::

        sage: from sage.structure.element import is_Element
        sage: is_Element(2/3)
        True
        sage: is_Element(QQ^3)
        False
    """
    return isinstance(x, Element)


cdef class Element(SageObject):
    """
    Generic element of a structure. All other types of elements
    (:class:`RingElement`, :class:`ModuleElement`, etc)
    derive from this type.

    Subtypes must either call ``__init__()`` to set ``_parent``, or may
    set ``_parent`` themselves if that would be more efficient.

    .. automethod:: _cmp_
    .. automethod:: _richcmp_
    .. automethod:: __add__
    .. automethod:: __sub__
    .. automethod:: __neg__
    .. automethod:: __mul__
    .. automethod:: __div__
    .. automethod:: __truediv__
    .. automethod:: __floordiv__
    .. automethod:: __mod__
    """
    def __getmetaclass__(_):
        from sage.misc.inherit_comparison import InheritComparisonMetaclass
        return InheritComparisonMetaclass

    def __init__(self, parent):
        r"""
        INPUT:

        - ``parent`` - a SageObject
        """
        self._parent = parent

    def _set_parent(self, parent):
        r"""
        INPUT:

        - ``parent`` - a SageObject
        """
        self._parent = parent

    def __getattr__(self, name):
        """
        Lookup a method or attribute from the category abstract classes.

        Let ``P`` be a parent in a category ``C``. Usually the methods
        of ``C.element_class`` are made directly available to elements
        of ``P`` via standard class inheritance. This is not the case
        any more if the elements of ``P`` are instances of an
        extension type. See :class:`Category` for details.

        The purpose of this method is to emulate this inheritance: for
        ``e`` and element of ``P``, if an attribute or method
        ``e.foo`` is not found in the super classes of ``e``, it's
        looked up manually in ``C.element_class`` and bound to ``e``.

        .. NOTES::

            - The attribute or method is actually looked up in
              ``P._abstract_element_class``. In most cases this is
              just an alias for ``C.element_class``, but some parents,
              notably homsets, customizes this to let elements also
              inherit from other abstract classes. See
              :meth:`Parent._abstract_element_class` and
              :meth:`Homset._abstract_element_class` for details.

            - This mechanism may also enter into action when the
              category of `P` is refined on the fly, leaving
              previously constructed elements in an outdated element
              class.

              See :class:`~sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing_generic`
              for an example.

        EXAMPLES:

        We test that ``1`` (an instance of the extension type
        ``Integer``) inherits the methods from the categories of
        ``ZZ``, that is from ``CommutativeRings().element_class``::

            sage: 1.is_idempotent(), 2.is_idempotent()
            (True, False)

        This method is actually provided by the ``Magmas()`` super
        category of ``CommutativeRings()``::

            sage: 1.is_idempotent
            <bound method JoinCategory.element_class.is_idempotent of 1>
            sage: 1.is_idempotent.__module__
            'sage.categories.magmas'

        TESTS::

            sage: 1.blah_blah
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.rings.integer.Integer' object has no attribute 'blah_blah'
            sage: Semigroups().example().an_element().is_idempotent
            <bound method LeftZeroSemigroup_with_category.element_class.is_idempotent of 42>
            sage: Semigroups().example().an_element().blah_blah
            Traceback (most recent call last):
            ...
            AttributeError: 'LeftZeroSemigroup_with_category.element_class' object has no attribute 'blah_blah'
        """
        return self.getattr_from_category(name)

    cdef getattr_from_category(self, name):
        # Lookup a method or attribute from the category abstract classes.
        # See __getattr__ above for documentation.
        cdef Parent P = self._parent
        if P is None:
            # This is highly unlikely but we deal with it anyway...
            # Usually, this will just raise AttributeError in
            # getattr_from_other_class().
            cls = type
        else:
            cls = P._abstract_element_class
        return getattr_from_other_class(self, cls, name)

    def __dir__(self):
        """
        Let cat be the category of the parent of ``self``. This method
        emulates ``self`` being an instance of both ``Element`` and
        ``cat.element_class``, in that order, for attribute directory.

        EXAMPLES::

            sage: dir(1/2)
            ['N', ..., 'is_idempotent', 'is_integer', 'is_integral', ...]

        Caveat: dir on Integer's and some other extension types seem to ignore __dir__::

            sage: 1.__dir__()
            ['N', ..., 'is_idempotent', 'is_integer', 'is_integral', ...]
            sage: dir(1)         # todo: not implemented
            ['N', ..., 'is_idempotent', 'is_integer', 'is_integral', ...]
        """
        from .misc import dir_with_other_class
        return dir_with_other_class(self, self.parent().category().element_class)

    def _repr_(self):
        return "Generic element of a structure"

    def __getstate__(self):
        """
        Return a tuple describing the state of your object.

        This should return all information that will be required to unpickle
        the object. The functionality for unpickling is implemented in
        __setstate__().

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: i = ideal(x^2 - y^2 + 1)
            sage: i.__getstate__()
            (Monoid of ideals of Multivariate Polynomial Ring in x, y over Rational Field,
             {'_Ideal_generic__gens': (x^2 - y^2 + 1,),
              '_Ideal_generic__ring': Multivariate Polynomial Ring in x, y over Rational Field,
              '_gb_by_ordering': {}})
        """
        return (self._parent, self.__dict__)

    def __setstate__(self, state):
        """
        Initializes the state of the object from data saved in a pickle.

        During unpickling __init__ methods of classes are not called, the saved
        data is passed to the class via this function instead.

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: i = ideal(x); i
            Ideal (x) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: S.<x,y,z> = ZZ[]
            sage: i.__setstate__((R,{'_Ideal_generic__ring':S,'_Ideal_generic__gens': (S(x^2 - y^2 + 1),)}))
            sage: i
            Ideal (x^2 - y^2 + 1) of Multivariate Polynomial Ring in x, y, z over Integer Ring
        """
        self._set_parent(state[0])
        self.__dict__ = state[1]

    def __copy__(self):
        """
        Return a copy of ``self``.

        OUTPUT:

          - a new object which is a copy of ``self``.

        This implementation ensures that ``self.__dict__`` is properly copied
        when it exists (typically for instances of classes deriving from
        :class:`Element`).

        TESTS::

            sage: from sage.structure.element import Element
            sage: el = Element(parent = ZZ)
            sage: el1 = copy(el)
            sage: el1 is el
            False

            sage: class Demo(Element): pass
            sage: el = Demo(parent = ZZ)
            sage: el.x = [1,2,3]
            sage: el1 = copy(el)
            sage: el1 is el
            False
            sage: el1.__dict__ is el.__dict__
            False
        """
        cls = self.__class__
        cdef Element res = cls.__new__(cls)
        res._parent = self._parent
        try:
            D = self.__dict__
        except AttributeError:
            return res
        for k,v in D.iteritems():
            try:
                setattr(res, k, v)
            except AttributeError:
                pass
        return res

    def _im_gens_(self, codomain, im_gens):
        """
        Return the image of ``self`` in codomain under the map that sends
        the images of the generators of the parent of ``self`` to the
        tuple of elements of im_gens.
        """
        raise NotImplementedError

    cpdef base_extend(self, R):
        cdef Parent V
        V = self._parent.base_extend(R)
        return V.coerce(self)

    def base_ring(self):
        """
        Return the base ring of this element's parent (if that makes sense).

        TESTS::

            sage: QQ.base_ring()
            Rational Field
            sage: identity_matrix(3).base_ring()
            Integer Ring
        """
        return self._parent.base_ring()

    def category(self):
        from sage.categories.all import Elements
        return Elements(self._parent)


    def _test_category(self, **options):
        """
        Run generic tests on the method :meth:`.category`.

        See also: :class:`TestSuite`.

        EXAMPLES::

            sage: 3._test_category()

        Let us now write a broken :meth:`.category` method::

            sage: from sage.categories.examples.sets_cat import PrimeNumbers
            sage: class CCls(PrimeNumbers):
            ....:     def an_element(self):
            ....:         return 18
            sage: CC = CCls()
            sage: CC._test_an_element()
            Traceback (most recent call last):
            ...
            AssertionError: self.an_element() is not in self
        """
        from sage.categories.objects    import Objects
        tester = self._tester(**options)
        SageObject._test_category(self, tester = tester)
        category = self.category()
        # Tests that self inherits methods from the categories
        if not is_extension_type(self.__class__):
            # For usual Python classes, that should be done with
            # standard inheritance
            tester.assert_(isinstance(self, self.parent().category().element_class))
        else:
            # For extension types we just check that inheritance
            # occurs on a dummy attribute of Sets().ElementMethods
            tester.assert_(hasattr(self, "_dummy_attribute"))

    def _test_eq(self, **options):
        """
        Test that ``self`` is equal to ``self`` and different to ``None``.

        See also: :class:`TestSuite`.

        TESTS::

            sage: from sage.structure.element import Element
            sage: O = Element(Parent())
            sage: O._test_eq()

        Let us now write a broken class method::

            sage: class CCls(Element):
            ....:     def __eq__(self, other):
            ....:         return True
            sage: CCls(Parent())._test_eq()
            Traceback (most recent call last):
            ...
            AssertionError: broken equality: Generic element of a structure == None

        Let us now break inequality::

            sage: class CCls(Element):
            ....:     def __ne__(self, other):
            ....:         return True
            sage: CCls(Parent())._test_eq()
            Traceback (most recent call last):
            ...
            AssertionError: broken non-equality: Generic element of a structure != itself
        """
        tester = self._tester(**options)
        # We don't use assertEqual / assertNonEqual in order to be
        # 100% sure we indeed call the operators == and !=, whatever
        # the version of Python is (see #11236)
        tester.assertTrue(self == self,
                   LazyFormat("broken equality: %s == itself is False")%self)
        tester.assertFalse(self == None,
                   LazyFormat("broken equality: %s == None")%self)
        tester.assertFalse(self != self,
                           LazyFormat("broken non-equality: %s != itself")%self)
        tester.assertTrue(self != None,
                          LazyFormat("broken non-equality: %s is not != None")%self)

    def parent(self, x=None):
        """
        Return the parent of this element; or, if the optional argument x is
        supplied, the result of coercing x into the parent of this element.
        """
        if x is None:
            return self._parent
        else:
            return self._parent(x)


    def subs(self, in_dict=None, **kwds):
        """
        Substitutes given generators with given values while not touching
        other generators. This is a generic wrapper around ``__call__``.
        The syntax is meant to be compatible with the corresponding method
        for symbolic expressions.

        INPUT:

        - ``in_dict`` - (optional) dictionary of inputs

        - ``**kwds`` - named parameters

        OUTPUT:

        - new object if substitution is possible, otherwise self.

        EXAMPLES::

            sage: x, y = PolynomialRing(ZZ,2,'xy').gens()
            sage: f = x^2 + y + x^2*y^2 + 5
            sage: f((5,y))
            25*y^2 + y + 30
            sage: f.subs({x:5})
            25*y^2 + y + 30
            sage: f.subs(x=5)
            25*y^2 + y + 30
            sage: (1/f).subs(x=5)
            1/(25*y^2 + y + 30)
            sage: Integer(5).subs(x=4)
            5
        """
        if not hasattr(self, '__call__'):
            return self
        parent = self._parent
        try:
            ngens = parent.ngens()
        except (AttributeError, NotImplementedError, TypeError):
            return self
        variables=[]
        # use "gen" instead of "gens" as a ParentWithGens is not
        # required to have the latter
        for i in xrange(0,ngens):
            gen=parent.gen(i)
            if str(gen) in kwds:
                variables.append(kwds[str(gen)])
            elif in_dict and gen in in_dict:
                variables.append(in_dict[gen])
            else:
                variables.append(gen)
        return self(*variables)

    def numerical_approx(self, prec=None, digits=None, algorithm=None):
        """
        Return a numerical approximation of ``self`` with ``prec`` bits
        (or decimal ``digits``) of precision.

        No guarantee is made about the accuracy of the result.

        INPUT:

        - ``prec`` -- precision in bits

        - ``digits`` -- precision in decimal digits (only used if
          ``prec`` is not given)

        - ``algorithm`` -- which algorithm to use to compute this
          approximation (the accepted algorithms depend on the object)

        If neither ``prec`` nor ``digits`` is given, the default
        precision is 53 bits (roughly 16 digits).

        EXAMPLES::

            sage: (2/3).numerical_approx()
            0.666666666666667
            sage: pi.n(digits=10)
            3.141592654
            sage: pi.n(prec=20)
            3.1416

        TESTS:

        Check that :trac:`14778` is fixed::

            sage: (0).n(algorithm='foo')
            0.000000000000000

        The ``.N`` method is a deprecated alias::

            sage: 0.N()
            doctest:...: DeprecationWarning: N is deprecated. Please use n instead.
            See http://trac.sagemath.org/13055 for details.
            0.000000000000000
        """
        from sage.arith.numerical_approx import numerical_approx_generic
        if prec is None:
            prec = digits_to_bits(digits)
        return numerical_approx_generic(self, prec)

    def n(self, prec=None, digits=None, algorithm=None):
        """
        Alias for :meth:`numerical_approx`.

        EXAMPLES::

            sage: (2/3).n()
            0.666666666666667
        """
        return self.numerical_approx(prec, digits, algorithm)

    N = deprecated_function_alias(13055, n)

    def _mpmath_(self, prec=53, rounding=None):
        """
        Evaluates numerically and returns an mpmath number.
        Used as fallback for conversion by mpmath.mpmathify().

        .. NOTE::

            Currently, the rounding mode is ignored.

        EXAMPLES::

            sage: from sage.libs.mpmath.all import mp, mpmathify
            sage: mp.dps = 30
            sage: 25._mpmath_(53)
            mpf('25.0')
            sage: mpmathify(3+4*I)
            mpc(real='3.0', imag='4.0')
            sage: mpmathify(1+pi)
            mpf('4.14159265358979323846264338327933')
            sage: (1+pi)._mpmath_(10)
            mpf('4.140625')
            sage: (1+pi)._mpmath_(mp.prec)
            mpf('4.14159265358979323846264338327933')
        """
        return self.n(prec)._mpmath_(prec=prec)

    def substitute(self,in_dict=None,**kwds):
        """
        This is an alias for self.subs().

        INPUT:

        - ``in_dict`` - (optional) dictionary of inputs

        - ``**kwds``  - named parameters

        OUTPUT:

        - new object if substitution is possible, otherwise self.

        EXAMPLES::

            sage: x, y = PolynomialRing(ZZ,2,'xy').gens()
            sage: f = x^2 + y + x^2*y^2 + 5
            sage: f((5,y))
            25*y^2 + y + 30
            sage: f.substitute({x:5})
            25*y^2 + y + 30
            sage: f.substitute(x=5)
            25*y^2 + y + 30
            sage: (1/f).substitute(x=5)
            1/(25*y^2 + y + 30)
            sage: Integer(5).substitute(x=4)
            5
         """
        return self.subs(in_dict,**kwds)

    cpdef _act_on_(self, x, bint self_on_left):
        """
        Use this method to implement ``self`` acting on ``x``.

        Return None or raise a CoercionException if no
        such action is defined here.
        """
        return None

    cpdef _acted_upon_(self, x, bint self_on_left):
        """
        Use this method to implement ``self`` acted on by x.

        Return None or raise a CoercionException if no
        such action is defined here.
        """
        return None


    def __xor__(self, right):
        raise RuntimeError("Use ** for exponentiation, not '^', which means xor\n"+\
              "in Python, and has the wrong precedence.")

    def __pos__(self):
        return self

    def _coeff_repr(self, no_space=True):
        if self._is_atomic():
            s = repr(self)
        else:
            s = "(%s)"%repr(self)
        if no_space:
            return s.replace(' ','')
        return s

    def _latex_coeff_repr(self):
        try:
            s = self._latex_()
        except AttributeError:
            s = str(self)
        if self._is_atomic():
            return s
        else:
            return "\\left(%s\\right)"%s

    def _is_atomic(self):
        """
        Return ``True`` if and only if parenthesis are not required when
        *printing* out any of `x - s`, `x + s`, `x^s` and `x/s`.

        EXAMPLES::

            sage: n = 5; n._is_atomic()
            True
            sage: n = x+1; n._is_atomic()
            False
        """
        if self._parent._repr_option('element_is_atomic'):
            return True
        s = str(self)
        return s.find("+") == -1 and s.find("-") == -1 and s.find(" ") == -1

    def __nonzero__(self):
        r"""
        Return ``True`` if ``self`` does not equal self.parent()(0).

        Note that this is automatically called when converting to
        boolean, as in the conditional of an if or while statement.

        TESTS:
        Verify that :trac:`5185` is fixed.

        ::

            sage: v = vector({1: 1, 3: -1})
            sage: w = vector({1: -1, 3: 1})
            sage: v + w
            (0, 0, 0, 0)
            sage: (v+w).is_zero()
            True
            sage: bool(v+w)
            False
            sage: (v+w).__nonzero__()
            False
        """
        return self != self._parent.zero()

    def is_zero(self):
        """
        Return ``True`` if ``self`` equals ``self.parent()(0)``.

        The default implementation is to fall back to ``not
        self.__nonzero__``.

        .. WARNING::

            Do not re-implement this method in your subclass but
            implement ``__nonzero__`` instead.
        """
        return not self

    def __cmp__(self, other):
        """
        Compare ``left`` and ``right`` using the coercion framework.

        ``self`` and ``other`` are coerced to a common parent and then
        ``_cmp_`` and ``_richcmp_`` are tried.

        EXAMPLES:

        We create an ``Element`` class where we define ``_richcmp_``
        and check that comparison works::

            sage: cython('''
            ....: from sage.structure.sage_object cimport rich_to_bool
            ....: from sage.structure.element cimport Element
            ....: cdef class FloatCmp(Element):
            ....:     cdef float x
            ....:     def __init__(self, float v):
            ....:         self.x = v
            ....:     cpdef _richcmp_(self, other, int op):
            ....:         cdef float x1 = (<FloatCmp>self).x
            ....:         cdef float x2 = (<FloatCmp>other).x
            ....:         return rich_to_bool(op, (x1 > x2) - (x1 < x2) )
            ....: ''')
            sage: a = FloatCmp(1)
            sage: b = FloatCmp(2)
            sage: cmp(a, b)
            -1
            sage: b.__cmp__(a)
            1
            sage: a <= b, b <= a
            (True, False)

        This works despite ``_cmp_`` not being implemented::

            sage: a._cmp_(b)
            Traceback (most recent call last):
            ...
            NotImplementedError: comparison not implemented for <type '...FloatCmp'>
        """
        if have_same_parent(self, other):
            left = self
            right = other
        else:
            try:
                left, right = coercion_model.canonical_coercion(self, other)
            except TypeError:
                # Compare by id()
                if (<unsigned long><PyObject*>self) >= (<unsigned long><PyObject*>other):
                    # It cannot happen that self is other, since they don't
                    # have the same parent.
                    return 1
                else:
                    return -1

            if not isinstance(left, Element):
                assert type(left) is type(right)
                return cmp(left, right)

        # Now we have two Sage Elements with the same parent
        try:
            # First attempt: use _cmp_()
            return (<Element>left)._cmp_(<Element>right)
        except NotImplementedError:
            # Second attempt: use _richcmp_()
            if (<Element>left)._richcmp_(right, Py_EQ):
                return 0
            if (<Element>left)._richcmp_(right, Py_LT):
                return -1
            if (<Element>left)._richcmp_(right, Py_GT):
                return 1
            raise

    def _cache_key(self):
        """
        Provide a hashable key for an element if it is not hashable

        EXAMPLES::

            sage: a=sage.structure.element.Element(ZZ)
            sage: a._cache_key()
            (Integer Ring, 'Generic element of a structure')
        """
        return(self.parent(),str(self))

    ####################################################################
    # In a Cython or a Python class, you must define either _cmp_
    # (if your subclass is totally ordered), _richcmp_ (if your subclass
    # is partially ordered), or both (if your class has both a total order
    # and a partial order, or if that gives better performance).
    #
    # Rich comparisons (like a < b) will default to using _richcmp_,
    # three-way comparisons (like cmp(a,b)) will default to using
    # _cmp_. But if you define just one of _richcmp_ and _cmp_, it will
    # be used for all kinds of comparisons.
    #
    # In the _cmp_ and _richcmp_ methods, you can assume that both
    # arguments have identical parents.
    ####################################################################
    def __richcmp__(self, other, int op):
        """
        Compare ``self`` and ``other`` using the coercion framework,
        comparing according to the comparison operator ``op``.

        Normally, a class will not redefine ``__richcmp__`` but rely on
        this ``Element.__richcmp__`` method which uses coercion if
        needed to compare elements. After coercion (or if no coercion
        is needed), ``_richcmp_`` is called.

        If a class wants to implement rich comparison without coercion,
        then ``__richcmp__`` should be defined.
        See :class:`sage.numerical.linear_functions.LinearConstraint`
        for such an example.

        For efficiency reasons, a class can do certain "manual"
        coercions directly in ``__richcmp__``, using
        ``coercion_model.richcmp()`` for the remaining cases.
        This is done for example in :class:`Integer`.
        """
        if have_same_parent(self, other):
            # Same parents, in particular self and other must both be
            # an instance of Element. The explicit casts below make
            # Cython generate optimized code for this call.
            return (<Element>self)._richcmp_(other, op)
        else:
            return coercion_model.richcmp(self, other, op)

    cpdef _richcmp_(left, right, int op):
        r"""
        Default implementation of rich comparisons for elements with
        equal parents.

        It tries to see if ``_cmp_`` is implemented. Otherwise it does a
        comparison by id for ``==`` and ``!=``. Calling this default method
        with ``<``, ``<=``, ``>`` or ``>=`` will raise a
        ``NotImplementedError``.

        EXAMPLES::

            sage: from sage.structure.parent import Parent
            sage: from sage.structure.element import Element
            sage: P = Parent()
            sage: e1 = Element(P); e2 = Element(P)
            sage: e1 == e1    # indirect doctest
            True
            sage: e1 == e2    # indirect doctest
            False
            sage: e1 < e2     # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: comparison not implemented for <type 'sage.structure.element.Element'>
        """
        # Obvious case
        if left is right:
            return rich_to_bool(op, 0)

        cdef int c
        try:
            c = left._cmp_(right)
        except NotImplementedError:
            # Check equality by id(), knowing that left is not right
            if op == Py_EQ: return False
            if op == Py_NE: return True
            raise
        assert -1 <= c <= 1
        return rich_to_bool(op, c)

    cpdef int _cmp_(left, right) except -2:
        """
        Default three-way comparison method which only checks for a
        Python class defining ``__cmp__``.
        """
        left_cmp = left.__cmp__
        if isinstance(left_cmp, MethodType):
            return left_cmp(right)
        msg = LazyFormat("comparison not implemented for %r")%type(left)
        raise NotImplementedError(msg)

    ##################################################
    # Arithmetic using the coercion model
    ##################################################

    def __add__(left, right):
        """
        Top-level addition operator for :class:`Element` invoking
        the coercion model.

        See :ref:`element_arithmetic`.

        EXAMPLES::

            sage: from sage.structure.element import Element
            sage: class MyElement(Element):
            ....:     def _add_(self, other):
            ....:         return 42
            sage: e = MyElement(Parent())
            sage: e + e
            42

        TESTS::

            sage: e = Element(Parent())
            sage: e + e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '+': '<type 'sage.structure.parent.Parent'>' and '<type 'sage.structure.parent.Parent'>'
            sage: 1 + e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '+': 'Integer Ring' and '<type 'sage.structure.parent.Parent'>'
            sage: e + 1
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '+': '<type 'sage.structure.parent.Parent'>' and 'Integer Ring'
            sage: int(1) + e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for +: 'int' and 'sage.structure.element.Element'
            sage: e + int(1)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for +: 'sage.structure.element.Element' and 'int'
            sage: None + e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for +: 'NoneType' and 'sage.structure.element.Element'
            sage: e + None
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for +: 'sage.structure.element.Element' and 'NoneType'
        """
        cdef int cl = classify_elements(left, right)
        if HAVE_SAME_PARENT(cl):
            return (<Element>left)._add_(right)
        # Left and right are Sage elements => use coercion model
        if BOTH_ARE_ELEMENT(cl):
            return coercion_model.bin_op(left, right, add)

        try:
            # Special case addition with Python int
            if isinstance(right, int):
                return (<Element>left)._add_long(PyInt_AS_LONG(right))
            if isinstance(left, int):
                return (<Element>right)._add_long(PyInt_AS_LONG(left))
            return coercion_model.bin_op(left, right, add)
        except TypeError:
            # Either coercion failed or arithmetic is not defined.
            #
            # According to the Python convention, we should return
            # NotImplemented now. This will cause Python to try the
            # reversed addition (__radd__).
            return NotImplemented

    cdef _add_(self, other):
        """
        Virtual addition method for elements with identical parents.

        This default Cython implementation of ``_add_`` calls the
        Python method ``self._add_`` if it exists. This method may be
        defined in the ``ElementMethods`` of the category of the parent.
        If the method is not found, a ``TypeError`` is raised
        indicating that the operation is not supported.

        See :ref:`element_arithmetic`.

        EXAMPLES:

        This method is not visible from Python::

            sage: from sage.structure.element import Element
            sage: e = Element(Parent())
            sage: e._add_(e)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.structure.element.Element' object has no attribute '_add_'
        """
        try:
            python_op = (<object>self)._add_
        except AttributeError:
            raise bin_op_exception('+', self, other)
        else:
            return python_op(other)

    cdef _add_long(self, long n):
        """
        Generic path for adding a C long, assumed to commute.

        EXAMPLES::

            sage: cython(  # long time
            ....: '''
            ....: from sage.structure.element cimport Element
            ....: cdef class MyElement(Element):
            ....:     cdef _add_long(self, long n):
            ....:         return n
            ....: ''')
            sage: e = MyElement(Parent())  # long time
            sage: i = int(42)
            sage: i + e, e + i  # long time
            (42, 42)
        """
        return coercion_model.bin_op(self, n, add)

    def __sub__(left, right):
        """
        Top-level subtraction operator for :class:`Element` invoking
        the coercion model.

        See :ref:`element_arithmetic`.

        EXAMPLES::

            sage: from sage.structure.element import Element
            sage: class MyElement(Element):
            ....:     def _sub_(self, other):
            ....:         return 42
            sage: e = MyElement(Parent())
            sage: e - e
            42

        TESTS::

            sage: e = Element(Parent())
            sage: e - e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '-': '<type 'sage.structure.parent.Parent'>' and '<type 'sage.structure.parent.Parent'>'
            sage: 1 - e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '-': 'Integer Ring' and '<type 'sage.structure.parent.Parent'>'
            sage: e - 1
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '-': '<type 'sage.structure.parent.Parent'>' and 'Integer Ring'
            sage: int(1) - e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for -: 'int' and 'sage.structure.element.Element'
            sage: e - int(1)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for -: 'sage.structure.element.Element' and 'int'
            sage: None - e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for -: 'NoneType' and 'sage.structure.element.Element'
            sage: e - None
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for -: 'sage.structure.element.Element' and 'NoneType'
        """
        # See __add__ for comments
        cdef int cl = classify_elements(left, right)
        if HAVE_SAME_PARENT(cl):
            return (<Element>left)._sub_(right)
        if BOTH_ARE_ELEMENT(cl):
            return coercion_model.bin_op(left, right, sub)

        try:
            return coercion_model.bin_op(left, right, sub)
        except TypeError:
            return NotImplemented

    cdef _sub_(self, other):
        """
        Virtual subtraction method for elements with identical parents.

        This default Cython implementation of ``_sub_`` calls the
        Python method ``self._sub_`` if it exists. This method may be
        defined in the ``ElementMethods`` of the category of the parent.
        If the method is not found, a ``TypeError`` is raised
        indicating that the operation is not supported.

        See :ref:`element_arithmetic`.

        EXAMPLES:

        This method is not visible from Python::

            sage: from sage.structure.element import Element
            sage: e = Element(Parent())
            sage: e._sub_(e)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.structure.element.Element' object has no attribute '_sub_'
        """
        try:
            python_op = (<object>self)._sub_
        except AttributeError:
            raise bin_op_exception('-', self, other)
        else:
            return python_op(other)

    def __neg__(self):
        """
        Top-level negation operator for :class:`Element`.

        EXAMPLES::

            sage: from sage.structure.element import Element
            sage: class MyElement(Element):
            ....:     def _neg_(self):
            ....:         return 42
            sage: e = MyElement(Parent())
            sage: -e
            42

        TESTS::

            sage: e = Element(Parent())
            sage: -e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent for unary -: '<type 'sage.structure.parent.Parent'>'
        """
        return self._neg_()

    cdef _neg_(self):
        """
        Virtual unary negation method for elements.

        This default Cython implementation of ``_neg_`` calls the
        Python method ``self._neg_`` if it exists. This method may be
        defined in the ``ElementMethods`` of the category of the parent.
        If the method is not found, a ``TypeError`` is raised
        indicating that the operation is not supported.

        See :ref:`element_arithmetic`.

        EXAMPLES:

        This method is not visible from Python::

            sage: from sage.structure.element import Element
            sage: e = Element(Parent())
            sage: e._neg_()
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.structure.element.Element' object has no attribute '_neg_'
        """
        try:
            python_op = (<object>self)._neg_
        except AttributeError:
            raise unary_op_exception('unary -', self)
        else:
            return python_op()

    def __mul__(left, right):
        """
        Top-level multiplication operator for :class:`Element` invoking
        the coercion model.

        See :ref:`element_arithmetic`.

        EXAMPLES::

            sage: from sage.structure.element import Element
            sage: class MyElement(Element):
            ....:     def _mul_(self, other):
            ....:         return 42
            sage: e = MyElement(Parent())
            sage: e * e
            42

        TESTS::

            sage: e = Element(Parent())
            sage: e * e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': '<type 'sage.structure.parent.Parent'>' and '<type 'sage.structure.parent.Parent'>'
            sage: 1 * e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Integer Ring' and '<type 'sage.structure.parent.Parent'>'
            sage: e * 1
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': '<type 'sage.structure.parent.Parent'>' and 'Integer Ring'
            sage: int(1) * e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for *: 'int' and 'sage.structure.element.Element'
            sage: e * int(1)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for *: 'sage.structure.element.Element' and 'int'
            sage: None * e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for *: 'NoneType' and 'sage.structure.element.Element'
            sage: e * None
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for *: 'sage.structure.element.Element' and 'NoneType'

        ::

            sage: A = AlgebrasWithBasis(QQ).example(); A
            An example of an algebra with basis: the free algebra
            on the generators ('a', 'b', 'c') over Rational Field
            sage: x = A.an_element()
            sage: x
            B[word: ] + 2*B[word: a] + 3*B[word: b] + B[word: bab]
            sage: x.__mul__(x)
            B[word: ] + 4*B[word: a] + 4*B[word: aa] + 6*B[word: ab]
            + 2*B[word: abab] + 6*B[word: b] + 6*B[word: ba]
            + 2*B[word: bab] + 2*B[word: baba] + 3*B[word: babb]
            + B[word: babbab] + 9*B[word: bb] + 3*B[word: bbab]
        """
        cdef int cl = classify_elements(left, right)
        if HAVE_SAME_PARENT(cl):
            return (<Element>left)._mul_(right)
        if BOTH_ARE_ELEMENT(cl):
            return coercion_model.bin_op(left, right, mul)

        try:
            if isinstance(right, int):
                return (<Element>left)._mul_long(PyInt_AS_LONG(right))
            if isinstance(left, int):
                return (<Element>right)._mul_long(PyInt_AS_LONG(left))
            return coercion_model.bin_op(left, right, mul)
        except TypeError:
            return NotImplemented

    cdef _mul_(self, other):
        """
        Virtual multiplication method for elements with identical parents.

        This default Cython implementation of ``_mul_`` calls the
        Python method ``self._mul_`` if it exists. This method may be
        defined in the ``ElementMethods`` of the category of the parent.
        If the method is not found, a ``TypeError`` is raised
        indicating that the operation is not supported.

        See :ref:`element_arithmetic`.

        EXAMPLES:

        This method is not visible from Python::

            sage: from sage.structure.element import Element
            sage: e = Element(Parent())
            sage: e._mul_(e)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.structure.element.Element' object has no attribute '_mul_'
        """
        try:
            python_op = (<object>self)._mul_
        except AttributeError:
            raise bin_op_exception('*', self, other)
        else:
            return python_op(other)

    cdef _mul_long(self, long n):
        """
        Generic path for multiplying by a C long, assumed to commute.

        EXAMPLES::

            sage: cython(  # long time
            ....: '''
            ....: from sage.structure.element cimport Element
            ....: cdef class MyElement(Element):
            ....:     cdef _mul_long(self, long n):
            ....:         return n
            ....: ''')
            sage: e = MyElement(Parent())  # long time
            sage: i = int(42)
            sage: i * e, e * i  # long time
            (42, 42)
        """
        return coercion_model.bin_op(self, n, mul)

    def __div__(left, right):
        """
        Top-level division operator for :class:`Element` invoking
        the coercion model.

        See :ref:`element_arithmetic`.

        EXAMPLES::

            sage: 2 / 3
            2/3
            sage: pi / 3
            1/3*pi
            sage: K.<i> = NumberField(x^2+1)
            sage: 2 / K.ideal(i+1)
            Fractional ideal (-i + 1)

        ::

            sage: from sage.structure.element import Element
            sage: class MyElement(Element):
            ....:     def _div_(self, other):
            ....:         return 42
            sage: e = MyElement(Parent())
            sage: e / e
            42

        TESTS::

            sage: e = Element(Parent())
            sage: e / e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '/': '<type 'sage.structure.parent.Parent'>' and '<type 'sage.structure.parent.Parent'>'
            sage: 1 / e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '/': 'Integer Ring' and '<type 'sage.structure.parent.Parent'>'
            sage: e / 1
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '/': '<type 'sage.structure.parent.Parent'>' and 'Integer Ring'
            sage: int(1) / e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for /: 'int' and 'sage.structure.element.Element'
            sage: e / int(1)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for /: 'sage.structure.element.Element' and 'int'
            sage: None / e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for /: 'NoneType' and 'sage.structure.element.Element'
            sage: e / None
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for /: 'sage.structure.element.Element' and 'NoneType'
        """
        # See __add__ for comments
        cdef int cl = classify_elements(left, right)
        if HAVE_SAME_PARENT(cl):
            return (<Element>left)._div_(right)
        if BOTH_ARE_ELEMENT(cl):
            return coercion_model.bin_op(left, right, div)

        try:
            return coercion_model.bin_op(left, right, div)
        except TypeError:
            return NotImplemented

    def __truediv__(left, right):
        """
        Top-level true division operator for :class:`Element` invoking
        the coercion model.

        See :ref:`element_arithmetic`.

        EXAMPLES::

            sage: operator.truediv(2, 3)
            2/3
            sage: operator.truediv(pi, 3)
            1/3*pi
            sage: K.<i> = NumberField(x^2+1)
            sage: operator.truediv(2, K.ideal(i+1))
            Fractional ideal (-i + 1)

        ::

            sage: from sage.structure.element import Element
            sage: class MyElement(Element):
            ....:     def _div_(self, other):
            ....:         return 42
            sage: e = MyElement(Parent())
            sage: operator.truediv(e, e)
            42

        TESTS::

            sage: e = Element(Parent())
            sage: operator.truediv(e, e)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '/': '<type 'sage.structure.parent.Parent'>' and '<type 'sage.structure.parent.Parent'>'
            sage: operator.truediv(1, e)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '/': 'Integer Ring' and '<type 'sage.structure.parent.Parent'>'
            sage: operator.truediv(e, 1)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '/': '<type 'sage.structure.parent.Parent'>' and 'Integer Ring'
            sage: operator.truediv(int(1), e)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for /: 'int' and 'sage.structure.element.Element'
            sage: operator.truediv(e, int(1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for /: 'sage.structure.element.Element' and 'int'
            sage: operator.truediv(None, e)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for /: 'NoneType' and 'sage.structure.element.Element'
            sage: operator.truediv(e, None)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for /: 'sage.structure.element.Element' and 'NoneType'
        """
        # See __add__ for comments
        cdef int cl = classify_elements(left, right)
        if HAVE_SAME_PARENT(cl):
            return (<Element>left)._div_(right)
        if BOTH_ARE_ELEMENT(cl):
            return coercion_model.bin_op(left, right, truediv)

        try:
            return coercion_model.bin_op(left, right, truediv)
        except TypeError:
            return NotImplemented

    cdef _div_(self, other):
        """
        Virtual division method for elements with identical parents.
        This is called for Python 2 division as well as true division.

        This default Cython implementation of ``_div_`` calls the
        Python method ``self._div_`` if it exists. This method may be
        defined in the ``ElementMethods`` of the category of the parent.
        If the method is not found, a ``TypeError`` is raised
        indicating that the operation is not supported.

        See :ref:`element_arithmetic`.

        EXAMPLES:

        This method is not visible from Python::

            sage: from sage.structure.element import Element
            sage: e = Element(Parent())
            sage: e._div_(e)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.structure.element.Element' object has no attribute '_div_'
        """
        try:
            python_op = (<object>self)._div_
        except AttributeError:
            raise bin_op_exception('/', self, other)
        else:
            return python_op(other)

    def __floordiv__(left, right):
        """
        Top-level floor division operator for :class:`Element` invoking
        the coercion model.

        See :ref:`element_arithmetic`.

        EXAMPLES::

            sage: 7 // 3
            2
            sage: 7 // int(3)
            2
            sage: int(7) // 3
            2

        ::

            sage: from sage.structure.element import Element
            sage: class MyElement(Element):
            ....:     def _floordiv_(self, other):
            ....:         return 42
            sage: e = MyElement(Parent())
            sage: e // e
            42

        TESTS::

            sage: e = Element(Parent())
            sage: e // e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '//': '<type 'sage.structure.parent.Parent'>' and '<type 'sage.structure.parent.Parent'>'
            sage: 1 // e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '//': 'Integer Ring' and '<type 'sage.structure.parent.Parent'>'
            sage: e // 1
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '//': '<type 'sage.structure.parent.Parent'>' and 'Integer Ring'
            sage: int(1) // e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for //: 'int' and 'sage.structure.element.Element'
            sage: e // int(1)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for //: 'sage.structure.element.Element' and 'int'
            sage: None // e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for //: 'NoneType' and 'sage.structure.element.Element'
            sage: e // None
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for //: 'sage.structure.element.Element' and 'NoneType'
        """
        # See __add__ for comments
        cdef int cl = classify_elements(left, right)
        if HAVE_SAME_PARENT(cl):
            return (<Element>left)._floordiv_(right)
        if BOTH_ARE_ELEMENT(cl):
            return coercion_model.bin_op(left, right, floordiv)

        try:
            return coercion_model.bin_op(left, right, floordiv)
        except TypeError:
            return NotImplemented

    cdef _floordiv_(self, other):
        """
        Virtual floor division method for elements with identical parents.

        This default Cython implementation of ``_floordiv_`` calls the
        Python method ``self._floordiv_`` if it exists. This method may be
        defined in the ``ElementMethods`` of the category of the parent.
        If the method is not found, a ``TypeError`` is raised
        indicating that the operation is not supported.

        See :ref:`element_arithmetic`.

        EXAMPLES:

        This method is not visible from Python::

            sage: from sage.structure.element import Element
            sage: e = Element(Parent())
            sage: e._floordiv_(e)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.structure.element.Element' object has no attribute '_floordiv_'
        """
        try:
            python_op = (<object>self)._floordiv_
        except AttributeError:
            raise bin_op_exception('//', self, other)
        else:
            return python_op(other)

    def __mod__(left, right):
        """
        Top-level modulo operator for :class:`Element` invoking
        the coercion model.

        See :ref:`element_arithmetic`.

        EXAMPLES::

            sage: 7 % 3
            1
            sage: 7 % int(3)
            1
            sage: int(7) % 3
            1

        ::

            sage: from sage.structure.element import Element
            sage: class MyElement(Element):
            ....:     def _mod_(self, other):
            ....:         return 42
            sage: e = MyElement(Parent())
            sage: e % e
            42

        TESTS::

            sage: e = Element(Parent())
            sage: e % e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '%': '<type 'sage.structure.parent.Parent'>' and '<type 'sage.structure.parent.Parent'>'
            sage: 1 % e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '%': 'Integer Ring' and '<type 'sage.structure.parent.Parent'>'
            sage: e % 1
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '%': '<type 'sage.structure.parent.Parent'>' and 'Integer Ring'
            sage: int(1) % e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for %: 'int' and 'sage.structure.element.Element'
            sage: e % int(1)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for %: 'sage.structure.element.Element' and 'int'
            sage: None % e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for %: 'NoneType' and 'sage.structure.element.Element'
            sage: e % None
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for %: 'sage.structure.element.Element' and 'NoneType'
        """
        # See __add__ for comments
        cdef int cl = classify_elements(left, right)
        if HAVE_SAME_PARENT(cl):
            return (<Element>left)._mod_(right)
        if BOTH_ARE_ELEMENT(cl):
            return coercion_model.bin_op(left, right, mod)

        try:
            return coercion_model.bin_op(left, right, mod)
        except TypeError:
            return NotImplemented

    cdef _mod_(self, other):
        """
        Virtual modulo method for elements with identical parents.

        This default Cython implementation of ``_mod_`` calls the
        Python method ``self._mod_`` if it exists. This method may be
        defined in the ``ElementMethods`` of the category of the parent.
        If the method is not found, a ``TypeError`` is raised
        indicating that the operation is not supported.

        See :ref:`element_arithmetic`.

        EXAMPLES:

        This method is not visible from Python::

            sage: from sage.structure.element import Element
            sage: e = Element(Parent())
            sage: e._mod_(e)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.structure.element.Element' object has no attribute '_mod_'
        """
        try:
            python_op = (<object>self)._mod_
        except AttributeError:
            raise bin_op_exception('%', self, other)
        else:
            return python_op(other)


def is_ModuleElement(x):
    """
    Return ``True`` if x is of type ModuleElement.

    This is even faster than using isinstance inline.

    EXAMPLES::

        sage: from sage.structure.element import is_ModuleElement
        sage: is_ModuleElement(2/3)
        True
        sage: is_ModuleElement((QQ^3).0)
        True
        sage: is_ModuleElement('a')
        False
    """
    return isinstance(x, ModuleElement)

cdef class ElementWithCachedMethod(Element):
    r"""
    An element class that fully supports cached methods.

    NOTE:

    The :class:`~sage.misc.cachefunc.cached_method` decorator provides
    a convenient way to automatically cache the result of a computation.
    Since :trac:`11115`, the cached method decorator applied to a
    method without optional arguments is faster than a hand-written cache
    in Python, and a cached method without any arguments (except ``self``)
    is actually faster than a Python method that does nothing more but
    to return ``1``. A cached method can also be inherited from the parent
    or element class of a category.

    However, this holds true only if attribute assignment is supported.
    If you write an extension class in Cython that does not accept attribute
    assignment then a cached method inherited from the category will be
    slower (for :class:`~sage.structure.parent.Parent`) or the cache would
    even break (for :class:`Element`).

    This class should be used if you write an element class, can not provide
    it with attribute assignment, but want that it inherits a cached method
    from the category. Under these conditions, your class should inherit
    from this class rather than :class:`Element`. Then, the cache will work,
    but certainly slower than with attribute assignment. Lazy attributes
    work as well.

    EXAMPLE:

    We define three element extension classes. The first inherits from
    :class:`Element`, the second from this class, and the third simply
    is a Python class. We also define a parent class and, in Python, a
    category whose element and parent classes define cached methods.
    ::

        sage: cython_code = ["from sage.structure.element cimport Element, ElementWithCachedMethod",
        ....:     "cdef class MyBrokenElement(Element):",
        ....:     "    cdef public object x",
        ....:     "    def __init__(self,P,x):",
        ....:     "        self.x=x",
        ....:     "        Element.__init__(self,P)",
        ....:     "    def __neg__(self):",
        ....:     "        return MyBrokenElement(self.parent(),-self.x)",
        ....:     "    def _repr_(self):",
        ....:     "        return '<%s>'%self.x",
        ....:     "    def __hash__(self):",
        ....:     "        return hash(self.x)",
        ....:     "    cpdef int _cmp_(left, right) except -2:",
        ....:     "        return cmp(left.x,right.x)",
        ....:     "    def raw_test(self):",
        ....:     "        return -self",
        ....:     "cdef class MyElement(ElementWithCachedMethod):",
        ....:     "    cdef public object x",
        ....:     "    def __init__(self,P,x):",
        ....:     "        self.x=x",
        ....:     "        Element.__init__(self,P)",
        ....:     "    def __neg__(self):",
        ....:     "        return MyElement(self.parent(),-self.x)",
        ....:     "    def _repr_(self):",
        ....:     "        return '<%s>'%self.x",
        ....:     "    def __hash__(self):",
        ....:     "        return hash(self.x)",
        ....:     "    cpdef int _cmp_(left, right) except -2:",
        ....:     "        return cmp(left.x,right.x)",
        ....:     "    def raw_test(self):",
        ....:     "        return -self",
        ....:     "class MyPythonElement(MyBrokenElement): pass",
        ....:     "from sage.structure.parent cimport Parent",
        ....:     "cdef class MyParent(Parent):",
        ....:     "    Element = MyElement"]
        sage: cython('\n'.join(cython_code))
        sage: cython_code = ["from sage.all import cached_method, cached_in_parent_method, Category, Objects",
        ....:     "class MyCategory(Category):",
        ....:     "    @cached_method",
        ....:     "    def super_categories(self):",
        ....:     "        return [Objects()]",
        ....:     "    class ElementMethods:",
        ....:     "        @cached_method",
        ....:     "        def element_cache_test(self):",
        ....:     "            return -self",
        ....:     "        @cached_in_parent_method",
        ....:     "        def element_via_parent_test(self):",
        ....:     "            return -self",
        ....:     "    class ParentMethods:",
        ....:     "        @cached_method",
        ....:     "        def one(self):",
        ....:     "            return self.element_class(self,1)",
        ....:     "        @cached_method",
        ....:     "        def invert(self, x):",
        ....:     "            return -x"]
        sage: cython('\n'.join(cython_code))
        sage: C = MyCategory()
        sage: P = MyParent(category=C)
        sage: ebroken = MyBrokenElement(P,5)
        sage: e = MyElement(P,5)

    The cached methods inherited by ``MyElement`` works::

        sage: e.element_cache_test()
        <-5>
        sage: e.element_cache_test() is e.element_cache_test()
        True
        sage: e.element_via_parent_test()
        <-5>
        sage: e.element_via_parent_test() is e.element_via_parent_test()
        True

    The other element class can only inherit a
    ``cached_in_parent_method``, since the cache is stored in the
    parent. In fact, equal elements share the cache, even if they are
    of different types::

        sage: e == ebroken
        True
        sage: type(e) == type(ebroken)
        False
        sage: ebroken.element_via_parent_test() is e.element_via_parent_test()
        True

    However, the cache of the other inherited method breaks, although the method
    as such works::

        sage: ebroken.element_cache_test()
        <-5>
        sage: ebroken.element_cache_test() is ebroken.element_cache_test()
        False

    Since ``e`` and ``ebroken`` share the cache, when we empty it for one element
    it is empty for the other as well::

        sage: b = ebroken.element_via_parent_test()
        sage: e.element_via_parent_test.clear_cache()
        sage: b is ebroken.element_via_parent_test()
        False

    Note that the cache only breaks for elements that do no allow attribute assignment.
    A Python version of ``MyBrokenElement`` therefore allows for cached methods::

        sage: epython = MyPythonElement(P,5)
        sage: epython.element_cache_test()
        <-5>
        sage: epython.element_cache_test() is epython.element_cache_test()
        True

    """
    cdef getattr_from_category(self, name):
        """
        This getattr method ensures that cached methods and lazy attributes
        can be inherited from the element class of a category.

        .. NOTE::

            The use of cached methods is demonstrated in the main doc
            string of this class. Here, we demonstrate lazy
            attributes.

        EXAMPLE::

            sage: cython('''
            ....: from sage.structure.element cimport ElementWithCachedMethod
            ....: cdef class MyElement(ElementWithCachedMethod):
            ....:     cdef public object x
            ....:     def __init__(self, P, x):
            ....:         self.x = x
            ....:         ElementWithCachedMethod.__init__(self,P)
            ....:     def _repr_(self):
            ....:         return '<%s>'%self.x
            ....: from sage.structure.parent cimport Parent
            ....: cdef class MyParent(Parent):
            ....:     Element = MyElement
            ....: from sage.all import cached_method, lazy_attribute, Category, Objects
            ....: class MyCategory(Category):
            ....:     @cached_method
            ....:     def super_categories(self):
            ....:         return [Objects()]
            ....:     class ElementMethods:
            ....:         @lazy_attribute
            ....:         def my_lazy_attr(self):
            ....:             return 'lazy attribute of <%s>'%self.x
            ....: ''')
            sage: C = MyCategory()
            sage: P = MyParent(category=C)
            sage: e = MyElement(P,5)
            sage: e.my_lazy_attr
            'lazy attribute of <5>'
            sage: e.my_lazy_attr is e.my_lazy_attr
            True
        """
        try:
            return self.__cached_methods[name]
        except KeyError:
            attr = getattr_from_other_class(self,
                                        self._parent.category().element_class,
                                        name)
            self.__cached_methods[name] = attr
            return attr
        except TypeError:
            attr = getattr_from_other_class(self,
                                        self._parent.category().element_class,
                                        name)
            self.__cached_methods = {name : attr}
            return attr

cdef class ModuleElement(Element):
    """
    Generic element of a module.
    """
    cdef _add_long(self, long n):
        """
        Generic path for adding a C long, assumed to commute.
        """
        if n == 0:
            return self
        return coercion_model.bin_op(self, n, add)

    cpdef _sub_(self, other):
        """
        Default implementation of subtraction using addition and
        negation.
        """
        return self + (-other)

    cpdef _neg_(self):
        """
        Default implementation of negation using multiplication
        with -1.
        """
        return self._mul_long(-1)

    cdef _mul_long(self, long n):
        """
        Generic path for multiplying by a C long, assumed to commute.
        """
        if n == 1:
            return self
        return coercion_model.bin_op(self, n, mul)

    # rmul -- left * self
    cpdef _rmul_(self, RingElement left):
        """
        Reversed scalar multiplication for module elements with the
        module element on the right and the scalar on the left.

        By default, we assume commutativity and reverse the arguments.
        """
        return self._lmul_(left)

    # lmul -- self * right
    cpdef _lmul_(self, RingElement right):
        """
        Scalar multiplication for module elements with the module
        element on the left and the scalar on the right.

        Returning None indicates that this action is not implemented here.
        """
        return None

    ##################################################
    # Other properties
    ##################################################
    def order(self):              ### DO NOT OVERRIDE THIS!!! Instead, override additive_order.
        """
        Return the additive order of self.
        """
        return self.additive_order()

    def additive_order(self):
        """
        Return the additive order of self.
        """
        raise NotImplementedError

########################################################################
# Monoid
########################################################################

def is_MonoidElement(x):
    """
    Return ``True`` if x is of type MonoidElement.
    """
    return isinstance(x, MonoidElement)

cdef class MonoidElement(Element):
    """
    Generic element of a monoid.
    """

    #############################################################
    # Other generic functions that should be available to
    # any monoid element.
    #############################################################
    def order(self):
        """
        Return the multiplicative order of self.
        """
        return self.multiplicative_order()

    def multiplicative_order(self):
        """
        Return the multiplicative order of self.
        """
        raise NotImplementedError

    def __pow__(self, n, dummy):
        """
        Return the (integral) power of self.
        """
        if dummy is not None:
            raise RuntimeError("__pow__ dummy argument not used")
        return generic_power_c(self,n,None)

    def powers(self, n):
        r"""
        Return the list `[x^0, x^1, \ldots, x^{n-1}]`.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: g = G([2, 3, 4, 1])
            sage: g.powers(4)
            [(), (1,2,3,4), (1,3)(2,4), (1,4,3,2)]
        """
        if n < 0:
            raise ValueError("negative number of powers requested")
        elif n == 0:
            return []
        x = self._parent.one()
        l = [x]
        for i in xrange(n - 1):
            x = x * self
            l.append(x)
        return l

    def __nonzero__(self):
        return True

def is_AdditiveGroupElement(x):
    """
    Return ``True`` if x is of type AdditiveGroupElement.
    """
    return isinstance(x, AdditiveGroupElement)

cdef class AdditiveGroupElement(ModuleElement):
    """
    Generic element of an additive group.
    """
    def order(self):
        """
        Return additive order of element
        """
        return self.additive_order()

    def __invert__(self):
        raise NotImplementedError("multiplicative inverse not defined for additive group elements")


def is_MultiplicativeGroupElement(x):
    """
    Return ``True`` if x is of type MultiplicativeGroupElement.
    """
    return isinstance(x, MultiplicativeGroupElement)

cdef class MultiplicativeGroupElement(MonoidElement):
    """
    Generic element of a multiplicative group.
    """
    def order(self):
        """
        Return the multiplicative order of self.
        """
        return self.multiplicative_order()

    cpdef _div_(self, right):
        """
        Default implementation of division using multiplication by
        the inverse.
        """
        return self * ~right

    def __invert__(self):
        r"""
        Return the inverse of ``self``.
        """
        return self._parent.one() / self


def is_RingElement(x):
    """
    Return ``True`` if x is of type RingElement.
    """
    return isinstance(x, RingElement)

cdef class RingElement(ModuleElement):
    def is_one(self):
        return self == self._parent.one()

    def __pow__(self, n, dummy):
        """
        Return the (integral) power of self.

        EXAMPLE::

            sage: a = Integers(389)['x']['y'](37)
            sage: p = sage.structure.element.RingElement.__pow__
            sage: p(a,2)
            202
            sage: p(a,2,1)
            Traceback (most recent call last):
            ...
            RuntimeError: __pow__ dummy argument not used
            sage: p(a,388)
            1
            sage: p(a,2^120)
            81
            sage: p(a,0)
            1
            sage: p(a,1) == a
            True
            sage: p(a,2) * p(a,3) == p(a,5)
            True
            sage: p(a,3)^2 == p(a,6)
            True
            sage: p(a,57) * p(a,43) == p(a,100)
            True
            sage: p(a,-1) == 1/a
            True
            sage: p(a,200) * p(a,-64) == p(a,136)
            True
            sage: p(2, 1/2)
            Traceback (most recent call last):
            ...
            NotImplementedError: non-integral exponents not supported

        TESTS::

        These aren't testing this code, but they are probably good to have around::

            sage: 2r**(SR(2)-1-1r)
            1
            sage: 2r^(1/2)
            sqrt(2)

        Exponent overflow should throw an OverflowError (:trac:`2956`)::

            sage: K.<x,y> = AA[]
            sage: x^(2^64 + 12345)
            Traceback (most recent call last):
            ...
            OverflowError: Exponent overflow (2147483648).

        Another example from :trac:`2956` which always overflows
        with Singular 4::

            sage: K.<x,y> = ZZ[]
            sage: (x^12345)^54321
            Traceback (most recent call last):
            ...
            OverflowError: exponent overflow (670592745)

        """
        if dummy is not None:
            raise RuntimeError("__pow__ dummy argument not used")
        return generic_power_c(self,n,None)

    def powers(self, n):
        r"""
        Return the list `[x^0, x^1, \ldots, x^{n-1}]`.

        EXAMPLES::

            sage: 5.powers(3)
            [1, 5, 25]
        """
        if n < 0:
            raise ValueError("negative number of powers requested")
        elif n == 0:
            return []
        x = self._parent.one()
        l = [x]
        for i in xrange(n - 1):
            x = x * self
            l.append(x)
        return l

    cpdef _div_(self, other):
        """
        Default implementation of division using the fraction field.
        """
        try:
            frac = self._parent.fraction_field()
        except AttributeError:
            raise bin_op_exception('/', self, other)
        return frac(self, other)

    def __invert__(self):
        return self._parent.one() / self

    def additive_order(self):
        """
        Return the additive order of ``self``.
        """
        raise NotImplementedError

    def multiplicative_order(self):
        r"""
        Return the multiplicative order of ``self``, if ``self`` is a unit,
        or raise ``ArithmeticError`` otherwise.
        """
        if not self.is_unit():
            raise ArithmeticError("self (=%s) must be a unit to have a multiplicative order.")
        raise NotImplementedError

    def is_nilpotent(self):
        """
        Return ``True`` if ``self`` is nilpotent, i.e., some power of ``self``
        is 0.

        TESTS::

            sage: a = QQ(2)
            sage: a.is_nilpotent()
            False
            sage: a = QQ(0)
            sage: a.is_nilpotent()
            True
            sage: m = matrix(QQ,3,[[3,2,3],[9,0,3],[-9,0,-3]])
            sage: m.is_nilpotent()
            Traceback (most recent call last):
            ...
            AttributeError: ... object has no attribute 'is_nilpotent'
        """
        if self.is_unit():
            return False
        if self.is_zero():
            return True
        raise NotImplementedError

    def abs(self):
        """
        Return the absolute value of ``self``.  (This just calls the ``__abs__``
        method, so it is equivalent to the ``abs()`` built-in function.)

        EXAMPLES::

            sage: RR(-1).abs()
            1.00000000000000
            sage: ZZ(-1).abs()
            1
            sage: CC(I).abs()
            1.00000000000000
            sage: Mod(-15, 37).abs()
            Traceback (most recent call last):
            ...
            ArithmeticError: absolute valued not defined on integers modulo n.
        """
        return abs(self)

    def is_prime(self):
        """
        Is ``self`` a prime element?

        A *prime* element is a non-zero, non-unit element `p` such that,
        whenever `p` divides `ab` for some `a` and `b`, then `p`
        divides `a` or `p` divides `b`.

        EXAMPLES:

        For polynomial rings, prime is the same as irreducible::

            sage: R.<x,y> = QQ[]
            sage: x.is_prime()
            True
            sage: (x^2 + y^3).is_prime()
            True
            sage: (x^2 - y^2).is_prime()
            False
            sage: R(0).is_prime()
            False
            sage: R(2).is_prime()
            False

        For the Gaussian integers::

            sage: K.<i> = QuadraticField(-1)
            sage: ZI = K.ring_of_integers()
            sage: ZI(3).is_prime()
            True
            sage: ZI(5).is_prime()
            False
            sage: ZI(2+i).is_prime()
            True
            sage: ZI(0).is_prime()
            False
            sage: ZI(1).is_prime()
            False

        In fields, an element is never prime::

            sage: RR(0).is_prime()
            False
            sage: RR(2).is_prime()
            False

        For integers, prime numbers are redefined to be positive::

            sage: RingElement.is_prime(-2)
            True
            sage: Integer.is_prime(-2)
            False
        """
        if not self:  # We exclude the 0 element
            return False
        return self._parent.ideal(self).is_prime()


def is_CommutativeRingElement(x):
    """
    Return ``True`` if x is of type CommutativeRingElement.

    TESTS::

        sage: from sage.structure.element import is_CommutativeRingElement
        sage: is_CommutativeRingElement(oo)
        False

        sage: is_CommutativeRingElement(1)
        True
    """
    return isinstance(x, CommutativeRingElement)

cdef class CommutativeRingElement(RingElement):
    """
    Base class for elements of commutative rings.
    """
    def inverse_mod(self, I):
        r"""
        Return an inverse of ``self`` modulo the ideal `I`, if defined,
        i.e., if `I` and ``self`` together generate the unit ideal.
        """
        raise NotImplementedError

    def divides(self, x):
        """
        Return ``True`` if ``self`` divides x.

        EXAMPLES::

            sage: P.<x> = PolynomialRing(QQ)
            sage: x.divides(x^2)
            True
            sage: x.divides(x^2+2)
            False
            sage: (x^2+2).divides(x)
            False
            sage: P.<x> = PolynomialRing(ZZ)
            sage: x.divides(x^2)
            True
            sage: x.divides(x^2+2)
            False
            sage: (x^2+2).divides(x)
            False

        :trac:`5347` has been fixed::

            sage: K = GF(7)
            sage: K(3).divides(1)
            True
            sage: K(3).divides(K(1))
            True

        ::

            sage: R = Integers(128)
            sage: R(0).divides(1)
            False
            sage: R(0).divides(0)
            True
            sage: R(0).divides(R(0))
            True
            sage: R(1).divides(0)
            True
            sage: R(121).divides(R(120))
            True
            sage: R(120).divides(R(121))
            Traceback (most recent call last):
            ...
            ArithmeticError: reduction modulo 120 not defined

        If ``x`` has different parent than ``self``, they are first coerced to a
        common parent if possible. If this coercion fails, it returns a
        TypeError. This fixes :trac:`5759`. ::

            sage: Zmod(2)(0).divides(Zmod(2)(0))
            True
            sage: Zmod(2)(0).divides(Zmod(2)(1))
            False
            sage: Zmod(5)(1).divides(Zmod(2)(1))
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: 'Ring of integers modulo 5' and 'Ring of integers modulo 2'
            sage: Zmod(35)(4).divides(Zmod(7)(1))
            True
            sage: Zmod(35)(7).divides(Zmod(7)(1))
            False
        """
        if have_same_parent(self, x):
            # First we test some generic conditions:
            try:
                if x.is_zero():
                    return True # everything divides 0
            except (AttributeError, NotImplementedError):
                pass

            try:
                if self.is_zero():
                    return False # 0 divides nothing else
            except (AttributeError, NotImplementedError):
                pass

            try:
                if self.is_unit():
                    return True # units divide everything
            except (AttributeError, NotImplementedError):
                pass

            try:
                if self.is_one():
                    return True # 1 divides everything
                                # (is_unit() may not be implemented)
            except (AttributeError, NotImplementedError):
                pass

            return (x % self) == 0

        else:
            #Different parents, use coercion
            a, b = coercion_model.canonical_coercion(self, x)
            return a.divides(b)

    def mod(self, I):
        r"""
        Return a representative for ``self`` modulo the ideal I (or the ideal
        generated by the elements of I if I is not an ideal.)

        EXAMPLE:  Integers
        Reduction of 5 modulo an ideal::

            sage: n = 5
            sage: n.mod(3*ZZ)
            2

        Reduction of 5 modulo the ideal generated by 3::

            sage: n.mod(3)
            2

        Reduction of 5 modulo the ideal generated by 15 and 6, which is `(3)`.

        ::

            sage: n.mod([15,6])
            2


        EXAMPLE: Univariate polynomials

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^3 + x + 1
            sage: f.mod(x + 1)
            -1

        Reduction for `\ZZ[x]`::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = x^3 + x + 1
            sage: f.mod(x + 1)
            -1

        When little is implemented about a given ring, then mod may
        return simply return `f`.

        EXAMPLE: Multivariate polynomials
        We reduce a polynomial in two variables modulo a polynomial
        and an ideal::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3)
            sage: (x^2 + y^2 + z^2).mod(x+y+z)
            2*y^2 + 2*y*z + 2*z^2

        Notice above that `x` is eliminated.  In the next example,
        both `y` and `z` are eliminated::

            sage: (x^2 + y^2 + z^2).mod( (x - y, y - z) )
            3*z^2
            sage: f = (x^2 + y^2 + z^2)^2; f
            x^4 + 2*x^2*y^2 + y^4 + 2*x^2*z^2 + 2*y^2*z^2 + z^4
            sage: f.mod( (x - y, y - z) )
            9*z^4

        In this example `y` is eliminated::

            sage: (x^2 + y^2 + z^2).mod( (x^3, y - z) )
            x^2 + 2*z^2
        """
        from sage.rings.ideal import is_Ideal
        if not is_Ideal(I) or not I.ring() is self._parent:
            I = self._parent.ideal(I)
            #raise TypeError, "I = %s must be an ideal in %s"%(I, self.parent())
        return I.reduce(self)

    ##################################################
    # square roots
    ##################################################

    def is_square(self, root=False):
        """
        Return whether or not the ring element ``self`` is a square.

        If the optional argument root is ``True``, then also return
        the square root (or ``None``, if it is not a square).

        INPUT:

        - ``root`` - whether or not to also return a square
          root (default: ``False``)

        OUTPUT:

        - ``bool`` -- whether or not a square

        - ``object`` -- (optional) an actual square root if
          found, and ``None`` otherwise.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = 12*(x+1)^2 * (x+3)^2
            sage: f.is_square()
            False
            sage: f.is_square(root=True)
            (False, None)
            sage: h = f/3
            sage: h.is_square()
            True
            sage: h.is_square(root=True)
            (True, 2*x^2 + 8*x + 6)

        .. NOTE::

            This is the is_square implementation for general
            commutative ring elements. It's implementation is to raise
            a NotImplementedError. The function definition is here to
            show what functionality is expected and provide a general
            framework.
        """
        raise NotImplementedError("is_square() not implemented for elements of %s" % self.parent())

    def sqrt(self, extend=True, all=False, name=None):
        """
        It computes the square root.

        INPUT:

        -  ``extend`` - Whether to make a ring extension containing a square root if ``self`` is not a square (default: ``True``)

        -  ``all`` - Whether to return a list of all square roots or just a square root (default: False)

        -  ``name`` - Required when ``extend=True`` and ``self`` is not a square. This will be the name of the generator extension.

        OUTPUT:

        - if ``all=False`` it returns a square root. (throws an error if ``extend=False`` and ``self`` is not a square)

        - if ``all=True`` it returns a list of all the square roots (could be empty if ``extend=False`` and ``self`` is not a square)

        ALGORITHM:

        It uses ``is_square(root=true)`` for the hard part of the work, the rest is just wrapper code.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: (x^2).sqrt()
            x
            sage: f=x^2-4*x+4; f.sqrt(all=True)
            [x - 2, -x + 2]
            sage: sqrtx=x.sqrt(name="y"); sqrtx
            y
            sage: sqrtx^2
            x
            sage: x.sqrt(all=true,name="y")
            [y, -y]
            sage: x.sqrt(extend=False,all=True)
            []
            sage: x.sqrt()
            Traceback (most recent call last):
            ...
            TypeError: Polynomial is not a square. You must specify the name of the square root when using the default extend = True
            sage: x.sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: trying to take square root of non-square x with extend = False

        TESTS::

            sage: f = (x+3)^2; f.sqrt()
            x + 3
            sage: f = (x+3)^2; f.sqrt(all=True)
            [x + 3, -x - 3]
            sage: f = (x^2 - x + 3)^2; f.sqrt()
            x^2 - x + 3
            sage: f = (x^2 - x + 3)^6; f.sqrt()
            x^6 - 3*x^5 + 12*x^4 - 19*x^3 + 36*x^2 - 27*x + 27
            sage: g = (R.random_element(15))^2
            sage: g.sqrt()^2 == g
            True

            sage: R.<x> = GF(250037)[]
            sage: f = x^2/(x+1)^2; f.sqrt()
            x/(x + 1)
            sage: f = 9 * x^4 / (x+1)^2; f.sqrt()
            3*x^2/(x + 1)
            sage: f = 9 * x^4 / (x+1)^2; f.sqrt(all=True)
            [3*x^2/(x + 1), 250034*x^2/(x + 1)]

            sage: R.<x> = QQ[]
            sage: a = 2*(x+1)^2 / (2*(x-1)^2); a.sqrt()
            (2*x + 2)/(2*x - 2)
            sage: sqrtx=(1/x).sqrt(name="y"); sqrtx
            y
            sage: sqrtx^2
            1/x
            sage: (1/x).sqrt(all=true,name="y")
            [y, -y]
            sage: (1/x).sqrt(extend=False,all=True)
            []
            sage: (1/(x^2-1)).sqrt()
            Traceback (most recent call last):
            ...
            TypeError: Polynomial is not a square. You must specify the name of the square root when using the default extend = True
            sage: (1/(x^2-3)).sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: trying to take square root of non-square 1/(x^2 - 3) with extend = False
        """
        #This code is very general, it works for all integral domains that have the
        #is_square(root = True) option

        from sage.rings.ring import IntegralDomain
        P = self._parent
        is_sqr, sq_rt = self.is_square(root=True)
        if is_sqr:
            if all:
                if not isinstance(P, IntegralDomain):
                    raise NotImplementedError('sqrt() with all=True is only implemented for integral domains, not for %s' % P)
                if P.characteristic()==2 or sq_rt==0:
                    #0 has only one square root, and in charasteristic 2 everything also has only 1 root
                    return [ sq_rt ]
                return [ sq_rt, -sq_rt ]
            return sq_rt
        #from now on we know that self is not a square
        if not isinstance(P, IntegralDomain):
            raise NotImplementedError('sqrt() of non squares is only implemented for integral domains, not for %s' % P)
        if not extend:
            #all square roots of a non-square should be an empty list
            if all:
                return []
            raise ValueError('trying to take square root of non-square %s with extend = False' % self)

        if name is None:
            raise TypeError("Polynomial is not a square. You must specify the name of the square root when using the default extend = True")
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        PY = PolynomialRing(P, 'y')
        y = PY.gen()
        sq_rt = PY.quotient(y**2-self, names = name)(y)
        if all:
            if P.characteristic() == 2:
                return [ sq_rt ]
            return [ sq_rt, -sq_rt ]
        return sq_rt

    ##############################################

cdef class Vector(ModuleElement):
    cdef bint is_sparse_c(self):
        raise NotImplementedError

    cdef bint is_dense_c(self):
        raise NotImplementedError

    def __mul__(left, right):
        """
        Multiplication of vector by vector, matrix, or scalar

        AUTHOR:

        - Gonzalo Tornaria (2007-06-21) - write test cases and fix them

        .. NOTE::

            scalar * vector is implemented (and tested) in class RingElement
            matrix * vector is implemented (and tested) in class Matrix

        TESTS:

        Here we test (vector * vector) multiplication::

            sage: parent(vector(ZZ,[1,2])*vector(ZZ,[1,2]))
            Integer Ring
            sage: parent(vector(ZZ,[1,2])*vector(QQ,[1,2]))
            Rational Field
            sage: parent(vector(QQ,[1,2])*vector(ZZ,[1,2]))
            Rational Field
            sage: parent(vector(QQ,[1,2])*vector(QQ,[1,2]))
            Rational Field

            sage: parent(vector(QQ,[1,2,3,4])*vector(ZZ['x'],[1,2,3,4]))
            Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x'],[1,2,3,4])*vector(QQ,[1,2,3,4]))
            Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ,[1,2,3,4])*vector(ZZ['x']['y'],[1,2,3,4]))
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x']['y'],[1,2,3,4])*vector(QQ,[1,2,3,4]))
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ['x'],[1,2,3,4])*vector(ZZ['x']['y'],[1,2,3,4]))
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x']['y'],[1,2,3,4])*vector(QQ['x'],[1,2,3,4]))
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ['y'],[1,2,3,4])*vector(ZZ['x']['y'],[1,2,3,4]))
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x']['y'],[1,2,3,4])*vector(QQ['y'],[1,2,3,4]))
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(ZZ['x'],[1,2,3,4])*vector(ZZ['y'],[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(ZZ['x'],[1,2,3,4])*vector(QQ['y'],[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'
            sage: parent(vector(QQ['x'],[1,2,3,4])*vector(ZZ['y'],[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(QQ['x'],[1,2,3,4])*vector(QQ['y'],[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'

        Here we test (vector * matrix) multiplication::

            sage: parent(vector(ZZ,[1,2])*matrix(ZZ,2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: parent(vector(QQ,[1,2])*matrix(ZZ,2,2,[1,2,3,4]))
            Vector space of dimension 2 over Rational Field
            sage: parent(vector(ZZ,[1,2])*matrix(QQ,2,2,[1,2,3,4]))
            Vector space of dimension 2 over Rational Field
            sage: parent(vector(QQ,[1,2])*matrix(QQ,2,2,[1,2,3,4]))
            Vector space of dimension 2 over Rational Field

            sage: parent(vector(QQ,[1,2])*matrix(ZZ['x'],2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x'],[1,2])*matrix(QQ,2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ,[1,2])*matrix(ZZ['x']['y'],2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x']['y'],[1,2])*matrix(QQ,2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ['x'],[1,2])*matrix(ZZ['x']['y'],2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x']['y'],[1,2])*matrix(QQ['x'],2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ['y'],[1,2])*matrix(ZZ['x']['y'],2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x']['y'],[1,2])*matrix(QQ['y'],2,2,[1,2,3,4]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(ZZ['x'],[1,2])*matrix(ZZ['y'],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(ZZ['x'],[1,2])*matrix(QQ['y'],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field'
            sage: parent(vector(QQ['x'],[1,2])*matrix(ZZ['y'],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(QQ['x'],[1,2])*matrix(QQ['y'],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field'

        Here we test (vector * scalar) multiplication::

            sage: parent(vector(ZZ,[1,2])*ZZ(1))
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: parent(vector(QQ,[1,2])*ZZ(1))
            Vector space of dimension 2 over Rational Field
            sage: parent(vector(ZZ,[1,2])*QQ(1))
            Vector space of dimension 2 over Rational Field
            sage: parent(vector(QQ,[1,2])*QQ(1))
            Vector space of dimension 2 over Rational Field

            sage: parent(vector(QQ,[1,2])*ZZ['x'](1))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x'],[1,2])*QQ(1))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ,[1,2])*ZZ['x']['y'](1))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x']['y'],[1,2])*QQ(1))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ['x'],[1,2])*ZZ['x']['y'](1))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x']['y'],[1,2])*QQ['x'](1))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ['y'],[1,2])*ZZ['x']['y'](1))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ['x']['y'],[1,2])*QQ['y'](1))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(ZZ['x'],[1,2])*ZZ['y'](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(ZZ['x'],[1,2])*QQ['y'](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Univariate Polynomial Ring in y over Rational Field'
            sage: parent(vector(QQ['x'],[1,2])*ZZ['y'](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(QQ['x'],[1,2])*QQ['y'](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in y over Rational Field'

        Here we test (scalar * vector) multiplication::

            sage: parent(ZZ(1)*vector(ZZ,[1,2]))
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: parent(QQ(1)*vector(ZZ,[1,2]))
            Vector space of dimension 2 over Rational Field
            sage: parent(ZZ(1)*vector(QQ,[1,2]))
            Vector space of dimension 2 over Rational Field
            sage: parent(QQ(1)*vector(QQ,[1,2]))
            Vector space of dimension 2 over Rational Field

            sage: parent(QQ(1)*vector(ZZ['x'],[1,2]))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ['x'](1)*vector(QQ,[1,2]))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field

            sage: parent(QQ(1)*vector(ZZ['x']['y'],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ['x']['y'](1)*vector(QQ,[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(QQ['x'](1)*vector(ZZ['x']['y'],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ['x']['y'](1)*vector(QQ['x'],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(QQ['y'](1)*vector(ZZ['x']['y'],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ['x']['y'](1)*vector(QQ['y'],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(ZZ['x'](1)*vector(ZZ['y'],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(ZZ['x'](1)*vector(QQ['y'],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'
            sage: parent(QQ['x'](1)*vector(ZZ['y'],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(QQ['x'](1)*vector(QQ['y'],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'
        """
        if have_same_parent(left, right):
            return (<Vector>left)._dot_product_(<Vector>right)
        return coercion_model.bin_op(left, right, mul)

    cpdef _dot_product_(Vector left, Vector right):
        return left._dot_product_coerce_(right)

    cpdef _dot_product_coerce_(Vector left, Vector right):
        raise bin_op_exception('*', left, right)

    cpdef _pairwise_product_(Vector left, Vector right):
        raise TypeError("unsupported operation for '%s' and '%s'"%(parent(left), parent(right)))

    def __div__(self, other):
        return self / other

    def __truediv__(self, right):
        right = py_scalar_to_element(right)
        if isinstance(right, RingElement):
            # Let __mul__ do the job
            return self * ~right
        if isinstance(right, Vector):
            try:
                W = (<Vector>right)._parent.submodule([right])
                return W.coordinates(self)[0] / W.coordinates(right)[0]
            except ArithmeticError:
                if right.is_zero():
                    raise ZeroDivisionError("division by zero vector")
                else:
                    raise ArithmeticError("vector is not in free module")
        raise bin_op_exception('/', self, right)

    def _magma_init_(self, magma):
        """
        Return string that evaluates in Magma to something equivalent
        to this vector.

        EXAMPLES::

            sage: v = vector([1,2,3])
            sage: v._magma_init_(magma)                 # optional - magma
            '_sage_[...]![1,2,3]'
            sage: mv = magma(v); mv                     # optional - magma
            (1 2 3)
            sage: mv.Type()                             # optional - magma
            ModTupRngElt
            sage: mv.Parent()                           # optional - magma
            Full RSpace of degree 3 over Integer Ring

            sage: v = vector(QQ, [1/2, 3/4, 5/6])
            sage: mv = magma(v); mv                     # optional - magma
            (1/2 3/4 5/6)
            sage: mv.Type()                             # optional - magma
            ModTupFldElt
            sage: mv.Parent()                           # optional - magma
            Full Vector space of degree 3 over Rational Field

        A more demanding example::

            sage: R.<x,y,z> = QQ[]
            sage: v = vector([x^3, y, 2/3*z + x/y])
            sage: magma(v)                              # optional - magma
            (            x^3               y (2/3*y*z + x)/y)
            sage: magma(v).Parent()                     # optional - magma
            Full Vector space of degree 3 over Multivariate rational function field of rank 3 over Rational Field
        """
        V = magma(self._parent)
        v = [x._magma_init_(magma) for x in self.list()]
        return '%s![%s]'%(V.name(), ','.join(v))

def is_Vector(x):
    return isinstance(x, Vector)

cdef class Matrix(ModuleElement):

    cdef bint is_sparse_c(self):
        raise NotImplementedError

    cdef bint is_dense_c(self):
        raise NotImplementedError

    def __mul__(left, right):
        """
        Multiplication of matrix by matrix, vector, or scalar

        AUTHOR:

        - Gonzalo Tornaria (2007-06-25) - write test cases and fix them

        .. NOTE::

            scalar * matrix is implemented (and tested) in class RingElement
            vector * matrix is implemented (and tested) in class Vector

        TESTS:

        Here we test (matrix * matrix) multiplication::

            sage: parent(matrix(ZZ,2,2,[1,2,3,4])*matrix(ZZ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: parent(matrix(QQ,2,2,[1,2,3,4])*matrix(ZZ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: parent(matrix(ZZ,2,2,[1,2,3,4])*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: parent(matrix(QQ,2,2,[1,2,3,4])*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            sage: parent(matrix(QQ,2,2,[1,2,3,4])*matrix(ZZ['x'],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ['x'],2,2,[1,2,3,4])*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ,2,2,[1,2,3,4])*matrix(ZZ['x']['y'],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ['x']['y'],2,2,[1,2,3,4])*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ['x'],2,2,[1,2,3,4])*matrix(ZZ['x']['y'],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ['x']['y'],2,2,[1,2,3,4])*matrix(QQ['x'],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ['y'],2,2,[1,2,3,4])*matrix(ZZ['x']['y'],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ['x']['y'],2,2,[1,2,3,4])*matrix(QQ['y'],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(ZZ['x'],2,2,[1,2,3,4])*matrix(ZZ['y'],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(matrix(ZZ['x'],2,2,[1,2,3,4])*matrix(QQ['y'],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field'
            sage: parent(matrix(QQ['x'],2,2,[1,2,3,4])*matrix(ZZ['y'],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(matrix(QQ['x'],2,2,[1,2,3,4])*matrix(QQ['y'],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field'

        Here we test (matrix * vector) multiplication::

            sage: parent(matrix(ZZ,2,2,[1,2,3,4])*vector(ZZ,[1,2]))
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: parent(matrix(QQ,2,2,[1,2,3,4])*vector(ZZ,[1,2]))
            Vector space of dimension 2 over Rational Field
            sage: parent(matrix(ZZ,2,2,[1,2,3,4])*vector(QQ,[1,2]))
            Vector space of dimension 2 over Rational Field
            sage: parent(matrix(QQ,2,2,[1,2,3,4])*vector(QQ,[1,2]))
            Vector space of dimension 2 over Rational Field

            sage: parent(matrix(QQ,2,2,[1,2,3,4])*vector(ZZ['x'],[1,2]))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ['x'],2,2,[1,2,3,4])*vector(QQ,[1,2]))
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ,2,2,[1,2,3,4])*vector(ZZ['x']['y'],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ['x']['y'],2,2,[1,2,3,4])*vector(QQ,[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ['x'],2,2,[1,2,3,4])*vector(ZZ['x']['y'],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ['x']['y'],2,2,[1,2,3,4])*vector(QQ['x'],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ['y'],2,2,[1,2,3,4])*vector(ZZ['x']['y'],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ['x']['y'],2,2,[1,2,3,4])*vector(QQ['y'],[1,2]))
            Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(ZZ['x'],2,2,[1,2,3,4])*vector(ZZ['y'],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(matrix(ZZ['x'],2,2,[1,2,3,4])*vector(QQ['y'],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'
            sage: parent(matrix(QQ['x'],2,2,[1,2,3,4])*vector(ZZ['y'],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(matrix(QQ['x'],2,2,[1,2,3,4])*vector(QQ['y'],[1,2]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'

        Here we test (matrix * scalar) multiplication::

            sage: parent(matrix(ZZ,2,2,[1,2,3,4])*ZZ(1))
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: parent(matrix(QQ,2,2,[1,2,3,4])*ZZ(1))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: parent(matrix(ZZ,2,2,[1,2,3,4])*QQ(1))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: parent(matrix(QQ,2,2,[1,2,3,4])*QQ(1))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            sage: parent(matrix(QQ,2,2,[1,2,3,4])*ZZ['x'](1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ['x'],2,2,[1,2,3,4])*QQ(1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ,2,2,[1,2,3,4])*ZZ['x']['y'](1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ['x']['y'],2,2,[1,2,3,4])*QQ(1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ['x'],2,2,[1,2,3,4])*ZZ['x']['y'](1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ['x']['y'],2,2,[1,2,3,4])*QQ['x'](1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(QQ['y'],2,2,[1,2,3,4])*ZZ['x']['y'](1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(matrix(ZZ['x']['y'],2,2,[1,2,3,4])*QQ['y'](1))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(matrix(ZZ['x'],2,2,[1,2,3,4])*ZZ['y'](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring' and 'Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(matrix(ZZ['x'],2,2,[1,2,3,4])*QQ['y'](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring' and 'Univariate Polynomial Ring in y over Rational Field'
            sage: parent(matrix(QQ['x'],2,2,[1,2,3,4])*ZZ['y'](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(matrix(QQ['x'],2,2,[1,2,3,4])*QQ['y'](1))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in y over Rational Field'

        Here we test (scalar * matrix) multiplication::

            sage: parent(ZZ(1)*matrix(ZZ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: parent(QQ(1)*matrix(ZZ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: parent(ZZ(1)*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: parent(QQ(1)*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            sage: parent(QQ(1)*matrix(ZZ['x'],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ['x'](1)*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field

            sage: parent(QQ(1)*matrix(ZZ['x']['y'],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ['x']['y'](1)*matrix(QQ,2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(QQ['x'](1)*matrix(ZZ['x']['y'],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ['x']['y'](1)*matrix(QQ['x'],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(QQ['y'](1)*matrix(ZZ['x']['y'],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(ZZ['x']['y'](1)*matrix(QQ['y'],2,2,[1,2,3,4]))
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(ZZ['x'](1)*matrix(ZZ['y'],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Integer Ring' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(ZZ['x'](1)*matrix(QQ['y'],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Integer Ring' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field'
            sage: parent(QQ['x'](1)*matrix(ZZ['y'],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Rational Field' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(QQ['x'](1)*matrix(QQ['y'],2,2,[1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Univariate Polynomial Ring in x over Rational Field' and 'Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field'

        Examples with matrices having matrix coefficients::

            sage: m = matrix
            sage: a = m([[m([[1,2],[3,4]]),m([[5,6],[7,8]])],[m([[9,10],[11,12]]),m([[13,14],[15,16]])]])
            sage: 3*a
            [[ 3  6]
            [ 9 12] [15 18]
            [21 24]]
            [[27 30]
            [33 36] [39 42]
            [45 48]]

            sage: m = matrix
            sage: a = m([[m([[1,2],[3,4]]),m([[5,6],[7,8]])],[m([[9,10],[11,12]]),m([[13,14],[15,16]])]])
            sage: a*3
            [[ 3  6]
            [ 9 12] [15 18]
            [21 24]]
            [[27 30]
            [33 36] [39 42]
            [45 48]]
        """
        if have_same_parent(left, right):
            return (<Matrix>left)._matrix_times_matrix_(<Matrix>right)
        return coercion_model.bin_op(left, right, mul)

    def __truediv__(left, right):
        """
        Division of the matrix ``left`` by the matrix or scalar
        ``right``.

        EXAMPLES::

            sage: a = matrix(ZZ, 2, range(4))
            sage: operator.truediv(a, 5)
            [ 0 1/5]
            [2/5 3/5]
            sage: a = matrix(ZZ, 2, range(4))
            sage: b = matrix(ZZ, 2, [1,1,0,5])
            sage: operator.truediv(a, b)
            [  0 1/5]
            [  2 1/5]
            sage: c = matrix(QQ, 2, [3,2,5,7])
            sage: operator.truediv(c, a)
            [-5/2  3/2]
            [-1/2  5/2]
        """
        if have_same_parent(left, right):
            return left * ~right
        return coercion_model.bin_op(left, right, truediv)

    def __div__(left, right):
        """
        Division of the matrix ``left`` by the matrix or scalar ``right``.

        EXAMPLES::

            sage: a = matrix(ZZ, 2, range(4))
            sage: a / 5
            [ 0 1/5]
            [2/5 3/5]
            sage: a = matrix(ZZ, 2, range(4))
            sage: b = matrix(ZZ, 2, [1,1,0,5])
            sage: a / b
            [  0 1/5]
            [  2 1/5]
            sage: c = matrix(QQ, 2, [3,2,5,7])
            sage: c / a
            [-5/2  3/2]
            [-1/2  5/2]
            sage: a / c
            [-5/11  3/11]
            [-1/11  5/11]
            sage: a / 7
            [  0 1/7]
            [2/7 3/7]

        Other rings work just as well::

            sage: a = matrix(GF(3),2,2,[0,1,2,0])
            sage: b = matrix(ZZ,2,2,[4,6,1,2])
            sage: a / b
            [1 2]
            [2 0]
            sage: c = matrix(GF(3),2,2,[1,2,1,1])
            sage: a / c
            [1 2]
            [1 1]
            sage: a = matrix(RDF,2,2,[.1,-.4,1.2,-.6])
            sage: b = matrix(RDF,2,2,[.3,.1,-.5,1.3])
            sage: a / b # rel tol 1e-10
            [-0.15909090909090906 -0.29545454545454547]
            [   2.863636363636364  -0.6818181818181817]
            sage: R.<t> = ZZ['t']
            sage: a = matrix(R,2,2,[t^2,t+1,-t,t+2])
            sage: b = matrix(R,2,2,[t^3-1,t,-t+3,t^2])
            sage: a / b
            [      (t^4 + t^2 - 2*t - 3)/(t^5 - 3*t)               (t^4 - t - 1)/(t^5 - 3*t)]
            [       (-t^3 + t^2 - t - 6)/(t^5 - 3*t) (t^4 + 2*t^3 + t^2 - t - 2)/(t^5 - 3*t)]
        """
        if have_same_parent(left, right):
            return left * ~right
        return coercion_model.bin_op(left, right, div)

    cdef _vector_times_matrix_(matrix_right, Vector vector_left):
        raise TypeError

    cdef _matrix_times_vector_(matrix_left, Vector vector_right):
        raise TypeError

    cdef _matrix_times_matrix_(left, Matrix right):
        raise TypeError



def is_Matrix(x):
    return isinstance(x, Matrix)

def is_IntegralDomainElement(x):
    """
    Return ``True`` if x is of type IntegralDomainElement.
    """
    return isinstance(x, IntegralDomainElement)

cdef class IntegralDomainElement(CommutativeRingElement):
    def is_nilpotent(self):
        return self.is_zero()


def is_DedekindDomainElement(x):
    """
    Return ``True`` if x is of type DedekindDomainElement.
    """
    return isinstance(x, DedekindDomainElement)

cdef class DedekindDomainElement(IntegralDomainElement):
    pass

def is_PrincipalIdealDomainElement(x):
    """
    Return ``True`` if x is of type PrincipalIdealDomainElement.
    """
    return isinstance(x, PrincipalIdealDomainElement)

cdef class PrincipalIdealDomainElement(DedekindDomainElement):
    def lcm(self, right):
        """
        Return the least common multiple of ``self`` and ``right``.
        """
        if not isinstance(right, Element) or not ((<Element>right)._parent is self._parent):
            from sage.arith.all import lcm
            return coercion_model.bin_op(self, right, lcm)
        return self._lcm(right)


# This is pretty nasty low level stuff. The idea is to speed up construction
# of EuclideanDomainElements (in particular Integers) by skipping some tp_new
# calls up the inheritance tree.
PY_SET_TP_NEW(EuclideanDomainElement, Element)

def is_EuclideanDomainElement(x):
    """
    Return ``True`` if x is of type EuclideanDomainElement.
    """
    return isinstance(x, EuclideanDomainElement)

cdef class EuclideanDomainElement(PrincipalIdealDomainElement):

    def degree(self):
        raise NotImplementedError

    def leading_coefficient(self):
        raise NotImplementedError

    def quo_rem(self, other):
        raise NotImplementedError

    def __divmod__(self, other):
        """
        Return the quotient and remainder of ``self`` divided by ``other``.

        EXAMPLES::

            sage: divmod(5,3)
            (1, 2)
            sage: divmod(25r,12)
            (2, 1)
            sage: divmod(25,12r)
            (2, 1)

        """
        if isinstance(self, Element):
            return self.quo_rem(other)
        else:
            x, y = canonical_coercion(self, other)
            return x.quo_rem(y)

    cpdef _floordiv_(self, right):
        """
        Quotient of division of ``self`` by other.  This is denoted //.

        This default implementation assumes that ``quo_rem`` has been
        implemented.

        EXAMPLES::

            sage: cython('''
            ....: from sage.structure.element cimport EuclideanDomainElement
            ....: cdef class MyElt(EuclideanDomainElement):
            ....:     def quo_rem(self, other):
            ....:         return self._parent.var('quo,rem')
            ....: ''')
            sage: e = MyElt(SR)
            sage: e // e
            quo
        """
        Q, _ = self.quo_rem(right)
        return Q

    cpdef _mod_(self, other):
        """
        Remainder of division of ``self`` by other.

        This default implementation assumes that ``quo_rem`` has been
        implemented.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: x % (x+1)
            -1
            sage: (x^3 + x - 1) % (x^2 - 1)
            2*x - 1

        ::

            sage: cython('''
            ....: from sage.structure.element cimport EuclideanDomainElement
            ....: cdef class MyElt(EuclideanDomainElement):
            ....:     def quo_rem(self, other):
            ....:         return self._parent.var('quo,rem')
            ....: ''')
            sage: e = MyElt(SR)
            sage: e % e
            rem
        """
        _, R = self.quo_rem(other)
        return R


def is_FieldElement(x):
    """
    Return ``True`` if x is of type FieldElement.
    """
    return isinstance(x, FieldElement)

cdef class FieldElement(CommutativeRingElement):
    cpdef _floordiv_(self, right):
        """
        Return the quotient of self and other. Since these are field
        elements, the floor division is exactly the same as usual division.

        EXAMPLES::

            sage: K.<b> = NumberField(x^4 + x^2 + 2/3)
            sage: c = (1+b) // (1-b); c
            3/4*b^3 + 3/4*b^2 + 3/2*b + 1/2
            sage: (1+b) / (1-b) == c
            True
            sage: c * (1-b)
            b + 1
        """
        return self._div_(right)

    def is_unit(self):
        r"""
        Return ``True`` if ``self`` is a unit in its parent ring.

        EXAMPLES::

            sage: a = 2/3; a.is_unit()
            True

        On the other hand, 2 is not a unit, since its parent is `\ZZ`.

        ::

            sage: a = 2; a.is_unit()
            False
            sage: parent(a)
            Integer Ring

        However, a is a unit when viewed as an element of QQ::

            sage: a = QQ(2); a.is_unit()
            True
        """
        return not not self

    def _lcm(self, FieldElement other):
        """
        Return the least common multiple of ``self`` and other.
        """
        if self.is_zero() and other.is_zero():
            return self
        else:
            return self._parent(1)

    def quo_rem(self, right):
        r"""
        Return the quotient and remainder obtained by dividing ``self`` by
        ``right``. Since this element lives in a field, the remainder is always
        zero and the quotient is ``self/right``.

        TESTS:

        Test if :trac:`8671` is fixed::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = R.quo(y^2 + 1)
            sage: S.is_field = lambda : False
            sage: F = Frac(S); u = F.one()
            sage: u.quo_rem(u)
            (1, 0)
        """
        if not isinstance(right, FieldElement) or not (parent(right) is self._parent):
            right = self.parent()(right)
        return self/right, 0

    def divides(self, FieldElement other):
        r"""
        Check whether ``self`` divides other, for field elements.

        Since this is a field, all values divide all other values,
        except that zero does not divide any non-zero values.

        EXAMPLES::

            sage: K.<rt3> = QQ[sqrt(3)]
            sage: K(0).divides(rt3)
            False
            sage: rt3.divides(K(17))
            True
            sage: K(0).divides(K(0))
            True
            sage: rt3.divides(K(0))
            True
        """
        if not (other._parent is self._parent):
            other = self.parent()(other)
        return bool(self) or other.is_zero()

def is_AlgebraElement(x):
    """
    Return ``True`` if x is of type AlgebraElement.

    TESTS::

        sage: from sage.structure.element import is_AlgebraElement
        sage: R.<x,y> = FreeAlgebra(QQ,2)
        sage: is_AlgebraElement(x*y)
        True

        sage: is_AlgebraElement(1)
        False
    """
    return isinstance(x, AlgebraElement)

cdef class AlgebraElement(RingElement):
    pass

def is_CommutativeAlgebraElement(x):
    """
    Return ``True`` if x is of type CommutativeAlgebraElement.
    """
    return isinstance(x, CommutativeAlgebraElement)

cdef class CommutativeAlgebraElement(CommutativeRingElement):
    pass

def is_InfinityElement(x):
    """
    Return ``True`` if x is of type InfinityElement.

    TESTS::

        sage: from sage.structure.element import is_InfinityElement
        sage: is_InfinityElement(1)
        False

        sage: is_InfinityElement(oo)
        True
    """
    return isinstance(x, InfinityElement)

cdef class InfinityElement(RingElement):
    def __invert__(self):
        from sage.rings.all import ZZ
        return ZZ(0)


#################################################################################
#
#  Coercion of elements
#
#################################################################################

cpdef canonical_coercion(x, y):
    """
    ``canonical_coercion(x,y)`` is what is called before doing an
    arithmetic operation between ``x`` and ``y``.  It returns a pair ``(z,w)``
    such that ``z`` is got from ``x`` and ``w`` from ``y`` via canonical coercion and
    the parents of ``z`` and ``w`` are identical.

    EXAMPLES::

        sage: A = Matrix([[0, 1], [1, 0]])
        sage: canonical_coercion(A, 1)
        (
        [0 1]  [1 0]
        [1 0], [0 1]
        )
    """
    return coercion_model.canonical_coercion(x,y)

cpdef bin_op(x, y, op):
    return coercion_model.bin_op(x,y,op)


def coerce(Parent p, x):
    try:
        return p._coerce_c(x)
    except AttributeError:
        return p(x)

# We define this base class here to avoid circular cimports.
cdef class CoercionModel:
    """
    Most basic coercion scheme. If it doesn't already match, throw an error.
    """
    cpdef canonical_coercion(self, x, y):
        if parent(x) is parent(y):
            return x,y
        raise TypeError("no common canonical parent for objects with parents: '%s' and '%s'"%(parent(x), parent(y)))

    cpdef bin_op(self, x, y, op):
        if parent(x) is parent(y):
            return op(x,y)
        raise bin_op_exception(op, x, y)

    cpdef richcmp(self, x, y, int op):
        x, y = self.canonical_coercion(x, y)
        return PyObject_RichCompare(x, y, op)


from . import coerce
cdef CoercionModel coercion_model = coerce.CoercionModel_cache_maps()


def get_coercion_model():
    """
    Return the global coercion model.

    EXAMPLES::

       sage: import sage.structure.element as e
       sage: cm = e.get_coercion_model()
       sage: cm
       <sage.structure.coerce.CoercionModel_cache_maps object at ...>
    """
    return coercion_model

def set_coercion_model(cm):
    global coercion_model
    coercion_model = cm

def coercion_traceback(dump=True):
    r"""
    This function is very helpful in debugging coercion errors. It prints
    the tracebacks of all the errors caught in the coercion detection. Note
    that failure is cached, so some errors may be omitted the second time
    around (as it remembers not to retry failed paths for speed reasons.

    For performance and caching reasons, exception recording must be
    explicitly enabled before using this function.

    EXAMPLES::

        sage: cm = sage.structure.element.get_coercion_model()
        sage: cm.record_exceptions()
        sage: 1 + 1/5
        6/5
        sage: coercion_traceback()  # Should be empty, as all went well.
        sage: 1/5 + GF(5).gen()
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand parent(s) for '+': 'Rational Field' and 'Finite Field of size 5'
        sage: coercion_traceback()
        Traceback (most recent call last):
        ...
        TypeError: no common canonical parent for objects with parents: 'Rational Field' and 'Finite Field of size 5'
    """
    if dump:
        for traceback in coercion_model.exception_stack():
            print(traceback)
    else:
        return coercion_model.exception_stack()


def coerce_binop(method):
    r"""
    Decorator for a binary operator method for applying coercion to the
    arguments before calling the method.

    Consider a parent class in the category framework, `S`, whose element class
    expose a method `binop`. If `a` and `b` are elements of `S`, then
    `a.binop(b)` behaves as expected. If `a` and `b` are not elements of `S`,
    but rather have a common parent `T` whose element class also exposes
    `binop`, we would rather expect `a.binop(b)` to compute `aa.binop(bb)`,
    where `aa = T(a)` and `bb = T(b)`. This decorator ensures that behaviour
    without having to otherwise modify the implementation of `binop` on the
    element class of `A`.

    Since coercion will be attempted on the arguments of the decorated method, a
    `TypeError` will be thrown if there is no common parent between the
    elements. An `AttributeError` or `NotImplementedError` or similar will be
    thrown if there is a common parent of the arguments, but its element class
    does not implement a method of the same name as the decorated method.

    EXAMPLES:

    Sparse polynomial rings uses `@coerce_binop` on `gcd`::

        sage: S.<x> = PolynomialRing(ZZ,sparse=True)
        sage: f = x^2
        sage: g = x
        sage: f.gcd(g)  #indirect doctest
        x
        sage: T = PolynomialRing(QQ, name='x', sparse=True)
        sage: h = 1/2*T(x)
        sage: u = f.gcd(h); u  #indirect doctest
        x
        sage: u.parent() == T
        True

    Another real example::

        sage: R1=QQ['x,y']
        sage: R2=QQ['x,y,z']
        sage: f=R1(1)
        sage: g=R1(2)
        sage: h=R2(1)
        sage: f.gcd(g)
        1
        sage: f.gcd(g,algorithm='modular')
        1
        sage: f.gcd(h)
        1
        sage: f.gcd(h,algorithm='modular')
        1
        sage: h.gcd(f)
        1
        sage: h.gcd(f,'modular')
        1

    We demonstrate a small class using `@coerce_binop` on a method::

        sage: from sage.structure.element import coerce_binop
        sage: class MyRational(Rational):
        ....:     def __init__(self,value):
        ....:         self.v = value
        ....:     @coerce_binop
        ....:     def test_add(self, other, keyword='z'):
        ....:         return (self.v, other, keyword)

    Calls func directly if the two arguments have the same parent::

        sage: x = MyRational(1)
        sage: x.test_add(1/2)
        (1, 1/2, 'z')
        sage: x.test_add(1/2, keyword=3)
        (1, 1/2, 3)

    Passes through coercion and does a method lookup if the left operand is not
    the same. If the common parent's element class does not have a method of the
    same name, an exception is raised::

        sage: x.test_add(2)
        (1, 2, 'z')
        sage: x.test_add(2, keyword=3)
        (1, 2, 3)
        sage: x.test_add(CC(2))
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.complex_number.ComplexNumber' object has no attribute 'test_add'

    TESTS:

    Test that additional arguments given to the method do not override
    the ``self`` argument, see :trac:`21322`::

        sage: f.gcd(g, 1)
        Traceback (most recent call last):
        ...
        TypeError: algorithm 1 not supported
    """
    @sage_wraps(method)
    def new_method(self, other, *args, **kwargs):
        if have_same_parent(self, other):
            return method(self, other, *args, **kwargs)
        else:
            a, b = coercion_model.canonical_coercion(self, other)
            if a is self:
                return method(a, b, *args, **kwargs)
            else:
                return getattr(a, method.__name__)(b, *args, **kwargs)
    return new_method


###############################################################################

def generic_power(a, n, one=None):
    """
    Computes `a^n`, where `n` is an integer, and `a` is an object which
    supports multiplication.  Optionally an additional argument,
    which is used in the case that ``n == 0``:

    - ``one`` - the "unit" element, returned directly (can be anything)

    If this is not supplied, ``int(1)`` is returned.

    EXAMPLES::

        sage: from sage.structure.element import generic_power
        sage: generic_power(int(12),int(0))
        1
        sage: generic_power(int(0),int(100))
        0
        sage: generic_power(Integer(10),Integer(0))
        1
        sage: generic_power(Integer(0),Integer(23))
        0
        sage: sum([generic_power(2,i) for i in range(17)]) #test all 4-bit combinations
        131071
        sage: F = Zmod(5)
        sage: a = generic_power(F(2), 5); a
        2
        sage: a.parent() is F
        True
        sage: a = generic_power(F(1), 2)
        sage: a.parent() is F
        True

        sage: generic_power(int(5), 0)
        1
    """

    return generic_power_c(a,n,one)

cdef generic_power_c(a, nn, one):
    try:
        n = PyNumber_Index(nn)
    except TypeError:
        try:
            # Try harder, since many things coerce to Integer.
            from sage.rings.integer import Integer
            n = int(Integer(nn))
        except TypeError:
            raise NotImplementedError("non-integral exponents not supported")

    if not n:
        if one is None:
            if isinstance(a, Element):
                return (<Element>a)._parent.one()
            try:
                try:
                    return a.parent().one()
                except AttributeError:
                    return type(a)(1)
            except Exception:
                return 1 #oops, the one sucks
        else:
            return one
    elif n < 0:
        # I don't think raising division by zero is really my job. It should
        # be the one of ~a. Moreover, this does not handle the case of monoids
        # with partially defined division (e.g. the multiplicative monoid of a
        # ring such as ZZ/12ZZ)
        #        if not a:
        #            raise ZeroDivisionError
        a = ~a
        n = -n

    if n < 4:
        # These cases will probably be called often
        # and don't benefit from the code below
        if n == 1:
            return a
        elif n == 2:
            return a*a
        elif n == 3:
            return a*a*a

    # check for idempotence, and store the result otherwise
    aa = a*a
    if aa == a:
        return a

    # since we've computed a^2, let's start squaring there
    # so, let's keep the least-significant bit around, just
    # in case.
    m = n & 1
    n = n >> 1

    # One multiplication can be saved by starting with
    # the second-smallest power needed rather than with 1
    # we've already squared a, so let's start there.
    apow = aa
    while n&1 == 0:
        apow = apow*apow
        n = n >> 1
    power = apow
    n = n >> 1

    # now multiply that least-significant bit in...
    if m:
        power = power * a

    # and this is straight from the book.
    while n != 0:
        apow = apow*apow
        if n&1 != 0:
            power = power*apow
        n = n >> 1

    return power
