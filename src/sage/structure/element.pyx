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

The Abstract Element Class Hierarchy
------------------------------------

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
                            FiniteFieldElement
                        CommutativeAlgebraElement
                    AlgebraElement   (note -- can't derive from module, since no multiple inheritance)
                        CommutativeAlgebra ??? (should be removed from element.pxd)
                        Matrix
                    InfinityElement
                        PlusInfinityElement
                        MinusInfinityElement
                AdditiveGroupElement
                Vector

            MonoidElement
                MultiplicativeGroupElement
        ElementWithCachedMethod


How to Define a New Element Class
---------------------------------

Elements typically define a method ``_new_c``, e.g.,

::

    cdef _new_c(self, defining data):
        cdef FreeModuleElement_generic_dense x
        x = FreeModuleElement_generic_dense.__new__(FreeModuleElement_generic_dense)
        x._parent = self._parent
        x._entries = v

that creates a new sibling very quickly from defining data
with assumed properties.

Sage has a special system in place for handling arithmetic operations
for all Element subclasses. There are various rules that must be
followed by both arithmetic implementers and callers.

A quick summary for the impatient:

- To implement addition for any Element class, override ``def _add_()``.
- If you want to add ``x`` and ``y``, whose parents you know are **identical**,
  you may call ``_add_(x, y)``. This will be the fastest way to guarantee
  that the correct implementation gets called. Of course you can still
  always use ``x + y``.

Now in more detail. The aims of this system are to provide (1) an efficient
calling protocol from both Python and Cython, (2) uniform coercion semantics
across Sage, (3) ease of use, (4) readability of code.

We will take addition of RingElements as an example; all other
operators and classes are similar. There are three relevant functions,
with subtly differing names (``add`` vs.  ``iadd``, single vs. double
underscores).

-  **def RingElement.__add__**

   This function is called by Python or Cython when the binary "+" operator
   is encountered. It **assumes** that at least one of its arguments is a
   RingElement; only a really twisted programmer would violate this
   condition. It has a fast pathway to deal with the most common case
   where the arguments have the same parent. Otherwise, it uses the coercion
   module to work out how to make them have the same parent. After any
   necessary coercions have been performed, it calls ``_add_`` to dispatch to
   the correct underlying addition implementation.

   Note that although this function is declared as ``def``, it doesn't have the
   usual overheads associated with Python functions (either for the caller or
   for ``__add__`` itself). This is because Python has optimised calling
   protocols for such special functions.

-  **def RingElement._add_**

   This is the function you should override to implement addition in a
   subclass of ``RingElement``.

   The two arguments to this function are guaranteed to have the **same
   parent**. Its return value **must** have the **same parent** as its
   arguments.

   If you want to add two objects and you know that their parents are
   the same object, you are encouraged to call this function directly,
   instead of using ``x + y``.

   When implementing ``_add_`` in a Cython extension class, use
   ``cpdef _add_`` instead of ``def _add_``.

For speed, there are also *inplace* versions of the arithmetic commands.
**Do not** call them directly, they may mutate the object and will be called
when and only when it has been determined that the old object will no longer
be accessible from the calling function after this operation.

-  **def RingElement._iadd_**

   This is the function you should override to implement inplace
   addition in a Python subclass of RingElement.

   The two arguments to this function are guaranteed to have the **same
   parent**. Its return value **must** have the **same parent** as its
   arguments.

   The default implementation of this function is to call ``_add_``,
   so if no one has defined a Python implementation, the correct Cython
   implementation will get called.
"""

##################################################################
# Generic element, so all this functionality must be defined
# by any element.  Derived class must call __init__
##################################################################

include "sage/ext/cdefs.pxi"
include "sage/ext/python.pxi"
include "coerce.pxi"
from sage.ext.stdsage cimport *

from cpython.ref cimport PyObject

import types
cdef add, sub, mul, div, iadd, isub, imul, idiv
from operator import add, sub, mul, div, iadd, isub, imul, idiv

cdef MethodType
from types import MethodType

from sage.structure.sage_object cimport rich_to_bool
from sage.structure.coerce cimport py_scalar_to_element
from sage.structure.parent cimport Parent
from sage.structure.misc import is_extension_type, getattr_from_other_class
from sage.misc.lazy_format import LazyFormat
from sage.misc import sageinspect

# Create a dummy attribute error, using some kind of lazy error message,
# so that neither the error itself not the message need to be created
# repeatedly, which would cost time.
from sage.structure.misc cimport AttributeErrorMessage
cdef AttributeErrorMessage dummy_error_message = AttributeErrorMessage(None, '')
dummy_attribute_error = AttributeError(dummy_error_message)


def make_element(_class, _dict, parent):
    """
    This function is only here to support old pickles.

    Pickling functionality is moved to Element.{__getstate__,__setstate__}
    functions.
    """
    from sage.misc.pickle_old import make_element_old
    return make_element_old(_class, _dict, parent)


def parent(x):
    """
    Return the parent of the element ``x``.

    Usually, this means the mathematical object of which ``x`` is an
    element.

    INPUT:

    - ``x`` -- an element

    OUTPUT:

    - if ``x`` is a Sage :class:`Element`, return ``x.parent()``.

    - if ``x`` has a ``parent`` method and ``x`` does not have an
      ``__int__`` or ``__float__`` method, return ``x.parent()``.

    - otherwise, return ``type(x)``.

    .. SEEALSO::

        `Parents, Conversion and Coercion <http://www.sagemath.org/doc/tutorial/tour_coercion.html>`_
            Section in the Sage Tutorial

    EXAMPLES::

        sage: a = 42
        sage: parent(a)
        Integer Ring
        sage: b = 42/1
        sage: parent(b)
        Rational Field
        sage: c = 42.0
        sage: parent(c)
        Real Field with 53 bits of precision

    Some more complicated examples::

        sage: x = Partition([3,2,1,1,1])
        sage: parent(x)
        Partitions
        sage: v = vector(RDF, [1,2,3])
        sage: parent(v)
        Vector space of dimension 3 over Real Double Field

    The following are not considered to be elements, so the type is
    returned::

        sage: d = int(42)  # Python int
        sage: parent(d)
        <type 'int'>
        sage: L = range(10)
        sage: parent(L)
        <type 'list'>
    """
    return parent_c(x)

def have_same_parent(left, right):
    """
    Return ``True`` if and only if ``left`` and ``right`` have the
    same parent.

    .. WARNING::

        This function assumes that at least one of the arguments is a
        Sage :class:`Element`. When in doubt, use the slower
        ``parent(left) is parent(right)`` instead.

    EXAMPLES::

        sage: from sage.structure.element import have_same_parent
        sage: have_same_parent(1, 3)
        True
        sage: have_same_parent(1, 1/2)
        False
        sage: have_same_parent(gap(1), gap(1/2))
        True

    These have different types but the same parent::

        sage: a = RLF(2)
        sage: b = exp(a)
        sage: type(a)
        <type 'sage.rings.real_lazy.LazyWrapper'>
        sage: type(b)
        <type 'sage.rings.real_lazy.LazyNamedUnop'>
        sage: have_same_parent(a, b)
        True
    """
    return have_same_parent_c(left, right)


cdef dict _coerce_op_symbols = {'mul':'*', 'add':'+', 'sub':'-', 'div':'/', 'imul': '*', 'iadd': '+', 'isub':'-', 'idiv':'/'}

cdef str arith_error_message(x, y, op):
    name = op.__name__
    try:
        name = _coerce_op_symbols[name]
    except KeyError:
        pass
    return "unsupported operand parent(s) for '%s': '%s' and '%s'"%(name, parent_c(x), parent_c(y))


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
    (RingElement, ModuleElement, etc) derive from this type.

    Subtypes must either call ``__init__()`` to set ``_parent``, or may
    set ``_parent`` themselves if that would be more efficient.

    .. automethod:: _cmp_
    .. automethod:: _richcmp_
    """
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

    cdef _set_parent_c(self, Parent parent):
        self._parent = parent

    def _make_new_with_parent_c(self, Parent parent):
        self._parent = parent
        return self


    def __getattr__(self, str name):
        """
        Lookup a method or attribute from the category abstract classes.

        Let ``P`` be a parent in a category ``C``. Usually the methods
        of ``C.element_class`` are made directly available to elements
        of ``P`` via standard class inheritance. This is not the case
        any more if the elements of ``P`` are instances of an
        extension type. See :class:`Category`. for details.

        The purpose of this method is to emulate this inheritance: for
        ``e`` and element of ``P``, if an attribute or method
        ``e.foo`` is not found in the super classes of ``e``, it's
        looked up manually in ``C.element_class`` and bound to ``e``.

        .. NOTES::

            - Attributes beginning with two underscores but not ending
              with an unnderscore are considered private and are thus
              exempted from the lookup in ``cat.element_class``.

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

        We test that "private" attributes are not requested from the element class
        of the category (:trac:`10467`)::

            sage: C = EuclideanDomains()
            sage: P.<x> = QQ[]
            sage: C.element_class.__foo = 'bar'
            sage: x.parent() in C
            True
            sage: x.__foo
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint' object has no attribute '__foo'
        """
        if (name.startswith('__') and not name.endswith('_')):
            dummy_error_message.cls = type(self)
            dummy_error_message.name = name
            raise dummy_attribute_error
        cdef Parent P = self._parent or self.parent()
        if P is None or P._category is None:
            dummy_error_message.cls = type(self)
            dummy_error_message.name = name
            raise dummy_attribute_error
        return getattr_from_other_class(self, P._abstract_element_class, name)

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
        from sage.structure.parent import dir_with_other_class
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

    def __hash__(self):
        return hash(str(self))

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
        if not hasattr(self,'__call__'):
            return self
        parent=self._parent
        # We should better not test for ParentWIthGens,
        # as this is essentially deprecated.
        #from sage.structure.parent_gens import is_ParentWithGens
        #if not is_ParentWithGens(parent):
        #    return self
        # Better: Duck typing!
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
        Return a numerical approximation of x with at least prec bits of
        precision.

        EXAMPLES::

            sage: (2/3).n()
            0.666666666666667
            sage: pi.n(digits=10)  # indirect doctest
            3.141592654
            sage: pi.n(prec=20)   # indirect doctest
            3.1416

        TESTS:

        Check that :trac:`14778` is fixed::

            sage: (0).n(algorithm='foo')
            0.000000000000000
        """
        from sage.misc.functional import numerical_approx
        return numerical_approx(self, prec=prec, digits=digits,
                                algorithm=algorithm)
    n = numerical_approx
    N = n

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

    cdef int _cmp(self, other) except -2:
        """
        Compare ``left`` and ``right`` using the coercion framework.

        ``self`` and ``other`` are coerced to a common parent and then
        ``_cmp_`` and ``_richcmp_`` are tried.

        EXAMPLES:

        We create an ``Element`` class where we define ``__cmp__``
        and ``_richcmp_`` and check that comparison works::

            sage: cython('''
            ....: from sage.structure.sage_object cimport rich_to_bool
            ....: from sage.structure.element cimport Element
            ....: cdef class FloatCmp(Element):
            ....:     cdef float x
            ....:     def __init__(self, float v):
            ....:         self.x = v
            ....:     def __cmp__(self, other):
            ....:         return (<Element>self)._cmp(other)
            ....:     cpdef _richcmp_(self, Element other, int op):
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
        if have_same_parent_c(self, other):
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
            if (<Element>left)._richcmp_(<Element>right, Py_EQ):
                return 0
            if (<Element>left)._richcmp_(<Element>right, Py_LT):
                return -1
            if (<Element>left)._richcmp_(<Element>right, Py_GT):
                return 1
            raise

    cdef _richcmp(self, other, int op):
        """
        Compare ``self`` and ``other`` using the coercion framework,
        comparing according to the comparison operator ``op``.

        This method exists only because of the strange way that Python
        handles inheritance of ``__richcmp__``. A Cython class should
        always define ``__richcmp__`` as calling ``_richcmp``.

        Normally, a class will not redefine ``_richcmp`` but rely on
        this ``Element._richcmp`` method which uses coercion to
        compare elements. Then ``_richcmp_`` is called on the coerced
        elements.

        If a class wants to implement rich comparison without coercion,
        then ``_richcmp`` should be defined (as well as ``__richcmp__``
        as usual).
        See :class:`sage.numerical.linear_functions.LinearConstraint`
        for such an example.

        For efficiency reasons, a class can do certain "manual"
        coercions directly in ``__richcmp__``, using ``_richcmp``
        for the remaining cases. This is done for example in
        :class:`Integer`.
        """
        if have_same_parent_c(self, other):
            # Same parents, in particular other must be an Element too.
            # The explicit cast <Element>other tells Cython to omit the
            # check isinstance(other, Element) when calling _richcmp_
            return self._richcmp_(<Element>other, op)

        # Some very special cases
        if self is None or self is Ellipsis:
            return rich_to_bool(op, -1)
        if other is None or other is Ellipsis:
            return rich_to_bool(op, 1)

        # Different parents => coerce
        try:
            left, right = coercion_model.canonical_coercion(self, other)
        except (TypeError, NotImplementedError):
            pass
        else:
            if isinstance(left, Element):
                return (<Element>left)._richcmp(<Element>right, op)
            # left and right are the same non-Element type:
            # use a plain cmp()
            return rich_to_bool(op, cmp(left, right))

        # Comparing with coercion didn't work, try something else.

        # Often things are compared against 0, even when there is no
        # canonical coercion from ZZ. For ModuleElements, we manually
        # convert zero to handle this case.
        from sage.rings.integer import Integer
        cdef Element zero
        try:
            if isinstance(self, ModuleElement) and isinstance(other, (int, float, Integer)) and not other:
                zero = (<Element>self)._parent(0)
                return self._richcmp_(zero, op)
            elif isinstance(other, ModuleElement) and isinstance(self, (int, float, Integer)) and not self:
                zero = (<Element>other)._parent(0)
                return zero._richcmp_(other, op)
        except (TypeError, AttributeError):
            pass

        # If types are not equal: compare types
        cdef int r = cmp(type(self), type(other))
        if r:
            return rich_to_bool(op, r)

        # Final attempt: compare by id()
        if (<unsigned long><PyObject*>self) >= (<unsigned long><PyObject*>other):
            # It cannot happen that self is other, since they don't
            # have the same parent.
            return rich_to_bool(op, 1)
        else:
            return rich_to_bool(op, -1)

    ####################################################################
    # For a derived Cython class, you **must** put the __richcmp__
    # method below in your subclasses, in order for it to take
    # advantage of the above generic comparison code.
    #
    # You must also define either _cmp_ (if your subclass is totally
    # ordered), _richcmp_ (if your subclass is partially ordered), or
    # both (if your class has both a total order and a partial order,
    # or if that gives better performance).
    #
    # Rich comparisons (like a < b) will default to using _richcmp_,
    # three-way comparisons (like cmp(a,b)) will default to using
    # _cmp_. But if you define just one of _richcmp_ and _cmp_, it will
    # be used for all kinds of comparisons.
    #
    # In the _cmp_ and _richcmp_ methods, you can assume that both
    # arguments have identical parents.
    ####################################################################
    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    def __cmp__(left, right):
        return (<Element>left)._cmp(right)

    cpdef _richcmp_(left, Element right, int op):
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

    cpdef int _cmp_(left, Element right) except -2:
        """
        Default three-way comparison method which only checks for a
        Python class defining ``__cmp__``.
        """
        left_cmp = left.__cmp__
        if isinstance(left_cmp, MethodType):
            return left_cmp(right)
        msg = LazyFormat("comparison not implemented for %r")%type(left)
        raise NotImplementedError(msg)


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
        ....:     "    def __cmp__(left, right):",
        ....:     "        return (<Element>left)._cmp(right)",
        ....:     "    def __richcmp__(left, right, op):",
        ....:     "        return (<Element>left)._richcmp(right,op)",
        ....:     "    cpdef int _cmp_(left, Element right) except -2:",
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
        ....:     "    def __cmp__(left, right):",
        ....:     "        return (<Element>left)._cmp(right)",
        ....:     "    def __richcmp__(left, right, op):",
        ....:     "        return (<Element>left)._richcmp(right,op)",
        ....:     "    cpdef int _cmp_(left, Element right) except -2:",
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
    def __getattr__(self, name):
        """
        This getattr method ensures that cached methods and lazy attributes
        can be inherited from the element class of a category.

        .. NOTE::

            The use of cached methods is demonstrated in the main doc
            string of this class. Here, we demonstrate lazy
            attributes.

        EXAMPLE::

            sage: cython_code = ["from sage.structure.element cimport ElementWithCachedMethod",
            ... "cdef class MyElement(ElementWithCachedMethod):",
            ... "    cdef public object x",
            ... "    def __init__(self,P,x):",
            ... "        self.x=x",
            ... "        ElementWithCachedMethod.__init__(self,P)",
            ... "    def _repr_(self):",
            ... "        return '<%s>'%self.x",
            ... "from sage.structure.parent cimport Parent",
            ... "cdef class MyParent(Parent):",
            ... "    Element = MyElement",
            ... "from sage.all import cached_method, lazy_attribute, Category, Objects",
            ... "class MyCategory(Category):",
            ... "    @cached_method",
            ... "    def super_categories(self):",
            ... "        return [Objects()]",
            ... "    class ElementMethods:",
            ... "        @lazy_attribute",
            ... "        def my_lazy_attr(self):",
            ... "            return 'lazy attribute of <%d>'%self.x"]
            sage: cython('\n'.join(cython_code))
            sage: C = MyCategory()
            sage: P = MyParent(category=C)
            sage: e = MyElement(P,5)
            sage: e.my_lazy_attr
            'lazy attribute of <5>'
            sage: e.my_lazy_attr is e.my_lazy_attr
            True

        """
        if name.startswith('__') and not name.endswith('_'):
            dummy_error_message.cls = type(self)
            dummy_error_message.name = name
            raise dummy_attribute_error
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

    ##################################################
    # Addition
    ##################################################
    def __add__(left, right):
        """
        Top-level addition operator for ModuleElements.

        See extensive documentation at the top of element.pyx.
        """
        # Try fast pathway if they are both ModuleElements and the parents
        # match.

        # (We know at least one of the arguments is a ModuleElement. So if
        # their types are *equal* (fast to check) then they are both
        # ModuleElements. Otherwise use the slower test via isinstance.)
        if have_same_parent_c(left, right):
            # If we hold the only references to this object, we can
            # safely mutate it. NOTE the threshold is different by one
            # for __add__ and __iadd__.
            if (<PyObject *>left).ob_refcnt < inplace_threshold:
                return (<ModuleElement>left)._iadd_(<ModuleElement>right)
            else:
                return (<ModuleElement>left)._add_(<ModuleElement>right)

        global coercion_model
        return coercion_model.bin_op(left, right, add)

    cpdef ModuleElement _add_(left, ModuleElement right):
        raise TypeError(arith_error_message(left, right, add))

    def __iadd__(ModuleElement self, right):
        if have_same_parent_c(self, right):
            if (<PyObject *>self).ob_refcnt <= inplace_threshold:
                return self._iadd_(<ModuleElement>right)
            else:
                return self._add_(<ModuleElement>right)
        else:
            global coercion_model
            return coercion_model.bin_op(self, right, iadd)

    cpdef ModuleElement _iadd_(self, ModuleElement right):
        return self._add_(right)

    ##################################################
    # Subtraction
    ##################################################

    def __sub__(left, right):
        """
        Top-level subtraction operator for ModuleElements.
        See extensive documentation at the top of element.pyx.
        """
        if have_same_parent_c(left, right):
            if (<PyObject *>left).ob_refcnt < inplace_threshold:
                return (<ModuleElement>left)._isub_(<ModuleElement>right)
            else:
                return (<ModuleElement>left)._sub_(<ModuleElement>right)
        global coercion_model
        return coercion_model.bin_op(left, right, sub)

    cpdef ModuleElement _sub_(left, ModuleElement right):
        # default implementation is to use the negation and addition
        # dispatchers:
        return left._add_(-right)

    def __isub__(ModuleElement self, right):
        if have_same_parent_c(self, right):
            if (<PyObject *>self).ob_refcnt <= inplace_threshold:
                return self._isub_(<ModuleElement>right)
            else:
                return self._sub_(<ModuleElement>right)
        else:
            global coercion_model
            return coercion_model.bin_op(self, right, isub)

    cpdef ModuleElement _isub_(self, ModuleElement right):
        return self._sub_(right)

    ##################################################
    # Negation
    ##################################################

    def __neg__(self):
        """
        Top-level negation operator for ModuleElements, which
        may choose to implement _neg_ rather than __neg__ for
        consistency.
        """
        return self._neg_()

    cpdef ModuleElement _neg_(self):
        # default implementation is to try multiplying by -1.
        global coercion_model
        if self._parent._base is None:
            return coercion_model.bin_op(-1, self, mul)
        else:
            return coercion_model.bin_op(self._parent._base(-1), self, mul)

    ##################################################
    # Module element multiplication (scalars, etc.)
    ##################################################
    def __mul__(left, right):
        if PyInt_CheckExact(right):
            return (<ModuleElement>left)._mul_long(PyInt_AS_LONG(right))
        if PyInt_CheckExact(left):
            return (<ModuleElement>right)._mul_long(PyInt_AS_LONG(left))
        if have_same_parent_c(left, right):
            raise TypeError(arith_error_message(left, right, mul))
        # Always do this
        global coercion_model
        return coercion_model.bin_op(left, right, mul)

    def __imul__(left, right):
        if have_same_parent_c(left, right):
             raise TypeError
        # Always do this
        global coercion_model
        return coercion_model.bin_op(left, right, imul)

    # rmul -- left * self
    cpdef ModuleElement _rmul_(self, RingElement left):
        """
        Default module left scalar multiplication, which is to try to
        canonically coerce the scalar to the integers and do that
        multiplication, which is always defined.

        Returning None indicates that this action is not implemented here.
        """
        return None

    # lmul -- self * right

    cpdef ModuleElement _lmul_(self, RingElement right):
        """
        Default module left scalar multiplication, which is to try to
        canonically coerce the scalar to the integers and do that
        multiplication, which is always defined.

        Returning None indicates that this action is not implemented here.
        """
        return None

    cdef ModuleElement _mul_long(self, long n):
        """
        Generic path for multiplying by a C long, assumed to commute.
        """
        return coercion_model.bin_op(self, n, mul)

    cpdef ModuleElement _ilmul_(self, RingElement right):
        return self._lmul_(right)

    cdef RingElement coerce_to_base_ring(self, x):
        if isinstance(x, Element) and (<Element>x)._parent is self._parent._base:
            return x
        try:
            return self._parent._base._coerce_c(x)
        except AttributeError:
            return self._parent._base(x)

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
    # Multiplication
    #############################################################
    def __mul__(left, right):
        """
        Top-level multiplication operator for monoid elements.
        See extensive documentation at the top of element.pyx.
        """
        global coercion_model
        if have_same_parent_c(left, right):
            return (<MonoidElement>left)._mul_(<MonoidElement>right)
        try:
            return coercion_model.bin_op(left, right, mul)
        except TypeError as msg:
            if isinstance(left, (int, long)) and left==1:
                return right
            elif isinstance(right, (int, long)) and right==1:
                return left
            raise


    cpdef MonoidElement _mul_(left, MonoidElement right):
        """
        Cython classes should override this function to implement multiplication.
        See extensive documentation at the top of element.pyx.
        """
        raise TypeError

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

    cpdef ModuleElement _rmul_(self, RingElement left):
        return self._lmul_(left)

    cpdef ModuleElement _lmul_(self, RingElement right):
        """
        Default module left scalar multiplication, which is to try to
        canonically coerce the scalar to the integers and do that
        multiplication, which is always defined.

        Returning None indicates this action is not implemented.
        """
        return None

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

    def _add_(self, x):
        raise ArithmeticError("addition not defined in a multiplicative group")

    def __div__(left, right):
        if have_same_parent_c(left, right):
            return left._div_(right)
        global coercion_model
        return coercion_model.bin_op(left, right, div)

    cpdef MultiplicativeGroupElement _div_(self, MultiplicativeGroupElement right):
        """
        Cython classes should override this function to implement division.
        See extensive documentation at the top of element.pyx.
        """
        return self * ~right

    def __invert__(self):
        r"""
        Return the inverse of ``self``.
        """
        if self.is_one():
            return self
        return self.parent().one()/self


def is_RingElement(x):
    """
    Return ``True`` if x is of type RingElement.
    """
    return isinstance(x, RingElement)

cdef class RingElement(ModuleElement):
    ##################################################
    def is_one(self):
        return self == self._parent.one()

    ##################################
    # Fast long add/sub path.
    ##################################

    def __add__(left, right):
        """
        Top-level addition operator for RingElements.

        See extensive documentation at the top of element.pyx.
        """
        if have_same_parent_c(left, right):
            if (<PyObject *>left).ob_refcnt < inplace_threshold:
                return (<ModuleElement>left)._iadd_(<ModuleElement>right)
            else:
                return (<ModuleElement>left)._add_(<ModuleElement>right)
        if PyInt_CheckExact(right):
            return (<RingElement>left)._add_long(PyInt_AS_LONG(right))
        elif PyInt_CheckExact(left):
            return (<RingElement>right)._add_long(PyInt_AS_LONG(left))
        return coercion_model.bin_op(left, right, add)

    cdef RingElement _add_long(self, long n):
        """
        Generic path for adding a C long, assumed to commute.
        """
        return coercion_model.bin_op(self, n, add)

    def __sub__(left, right):
        """
        Top-level subtraction operator for RingElements.

        See extensive documentation at the top of element.pyx.
        """
        cdef long n
        if have_same_parent_c(left, right):
            if (<PyObject *>left).ob_refcnt < inplace_threshold:
                return (<ModuleElement>left)._isub_(<ModuleElement>right)
            else:
                return (<ModuleElement>left)._sub_(<ModuleElement>right)
        if PyInt_CheckExact(right):
            n = PyInt_AS_LONG(right)
            # See UNARY_NEG_WOULD_OVERFLOW in Python's intobject.c
            if (n == 0) | (<unsigned long>n != 0 - <unsigned long>n):
                return (<RingElement>left)._add_long(-n)
        return coercion_model.bin_op(left, right, sub)

    ##################################
    # Multiplication
    ##################################

    cpdef ModuleElement _lmul_(self, RingElement right):
        # We return None to invoke the default action of coercing into self
        return None

    cpdef ModuleElement _rmul_(self, RingElement left):
        # We return None to invoke the default action of coercing into self
        return None

    def __mul__(left, right):
        """
        Top-level multiplication operator for ring elements.
        See extensive documentation at the top of element.pyx.

        AUTHOR:

        - Gonzalo Tornaria (2007-06-25) - write base-extending test cases and fix them

        TESTS:

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

        """
        # Try fast pathway if they are both RingElements and the parents match.
        # (We know at least one of the arguments is a RingElement. So if their
        # types are *equal* (fast to check) then they are both RingElements.
        # Otherwise use the slower test via isinstance.)
        if have_same_parent_c(left, right):
            if (<PyObject *>left).ob_refcnt < inplace_threshold:
                return (<RingElement>left)._imul_(<RingElement>right)
            else:
                return (<RingElement>left)._mul_(<RingElement>right)
        if PyInt_CheckExact(right):
            return (<ModuleElement>left)._mul_long(PyInt_AS_LONG(right))
        elif PyInt_CheckExact(left):
            return (<ModuleElement>right)._mul_long(PyInt_AS_LONG(left))
        return coercion_model.bin_op(left, right, mul)

    cpdef RingElement _mul_(self, RingElement right):
        """
        Cython classes should override this function to implement multiplication.
        See extensive documentation at the top of element.pyx.
        """
        raise TypeError(arith_error_message(self, right, mul))

    def __imul__(left, right):
        if have_same_parent_c(left, right):
            if (<PyObject *>left).ob_refcnt <= inplace_threshold:
                return (<RingElement>left)._imul_(<RingElement>right)
            else:
                return (<RingElement>left)._mul_(<RingElement>right)

        global coercion_model
        return coercion_model.bin_op(left, right, imul)

    cpdef RingElement _imul_(RingElement self, RingElement right):
        return self._mul_(right)

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

        These aren't testing this code, but they are probably good to have around.

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

        Another example from :trac:`2956`; this should overflow on x32
        and succeed on x64::

            sage: K.<x,y> = ZZ[]
            sage: (x^12345)^54321
            x^670592745                                   # 64-bit
            Traceback (most recent call last):            # 32-bit
            ...                                           # 32-bit
            OverflowError: Exponent overflow (670592745). # 32-bit

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

    ##################################
    # Division
    ##################################

    def __truediv__(self, right):
        # in sage all divs are true
        if not isinstance(self, Element):
            return coercion_model.bin_op(self, right, div)
        return self.__div__(right)

    def __div__(self, right):
        """
        Top-level multiplication operator for ring elements.
        See extensive documentation at the top of element.pyx.
        """
        if have_same_parent_c(self, right):
            if (<PyObject *>self).ob_refcnt < inplace_threshold:
                return (<RingElement>self)._idiv_(<RingElement>right)
            else:
                return (<RingElement>self)._div_(<RingElement>right)
        global coercion_model
        return coercion_model.bin_op(self, right, div)

    cpdef RingElement _div_(self, RingElement right):
        """
        Cython classes should override this function to implement division.
        See extensive documentation at the top of element.pyx.
        """
        try:
            return self._parent.fraction_field()(self, right)
        except AttributeError:
            if not right:
                raise ZeroDivisionError("Cannot divide by zero")
            else:
                raise TypeError(arith_error_message(self, right, div))

    def __idiv__(self, right):
        """
        Top-level division operator for ring elements.
        See extensive documentation at the top of element.pyx.
        """
        if have_same_parent_c(self, right):
            if (<PyObject *>self).ob_refcnt <= inplace_threshold:
                return (<RingElement>self)._idiv_(<RingElement>right)
            else:
                return (<RingElement>self)._div_(<RingElement>right)
        global coercion_model
        return coercion_model.bin_op(self, right, idiv)

    cpdef RingElement _idiv_(RingElement self, RingElement right):
        return self._div_(right)

    def __invert__(self):
        if self.is_one():
            return self
        return 1/self

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
        return self.__abs__()

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

        sage: from sage.rings.commutative_ring_element import is_CommutativeRingElement
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
            ZeroDivisionError: reduction modulo right not defined.

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
        #Check if the parents are the same:

        if have_same_parent_c(self, x):
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
            global coercion_model
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

        from sage.rings.integral_domain import is_IntegralDomain
        P=self._parent
        is_sqr, sq_rt = self.is_square( root = True )
        if is_sqr:
            if all:
                if not is_IntegralDomain(P):
                    raise NotImplementedError('sqrt() with all=True is only implemented for integral domains, not for %s' % P)
                if P.characteristic()==2 or sq_rt==0:
                    #0 has only one square root, and in charasteristic 2 everything also has only 1 root
                    return [ sq_rt ]
                return [ sq_rt, -sq_rt ]
            return sq_rt
        #from now on we know that self is not a square
        if not is_IntegralDomain(P):
            raise NotImplementedError('sqrt() of non squares is only implemented for integral domains, not for %s' % P)
        if not extend:
            #all square roots of a non-square should be an empty list
            if all:
                return []
            raise ValueError('trying to take square root of non-square %s with extend = False' % self)

        if name == None:
            raise TypeError("Polynomial is not a square. You must specify the name of the square root when using the default extend = True")
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        PY = PolynomialRing(P,'y')
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

    def __imul__(left, right):
        if have_same_parent_c(left, right):
            return (<Vector>left)._dot_product_(<Vector>right)
        # Always do this
        global coercion_model
        return coercion_model.bin_op(left, right, imul)


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

        """
        if have_same_parent_c(left, right):
            return (<Vector>left)._dot_product_(<Vector>right)
        # Always do this
        global coercion_model
        return coercion_model.bin_op(left, right, mul)

    cpdef Element _dot_product_(Vector left, Vector right):
        return left._dot_product_coerce_(right)

    cpdef Element _dot_product_coerce_(Vector left, Vector right):
        raise TypeError(arith_error_message(left, right, mul))

    cpdef Vector _pairwise_product_(Vector left, Vector right):
        raise TypeError("unsupported operation for '%s' and '%s'"%(parent_c(left), parent_c(right)))

    def __div__(self, right):
        right = py_scalar_to_element(right)
        if isinstance(right, RingElement):
            # Let __mul__ do the job
            return self.__mul__(~right)
        if isinstance(right, Vector):
            try:
                W = right.parent().submodule([right])
                return W.coordinates(self)[0] / W.coordinates(right)[0]
            except ArithmeticError:
                if right.is_zero():
                    raise ZeroDivisionError("division by zero vector")
                else:
                    raise ArithmeticError("vector is not in free module")
        raise TypeError(arith_error_message(self, right, div))

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

    def __imul__(left, right):
        if have_same_parent_c(left, right):
            return (<Matrix>left)._matrix_times_matrix_(<Matrix>right)
        else:
            global coercion_model
            return coercion_model.bin_op(left, right, imul)

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
        if have_same_parent_c(left, right):
            return (<Matrix>left)._matrix_times_matrix_(<Matrix>right)
        else:
            global coercion_model
            return coercion_model.bin_op(left, right, mul)

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
        cdef Matrix rightinv
        if have_same_parent_c(left, right):
            rightinv = ~<Matrix>right
            if have_same_parent_c(left,rightinv):
                return (<Matrix>left)._matrix_times_matrix_(rightinv)
            else:
                return (<Matrix>(left.change_ring(rightinv.parent().base_ring())))._matrix_times_matrix_(rightinv)
        else:
            global coercion_model
            return coercion_model.bin_op(left, right, div)

    cdef Vector _vector_times_matrix_(matrix_right, Vector vector_left):
        raise TypeError

    cdef Vector _matrix_times_vector_(matrix_left, Vector vector_right):
        raise TypeError

    cdef Matrix _matrix_times_matrix_(left, Matrix right):
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
            from sage.rings.arith import lcm
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

    def __floordiv__(self,right):
        """
        Quotient of division of ``self`` by other.  This is denoted //.
        """
        Q, _ = self.quo_rem(right)
        return Q

    def __mod__(self, other):
        """
        Remainder of division of ``self`` by other.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: x % (x+1)
            -1
            sage: (x**3 + x - 1) % (x**2 - 1)
            2*x - 1
        """
        _, R = self.quo_rem(other)
        return R

def is_FieldElement(x):
    """
    Return ``True`` if x is of type FieldElement.
    """
    return isinstance(x, FieldElement)

cdef class FieldElement(CommutativeRingElement):

    def __floordiv__(self, other):
        return self / other

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

cdef class PlusInfinityElement(InfinityElement):
    pass

cdef class MinusInfinityElement(InfinityElement):
    pass


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
    global coercion_model
    return coercion_model.canonical_coercion(x,y)

cpdef bin_op(x, y, op):
    global coercion_model
    return coercion_model.bin_op(x,y,op)


def coerce(Parent p, x):
    try:
        return p._coerce_c(x)
    except AttributeError:
        return p(x)

def coerce_cmp(x,y):
    from sage.misc.superseded import deprecation
    deprecation(18322, 'the global coerce_cmp() function is deprecated')
    cdef int c
    try:
        x, y = coercion_model.canonical_coercion(x, y)
        return cmp(x,y)
    except TypeError:
        c = cmp(type(x), type(y))
        if c == 0: c = -1
        return c

# We define this base class here to avoid circular cimports.
cdef class CoercionModel:
    """
    Most basic coercion scheme. If it doesn't already match, throw an error.
    """
    cpdef canonical_coercion(self, x, y):
        if parent_c(x) is parent_c(y):
            return x,y
        raise TypeError("no common canonical parent for objects with parents: '%s' and '%s'"%(parent_c(x), parent_c(y)))

    cpdef bin_op(self, x, y, op):
        if parent_c(x) is parent_c(y):
            return op(x,y)
        raise TypeError(arith_error_message(x,y,op))

import coerce
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
            print traceback
    else:
        return coercion_model.exception_stack()


cdef class NamedBinopMethod:
    """
    A decorator to be used on binary operation methods that should operate
    on elements of the same parent. If the parents of the arguments differ,
    coercion is performed, then the method is re-looked up by name on the
    first argument.

    In short, using the ``NamedBinopMethod`` (alias ``coerce_binop``) decorator
    on a method gives it the exact same semantics of the basic arithmetic
    operations like ``_add_``, ``_sub_``, etc. in that both operands are
    guaranteed to have exactly the same parent.
    """
    cdef _self
    cdef _func
    cdef _name

    def __init__(self, func, name=None, obj=None):
        """
        TESTS::

            sage: from sage.structure.element import NamedBinopMethod
            sage: NamedBinopMethod(gcd)(12, 15)
            3
        """
        self._func = func
        if name is None:
            if isinstance(func, types.FunctionType):
                name = func.__name__
            if isinstance(func, types.UnboundMethodType):
                name = func.__func__.__name__
            else:
                name = func.__name__
        self._name = name
        self._self = obj

    def __call__(self, x, y=None, **kwds):
        """
        TESTS::

            sage: from sage.structure.element import NamedBinopMethod
            sage: test_func = NamedBinopMethod(lambda x, y, **kwds: (x, y, kwds), '_add_')
            sage: class test_class(Rational):
            ....:     def __init__(self,value):
            ....:         self.v = value
            ....:     @NamedBinopMethod
            ....:     def test_add(self, other, keyword='z'):
            ....:         return (self.v, other, keyword)

        Calls func directly if the two arguments have the same parent::

            sage: test_func(1, 2)
            (1, 2, {})
            sage: x = test_class(1)
            sage: x.test_add(1/2)
            (1, 1/2, 'z')
            sage: x.test_add(1/2, keyword=3)
            (1, 1/2, 3)

        Passes through coercion and does a method lookup if the
        left operand is not the same::

            sage: test_func(0.5, 1)
            (0.500000000000000, 1.00000000000000, {})
            sage: test_func(1, 2, algorithm='fast')
            (1, 2, {'algorithm': 'fast'})
            sage: test_func(1, 1/2)
            3/2
            sage: x.test_add(2)
            (1, 2, 'z')
            sage: x.test_add(2, keyword=3)
            (1, 2, 3)

        A real example::

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
            sage: h.gcd(f,algorithm='modular')
            1
        """
        if y is None:
            if self._self is None:
                self._func(x, **kwds)
            else:
                x, y = self._self, x
        if not have_same_parent_c(x, y):
            old_x = x
            x,y = coercion_model.canonical_coercion(x, y)
            if old_x is x:
                return self._func(x,y, **kwds)
            else:
                return getattr(x, self._name)(y, **kwds)
        else:
            return self._func(x,y, **kwds)

    def __get__(self, obj, objtype):
        """
        Used to transform from an unbound to a bound method.

        TESTS::
            sage: from sage.structure.element import NamedBinopMethod
            sage: R.<x> = ZZ[]
            sage: isinstance(x.quo_rem, NamedBinopMethod)
            True
            sage: x.quo_rem(x)
            (1, 0)
            sage: type(x).quo_rem(x,x)
            (1, 0)
        """
        return NamedBinopMethod(self._func, self._name, obj)

    def _sage_doc_(self):
        """
        Return the docstring of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.structure.element import NamedBinopMethod
            sage: g = NamedBinopMethod(gcd)
            sage: from sage.misc.sageinspect import sage_getdoc
            sage: sage_getdoc(g) == sage_getdoc(gcd)
            True
        """
        return sageinspect._sage_getdoc_unformatted(self._func)

    def _sage_src_(self):
        """
        Return the source of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.structure.element import NamedBinopMethod
            sage: g = NamedBinopMethod(gcd)
            sage: 'def gcd(' in g._sage_src_()
            True
        """
        return sageinspect.sage_getsource(self._func)

    def _sage_argspec_(self):
        """
        Return the argspec of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.structure.element import NamedBinopMethod
            sage: g = NamedBinopMethod(gcd)
            sage: g._sage_argspec_()
            ArgSpec(args=['a', 'b'], varargs=None, keywords='kwargs', defaults=(None,))
        """
        return sageinspect.sage_getargspec(self._func)

coerce_binop = NamedBinopMethod

###############################################################################

from sage.misc.lazy_import import lazy_import
lazy_import('sage.rings.arith', ['gcd', 'xgcd', 'lcm'], deprecation=10779)


######################

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
