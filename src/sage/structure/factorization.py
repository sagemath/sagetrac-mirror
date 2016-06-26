r"""
Factorizations

The :class:`Factorization` class provides a structure for holding quite
general lists of objects with integer multiplicities.  These may hold
the results of an arithmetic or algebraic factorization, where the
objects may be primes or irreducible polynomials and the
multiplicities are the (non-zero) exponents in the factorization.  For
other types of examples, see below.

:class:`Factorization` class objects contain a ``list``, so can be
printed nicely and be manipulated like a list of prime-exponent pairs,
or easily turned into a plain list.  For example, we factor the
integer `-45`::

    sage: F = factor(-45)

This returns an object of type :class:`Factorization`::

    sage: type(F)
    <class 'sage.structure.factorization_integer.IntegerFactorization'>

It prints in a nice factored form::

    sage: F
    -1 * 3^2 * 5

There is an underlying list representation, which ignores the unit part::

    sage: list(F)
    [(3, 2), (5, 1)]

A :class:`Factorization` is not actually a list::

    sage: isinstance(F, list)
    False

However, we can access the :class:`Factorization` F itself as if it were a list::

    sage: F[0]
    (3, 2)
    sage: F[1]
    (5, 1)

To get at the unit part, use the :meth:`Factorization.unit` function::

    sage: F.unit()
    -1

All factorizations are immutable, up to ordering with ``sort()`` and
simplifying with ``simplify()``.  Thus if you write a function that
returns a cached version of a factorization, you do not have to return
a copy.

::

    sage: F = factor(-12); F
    -1 * 2^2 * 3
    sage: F[0] = (5,4)
    Traceback (most recent call last):
    ...
    TypeError: 'Factorization' object does not support item assignment

EXAMPLES:

This more complicated example involving polynomials also illustrates
that the unit part is not discarded from factorizations::

    sage: x = QQ['x'].0
    sage: f = -5*(x-2)*(x-3)
    sage: f
    -5*x^2 + 25*x - 30
    sage: F = f.factor(); F
    (-5) * (x - 3) * (x - 2)
    sage: F.unit()
    -5
    sage: F.value()
    -5*x^2 + 25*x - 30

The underlying list is the list of pairs `(p_i, e_i)`, where each
`p_i` is a 'prime' and each `e_i` is an integer. The unit part
is discarded by the list::

    sage: list(F)
    [(x - 3, 1), (x - 2, 1)]
    sage: len(F)
    2
    sage: F[1]
    (x - 2, 1)

In the ring `\ZZ[x]`, the integer `-5` is not a unit, so the
factorization has three factors::

    sage: x = ZZ['x'].0
    sage: f = -5*(x-2)*(x-3)
    sage: f
    -5*x^2 + 25*x - 30
    sage: F = f.factor(); F
    (-1) * 5 * (x - 3) * (x - 2)
    sage: F.universe()
    Univariate Polynomial Ring in x over Integer Ring
    sage: F.unit()
    -1
    sage: list(F)
    [(5, 1), (x - 3, 1), (x - 2, 1)]
    sage: F.value()
    -5*x^2 + 25*x - 30
    sage: len(F)
    3

On the other hand, -1 is a unit in `\ZZ`, so it is included in the unit::

    sage: x = ZZ['x'].0
    sage: f = -1*(x-2)*(x-3)
    sage: F = f.factor(); F
    (-1) * (x - 3) * (x - 2)
    sage: F.unit()
    -1
    sage: list(F)
    [(x - 3, 1), (x - 2, 1)]

Factorizations can involve fairly abstract mathematical objects::

    sage: F = ModularSymbols(11,4).factorization()
    sage: F
    (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 6 for Gamma_0(11) of weight 4 with sign 0 over Rational Field) *
    (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 6 for Gamma_0(11) of weight 4 with sign 0 over Rational Field) *
    (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 6 for Gamma_0(11) of weight 4 with sign 0 over Rational Field)
    sage: type(F)
    <class 'sage.structure.factorization.Factorization'>


    sage: K.<a> = NumberField(x^2 + 3); K
    Number Field in a with defining polynomial x^2 + 3
    sage: f = K.factor(15); f
    (Fractional ideal (-a))^2 * (Fractional ideal (5))
    sage: f.universe()
    Monoid of ideals of Number Field in a with defining polynomial x^2 + 3
    sage: f.unit()
    Fractional ideal (1)
    sage: g=K.factor(9); g
    (Fractional ideal (-a))^4
    sage: f.lcm(g)
    (Fractional ideal (-a))^4 * (Fractional ideal (5))
    sage: f.gcd(g)
    (Fractional ideal (-a))^2
    sage: f.is_integral()
    True

TESTS::

    sage: F = factor(-20); F
    -1 * 2^2 * 5
    sage: G = loads(dumps(F)); G
    -1 * 2^2 * 5
    sage: G == F
    True
    sage: G is F
    False

AUTHORS:

- William Stein (2006-01-22): added unit part as suggested by David Kohel.

- William Stein (2008-01-17): wrote much of the documentation and
  fixed a couple of bugs.

- Nick Alexander (2008-01-19): added support for non-commuting factors.

- John Cremona (2008-08-22): added division, lcm, gcd, is_integral and
  universe functions
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from six.moves import range
from six import iteritems, integer_types

from sage.structure.sage_object import SageObject
from sage.structure.element import Element, parent
from sage.structure.parent import Parent
from sage.structure.sequence import Sequence
from sage.structure.richcmp import richcmp_method, richcmp, richcmp_not_equal
from sage.rings.integer import Integer
from sage.misc.all import prod
from sage.misc.cachefunc import cached_method



@richcmp_method
class Factorization(SageObject):
    """
    A formal factorization of an object.

    EXAMPLES::

        sage: N = 2006
        sage: F = N.factor(); F
        2 * 17 * 59
        sage: F.unit()
        1
        sage: F = factor(-2006); F
        -1 * 2 * 17 * 59
        sage: F.unit()
        -1
        sage: loads(F.dumps()) == F
        True
        sage: F = Factorization([(x,1/3)])
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer
    """
    def __init__(self, x, unit=None, universe=None, cr=False, sort=True, simplify=True):
        """
        Create a :class:`Factorization` object.

        INPUT:

        - ``x`` - a list of pairs (p, e) with e an integer;
          otherwise a TypeError is raised

        - ``unit`` - (default: 1) the unit part of the factorization.

        - ``universe`` - (default: None) a common parent for unit and factors.

        - ``cr`` - (default: False) if True, print the factorization with
          carriage returns between factors.

        - ``sort`` - (default: True) if True, sort the factors by calling
          the sort function ``self.sort()`` after creating the factorization

        - ``simplify`` - (default: True) if True, remove duplicate
          factors from the factorization.  See the documentation for
          self.simplify.

        EXAMPLES:

        We create a factorization with all the default options::

            sage: Factorization([(2,3), (5, 1)])
            2^3 * 5

        We create a factorization with a specified unit part::

            sage: Factorization([(2,3), (5, 1)], unit=-1)
            -1 * 2^3 * 5

        We try to create a factorization but with a string an exponent, which
        results in a TypeError::

            sage: Factorization([(2,3), (5, 'x')])
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'x' to an integer

        We create a factorization that puts newlines after each multiply sign
        when printing.  This is mainly useful when the primes are large::

            sage: Factorization([(2,3), (5, 2)], cr=True)
            2^3 *
            5^2

        Another factorization with newlines and nontrivial unit part, which
        appears on a line by itself::

            sage: Factorization([(2,3), (5, 2)], cr=True, unit=-2)
            -2 *
            2^3 *
            5^2

        A factorization, but where we do not sort the factors::

            sage: Factorization([(5,3), (2, 3)], sort=False)
            5^3 * 2^3

        By default, in the commutative case, factorizations are sorted by the
        prime base::

            sage: Factorization([(2, 7), (5,2), (2, 5)])
            2^12 * 5^2
            sage: R.<a,b> = FreeAlgebra(QQ,2)
            sage: Factorization([(a,1),(b,1),(a,2)])
            a * b * a^2

        TESTS:

        Check that list are copied (:trac:`20214`)::

            sage: t0 = (313283749, 2)
            sage: t1 = (111231489, 5)
            sage: l = [t0 ,t1]
            sage: F = Factorization(l, cr=True)
            sage: assert list(F) == [t1,t0] and l[0] is t0 and l[1] is t1

        Check that :trac:`20607` is fixed::

            sage: L.<x,y> = LaurentPolynomialRing(ZZ)
            sage: f = (y + x/y).factor()
            sage: f
            (y^-1) * (y^2 + x)
            sage: f.universe()
            Multivariate Laurent Polynomial Ring in x, y over Integer Ring
            sage: f.unit().parent() is f.universe()
            True
        """
        if not isinstance(x, (tuple, list)):
            raise TypeError("x must be a list")

        xx = []
        for t in x:
            if not isinstance(t, (tuple, list)) or len(t) != 2:
                raise TypeError("x must be a list of pairs (p, e) with e an integer")
            xx.append((t[0], Integer(t[1])))
        x = xx

        if universe is None:
            if unit is None:
                try:
                    universe = Sequence(t[0] for t in x).universe()
                except TypeError:
                    pass
            else:
                try:
                    universe = Sequence([t[0] for t in x] + [unit]).universe()
                except TypeError:
                    pass

        if isinstance(universe, Parent):
            x = [(universe.coerce(i),j) for i,j in x]
            if unit is not None:
                unit = universe.coerce(unit)

        elif universe is not None:
            x = [(universe(i),j) for i,j in x]
            if unit is not None:
                unit = universe(unit)

        if unit is None and universe is not None:
            try:
                unit = universe.one()
            except AttributeError:
                try:
                    unit = universe(1)
                except (ValueError,TypeError):
                    pass

        self.__x = x
        self.__unit = unit
        self.__cr = cr
        self.__universe = universe

        if sort and self.is_commutative():
            self.sort()
        if simplify:
            self.simplify()

    def __getitem__(self, i):
        """
        Return `i^{th}` factor of self.

        EXAMPLES::

            sage: a = factor(-75); a
            -1 * 3 * 5^2
            sage: a[0]
            (3, 1)
            sage: a[1]
            (5, 2)
            sage: a[-1]
            (5, 2)
            sage: a[5]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self.__x[i]

    def __setitem__(self, i, v):
        """
        Set the `i^{th}` factor of self.

        .. warning::

           NOT ALLOWED -- Factorizations are immutable.

        EXAMPLES::

            sage: a = factor(-75); a
            -1 * 3 * 5^2
            sage: a[0] = (2,3)
            Traceback (most recent call last):
            ...
            TypeError: 'Factorization' object does not support item assignment
        """
        raise TypeError("'Factorization' object does not support item assignment")

    def __len__(self):
        """
        Return the number of prime factors of self, not counting
        the unit part.

        EXAMPLES::

            sage: len(factor(15))
            2

        Note that the unit part is not included in the count::

            sage: a = factor(-75); a
            -1 * 3 * 5^2
            sage: len(a)
            2
            sage: list(a)
            [(3, 1), (5, 2)]
            sage: len(list(a))
            2
        """
        return len(self.__x)

    def __richcmp__(self, other, op):
        """
        Compare ``self`` and ``other``.

        This first compares the values.

        If values are equal, this compares the units.

        If units are equal, this compares the underlying lists of
        ``self`` and ``other``.

        EXAMPLES:

        We compare two contrived formal factorizations::

            sage: a = Factorization([(2, 7), (5,2), (2, 5)])
            sage: b = Factorization([(2, 7), (5,10), (7, 3)])
            sage: a
            2^12 * 5^2
            sage: b
            2^7 * 5^10 * 7^3
            sage: a < b
            True
            sage: b < a
            False
            sage: a.value()
            102400
            sage: b.value()
            428750000000

        We compare factorizations of some polynomials::

            sage: x = polygen(QQ)
            sage: x^2 - 1 > x^2 - 4
            True
            sage: factor(x^2 - 1) > factor(x^2 - 4)
            True
        """
        if not isinstance(other, Factorization):
            return NotImplemented

        lx = self.value()
        rx = other.value()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = self.__unit
        rx = other.__unit
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        return richcmp(self.__x, other.__x, op)

    def __copy__(self):
        r"""
        Return a copy of self.

        This is *not* a deepcopy -- only references to the factors are
        returned, not copies of them.  Use ``deepcopy(self)`` if you need
        a deep copy of self.

        EXAMPLES:

        We create a factorization that has mutable primes::

            sage: F = Factorization([([1,2], 5), ([5,6], 10)]); F
            ([1, 2])^5 * ([5, 6])^10

        We make a copy of it::

            sage: G = copy(F); G
            ([1, 2])^5 * ([5, 6])^10
            sage: G is F
            False

        Note that if we change one of the mutable "primes" of F, this does
        not change G::

            sage: F[1][0][0] = 'hello'
            sage: G
            ([1, 2])^5 * ([5, 6])^10
        """
        # No need to sort, since the factorization is already sorted
        # in whatever order is desired.
        return Factorization(self.__x, unit=self.__unit,
                universe=self.__universe, cr=self.__cr,
                sort=False, simplify=False)

    def __deepcopy__(self, memo):
        r"""
        Return a deep copy of self.

        EXAMPLES:

        We make a factorization that has mutable entries::

            sage: F = Factorization([([1,2], 5), ([5,6], 10)]); F
            ([1, 2])^5 * ([5, 6])^10

        Now we make a copy of it and a deep copy::

            sage: K = copy(F)
            sage: G = deepcopy(F); G
            ([1, 2])^5 * ([5, 6])^10

        We change one of the mutable entries of F::

            sage: F[0][0][0] = 10

        This of course changes F::

            sage: F
            ([10, 2])^5 * ([5, 6])^10

        It does not change the copy K of F::

            sage: K
            ([1, 2])^5 * ([5, 6])^10

        It does *not* change the deep copy G::

            sage: G
            ([1, 2])^5 * ([5, 6])^10
        """
        import copy
        return Factorization(copy.deepcopy(list(self), memo),
                             cr=self.__cr, sort=False, simplify=False)

    def universe(self):
        r"""
        Return the parent structure of my factors.

        .. note::

           This used to be called ``base_ring``, but the universe
           of a factorization need not be a ring.

        EXAMPLES::

            sage: F = factor(2006)
            sage: F.universe()
            Integer Ring

            sage: F = Factorization([(1/3, 2)], 3)
            sage: (F*F^-1).universe()
            Rational Field

            sage: F = ModularSymbols(11,4).factorization()
            sage: F.universe()
        """
        return self.__universe

    def base_change(self, U):
        """
        Return the factorization self, with its factors (including the
        unit part) coerced into the universe `U`.

        EXAMPLES::

            sage: F = factor(2006)
            sage: F.universe()
            Integer Ring
            sage: P.<x> = ZZ[]
            sage: F.base_change(P).universe()
            Univariate Polynomial Ring in x over Integer Ring

        This method will return a TypeError if the coercion is not
        possible::

            sage: g = x^2 - 1
            sage: F = factor(g); F
            (x - 1) * (x + 1)
            sage: F.universe()
            Univariate Polynomial Ring in x over Integer Ring
            sage: F.base_change(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: Impossible to coerce the factors of (x - 1) * (x + 1) into Integer Ring

        TESTS:

        Check that the method properly converts units (:trac:`20214`)::

            sage: f = 1.factor().base_change(AA)
            sage: f.universe() is AA and f.unit().parent() is AA
            True
            sage: f = ZZ['x'].gen().factor().base_change(QQ['x'])
            sage: f.universe() is QQ['x'] and f.unit().parent() is QQ['x']
            True
        """
        try:
            return Factorization(self.__x, unit=self.__unit, universe=U)
        except TypeError:
            raise TypeError("Impossible to coerce the factors of %s into %s"%(self, U))

    def is_commutative(self):
        """
        Return True if my factors commute.

        EXAMPLES::

            sage: F = factor(2006)
            sage: F.is_commutative()
            True
            sage: K = QuadraticField(23, 'a')
            sage: F = K.factor(13)
            sage: F.is_commutative()
            True

            sage: R = GL(2, GF(3))
            sage: a,b = R.gens()
            sage: F = Factorization([(a, 2)])
            sage: F.is_commutative()
            False
            sage: (F*F^-1).is_commutative()
            False
        """
        try:
            return self.universe().is_commutative()
        except Exception:
            # This is not the mathematically correct default, but agrees with
            # history -- we've always assumed factored things commute
            return True

    def _set_cr(self, cr):
        """
        Change whether or not the factorization is printed with
        carriage returns after each factor.

        EXAMPLES::

            sage: x = polygen(QQ,'x')
            sage: F = factor(x^6 - 1); F
            (x - 1) * (x + 1) * (x^2 - x + 1) * (x^2 + x + 1)
            sage: F._set_cr(True); F
            (x - 1) *
            (x + 1) *
            (x^2 - x + 1) *
            (x^2 + x + 1)
            sage: F._set_cr(False); F
            (x - 1) * (x + 1) * (x^2 - x + 1) * (x^2 + x + 1)
        """
        self.__cr = bool(cr)

    def simplify(self):
        """
        Combine adjacent products as much as possible.

        TESTS::

            sage: R.<x,y> = FreeAlgebra(ZZ, 2)
            sage: F = Factorization([(x,3), (y, 2), (y,2)], simplify=False); F
            x^3 * y^2 * y^2
            sage: F.simplify(); F
            x^3 * y^4
            sage: F * Factorization([(y, -2)], 2)
            (2) * x^3 * y^2
        """
        repeat = False
        simp = []
        import itertools
        for obj, agroup in itertools.groupby(list(self), lambda x: x[0]):
            xs = list(agroup)
            if len(xs) > 1:
                repeat = True
            n = sum([x[1] for x in xs])
            if n != 0:
                simp.append((obj, n))
        self.__x[0:] = simp
        if repeat:
            self.simplify()

    def sort(self, key=None):
        r"""
        Sort the factors in this factorization.

        INPUT:

        - ``key`` - (default: ``None``) comparison key

        OUTPUT:

        - changes this factorization to be sorted (inplace)

        EXAMPLES:

        We create a factored polynomial::

            sage: x = polygen(QQ,'x')
            sage: F = factor(x^3 + 1); F
            (x + 1) * (x^2 - x + 1)

        We sort it by decreasing degree::

            sage: F.sort(key=lambda x:(-x[0].degree(), x))
            sage: F
            (x^2 - x + 1) * (x + 1)
        """
        if self.__x:
            if key is not None:
                self.__x.sort(key=key)
            else:
                self.__x.sort()

    def unit(self):
        r"""
        Return the unit part of this factorization.

        EXAMPLES:

        We create a polynomial over the real double field and factor it::

            sage: x = polygen(RDF, 'x')
            sage: F = factor(-2*x^2 - 1); F
            (-2.0) * (x^2 + 0.5000000000000001)

        Note that the unit part of the factorization is `-2.0`::

            sage: F.unit()
            -2.0

            sage: F = factor(-2006); F
            -1 * 2 * 17 * 59
            sage: F.unit()
            -1
        """
        return self.__unit

    def _cr(self):
        """
        Return whether or not factorizations are printed with carriage
        returns between factors.

        EXAMPLES:

        Our first example involves factoring an integer::

            sage: F = factor(-93930); F
            -1 * 2 * 3 * 5 * 31 * 101
            sage: F._cr()
            False
            sage: F._set_cr(True)
            sage: F._cr()
            True

        This of course looks funny::

            sage: F
            -1 *
            2 *
            3 *
            5 *
            31 *
            101

        Next we factor a modular symbols space::

            sage: F = ModularSymbols(11).factor(); F
            (Modular Symbols subspace of dimension 1 of ...) *
            (Modular Symbols subspace of dimension 1 of ...) *
            (Modular Symbols subspace of dimension 1 of ...)
        """
        try:
            return self.__cr
        except AttributeError:
            self.__cr = False
            return False

    def _repr_(self):
        """
        Return the string representation of this factorization.

        EXAMPLES::

            sage: f = factor(-100); f
            -1 * 2^2 * 5^2
            sage: f._repr_()
            '-1 * 2^2 * 5^2'

        Note that the default printing of a factorization can be overloaded
        using the rename method::

            sage: f.rename('factorization of -100')
            sage: f
            factorization of -100

        However _repr_ always prints normally::

            sage: f._repr_()
            '-1 * 2^2 * 5^2'

        EXAMPLES::

           sage: x = polygen(QQ)
           sage: Factorization([(x-1,1), (x-2,2)])
           (x - 2)^2 * (x - 1)
           sage: Factorization([(x + 1, -3)])
           (x + 1)^-3
        """
        cr = self._cr()
        unit = self.unit()

        if len(self) == 0:
            if unit is None:
                return ''
            else:
                return repr(unit)

        x = self.__x[0][0]
        try:
            atomic = (isinstance(x, integer_types) or
                      self.universe()._repr_option('element_is_atomic'))
        except AttributeError:
            atomic = False

        l = []
        if unit is not None:
            try:
                unit_is_one = unit.is_one()
            except AttributeError:
                try:
                    unit_is_one = x == self.universe()(1)
                except (ValueError,TypeError):
                    unit_is_one = False
            if not unit_is_one:
                l.append('{!r}'.format(unit) if atomic else '({!r})'.format(unit))
        else:
            unit_is_one = True

        for t,n in self:
            t = repr(t)
            if not atomic and (n != 1 or len(self) > 1 or not unit_is_one):
                if '+' in t or '-' in t or ' ' in t:
                    t = '(%s)'%t
            if n != 1:
                t += '^%s'%n
            l.append(t)

        mul =  ' * '
        if cr:
            mul += '\n'
        return mul.join(l)

    def _latex_(self):
        r"""
        Return the LaTeX representation of this factorization.

        EXAMPLES::

            sage: f = factor(-100); f
            -1 * 2^2 * 5^2
            sage: latex(f)
            -1 \cdot 2^{2} \cdot 5^{2}
            sage: f._latex_()
            '-1 \\cdot 2^{2} \\cdot 5^{2}'
            sage: x = AA['x'].0; factor(x^2 + x + 1)._latex_() # trac 12178
            '(x^{2} + x + 1.000000000000000?)'
        """
        if len(self) == 0:
            return self.__unit._latex_()
        try:
            atomic = (isinstance(self.__x[0][0], integer_types) or
                      self.universe()._repr_option('element_is_atomic'))
        except AttributeError:
            atomic = False
        s = ''
        for i in range(len(self)):
            t = self.__x[i][0]._latex_()
            if not atomic and ('+' in t or '-' in t or ' ' in t):
                t = '(%s)'%t
            n = self.__x[i][1]
            if n != 1:
                t += '^{%s}'%n
            s += t
            if i < len(self)-1:
                s += ' \\cdot '
        if self.__unit != 1:
            if atomic:
                u = self.__unit._latex_()
            else:
                u = '\\left(%s\\right)'%self.__unit._latex_()
            s =  u + ' \\cdot ' + s
        return s

    @cached_method
    def __pari__(self):
        """
        Return the PARI factorization matrix corresponding to ``self``.

        EXAMPLES::

            sage: f = factor(-24)
            sage: pari(f)
            [-1, 1; 2, 3; 3, 1]

            sage: R.<x> = QQ[]
            sage: g = factor(x^10 - 1)
            sage: pari(g)
            [x - 1, 1; x + 1, 1; x^4 - x^3 + x^2 - x + 1, 1; x^4 + x^3 + x^2 + x + 1, 1]

        """
        from sage.libs.pari.all import pari
        from itertools import chain

        n = len(self)
        if self.__unit == 1:
            init = ()
        else:
            init = (self.__unit, 1)
            n += 1
        # concatenate (p, e) tuples
        entries = init + tuple(chain.from_iterable(self))
        return pari.matrix(n, 2, entries)

    def __add__(self, other):
        """
        Return the (unfactored) sum of self and other.

        EXAMPLES::

            sage: factor(-10) + 16
            6
            sage: factor(10) - 16
            -6
            sage: factor(100) + factor(19)
            119
        """
        if isinstance(other, Factorization):
            other = other.value()
        return self.value() + other

    def __sub__(self, other):
        """
        Return the (unfactored) difference of self and other.

        EXAMPLES::

            sage: factor(-10) + 16
            6
            sage: factor(10) - 16
            -6
        """
        if isinstance(other, Factorization):
            other = other.value()
        return self.value() - other

    def __radd__(self, left):
        """
        Return the (unfactored) sum of self and left.

        EXAMPLES::

            sage: 16 + factor(-10)
            6
        """
        return self.value() + left


    def __rsub__(self, left):
        """
        Return the (unfactored) difference of left and self.

        EXAMPLES::

            sage: 16 - factor(10)
            6
        """
        return left - self.value()

    def __neg__(self):
        """
        Return negative of this factorization.

        EXAMPLES::

            sage: a = factor(-75); a
            -1 * 3 * 5^2
            sage: -a
            3 * 5^2
            sage: (-a).unit()
            1

        TESTS:

        Check error when there is no unit (:trac:`20214`)::

            sage: fac = Factorization([[(1,2,3),2], [(1,1),3]])
            sage: fac.universe()
            <type 'tuple'>
            sage: -fac
            Traceback (most recent call last):
            ...
            ValueError: can not take negative of factorization without unit
        """
        unit = self.unit()
        if unit is None:
            raise ValueError("can not take negative of factorization without unit")
        unit = -unit
        return Factorization(list(self), unit, self.__universe, self.__cr,
                             sort=False, simplify=False)

    def __rmul__(self, left):
        """
        Return the product left * self, where left is not a Factorization.

        EXAMPLES::

            sage: a = factor(15); a
            3 * 5
            sage: -2 * a
            -2 * 3 * 5
            sage: a * -2
            -2 * 3 * 5
            sage: R.<x,y> = FreeAlgebra(QQ,2)
            sage: f = Factorization([(x,2),(y,3)]); f
            x^2 * y^3
            sage: x * f
            x^3 * y^3
            sage: f * x
            x^2 * y^3 * x

        Note that this does not automatically factor ``left``::

            sage: F = Factorization([(5,3), (2,3)])
            sage: 46 * F
            2^3 * 5^3 * 46
        """
        return Factorization([(left, 1)]) * self

    def __mul__(self, other):
        r"""
        Return the product of two factorizations, which is obtained by
        combining together like factors.

        If the two factorizations have different universes, this
        method will attempt to find a common universe for the
        product.  A TypeError is raised if this is impossible.

        EXAMPLES::

            sage: factor(-10) * factor(-16)
            2^5 * 5
            sage: factor(-10) * factor(16)
            -1 * 2^5 * 5

            sage: R.<x,y> = FreeAlgebra(ZZ, 2)
            sage: F = Factorization([(x,3), (y, 2), (x,1)]); F
            x^3 * y^2 * x
            sage: F*F
            x^3 * y^2 * x^4 * y^2 * x
            sage: -1 * F
            (-1) * x^3 * y^2 * x

            sage: P.<x> = ZZ[]
            sage: f = 2*x + 2
            sage: c = f.content(); g = f//c
            sage: Fc = factor(c); Fc.universe()
            Integer Ring
            sage: Fg = factor(g); Fg.universe()
            Univariate Polynomial Ring in x over Integer Ring
            sage: F = Fc * Fg; F.universe()
            Univariate Polynomial Ring in x over Integer Ring
            sage: [type(a[0]) for a in F]
            [<... 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>,
             <... 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>]

        Check that it only depends on the universe (:trac:`20214`)::

            sage: (Factorization([], universe=ZZ) * Factorization([], universe=QQ)).universe()
            Rational Field
            sage: (Factorization([], universe=QQ) * Factorization([], universe=ZZ)).universe()
            Rational Field
        """
        if not isinstance(other, Factorization):
            return self * Factorization([(other, 1)])

        from sage.structure.element import get_coercion_model
        cm = get_coercion_model()
        su = self.universe()
        ou = other.universe()
        if su is not None and ou is not None:
            try:
                universe = cm.common_parent(su, ou)
            except TypeError:
                raise TypeError("Cannot multiply %s and %s because they cannot be coerced into a common universe"%(self,other))

        if self.is_commutative() and other.is_commutative():
            d1 = dict(self)
            d2 = dict(other)
            s = {}
            for a in set(d1).union(set(d2)):
                s[a] = d1.get(a,0) + d2.get(a,0)
            return Factorization(list(iteritems(s)),
                    unit=self.unit()*other.unit(),
                    universe=universe)
        else:
            return Factorization(list(self) + list(other),
                    unit=self.unit()*other.unit(),
                    universe=universe)

    def __pow__(self, n):
        """
        Return the `n^{th}` power of a factorization, which is got by
        combining together like factors.

        EXAMPLES::

            sage: f = factor(-100); f
            -1 * 2^2 * 5^2
            sage: f^3
            -1 * 2^6 * 5^6
            sage: f^4
            2^8 * 5^8

            sage: K.<a> = NumberField(x^3 - 39*x - 91)
            sage: F = K.factor(7); F
            (Fractional ideal (7, a)) * (Fractional ideal (7, a + 2)) * (Fractional ideal (7, a - 2))
            sage: F^9
            (Fractional ideal (7, a))^9 * (Fractional ideal (7, a + 2))^9 * (Fractional ideal (7, a - 2))^9

            sage: R.<x,y> = FreeAlgebra(ZZ, 2)
            sage: F = Factorization([(x,3), (y, 2), (x,1)]); F
            x^3 * y^2 * x
            sage: F**2
            x^3 * y^2 * x^4 * y^2 * x

            sage: (ZZ.one().factor()**-3).universe()
            Integer Ring
        """
        if not isinstance(n, Integer):
            try:
                n = Integer(n)
            except TypeError:
                raise TypeError("Exponent n (= %s) must be an integer." % n)
        universe = self.universe()
        if n.is_zero():
            return Factorization([], universe=universe)
        elif n.is_one():
            return self
        elif n < 0:
            self = ~self
            n = -n

        if self.is_commutative():
            unit = self.unit()**n
            terms = [(p, n*e) for p, e in self]
            try:
                unit = universe(unit)
                return Factorization(terms, unit=unit, universe=universe, cr=self.__cr, sort=False, simplify=False)
            except (ValueError,TypeError):
                return Factorization(terms, unit=unit, cr=self.__cr, sort=False, simplify=False)

        from sage.arith.power import generic_power
        return generic_power(self, n)

    def __invert__(self):
        r"""
        Return the formal inverse of the factors in the factorization.

        EXAMPLES::

            sage: F = factor(2006); F
            2 * 17 * 59
            sage: F^-1
            2^-1 * 17^-1 * 59^-1
            sage: (F^-1).universe()
            Integer Ring

            sage: (~ZZ.one().factor()).universe()
            Integer Ring
            sage: (ZZ.one().factor()^-1).universe()
            Integer Ring

        This is a fake example where the unit is actually not a unit. The
        universe gets changed::

            sage: F = Factorization([(2,3), (3, 2), (5,1)], 2); F
            2 * 2^3 * 3^2 * 5
            sage: ~F
            1/2 * 2^-3 * 3^-2 * 5^-1
            sage: F^-1
            1/2 * 2^-3 * 3^-2 * 5^-1
            sage: ~F == F^-1
            True
            sage: (F^-1).universe()
            Rational Field
            sage: (~F).universe()
            Rational Field
        """
        universe = self.universe()
        unit = self.unit()

        if unit is not None:
            unit_is_one = False
            try:
                unit_is_one = unit.is_one()
            except AttributeError:
                try:
                    one = universe.one()
                except AttributeError:
                    one = universe(1)
                else:
                    one = None

                if one is not None:
                    unit_is_one = unit == one

            if not unit_is_one:
                try:
                    unit = unit.inverse_of_unit()
                except AttributeError:
                    unit = unit ** (-1)
                    try:
                        unit = universe(unit)
                    except (ValueError, TypeError):
                        universe = parent(unit)
                except ArithmeticError:
                    unit = unit ** (-1)
                    universe = parent(unit)

        return Factorization([(p,-e) for p,e in reversed(self)],
                unit=unit, universe=universe, cr=self._cr())

    def __truediv__(self, other):
        r"""
        Return the quotient of two factorizations, which is obtained by
        multiplying the first by the inverse of the second.

        EXAMPLES::

            sage: factor(-10) / factor(-16)
            2^-3 * 5
            sage: factor(-10) / factor(16)
            -1 * 2^-3 * 5

            sage: R.<x,y> = FreeAlgebra(QQ, 2)
            sage: F = Factorization([(x,3), (y, 2), (x,1)]); F
            x^3 * y^2 * x
            sage: G = Factorization([(y, 1), (x,1)],1); G
            y * x
            sage: F / G
            x^3 * y
        """
        if not isinstance(other, Factorization):
            return self / Factorization([(other, 1)])
        return self * other**-1

    __div__ = __truediv__

    def value(self):
        """
        Return the product of the factors in the factorization, multiplied out.

        EXAMPLES::

            sage: F = factor(-2006); F
            -1 * 2 * 17 * 59
            sage: F.value()
            -2006

            sage: R.<x,y> = FreeAlgebra(ZZ, 2)
            sage: F = Factorization([(x,3), (y, 2), (x,1)]); F
            x^3 * y^2 * x
            sage: F.value()
            x^3*y^2*x
        """
        return prod([p**e for p, e in self.__x], self.__unit)

    # Two aliases for ``value(self)``.
    expand = value
    prod   = value

    def gcd(self, other):
        r"""
        Return the gcd of two factorizations.

        If the two factorizations have different universes, this
        method will attempt to find a common universe for the
        gcd.  A TypeError is raised if this is impossible.

        EXAMPLES::

            sage: factor(-30).gcd(factor(-160))
            2 * 5
            sage: factor(gcd(-30,160))
            2 * 5

            sage: R.<x> = ZZ[]
            sage: (factor(-20).gcd(factor(5*x+10))).universe()
            Univariate Polynomial Ring in x over Integer Ring
        """
        if not isinstance(other, Factorization):
            raise NotImplementedError("can't take gcd of factorization and non-factorization")

        if len(self) and len(other):
            try:
                # first get the two factorizations to have the same
                # universe
                U = Sequence([self[0][0], other[0][0]]).universe()
                self = self.base_change(U)
                other = other.base_change(U)
            except TypeError:
                raise TypeError("Cannot take the gcd of %s and %s because they cannot be coerced into a common universe"%(self,other))

        if self.is_commutative() and other.is_commutative():
            d1 = dict(self)
            d2 = dict(other)
            s = {}
            for a in set(d1).intersection(set(d2)):
                s[a] = min(d1[a],d2[a])
            return Factorization(list(iteritems(s)))
        else:
            raise NotImplementedError("gcd is not implemented for non-commutative factorizations")

    def lcm(self, other):
        r"""
        Return the lcm of two factorizations.

        If the two factorizations have different universes, this
        method will attempt to find a common universe for the
        lcm.  A TypeError is raised if this is impossible.

        EXAMPLES::

            sage: factor(-10).lcm(factor(-16))
            2^4 * 5
            sage: factor(lcm(-10,16))
            2^4 * 5

            sage: R.<x> = ZZ[]
            sage: (factor(-20).lcm(factor(5*x+10))).universe()
            Univariate Polynomial Ring in x over Integer Ring
        """
        if not isinstance(other, Factorization):
            raise NotImplementedError("can't take lcm of factorization and non-factorization")

        if len(self) and len(other):
            try:
                # first get the two factorizations to have the same
                # universe
                U = Sequence([self[0][0], other[0][0]]).universe()
                self = self.base_change(U)
                other = other.base_change(U)
            except TypeError:
                raise TypeError("Cannot take the lcm of %s and %s because they cannot be coerced into a common universe"%(self,other))

        if self.is_commutative() and other.is_commutative():
            d1 = dict(self)
            d2 = dict(other)
            s = {}
            for a in set(d1).union(set(d2)):
                s[a] = max(d1.get(a,0),d2.get(a,0))
            return Factorization(list(iteritems(s)))
        else:
            raise NotImplementedError("lcm is not implemented for non-commutative factorizations")

    def is_integral(self):
        r"""
        Return True iff all exponents of this Factorization are non-negative.

        EXAMPLES::

            sage: F = factor(-10); F
            -1 * 2 * 5
            sage: F.is_integral()
            True

            sage: F = factor(-10) / factor(16); F
            -1 * 2^-3 * 5
            sage: F.is_integral()
            False

        """
        return all([e >=0 for p,e in self.__x])

    def radical(self):
        """
        Return the factorization of the radical of the value of self.

        First, check that all exponents in the factorization are
        positive, raise ValueError otherwise.  If all exponents are
        positive, return self with all exponents set to 1 and with the
        unit set to 1.

        EXAMPLES::

            sage: F = factor(-100); F
            -1 * 2^2 * 5^2
            sage: F.radical()
            2 * 5
            sage: factor(1/2).radical()
            Traceback (most recent call last):
            ...
            ValueError: All exponents in the factorization must be positive.
        """
        if not all([e > 0 for p,e in self.__x]):
            raise ValueError("All exponents in the factorization must be positive.")
        return Factorization([(p,1) for p,e in self.__x], unit=self.unit().parent()(1), cr=self.__cr, sort=False, simplify=False)

    def radical_value(self):
        """
        Return the product of the prime factors in self.

        First, check that all exponents in the factorization are
        positive, raise ValueError otherwise.  If all exponents are
        positive, return the product of the prime factors in self.
        This should be functionally equivalent to
        self.radical().value()

        EXAMPLES::

            sage: F = factor(-100); F
            -1 * 2^2 * 5^2
            sage: F.radical_value()
            10
            sage: factor(1/2).radical_value()
            Traceback (most recent call last):
            ...
            ValueError: All exponents in the factorization must be positive.
        """
        if not all([e > 0 for p,e in self.__x]):
            raise ValueError("All exponents in the factorization must be positive.")
        return prod([p for p,e in self.__x])

