"""
Free algebras

AUTHORS:

- David Kohel (2005-09)

- William Stein (2006-11-01): add all doctests; implemented many
  things.

- Simon King (2011-04): Put free algebras into the category framework.
  Reimplement free algebra constructor, using a
  :class:`~sage.structure.factory.UniqueFactory` for handling
  different implementations of free algebras. Allow degree weights
  for free algebras in letterplace implementation.

EXAMPLES::

    sage: F = FreeAlgebra(ZZ,3,'x,y,z')
    sage: F.base_ring()
    Integer Ring
    sage: G = FreeAlgebra(F, 2, 'm,n'); G
    Free Algebra on 2 generators (m, n) over Free Algebra on 3 generators (x, y, z) over Integer Ring
    sage: G.base_ring()
    Free Algebra on 3 generators (x, y, z) over Integer Ring

The above free algebra is based on a generic implementation. By
:trac:`7797`, there is a different implementation
:class:`~sage.algebras.letterplace.free_algebra_letterplace.FreeAlgebra_letterplace`
based on Singular's letterplace rings. It is currently restricted to
weighted homogeneous elements and is therefore not the default. But the
arithmetic is much faster than in the generic implementation.
Moreover, we can compute Groebner bases with degree bound for its
two-sided ideals, and thus provide ideal containment tests::

    sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')
    sage: F
    Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
    sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F
    sage: I.groebner_basis(degbound=4)
    Twosided Ideal (y*z*y*y - y*z*y*z + y*z*z*y - y*z*z*z, y*z*y*x + y*z*y*z + y*z*z*x + y*z*z*z, y*y*z*y - y*y*z*z + y*z*z*y - y*z*z*z, y*y*z*x + y*y*z*z + y*z*z*x + y*z*z*z, y*y*y - y*y*z + y*z*y - y*z*z, y*y*x + y*y*z + y*z*x + y*z*z, x*y + y*z, x*x - y*x - y*y - y*z) of Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
    sage: y*z*y*y*z*z + 2*y*z*y*z*z*x + y*z*y*z*z*z - y*z*z*y*z*x + y*z*z*z*z*x in I
    True

Positive integral degree weights for the letterplace implementation
was introduced in trac ticket #...::

    sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace', degrees=[2,1,3])
    sage: x.degree()
    2
    sage: y.degree()
    1
    sage: z.degree()
    3
    sage: I = F*[x*y-y*x, x^2+2*y*z, (x*y)^2-z^2]*F
    sage: Q.<a,b,c> = F.quo(I)
    sage: TestSuite(Q).run()
    sage: a^2*b^2
    c*c

TESTS::

    sage: F = FreeAlgebra(GF(5),3,'x')
    sage: TestSuite(F).run()
    sage: F is loads(dumps(F))
    True
    sage: F = FreeAlgebra(GF(5),3,'x', implementation='letterplace')
    sage: TestSuite(F).run()
    sage: F is loads(dumps(F))
    True

::

    sage: F.<x,y,z> = FreeAlgebra(GF(5),3)
    sage: TestSuite(F).run()
    sage: F is loads(dumps(F))
    True
    sage: F.<x,y,z> = FreeAlgebra(GF(5),3, implementation='letterplace')
    sage: TestSuite(F).run()
    sage: F is loads(dumps(F))
    True

::

    sage: F = FreeAlgebra(GF(5),3, ['xx', 'zba', 'Y'])
    sage: TestSuite(F).run()
    sage: F is loads(dumps(F))
    True
    sage: F = FreeAlgebra(GF(5),3, ['xx', 'zba', 'Y'], implementation='letterplace')
    sage: TestSuite(F).run()
    sage: F is loads(dumps(F))
    True

::

    sage: F = FreeAlgebra(GF(5),3, 'abc')
    sage: TestSuite(F).run()
    sage: F is loads(dumps(F))
    True
    sage: F = FreeAlgebra(GF(5),3, 'abc', implementation='letterplace')
    sage: TestSuite(F).run()
    sage: F is loads(dumps(F))
    True

::

    sage: F = FreeAlgebra(FreeAlgebra(ZZ,2,'ab'), 2, 'x')
    sage: TestSuite(F).run()
    sage: F is loads(dumps(F))
    True

Note that the letterplace implementation can only be used if the corresponding
(multivariate) polynomial ring has an implementation in Singular::

    sage: FreeAlgebra(FreeAlgebra(ZZ,2,'ab'), 2, 'x', implementation='letterplace')
    Traceback (most recent call last):
    ...
    NotImplementedError: The letterplace implementation is not available for the free algebra you requested

"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#  Copyright (C) 2005,2006 William Stein <wstein@gmail.com>
#  Copyright (C) 2011 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.rings import Rings

from sage.monoids.free_monoid import FreeMonoid
from sage.monoids.free_monoid_element import FreeMonoidElement

from sage.algebras.free_algebra_element import FreeAlgebraElement
from sage.algebras.pbw_algebra import PBWBasisOfFreeAlgebra

import sage.structure.parent_gens

from sage.structure.factory import UniqueFactory
from sage.misc.cachefunc import cached_method
from sage.all import PolynomialRing
from sage.rings.ring import Algebra
from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
from sage.rings.noncommutative_ideals import Ideal_nc
from sage.rings.infinity import infinity
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from sage.combinat.words.word import Word

class FreeAlgebraFactory(UniqueFactory):
    """
    A constructor of free algebras.

    See :mod:`~sage.algebras.free_algebra` for examples and corner cases.

    EXAMPLES::

        sage: FreeAlgebra(GF(5),3,'x')
        Free Algebra on 3 generators (x0, x1, x2) over Finite Field of size 5
        sage: F.<x,y,z> = FreeAlgebra(GF(5),3)
        sage: (x+y+z)^2
        x^2 + x*y + x*z + y*x + y^2 + y*z + z*x + z*y + z^2
        sage: FreeAlgebra(GF(5),3, 'xx, zba, Y')
        Free Algebra on 3 generators (xx, zba, Y) over Finite Field of size 5
        sage: FreeAlgebra(GF(5),3, 'abc')
        Free Algebra on 3 generators (a, b, c) over Finite Field of size 5
        sage: FreeAlgebra(GF(5),1, 'z')
        Free Algebra on 1 generators (z,) over Finite Field of size 5
        sage: FreeAlgebra(GF(5),1, ['alpha'])
        Free Algebra on 1 generators (alpha,) over Finite Field of size 5
        sage: FreeAlgebra(FreeAlgebra(ZZ,1,'a'), 2, 'x')
        Free Algebra on 2 generators (x0, x1) over Free Algebra on 1 generators (a,) over Integer Ring

    Free algebras are globally unique::

        sage: F = FreeAlgebra(ZZ,3,'x,y,z')
        sage: G = FreeAlgebra(ZZ,3,'x,y,z')
        sage: F is G
        True
        sage: F.<x,y,z> = FreeAlgebra(GF(5),3)  # indirect doctest
        sage: F is loads(dumps(F))
        True
        sage: F is FreeAlgebra(GF(5),['x','y','z'])
        True
        sage: copy(F) is F is loads(dumps(F))
        True
        sage: TestSuite(F).run()

    By :trac:`7797`, we provide a different implementation of free
    algebras, based on Singular's "letterplace rings". Our letterplace
    wrapper allows for chosing positive integral degree weights for the
    generators of the free algebra. However, only (weighted) homogenous
    elements are supported. Of course, isomorphic algebras in different
    implementations are not identical::

        sage: G = FreeAlgebra(GF(5),['x','y','z'], implementation='letterplace')
        sage: F == G
        False
        sage: G is FreeAlgebra(GF(5),['x','y','z'], implementation='letterplace')
        True
        sage: copy(G) is G is loads(dumps(G))
        True
        sage: TestSuite(G).run()

    ::

        sage: H = FreeAlgebra(GF(5),['x','y','z'], implementation='letterplace', degrees=[1,2,3])
        sage: F != H != G
        True
        sage: H is FreeAlgebra(GF(5),['x','y','z'], implementation='letterplace', degrees=[1,2,3])
        True
        sage: copy(H) is H is loads(dumps(H))
        True
        sage: TestSuite(H).run()

    Free algebras commute with their base ring.
    ::

        sage: K.<a,b> = FreeAlgebra(QQ,2)
        sage: K.is_commutative()
        False
        sage: L.<c> = FreeAlgebra(K,1)
        sage: L.is_commutative()
        False
        sage: s = a*b^2 * c^3; s
        a*b^2*c^3
        sage: parent(s)
        Free Algebra on 1 generators (c,) over Free Algebra on 2 generators (a, b) over Rational Field
        sage: c^3 * a * b^2
        a*b^2*c^3
    """
    def create_key(self,base_ring, arg1=None, arg2=None,
                                      sparse=False, order='degrevlex',
                                      names=None, name=None,
                                      implementation=None, degrees=None):
        """
        Create the key under which a free algebra is stored.

        TESTS::

            sage: FreeAlgebra.create_key(GF(5),['x','y','z'])
            (Finite Field of size 5, ('x', 'y', 'z'))
            sage: FreeAlgebra.create_key(GF(5),['x','y','z'],3)
            (Finite Field of size 5, ('x', 'y', 'z'))
            sage: FreeAlgebra.create_key(GF(5),3,'xyz')
            (Finite Field of size 5, ('x', 'y', 'z'))
            sage: FreeAlgebra.create_key(GF(5),['x','y','z'], implementation='letterplace')
            (Multivariate Polynomial Ring in x, y, z over Finite Field of size 5,)
            sage: FreeAlgebra.create_key(GF(5),['x','y','z'],3, implementation='letterplace')
            (Multivariate Polynomial Ring in x, y, z over Finite Field of size 5,)
            sage: FreeAlgebra.create_key(GF(5),3,'xyz', implementation='letterplace')
            (Multivariate Polynomial Ring in x, y, z over Finite Field of size 5,)
            sage: FreeAlgebra.create_key(GF(5),3,'xyz', implementation='letterplace', degrees=[1,2,3])
            ((1, 2, 3), Multivariate Polynomial Ring in x, y, z, x_ over Finite Field of size 5)

        """
        if arg1 is None and arg2 is None and names is None:
            # this is used for pickling
            if degrees is None:
                return (base_ring,)
            return tuple(degrees),base_ring
        PolRing = None
        # test if we can use libSingular/letterplace
        if implementation is not None and implementation != 'generic':
            try:
                PolRing = PolynomialRing(base_ring, arg1, arg2,
                                   sparse=sparse, order=order,
                                   names=names, name=name,
                                   implementation=implementation if implementation != 'letterplace' else None)
                if not isinstance(PolRing, MPolynomialRing_libsingular):
                    if PolRing.ngens() == 1:
                        PolRing = PolynomialRing(base_ring, 1, PolRing.variable_names())
                        if not isinstance(PolRing, MPolynomialRing_libsingular):
                            raise TypeError
                    else:
                        raise TypeError
            except (TypeError, NotImplementedError) as msg:
                raise NotImplementedError("The letterplace implementation is not available for the free algebra you requested")
        if PolRing is not None:
            if degrees is None:
                return (PolRing,)
            from sage.all import TermOrder
            T = PolRing.term_order() + TermOrder('lex',1)
            varnames = list(PolRing.variable_names())
            newname = 'x'
            while newname in varnames:
                newname += '_'
            varnames.append(newname)
            return tuple(degrees),PolynomialRing(PolRing.base(), varnames,
                    sparse=sparse, order=T,
                    implementation=implementation if implementation != 'letterplace' else None)
        # normalise the generator names
        from sage.all import Integer
        if isinstance(arg1, (int, long, Integer)):
            arg1, arg2 = arg2, arg1
        if not names is None:
            arg1 = names
        elif not name is None:
            arg1 = name
        if arg2 is None:
            arg2 = len(arg1)
        names = sage.structure.parent_gens.normalize_names(arg2, arg1)
        return base_ring, names

    def create_object(self, version, key):
        """
        Construct the free algebra that belongs to a unique key.

        NOTE:

        Of course, that method should not be called directly,
        since it does not use the cache of free algebras.

        TESTS::

            sage: FreeAlgebra.create_object('4.7.1', (QQ['x','y'],))
            Free Associative Unital Algebra on 2 generators (x, y) over Rational Field
            sage: FreeAlgebra.create_object('4.7.1', (QQ['x','y'],)) is FreeAlgebra(QQ,['x','y'])
            False

        """
        if len(key) == 1:
            from sage.algebras.letterplace.free_algebra_letterplace import FreeAlgebra_letterplace
            return FreeAlgebra_letterplace(key[0])
        if isinstance(key[0], tuple):
            from sage.algebras.letterplace.free_algebra_letterplace import FreeAlgebra_letterplace
            return FreeAlgebra_letterplace(key[1], degrees=key[0])
        return FreeAlgebra_generic(key[0], len(key[1]), key[1])

FreeAlgebra = FreeAlgebraFactory('FreeAlgebra')


def is_FreeAlgebra(x):
    """
    Return True if x is a free algebra; otherwise, return False.

    EXAMPLES::

        sage: from sage.algebras.free_algebra import is_FreeAlgebra
        sage: is_FreeAlgebra(5)
        False
        sage: is_FreeAlgebra(ZZ)
        False
        sage: is_FreeAlgebra(FreeAlgebra(ZZ,100,'x'))
        True
        sage: is_FreeAlgebra(FreeAlgebra(ZZ,10,'x',implementation='letterplace'))
        True
        sage: is_FreeAlgebra(FreeAlgebra(ZZ,10,'x',implementation='letterplace', degrees=range(1,11)))
        True

    """
    from sage.algebras.letterplace.free_algebra_letterplace import FreeAlgebra_letterplace
    return isinstance(x, (FreeAlgebra_generic,FreeAlgebra_letterplace))


class FreeAlgebra_generic(CombinatorialFreeModule, Algebra):
    """
    The free algebra on `n` generators over a base ring.

    INPUT:

    - ``R`` -- a ring
    - ``n`` -- an integer
    - ``names`` -- the generator names

    EXAMPLES::

        sage: F.<x,y,z> = FreeAlgebra(QQ, 3); F
        Free Algebra on 3 generators (x, y, z) over Rational Field
        sage: mul(F.gens())
        x*y*z
        sage: mul([ F.gen(i%3) for i in range(12) ])
        x*y*z*x*y*z*x*y*z*x*y*z
        sage: mul([ F.gen(i%3) for i in range(12) ]) + mul([ F.gen(i%2) for i in range(12) ])
        x*y*x*y*x*y*x*y*x*y*x*y + x*y*z*x*y*z*x*y*z*x*y*z
        sage: (2 + x*z + x^2)^2 + (x - y)^2
        4 + 5*x^2 - x*y + 4*x*z - y*x + y^2 + x^4 + x^3*z + x*z*x^2 + x*z*x*z

    TESTS:

    Free algebras commute with their base ring.
    ::

        sage: K.<a,b> = FreeAlgebra(QQ)
        sage: K.is_commutative()
        False
        sage: L.<c,d> = FreeAlgebra(K)
        sage: L.is_commutative()
        False
        sage: s = a*b^2 * c^3; s
        a*b^2*c^3
        sage: parent(s)
        Free Algebra on 2 generators (c, d) over Free Algebra on 2 generators (a, b) over Rational Field
        sage: c^3 * a * b^2
        a*b^2*c^3

    """
    Element = FreeAlgebraElement
    def __init__(self, R, n, names):
        """
        The free algebra on `n` generators over a base ring.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, 3); F # indirect doctet
            Free Algebra on 3 generators (x, y, z) over Rational Field

        TEST:

        Note that the following is *not* the recommended way to create
        a free algebra::

            sage: from sage.algebras.free_algebra import FreeAlgebra_generic
            sage: FreeAlgebra_generic(ZZ, 3, 'abc')
            Free Algebra on 3 generators (a, b, c) over Integer Ring
        """
        if R not in Rings():
            raise TypeError("Argument R must be a ring.")
        self.__ngens = n
        indices = FreeMonoid(n, names=names)
        cat = AlgebrasWithBasis(R)
        CombinatorialFreeModule.__init__(self, R, indices, prefix='F',
                                         category=cat)
        self._assign_names(indices.variable_names())

    def one_basis(self):
        """
        Return the index of the basis element `1`.

        EXAMPLES::

            sage: F = FreeAlgebra(QQ, 2, 'x,y')
            sage: F.one_basis()
            1
            sage: F.one_basis().parent()
            Free monoid on 2 generators (x, y)
        """
        return self._basis_keys.one()

    def is_field(self, proof=True):
        """
        Return True if this Free Algebra is a field, which is only if the
        base ring is a field and there are no generators

        EXAMPLES::

            sage: A = FreeAlgebra(QQ,0,'')
            sage: A.is_field()
            True
            sage: A = FreeAlgebra(QQ,1,'x')
            sage: A.is_field()
            False
        """
        if self.__ngens == 0:
            return self.base_ring().is_field(proof)
        return False

    def is_commutative(self):
        """
        Return True if this free algebra is commutative.

        EXAMPLES::

            sage: R.<x> = FreeAlgebra(QQ,1)
            sage: R.is_commutative()
            True
            sage: R.<x,y> = FreeAlgebra(QQ,2)
            sage: R.is_commutative()
            False
        """
        return self.__ngens <= 1 and self.base_ring().is_commutative()

    def __cmp__(self, other):
        """
        Two free algebras are considered the same if they have the same
        base ring, number of generators and variable names, and the same
        implementation.

        EXAMPLES::

            sage: F = FreeAlgebra(QQ,3,'x')
            sage: F == FreeAlgebra(QQ,3,'x')
            True
            sage: F is FreeAlgebra(QQ,3,'x')
            True
            sage: F == FreeAlgebra(ZZ,3,'x')
            False
            sage: F == FreeAlgebra(QQ,4,'x')
            False
            sage: F == FreeAlgebra(QQ,3,'y')
            False

        Note that since :trac:`7797` there is a different
        implementation of free algebras. Two corresponding free
        algebras in different implementations are not equal, but there
        is a coercion::


        """
        if not isinstance(other, FreeAlgebra_generic):
            return -1
        c = cmp(self.base_ring(), other.base_ring())
        if c: return c
        c = cmp(self.__ngens, other.ngens())
        if c: return c
        c = cmp(self.variable_names(), other.variable_names())
        if c: return c
        return 0

    def _repr_(self):
        """
        Text representation of this free algebra.

        EXAMPLES::

            sage: F = FreeAlgebra(QQ,3,'x')
            sage: F  # indirect doctest
            Free Algebra on 3 generators (x0, x1, x2) over Rational Field
            sage: F.rename('QQ<<x0,x1,x2>>')
            sage: F #indirect doctest
            QQ<<x0,x1,x2>>
            sage: FreeAlgebra(ZZ,1,['a'])
            Free Algebra on 1 generators (a,) over Integer Ring
        """
        return "Free Algebra on {} generators {} over {}".format(
            self.__ngens, self.gens(), self.base_ring())

    def _element_constructor_(self, x):
        """
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: R.<x,y> = FreeAlgebra(QQ,2)
            sage: R(3) # indirect doctest
            3

        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(GF(5),3)
            sage: L.<x,y,z> = FreeAlgebra(ZZ,3,implementation='letterplace')
            sage: F(x)     # indirect doctest
            x
            sage: F.1*L.2
            y*z
            sage: (F.1*L.2).parent() is F
            True

       ::

            sage: K.<z> = GF(25)
            sage: F.<a,b,c> = FreeAlgebra(K,3)
            sage: L.<a,b,c> = FreeAlgebra(K,3, implementation='letterplace')
            sage: F.1+(z+1)*L.2
            b + (z+1)*c

        Check that :trac:`15169` is fixed::

            sage: A.<x> = FreeAlgebra(CC)
            sage: A(2)
            2.00000000000000

        We check that the string coercions work correctly over
        inexact fields::

            sage: F.<x,y> = FreeAlgebra(CC)
            sage: F('2')
            2.00000000000000
            sage: F('x')
            1.00000000000000*x

        Check that it also converts factorizations::

            sage: f = Factorization([(x,2),(y,3)]); f
            1.00000000000000*x^2 * 1.00000000000000*y^3
            sage: F(f)
            1.00000000000000*x^2*y^3
        """
        if isinstance(x, FreeAlgebraElement):
            P = x.parent()
            if P is self:
                return x
            if P is not self.base_ring():
                return self.element_class(self, x)
        elif hasattr(x,'letterplace_polynomial'):
            P = x.parent()
            if self.has_coerce_map_from(P): # letterplace versus generic
                ngens = P.ngens()
                M = self._basis_keys
                def exp_to_monomial(T):
                    out = []
                    for i in xrange(len(T)):
                        if T[i]:
                            out.append((i%ngens,T[i]))
                    return M(out)
                return self.element_class(self, dict([(exp_to_monomial(T),c) for T,c in x.letterplace_polynomial().dict().iteritems()]))
        # ok, not a free algebra element (or should not be viewed as one).
        if isinstance(x, basestring):
            from sage.all import sage_eval
            G = self.gens()
            d = {str(v): G[i] for i,v in enumerate(self.variable_names())}
            return self(sage_eval(x, locals=d))
        R = self.base_ring()
        # coercion from free monoid
        if isinstance(x, FreeMonoidElement) and x.parent() is self._basis_keys:
            return self.element_class(self, {x: R.one()})
        # coercion from the PBW basis
        if isinstance(x, PBWBasisOfFreeAlgebra.Element) \
                and self.has_coerce_map_from(x.parent()._alg):
            return self(x.parent().expansion(x))

        # Check if it's a factorization
        from sage.structure.factorization import Factorization
        if isinstance(x, Factorization):
            return self.prod(f**i for f,i in x)

        # coercion via base ring
        x = R(x)
        if x == 0:
            return self.element_class(self, {})
        return self.element_class(self, {self.one_basis(): x})

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coercion from ``R`` into ``self`` and
        ``False`` otherwise.  The things that coerce into ``self`` are:

        - This free algebra.

        - Anything with a coercion into ``self.monoid()``.

        - Free algebras in the same variables over a base with a coercion
          map into ``self.base_ring()``.

        - The underlying monoid.

        - The PBW basis of ``self``.

        - Anything with a coercion into ``self.base_ring()``.

        TESTS::

            sage: F = FreeAlgebra(ZZ, 3, 'x,y,z')
            sage: G = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: H = FreeAlgebra(ZZ, 1, 'y')
            sage: F._coerce_map_from_(G)
            False
            sage: G._coerce_map_from_(F)
            True
            sage: F._coerce_map_from_(H)
            False
            sage: F._coerce_map_from_(QQ)
            False
            sage: G._coerce_map_from_(QQ)
            True
            sage: F._coerce_map_from_(G.monoid())
            True
            sage: F._coerce_map_from_(F.pbw_basis())
            True
            sage: F.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False

            sage: K.<z> = GF(25)
            sage: F.<a,b,c> = FreeAlgebra(K,3)
            sage: F._coerce_map_from_(ZZ)
            True
            sage: F._coerce_map_from_(QQ)
            False
            sage: F._coerce_map_from_(F.monoid())
            True
            sage: F._coerce_map_from_(F.pbw_basis())
            True
            sage: G = FreeAlgebra(ZZ, 3, 'a,b,c')
            sage: F._coerce_map_from_(G)
            True
            sage: G._coerce_map_from_(F)
            False
            sage: L.<a,b,c> = FreeAlgebra(K,3, implementation='letterplace')
            sage: F.1 + (z+1) * L.2
            b + (z+1)*c
        """
        if self._basis_keys.has_coerce_map_from(R):
            return True

        # free algebras in the same variable over any base that coerces in:
        if is_FreeAlgebra(R):
            if R.variable_names() == self.variable_names():
                return self.base_ring().has_coerce_map_from(R.base_ring())
        if isinstance(R, PBWBasisOfFreeAlgebra):
            return self.has_coerce_map_from(R._alg)

        return self.base_ring().has_coerce_map_from(R)

    def gen(self, i):
        """
        The ``i``-th generator of the algebra.

        EXAMPLES::

            sage: F = FreeAlgebra(ZZ,3,'x,y,z')
            sage: F.gen(0)
            x
        """
        if i < 0 or not i < self.__ngens:
            raise IndexError("Argument i (= {}) must be between 0 and {}.".format(i, self.__ngens-1))
        R = self.base_ring()
        F = self._basis_keys
        return self.element_class(self, {F.gen(i): R.one()})

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: F = FreeAlgebra(ZZ,3,'x,y,z')
            sage: F.algebra_generators()
            Finite family {'y': y, 'x': x, 'z': z}
        """
        ret = {}
        for i in range(self.__ngens):
            x = self.gen(i)
            ret[str(x)] = x
        from sage.sets.family import Family
        return Family(ret)

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: F = FreeAlgebra(ZZ,3,'x,y,z')
            sage: F.gens()
            (x, y, z)
        """
        return tuple(self.gen(i) for i in range(self.__ngens))

    def product_on_basis(self, x, y):
        """
        Return the product of the basis elements indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: F = FreeAlgebra(ZZ,3,'x,y,z')
            sage: I = F.basis().keys()
            sage: x,y,z = I.gens()
            sage: F.product_on_basis(x*y, z*y)
            x*y*z*y
        """
        return self.monomial(x * y)

    def _ideal_class_(self, n=0):
        r"""
        Return the class for the ideals of ``self``.

        EXAMPLES::

            sage: F = FreeAlgebra(ZZ, 3, 'x,y,z')
            sage: F._ideal_class_()
            <class 'sage.algebras.free_algebra.FreeAlgebraIdeal'>
        """
        return FreeAlgebraIdeal

    def quotient(self, I, mats=None, names=None, category=None):
        """
        Return a quotient algebra.

        The quotient algebra can be defined in two ways:

        - Via an ideal `I`.

        - Via the action of a free algebra `A` on a (finitely generated) free
          module. The input for the quotient algebra is a list of monomials
          (in the underlying monoid for `A`) which form a free basis for the
          module of `A`, and a list of matrices, which give the action of the
          free generators of `A` on this monomial basis.

        EXAMPLES:

        Here is the quaternion algebra defined in terms of three generators::

            sage: n = 3
            sage: A = FreeAlgebra(QQ,n,'i')
            sage: F = A.monoid()
            sage: i, j, k = F.gens()
            sage: mons = [ F(1), i, j, k ]
            sage: M = MatrixSpace(QQ,4)
            sage: mats = [M([0,1,0,0, -1,0,0,0, 0,0,0,-1, 0,0,1,0]),  M([0,0,1,0, 0,0,0,1, -1,0,0,0, 0,-1,0,0]),  M([0,0,0,1, 0,0,-1,0, 0,1,0,0, -1,0,0,0]) ]
            sage: H.<i,j,k> = A.quotient(mons, mats); H
            Free algebra quotient on 3 generators ('i', 'j', 'k') and dimension 4 over Rational Field
        """
        if mats is not None:
            import free_algebra_quotient
            return free_algebra_quotient.FreeAlgebraQuotient(self, I, mats, names)
        #return super(FreeAlgebra_generic, self).quotient(mons, names)
        from sage.algebras.finitely_presented_algebra import FinitelyPresentedAlgebra, TwoSidedAlgebraIdeal
        if not isinstance(I, TwoSidedAlgebraIdeal):
            raise TypeError("must be a two sided algebra ideal")
        return FinitelyPresentedAlgebra(self, I, names)

    quo = quotient

    def ngens(self):
        """
        The number of generators of the algebra.

        EXAMPLES::

            sage: F = FreeAlgebra(ZZ,3,'x,y,z')
            sage: F.ngens()
            3
        """
        return self.__ngens

    def monoid(self):
        """
        The free monoid of generators of the algebra.

        EXAMPLES::

            sage: F = FreeAlgebra(ZZ,3,'x,y,z')
            sage: F.monoid()
            Free monoid on 3 generators (x, y, z)
        """
        return self._basis_keys

    def g_algebra(self, relations, names=None, order='degrevlex', check=True):
        """
        The `G`-Algebra derived from this algebra by relations.
        By default is assumed, that two variables commute.

        .. TODO::

            - Coercion doesn't work yet, there is some cheating about assumptions
            - The optional argument ``check`` controls checking the degeneracy
              conditions. Furthermore, the default values interfere with
              non-degeneracy conditions.

        EXAMPLES::

            sage: A.<x,y,z> = FreeAlgebra(QQ,3)
            sage: G = A.g_algebra({y*x: -x*y})
            sage: (x,y,z) = G.gens()
            sage: x*y
            x*y
            sage: y*x
            -x*y
            sage: z*x
            x*z
            sage: (x,y,z) = A.gens()
            sage: G = A.g_algebra({y*x: -x*y+1})
            sage: (x,y,z) = G.gens()
            sage: y*x
            -x*y + 1
            sage: (x,y,z) = A.gens()
            sage: G = A.g_algebra({y*x: -x*y+z})
            sage: (x,y,z) = G.gens()
            sage: y*x
            -x*y + z
        """
        from sage.matrix.constructor import Matrix

        base_ring = self.base_ring()
        n = self.__ngens
        cmat = Matrix(base_ring, n)
        dmat = Matrix(self, n)
        for i in xrange(n):
            for j in xrange(i+1,n):
                cmat[i,j] = 1
        for (to_commute,commuted) in relations.iteritems():
            #This is dirty, coercion is broken
            assert isinstance(to_commute, FreeAlgebraElement), to_commute.__class__
            assert isinstance(commuted, FreeAlgebraElement), commuted
            ((v1,e1),(v2,e2)) = list(list(to_commute)[0][0])
            assert e1 == 1
            assert e2 == 1
            assert v1 > v2
            c_coef = None
            d_poly = None
            for (m,c) in commuted:
                if list(m) == [(v2,1),(v1,1)]:
                    c_coef = c
                    #buggy coercion workaround
                    d_poly = commuted - self(c) * self(m)
                    break
            assert not c_coef is None,list(m)
            v2_ind = self.gens().index(v2)
            v1_ind = self.gens().index(v1)
            cmat[v2_ind,v1_ind] = c_coef
            if d_poly:
                dmat[v2_ind,v1_ind] = d_poly
        from sage.rings.polynomial.plural import g_Algebra
        return g_Algebra(base_ring, cmat, dmat, names = names or self.variable_names(),
                         order=order, check=check)

    def poincare_birkhoff_witt_basis(self):
        """
        Return the Poincare-Birkhoff-Witt (PBW) basis of ``self``.

        EXAMPLES::

            sage: F.<x,y> = FreeAlgebra(QQ, 2)
            sage: F.poincare_birkhoff_witt_basis()
            The Poincare-Birkhoff-Witt basis of Free Algebra on 2 generators (x, y) over Rational Field
        """
        return PBWBasisOfFreeAlgebra(self)

    pbw_basis = poincare_birkhoff_witt_basis

    def pbw_element(self, elt):
        """
        Return the element ``elt`` in the Poincare-Birkhoff-Witt basis.

        EXAMPLES::

            sage: F.<x,y> = FreeAlgebra(QQ, 2)
            sage: F.pbw_element(x*y - y*x + 2)
            2*PBW[1] + PBW[x*y]
            sage: F.pbw_element(F.one())
            PBW[1]
            sage: F.pbw_element(x*y*x + x^3*y)
            PBW[x*y]*PBW[x] + PBW[y]*PBW[x]^2 + PBW[x^3*y] + PBW[x^2*y]*PBW[x]
             + PBW[x*y]*PBW[x]^2 + PBW[y]*PBW[x]^3
        """
        PBW = self.pbw_basis()
        if elt == self.zero():
            return PBW.zero()

        l = {}
        while elt: # != 0
            lst = list(elt)
            min_elt, coeff = lst[0]
            min_word = min_elt.to_word()
            for item in lst[1:-1]:
                word = item[0].to_word()
                if min_word.lex_less(word):
                    min_elt, coeff = item
                    min_word = word
            l[min_elt] = l.get(min_elt, 0) + coeff
            elt = elt - coeff * self.lie_polynomial(min_elt)
        return PBW.sum_of_terms([(k, v) for k,v in l.items() if v != 0], distinct=True)

    @cached_method
    def lie_polynomial(self, w):
        """
        Return the Lie polynomial associated to the Lyndon word ``w``. If
        ``w`` is not Lyndon, then return the product of Lie polynomials of the
        Lyndon factorization of ``w``.

        INPUT:

        - ``w`` -- a word or an element of the free monoid

        EXAMPLES::

            sage: F = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: M.<x,y,z> = FreeMonoid(3)
            sage: F.lie_polynomial(x*y)
            x*y - y*x
            sage: F.lie_polynomial(y*x)
            y*x
            sage: F.lie_polynomial(x^2*y*x)
            x^2*y*x - x*y*x^2
            sage: F.lie_polynomial(y*z*x*z*x*z)
            y*z*x*z*x*z - y*z*x*z^2*x - y*z^2*x^2*z + y*z^2*x*z*x
             - z*y*x*z*x*z + z*y*x*z^2*x + z*y*z*x^2*z - z*y*z*x*z*x

        TESTS:

        We test some corner cases and alternative inputs::

            sage: F.lie_polynomial(Word('xy'))
            x*y - y*x
            sage: F.lie_polynomial('xy')
            x*y - y*x
            sage: F.lie_polynomial(M.one())
            1
            sage: F.lie_polynomial(Word([]))
            1
            sage: F.lie_polynomial('')
            1
        """
        if not w:
            return self.one()
        M = self._basis_keys

        if len(w) == 1:
            return self(M(w))

        ret = self.one()
        # We have to be careful about order here.
        # Since the Lyndon factors appear from left to right
        #   we must multiply from left to right as well.
        for factor in Word(w).lyndon_factorization():
            if len(factor) == 1:
                ret = ret * self(M(factor))
                continue
            x,y = factor.standard_factorization()
            x = M(x)
            y = M(y)
            ret = ret * (self(x * y) - self(y * x))
        return ret

from sage.misc.cache import Cache
cache = Cache(FreeAlgebra_generic)

###################
## Ideals

class FreeAlgebraIdeal(Ideal_nc):
    """
    An ideal of a free algebra.
    """
    def __init__(self, F, gens, coerce=True, side="twosided"):
        r"""
        Initialize ``self``.

        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: I = F.ideal(x^2 - y, x*y - y*z)

        Skip the category test because ideals aren't currently
        proper category objects::

            sage: TestSuite(I).run(skip="_test_category")
        """
        Ideal_nc.__init__(self, F, gens, coerce, side)

        gens = self.gens()
        if F.zero() in gens:
            self._gb = (F.zero(),)
            self._gb_todo = []
            return
        if F.one() in gens:
            self._gb = tuple(F.gens())
            self._gb_todo = []
            return

        gb = map(lambda x: x / x.leading_coefficient(FreeAlgebraIdeal._lead_cmp), gens)
        self._gb = gb
        self._gb_todo = [(g, h) for g in gb for h in gb]

    def free_algebra(self):
        r"""
        Return the ambient free algebra of ``self``.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: I = F.ideal(x^2 - y, x*y - y*z)
            sage: I.free_algebra() is F
            True
        """
        return self.ring()

    @staticmethod
    def _lead_cmp(x, y):
        r"""
        Compare ``x`` and ``y`` by degree then reverse lex so
        :meth:`leading_item()` returns the largest degree and smallest lex.
        """
        x = x.to_word()
        y = y.to_word()
        c = cmp(len(x), len(y))
        if c != 0:
            return c
        return cmp(y, x)

    def groebner_basis(self, max_steps=infinity):
        """
        Return a Groebner basis of ``self``.

        INPUT:

        - ``max_steps`` -- (default: infinity) the maximum number of steps to
          do before terminating

        .. WARNING::

            This will run forever if the Groebner basis is infinite and
            ``max_steps`` is not specified.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: I = F.ideal(x^2 - x*y)
            sage: I.groebner_basis(5)
            WARNING: returning an incomplete Groebner basis
            (x^2 - x*y,
             x*y*x - x*y^2,
             x*y^2*x - x*y^3,
             x*y^3*x - x*y^4,
             x*y^4*x - x*y^5)
            sage: I = F.ideal(z^2 - z*y)
            sage: I.groebner_basis()
            (z*y - z^2,)
            sage: I = F.ideal(x^2 - x*y, x*y*x - y*x*y)
            sage: I.groebner_basis()
            (x^2 - x*y, x*y*x - y*x*y, x*y^2 - y*x*y)
            sage: I = F.ideal(z^2 - z*y, z*y*z - y*z*y)
            sage: I.groebner_basis()
            (z*y - z^2, y*z*y - z*y*z, y*z^4 - z^5)
        """
        side = self.side()
        if side != "twosided":
            raise NotImplementedError

        # Setup variables and functions
        l_cmp = FreeAlgebraIdeal._lead_cmp
        F = self.ring()

        n = len(F.gens())
        FM = F._basis_keys

        # Run Groebner basis algorithm

        zero = F.zero()
        from sage.combinat.words.word import Word
        while len(self._gb_todo) > 0 and len(self._gb) < max_steps:
            p = self._gb_todo.pop(0)

            # Compute the essential common multiples of the leading terms
            f = map(lambda x: list(x.leading_support(l_cmp).to_word()), p)
            ell = min(len(f[0]), len(f[1]))
            for k in range(1, ell):
                if f[0][:k] == f[1][-k:]:
                    w0 = Word(f[1][:-k])
                    w1 = Word(f[0][k:])

                    # Compute the S-polynomial
                    h = F(FM(w0)) * p[0] - p[1] * F(FM(w1))
                    h = self._normal_form(h, self._gb)

                    if h != zero:
                        h = h / h.leading_coefficient(l_cmp)
                        self._gb_todo.extend([(g, h) for g in self._gb])
                        self._gb.append(h)

        if len(self._gb_todo) != 0:
            print("WARNING: returning an incomplete Groebner basis")

        return tuple(self._gb)

    def _normal_form(self, x, G):
        """
        Return ``x`` modulo ``G`` (i.e. the normal form with respect
        to ``G``).

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: I = F.ideal(x^2 - x*y)
            sage: GB = I.groebner_basis(5)
            WARNING: returning an incomplete Groebner basis
            sage: I._normal_form(x*y^4*x + x^2, GB)
            x*y + x*y^5

        The following reduces to ``x*y^6`` with a larger Groebner
        basis. However, since we've truncated it, we don't have any
        reduction steps::

            sage: I._normal_form(x*y^5*x, GB)
            x*y^5*x
        """
        F = self.ring()
        ret = F.zero()
        FM = F._basis_keys
        l_cmp = FreeAlgebraIdeal._lead_cmp
        side = self.side()

        while x != 0:
            u, la = x.leading_item(l_cmp)
            found = False
            for g in G:
                LM, mu = g.leading_item(l_cmp)

                # Factor u by LM if possible
                w = u.to_word()
                lmw = LM.to_word()
                if side == 'twosided':
                    f = lmw.first_pos_in(w)
                    if f is not None:
                        found = True
                        x -= la / mu * F(FM(w[:f])) * g * F(FM(w[f+len(LM):]))
                        break
                elif side == 'left':
                    if lmw.is_proper_suffix(w):
                        x -= la / mu * F(FM(w[len(lmw):])) * g
                        break
                elif side == 'right':
                    if lmw.is_proper_prefix(w):
                        x -= la / mu * g * F(FM(w[:-len(lmw)]))
                        break
            if not found:
                ret += la * F(u)
                x -= la * F(u)
        return ret

    def reduce(self, x, max_steps=infinity):
        """
        Return ``x`` modulo ``self``.

        INPUT:

        - ``max_steps`` -- (default: infinity) the maximum number of steps to
          compute in the Groebner basis before terminating

        .. WARNING::

            If this terminates after ``max_steps``, the resulting reduction
            may not fully reduce the element.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: I = F.ideal(x^2 - x*y, x*y*x - y*x*y)
            sage: I.reduce(x^2 + x*y + x*y*x + z^2)
            2*x*y + z^2 + y*x*y

        The following is with an infinite Groebner basis. The element
        `x y^5 x` reduces to `x y^6` with a large enough Groebner basis,
        but there are no reduction steps when our Groebner basis is too small.

            sage: I = F.ideal(x^2 - x*y)
            sage: I.reduce(x*y^5*x, 5)
            WARNING: returning an incomplete Groebner basis
            x*y^5*x
            sage: I.reduce(x*y^5*x, 6)
            WARNING: returning an incomplete Groebner basis
            x*y^6
        """
        if x == self.ring().zero():
            return x
        G = self.groebner_basis(max_steps)
        return self._normal_form(x, G)

