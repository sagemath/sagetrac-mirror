# -*- coding: utf-8 -*-
r"""
Stuffle algebras

AUTHORS:

- Matthieu Deneufchâtel (2013-07)
"""

#*****************************************************************************
#  Copyright (C) 2013 Matthieu Deneufchâtel <matthieudeneufch@yahoo.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.rings import Rings
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.coalgebras_with_basis import CoalgebrasWithBasis
from sage.categories.tensor import tensor
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.alphabet import build_alphabet
from sage.combinat.words.words import Words
from sage.combinat.words.word import Word
from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.rings.infinity import infinity

def stuffle(w1, w2, alphabet=None):
    r"""
    Return a list containing the words obtained by the quasi
    shuffle-product of words ``w1`` and ``w2``.

    If `Y = \left\{ y_1, \ldots, y_n\}`, then for all `y_i, y_j \in Y` and for
    all `u, v \in Y^*` the stuffle product `\ast` is defined recursively as

    .. MATH::

        \left\lbrace
        \begin{aligned}
        u \ast 1 & = 1 \ast u = u, \\
        y_i u \ast y_j v & = y_i ( u \ast y_j v )
        + y_j ( y_i u \ast v ) + y_{i+j} ( u \ast v).
        \end{aligned}
        \right.

    EXAMPLES::

        sage: from sage.algebras.stuffle_algebra import stuffle
        sage: A = build_alphabet(10, 'y')
        sage: W = Words(A)
        sage: w1 = W([W.alphabet()[0], W.alphabet()[0]])
        sage: w2 = W([W.alphabet()[1]])
        sage: stuffle(w1, w2)
        [word: y1,y1,y2, word: y1,y2,y1, word: y1,y3, word: y2,y1,y1, word: y3,y1]

        sage: A = build_alphabet(10,'g')
        sage: W = Words(A)
        sage: w1 = W([W.alphabet()[0],W.alphabet()[0]])
        sage: w2 = W([W.alphabet()[1]])
        sage: stuffle(w1, w2, A)
        [word: g1,g1,g2, word: g1,g2,g1, word: g1,g3, word: g2,g1,g1, word: g3,g1]
    """
    if alphabet is None:
        W = w1.parent()
    else:
        W = Words(alphabet)

    letters = W.alphabet()
    wrd = W._construct_word
    if len(w1) == 0:
        return [w2]
    if len(w2) == 0:
        return [w1]

    t1 = wrd([w1[0]])
    t2 = wrd([w2[0]])
    r1 = w1[1:]
    r2 = w2[1:]
    l1 = [wrd(t1 * i) for i in stuffle(r1, w2, W)]
    l2 = [wrd(t2 * i) for i in stuffle(w1, r2, W)]
    l = letters[letters.index(w1[0]) + letters.index(w2[0])+1]
    t3 = wrd([l])
    l3 = [wrd(t3 * i) for i in stuffle(r1, r2, W)]
    return l1 + l2 + l3

class StuffleAlgebra(CombinatorialFreeModule):
    r"""
    The stuffle algebra over a base ring.

    Stuffle algebras are commutative and associative algebras, with a
    basis indexed by words. The product of two words `w_1 \ast w_2`
    is given by a deformation of the shuffle product of `w_1` and
    `w_2` called stuffle product (see
    :func:`~sage.algebras.stuffle_algebra.stuffle()`).

    INPUT:

    - ``R`` -- a ring

    - ``names`` -- the name of the letters

    - ``n`` -- (optional) the number of variables with name ``names``

    EXAMPLES::

        sage: from sage.algebras.stuffle_algebra import StuffleAlgebra
        sage: A = StuffleAlgebra(QQ, 'y', 40); A
        Stuffle Algebra over Rational Field with variables y
        sage: mul(A.gens()[0:3])
        St[word: y1,y2,y3] + St[word: y1,y3,y2] + St[word: y1,y5] + St[word: y2,y1,y3]
         + St[word: y2,y3,y1] + St[word: y2,y4] + St[word: y3,y1,y2] + St[word: y3,y2,y1]
         + 2*St[word: y3,y3] + St[word: y4,y2] + St[word: y5,y1] + St[word: y6]

        sage: mul([A.gen(i) for i in range(2)])
        St[word: y1,y2] + St[word: y2,y1] + St[word: y3]

        sage: S = StuffleAlgebra(ZZ, 'x,y'); S
        Stuffle Algebra over Integer Ring with variables y
        sage: S.base_ring()
        Integer Ring

        sage: G = StuffleAlgebra(S, 'g'); G
        Stuffle Algebra over Stuffle Algebra over Integer Ring with variables y with variables g
        sage: G.base_ring()
        Stuffle Algebra over Integer Ring with variables y

    Stuffle algebras commute with their base ring::

        sage: K = StuffleAlgebra(QQ, 'a,b')
        sage: K.is_commutative()
        True
        sage: a,b = K.gens()
        sage: L = StuffleAlgebra(K, 'c,d')
        sage: L.is_commutative()
        True
        sage: c,d = L.gens()
        sage: s = a*b * c^2; s
        (2*St[word:a,b]+2*St[word:b,a]+2*St[word:y3])*St[word: g1,g1]
         + (St[word:a,b]+St[word:b,a]+St[word:y3])*St[word: g2]
        sage: parent(s)
        Stuffle Algebra over Stuffle Algebra over Rational Field with variables y with variables g
        sage: c^2 * a * b
        (2*St[word:y1,y2]+2*St[word:y2,y1]+2*St[word:y3])*St[word: g1,g1]
         + (St[word:y1,y2]+St[word:y2,y1]+St[word:y3])*St[word: g2]

    Stuffle algebras are commutative::

        sage: c^3 * b * a * b == c * a * c * b^2 * c
        True

    REFERENCES:

    .. [Hoffman2000] Michael Hoffman. *Quasi-shuffle products*.
       J. Algebraic Combin. **11** (2000). No. 1, 49-68.
       :arxiv:`math/9907173v1`.
    """
    @staticmethod
    def __classcall_private__(cls, R, names, n=None):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: St1 = StuffleAlgebra(QQ, ['x','y','z'])
            sage: St2 = StuffleAlgebra(QQ, ('x','y','z'))
            sage: St3 = StuffleAlgebra(QQ, 'x,y,z')
            sage: St1 is St2 and St2 is St3
            True
        """
        if n is None:
            if isinstance(names, str):
                names = names.split(',')
            names = tuple(names)
        return super(cls, StuffleAlgebra).__classcall__(cls, R, names, n)

    def __init__(self, R, names, n):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: F = StuffleAlgebra(QQ)
            sage: TestSuite(F).run()

        TESTS::

            sage: StuffleAlgebra(24)
            Traceback (most recent call last):
            ...
            TypeError: argument R must be a ring
        """
        if R not in Rings():
            raise TypeError("argument R must be a ring")

        self._n = n
        self._names = names
        if n is not None:
            alph = build_alphabet(n, names)
        else:
            alph = build_alphabet(names)
        self._alphabet = alph
        CombinatorialFreeModule.__init__(self, R, Words(self._alphabet),
            prefix="St", latex_prefix="",
            category=(AlgebrasWithBasis(R), CommutativeAlgebras(R),
                      CoalgebrasWithBasis(R)))

    def variable_names(self):
        r"""
        Return the names of the variables.

        EXAMPLES::

            sage: R = StuffleAlgebra(QQ)
            sage: R.variable_names()
            Lazy
        """
        return self._alphabet

    def is_commutative(self):
        r"""
        Return ``True`` as the stuffle algebra is commutative.

        EXAMPLES::

            sage: R = StuffleAlgebra(QQ)
            sage: R.is_commutative()
            True
            sage: R = StuffleAlgebra(QQ)
            sage: R.is_commutative()
            True
        """
        return True

    def _repr_(self):
        r"""
        Text representation of this stuffle algebra.

        EXAMPLES::

            sage: StuffleAlgebra(QQ)
            Stuffle Algebra with variables y over Rational Field

            sage: StuffleAlgebra(ZZ)
            Stuffle Algebra with variables y over Integer Ring
        """
        if self._n is None or self._n <= 10:
            names = self._alphabet
        elif self._n == infinity:
            names = "{na}0, {na}1, {na}2, ...".format(na=self._names)
        else:
            names = "{na}0, ..., {na}{n}".format(na=self._names, n=self._n-1)
        return "Stuffle Algebra with variables {} over {}".format(names, self.base_ring())

    @cached_method
    def one_basis(self):
        r"""
        Return the empty word, which index of `1` of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = StuffleAlgebra(QQ)
            sage: A.one_basis()
            word:
            sage: A.one()
            St[word: ]
        """
        return self.basis().keys()([])

    def product_on_basis(self, w1, w2):
        r"""
        Return the product of basis elements ``w1`` and ``w2``, as per
        :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis()`.

        INPUT:

        - ``w1``, ``w2`` -- basis elements

        EXAMPLES::

            sage: A = StuffleAlgebra(QQ)
            sage: W = A.basis().keys()
            sage: A.product_on_basis(W([W.alphabet()[0]]),W([W.alphabet()[1]]))
            St[word: y1,y2] + St[word: y2,y1] + St[word: y3]
        """
        return self.sum_of_monomials(u for u in stuffle(w1, w2))

    def gen(self, i):
        r"""
        The ``i``-th generator of the algebra.

        INPUT:

        - ``i`` -- an integer

        EXAMPLES::

            sage: A=StuffleAlgebra(ZZ)
            sage: A.gen(3)
            St[word: y3]
            sage: A.gen(43)
            St[word: y43]
        """
        if self._alphabet.cardinality() < infinity:
            return list(self.algebra_generators())[i]
        return self.algebra_generators()[i]

    def counit(self, x):
        """
        Return the counit of ``x``.

        INPUT:

        - ``x`` -- an element of ``self``

        EXAMPLES::

            sage: F = StuffleAlgebra(QQ)
            sage: f=F.an_element()
            sage: F = StuffleAlgebra(QQ)
            sage: f=F.an_element(); f
            St[word: ] + 2*St[word: y1] + 3*St[word: y2]
            sage: f.counit()
            1
            sage: g=F.basis()[Word(['y1','y2'])]
            sage: g.counit()
            0
        """
        dic = x.monomial_coefficients()
        if Word() not in dic:
            return 0
        return dic[Word()]

    def coproduct_on_basis(self, w):
        """
        Return the coproduct of the element of the basis indexed
        by the word ``w``.

        INPUT:

        - ``w`` -- a word

        EXAMPLES::

            sage: F = StuffleAlgebra(QQ)
            sage: F.coproduct_on_basis(Word(['y3']))
            St[word: ] # St[word: y3] + St[word: y1] # St[word: y2] + St[word: y2] # St[word: y1] + St[word: y3] # St[word: ]
            sage: F.coproduct_on_basis(Word(['y1','y2']))
            St[word: ] # St[word: y1,y2] + St[word: y1] # St[word: y1,y1] + St[word: y1] # St[word: y2] + St[word: y1,y1] # St[word: y1] + St[word: y1,y2] # St[word: ] + St[word: y2] # St[word: y1]
        """
        if len(w) == 0:
            return self.tensor_square().one()

        if len(w) == 1:
            S = self.tensor_square()
            A = self._alphabet
            W = self.basis().keys()._construct_word
            i = self._alphabet.index(w[0]) + 1
            l = [(W([A[j-1]]), W([A[i-j-1]])) for j in range(1, i)] + [(w, W([])), (W([]), w)]
            return S.sum_of_monomials(l)

        B = self.basis()
        result = self.coproduct_on_basis(Word([w[0]]))
        for i in w[1:]:
            temp1 = self.coproduct_on_basis(Word([i]))
            temp2 = 0
            for ((u1, u2), coeff1) in list(temp1):
                for ((v1, v2), coeff2) in list(result):
                    temp2 += coeff1 * coeff2 * tensor((St[Word(v1)*Word(u1)], St[Word(v2)*Word(u2)]))
            result = temp2
        return result

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators of this algebra.

        EXAMPLES::

            sage: A=StuffleAlgebra(QQ); A
            Stuffle Algebra over Rational Field with variables y
            sage: A.algebra_generators()
        """
        W = self.basis().keys()._construct_word
        if self._n == infinity:
            from sage.sets.non_negative_integers import NonNegativeIntegers
            return Family(NonNegativeIntegers(),
                          lambda i: self.monomial(W([self._alphabet[i]])),
                          name='generator')
        return Family(self._alphabet,
                      lambda a: self.monomial(W([a])),
                      name='generator')

    gens = algebra_generators

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: R = StuffleAlgebra(QQ)
            sage: x, y = R.gens()[0:2]
            sage: x,y
            (St[word: y1], St[word: y2])
            sage: R(3)
            3*St[word: ]
            sage: R(x)
            St[word: y1]
        """
        P = x.parent()
        if isinstance(P, StuffleAlgebra):
            if P is self:
                return x
            if P is not self.base_ring():
                return self.element_class(self, x.monomial_coefficients())

        # ok, not a stuffle algebra element (or should not be viewed as one).
        if isinstance(x, basestring):
            from sage.misc.sage_eval import sage_eval
            return sage_eval(x, locals=self.gens_dict())

        R = self.base_ring()
        # coercion via base ring
        x = R(x)
        if x == 0:
            return self.zero()
        return self.from_base_ring_from_one_basis(x)

    def _coerce_impl(self, x):
        r"""
        Canonical coercion of ``x`` into ``self``.

        Here is what canonically coerces to ``self``:

        - this stuffle algebra,

        - anything that coerces to the base ring of this stuffle algebra,

        - any stuffle algebra on the same variables, whose base ring
          coerces to the base ring of this stuffle algebra.

        EXAMPLES::

            sage: F = StuffleAlgebra(GF(7)); F
            Stuffle Algebra over Finite Field of size 7 with variables y

        Elements of the stuffle algebra canonically coerce in::

            sage: x, y, z = F.gens()[0:3]
            sage: F.coerce(x*y) # indirect doctest
            St[word: y1,y2] + St[word: y2,y1] + St[word: y3]

        Elements of the integers coerce in, since there is a coerce map
        from `\ZZ` to GF(7)::

            sage: F.coerce(1)       # indirect doctest
            St[word: ]

        There is no coerce map from `\QQ` to `\GF{7}`::

            sage: F.coerce(2/3)  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Stuffle Algebra over Finite Field of size 7 with variables y

        Elements of the base ring coerce in::

            sage: F.coerce(GF(7)(5))
            5*St[word: ]

        The stuffle algebra over `\ZZ` coerces in, since
        `\ZZ` coerces to `\GF{7}`::

            sage: G = StuffleAlgebra(ZZ)
            sage: Gx,Gy,Gz = G.gens()[0:3]
            sage: z = F.coerce(Gx**2 * Gy);z
            2*St[word: y1,y1,y2] + 2*St[word: y1,y2,y1] + 2*St[word: y1,y3] + 2*St[word: y2,y1,y1] + 2*St[word: y2,y2] + 2*St[word: y3,y1] + St[word: y4]
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the stuffle
        algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

            sage: G.coerce(x^3*y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Stuffle Algebra over Finite Field of size 7 with variables y to Stuffle Algebra over Integer Ring with variables y
        """
        try:
            R = x.parent()

            # stuffle algebras in the same variables over any base
            # that coerces in:
            if isinstance(R, StuffleAlgebra) and R.variable_names() == self.variable_names():
                if self.has_coerce_map_from(R.base_ring()):
                    return self(x)

                raise TypeError("no natural map between bases of stuffle algebras")

        except AttributeError:
            pass

        # any ring that coerces to the base ring of this stuffle algebra.
        return self._coerce_try(x, [self.base_ring()])

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are:

        - Stuffle Algebras in the same variables over a base with a coercion
          map into ``self.base_ring()``.

        - Anything with a coercion into ``self.base_ring()``.

        TESTS::

            sage: F = StuffleAlgebra(ZZ)
            sage: G = StuffleAlgebra(QQ)
            sage: H = StuffleAlgebra(ZZ, 'g')
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
        """
        # stuffle algebras in the same variable over any base that coerces in:
        if isinstance(R, StuffleAlgebra) and R.variable_names() == self.variable_names():
            return self.base_ring().has_coerce_map_from(R.base_ring())

        return self.base_ring().has_coerce_map_from(R)

