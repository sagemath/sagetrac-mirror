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
from sage.categories.tensor import TensorProductsCategory, tensor
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.alphabet import Alphabet
from sage.combinat.words.words import Words
from sage.combinat.words.word import Word
from sage.misc.cachefunc import cached_method
from sage.sets.family import Family

def indexed_letters(n,names='y'):
    """
    Return a n letter alphabet of the form y_1, \dots, y_n. The default name of the letters is y.

    Beware: the first letter of the generated alphabet is y1 and not y0!

    EXAMPLES::

    sage: from sage.algebras.stuffle_algebra import indexed_letters
    sage: indexed_letters(5)
    {'y1', 'y2', 'y3', 'y4', 'y5'}
    sage: indexed_letters(10,'let')
    {'let1', 'let2', 'let3', 'let4', 'let5', 'let6', 'let7', 'let8', 'let9', 'let10'}

    """
    return Alphabet([str(names)+str(i) for i in range(1,n+1)])

def stuffle(w1,w2,names='y'):
    """
    Return a list containing the words obtained by the quasi shuffle-product of words w1 and w2. 
    If $Y = \left\{ y_1 , \dots , y_n , \dots\right\}$, this product is defined recursively by:
    for all $y_i, \, y_j \in Y$ and for all $u , \, v \in Y^*$,
    \begin{equation}
      \left\lbrace 
    \begin{aligned}
    u \stuffle 1 & = 1 \stuffle u = u~; \\   
    y_i u \stuffle y_j v & = y_i ( u \stuffle y_j v ) + y_j ( y_i u \stuffle v ) + y_{i+j} ( u \stuffle v).
    \end{aligned}
    \right.
    \end{equation}

    By default, the letters are supposed to be y_j, j \geq 1 but one can give another name.

    EXAMPLES::

    sage: from sage.algebras.stuffle_algebra import indexed_letters
    sage: from sage.algebras.stuffle_algebra import stuffle
    sage: A=indexed_letters(10,'y')
    sage: W=Words(A)
    sage: w1=W([W.alphabet()[0],W.alphabet()[0]])
    sage: w2=W([W.alphabet()[1]])
    sage: stuffle(w1,w2)
    [word: y1,y1,y2, word: y1,y2,y1, word: y1,y3, word: y2,y1,y1, word: y3,y1]

    sage: from sage.algebras.stuffle_algebra import stuffle
    sage: from sage.algebras.stuffle_algebra import indexed_letters
    sage: A=indexed_letters(10,'g')
    sage: W=Words(A)
    sage: w1=W([W.alphabet()[0],W.alphabet()[0]])
    sage: w2=W([W.alphabet()[1]])
    sage: stuffle(w1,w2,names='g')
    [word: g1,g1,g2, word: g1,g2,g1, word: g1,g3, word: g2,g1,g1, word: g3,g1]
	    
    """
    W = Words(indexed_letters(40,names))
    letters = W.alphabet()
    if w1 == Word():
        return [w2]
    elif w2 == Word():
	return [w1]
    else:
        t1 = W(Word([w1[0]]))
	t2 = W(Word([w2[0]]))
	r1 = w1[1:]
	r2 = w2[1:]
	l1 = [W(t1 * i) for i in stuffle(r1,w2,names)]
	l2 = [W(t2 * i) for i in stuffle(w1,r2,names)]
        l = letters[letters.index(w1[0])+letters.index(w2[0])+1]
	t3 = W(Word([l]))
	l3 = [W(t3 * i) for i in stuffle(r1,r2,names)]
	return l1+l2+l3

class StuffleAlgebra(CombinatorialFreeModule):
    r"""
    The stuffle algebra over a base ring.

    Stuffle algebras are commutative and associative algebras, with a
    basis indexed by words. The product of two words `w_1 \cdot w_2` is given
    by a deformation of the shuffle product of `w_1` and `w_2` called stuffle product (see above).

    INPUT:

    -  ``R`` -- ring

    -  ``names`` -- name of the letters (default = 'y')

    - ``n`` -- number of letters in the alphabet (default = 40) Stuffle algebras are defined with an infinite alphabet but I do not know yet how to define this structure. Hence, I give a number of variables big enough for "finite time" computations.

    EXAMPLES::

	sage: from sage.algebras.stuffle_algebra import StuffleAlgebra
	sage: A=StuffleAlgebra(QQ); A
	Stuffle Algebra over Rational Field with variables y                                                           
	sage: mul(A.gens()[0:3])
	B[word: y1,y2,y3] + B[word: y1,y3,y2] + B[word: y1,y5] + B[word: y2,y1,y3] + B[word: y2,y3,y1] + B[word: y2,y4] + B[word: y3,y1,y2] + B[word: y3,y2,y1] + 2*B[word: y3,y3] + B[word: y4,y2] + B[word: y5,y1] + B[word: y6]

	sage: mul([ A.gen(i) for i in range(2) ])                                                                      B[word: y1,y2] + B[word: y2,y1] + B[word: y3]

	sage: S = StuffleAlgebra(ZZ); S                                                                             
	Stuffle Algebra over Integer Ring with variables y                                                             
	sage: S.base_ring()                                                                                
	Integer Ring

        sage: G = StuffleAlgebra(S,'g'); G
        Stuffle Algebra over Stuffle Algebra over Integer Ring with variables y with variables g
        sage: G.base_ring()
        Stuffle Algebra over Integer Ring with variables y

    Stuffle algebras commute with their base ring::

        sage: K = StuffleAlgebra(QQ)
	sage: K.is_commutative()
	True
	sage: a,b=K.gens()[0:2]
	sage: L = StuffleAlgebra(K,'g')
	sage: L.is_commutative()
	True
	sage: c,d = L.gens()[0:2]
	sage: s = a*b * c^2; s
	(2*B[word:y1,y2]+2*B[word:y2,y1]+2*B[word:y3])*B[word: g1,g1] + (B[word:y1,y2]+B[word:y2,y1]+B[word:y3])*B[word: g2]
        sage: parent(s)
        Stuffle Algebra over Stuffle Algebra over Rational Field with variables y with variables g
        sage: c^2 * a * b
        (2*B[word:y1,y2]+2*B[word:y2,y1]+2*B[word:y3])*B[word: g1,g1] + (B[word:y1,y2]+B[word:y2,y1]+B[word:y3])*B[word: g2]

    Stuffle algebras are commutative::

        sage: c^3 * b * a * b == c * a * c * b^2 * c
        True
    """
    def __init__(self, R, names='y', n=40):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: F = StuffleAlgebra(QQ); F
            Stuffle Algebra over Rational Field with variables y
            sage: TestSuite(F).run()
	    sage: F = StuffleAlgebra(QQ,'g',n=45); F
	    Stuffle Algebra over Rational Field with variables g
	    sage: F.gen(32)
	    B[word: g33]
	    sage: F = StuffleAlgebra(QQ,n=45); F
	    Stuffle Algebra over Rational Field with variables y
	    sage: F._ngens
	    45

        TESTS::

            sage: StuffleAlgebra(24)
            Traceback (most recent call last):
            ...
            TypeError: argument R must be a ring
        """
        if R not in Rings():
            raise TypeError("argument R must be a ring")
        self._ngens = n
        self._alphabet = indexed_letters(self._ngens,names)
	self._names = names
        self.__ngens = self._alphabet.cardinality()
        CombinatorialFreeModule.__init__(self, R, Words(self._alphabet),
            latex_prefix = "",
            category = (AlgebrasWithBasis(R), CommutativeAlgebras(R),CoalgebrasWithBasis(R)))

    def variable_names(self):
        r"""
        Return the names of the variables.

        EXAMPLES::

            sage: R = StuffleAlgebra(QQ)
            sage: R.variable_names()
	    {'y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8', 'y9', 'y10', 'y11', 'y12', 'y13', 'y14', 'y15', 'y16', 'y17', 'y18', 'y19', 'y20', 'y21', 'y22', 'y23', 'y24', 'y25', 'y26', 'y27', 'y28', 'y29', 'y30', 'y31', 'y32', 'y33', 'y34', 'y35', 'y36', 'y37', 'y38', 'y39', 'y40'}
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

            sage: F = StuffleAlgebra(QQ)
            sage: F  # indirect doctest
            Stuffle Algebra over Rational Field

            sage: StuffleAlgebra(ZZ)
            Stuffle Algebra over Integer Ring
        """
        return "Stuffle Algebra over %s with variables %s"%(self.base_ring(),self._names)

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
            B[word: ]
        """
        return self.basis().keys()([])

    def product_on_basis(self, w1, w2):
        r"""
        Return the product of basis elements ``w1`` and ``w2``, as per
        :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis()`.

        INPUT:

        - ``w1``, ``w2`` -- Basis elements

        EXAMPLES::

	    sage: A=StuffleAlgebra(QQ)
	    sage: W=A.basis().keys()
	    sage: A.product_on_basis(W([W.alphabet()[0]]),W([W.alphabet()[1]]))
	    B[word: y1,y2] + B[word: y2,y1] + B[word: y3]
        """
        return sum(self.basis()[u] for u in stuffle(w1,w2,self._names))

    def gen(self,i):
        r"""
        The ``i``-th generator of the algebra.

        INPUT:

        - ``i`` -- an integer

        EXAMPLES::

	    sage: A=StuffleAlgebra(ZZ)
	    sage: A.gen(2)
	    B[word: y3]

	    sage: A.gen(43)
	    IndexError: argument i (= 43) must be between 0 and 40

        """
        n = self.__ngens
        if i < 0 or not i < n:
            raise IndexError("argument i (= %s) must be between 0 and %s"%(i, n))
        return self.algebra_generators()[i]

    def counit(self,x):
	"""
	Return the counit of x.

	EXAMPLES::

	sage: F = StuffleAlgebra(QQ)
	sage: f=F.an_element()
	sage: F = StuffleAlgebra(QQ)
	sage: f=F.an_element(); f
	B[word: ] + 2*B[word: y1] + 3*B[word: y2]
	sage: f.counit()
	1
	sage: g=F.basis()[Word(['y1','y2'])]
	sage: g.counit()
	0

	INPUT:

	- ``x`` -- an element of self
	"""

        dic = x.monomial_coefficients()
	if dic.has_key(Word()):
            return dic[Word()]
        else:
            return 0

    def coproduct_on_basis(self, w):
	"""
	Return the coproduct of the element of the basis indexed by the word w.

	INPUT:

	- ``w`` -- a word

        EXAMPLES::

	sage: F = StuffleAlgebra(QQ)
	sage: F.coproduct_on_basis(Word(['y3']))
	B[word: ] # B[word: y3] + B[word: y1] # B[word: y2] + B[word: y2] # B[word: y1] + B[word: y3] # B[word: ]
        sage: F.coproduct_on_basis(Word(['y1','y2']))
	B[word: ] # B[word: y1,y2] + B[word: y1] # B[word: y1,y1] + B[word: y1] # B[word: y2] + B[word: y1,y1] # B[word: y1] + B[word: y1,y2] # B[word: ] + B[word: y2] # B[word: y1]

	"""
	if len(w) == 0:
            return tensor((self.one(),self.one()))
	if len(w) == 1:
	    m = self.basis()[w]
	    i = self._alphabet.index(w[0])+1
	    l = [tensor((self.basis()[Word([self._alphabet[j-1]])],self.basis()[Word([self._alphabet[i-j-1]])])) for j in range(1,i)]
	    return tensor((m,self.one()))+tensor((self.one(),m))+sum(l)
        else:
            B = self.basis()
	    result = self.coproduct_on_basis(Word([w[0]]))
	    for i in w[1:]:
		temp1 = self.coproduct_on_basis(Word([i]))
		temp2 = 0
		for ((u1,u2),coeff1) in list(temp1):
                    for ((v1,v2),coeff2) in list(result):
			temp2 += coeff1 * coeff2 * tensor((B[Word(v1)*Word(u1)],B[Word(v2)*Word(u2)]))
	        result = temp2
	    return result
	   
    def coproduct(self,S):
	"""
	Return the coproduct of the series S.

	INPUT:

	- ``S`` -- an element of self

	EXAMPLES::

	sage: F=StuffleAlgebra(QQ)
	sage: f=F.an_element(); f
	B[word: ] + 2*B[word: y1] + 3*B[word: y2]
	sage: F.coproduct(f)
	B[word: ] # B[word: ] + 2*B[word: ] # B[word: y1] + 3*B[word: ] # B[word: y2] + 2*B[word: y1] # B[word: ] + 3*B[word: y1] # B[word: y1] + 3*B[word: y2] # B[word: ]

	sage: F.coproduct(F.one())
	B[word: ] # B[word: ]

	"""
        dic = S.monomial_coefficients()
	return sum([dic[i] * self.coproduct_on_basis(i) for i in dic.keys()])

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators of this algebra.

        EXAMPLES::

	sage: A=StuffleAlgebra(QQ); A                                                                                
	Stuffle Algebra over Rational Field with variables y                                                           
	sage: A.algebra_generators()                                                                                  
	Family (B[word: y1], B[word: y2], B[word: y3], B[word: y4], B[word: y5], B[word: y6], B[word: y7], B[word: y8], B[word: y9], B[word: y10], B[word: y11], B[word: y12], B[word: y13], B[word: y14], B[word: y15], B[word: y16], B[word: y17], B[word: y18], B[word: y19], B[word: y20], B[word: y21], B[word: y22], B[word: y23], B[word: y24], B[word: y25], B[word: y26], B[word: y27], B[word: y28], B[word: y29], B[word: y30], B[word: y31], B[word: y32], B[word: y33], B[word: y34], B[word: y35], B[word: y36], B[word: y37], B[word: y38], B[word: y39], B[word: y40])
        """

        Words = self.basis().keys()
        return Family( [self.monomial(Words([a])) for a in self._alphabet] )
        # FIXME: use this once the keys argument of FiniteFamily will be honoured
        # for the specifying the order of the elements in the family
        #return Family(self._alphabet, lambda a: self.term(self.basis().keys()(a)))

    gens = algebra_generators

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        EXAMPLES::

	    sage: R = StuffleAlgebra(QQ)                                                                   
	    sage: x, y = R.gens()[0:2]                                                                   
	    sage: x,y                                                                                                      
	    (B[word: y1], B[word: y2])                                                                                     
	    sage: R(3)                                                                                   
	    3*B[word: ]                                                                                                    
	    sage: R(x)                                                                                                    
	    B[word: y1]
        """
        P = x.parent()
        if isinstance(P, StuffleAlgebra):
            if P is self:
                return x
            if not (P is self.base_ring()):
                return self.element_class(self, x.monomial_coefficients())
        # ok, not a stuffle algebra element (or should not be viewed as one).
        if isinstance(x, basestring):
            from sage.misc.sage_eval import sage_eval
            return sage_eval(x,locals=self.gens_dict())
        R = self.base_ring()
        # coercion via base ring
        x = R(x)
        if x == 0:
            return self.element_class(self,{})
        else:
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
	    B[word: y1,y2] + B[word: y2,y1] + B[word: y3]

        Elements of the integers coerce in, since there is a coerce map
        from `\ZZ` to GF(7)::

            sage: F.coerce(1)       # indirect doctest
            B[word: ]

        There is no coerce map from `\QQ` to `\GF{7}`::

            sage: F.coerce(2/3)  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Stuffle Algebra over Finite Field of size 7 with variables y

        Elements of the base ring coerce in::

            sage: F.coerce(GF(7)(5))
            5*B[word: ]

        The stuffle algebra over `\ZZ` coerces in, since
        `\ZZ` coerces to `\GF{7}`::

            sage: G = StuffleAlgebra(ZZ)
	    sage: Gx,Gy,Gz = G.gens()[0:3]
            sage: z = F.coerce(Gx**2 * Gy);z
	    2*B[word: y1,y1,y2] + 2*B[word: y1,y2,y1] + 2*B[word: y1,y3] + 2*B[word: y2,y1,y1] + 2*B[word: y2,y2] + 2*B[word: y3,y1] + B[word: y4]
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
            if isinstance(R,StuffleAlgebra):
                if R.variable_names() == self.variable_names():
                    if self.has_coerce_map_from(R.base_ring()):
                        return self(x)
                    else:
                        raise TypeError("no natural map between bases of stuffle algebras")

        except AttributeError:
            pass

        # any ring that coerces to the base ring of this stuffle algebra.
        return self._coerce_try(x, [self.base_ring()])

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

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
        if isinstance(R, StuffleAlgebra):
            if R.variable_names() == self.variable_names():
                if self.base_ring().has_coerce_map_from(R.base_ring()):
                    return True
                else:
                    return False

        return self.base_ring().has_coerce_map_from(R)

