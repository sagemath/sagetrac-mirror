r"""
Quantum Shuffle Algebras

Quantum shuffle algebras are noncommutative associative algebras
with bases generated from words.

This class is built using parent class ShuffleAlgebra() but uses a
customized shuffle product for calculating the appropriate Laurent
polynomial coefficients from the Cartan matrix.

The effective difference between a quantum shuffle algebra (QSA)
and a normal shuffle algebra (SA) is that the characters in the
words in a QSA are associated to a Cartan matrix that specifies
a connectivity between them. For example, if we have the generators
`a`,`b` and a Cartan matrix:

.. MATH::

    \begin{bmatrix}
    2 & -1 \\
    -1 & 2
    \end{bmatrix}

then the quantum shuffle product `a \cdot b` is given by `ab + q ba`,
powers of the Laurent monomial `q` are introduced when characters from
when the second word move past characters from the first word. These
powers are the negative of the corresponding entry in the Cartan matrix.

AUTHORS: 

- Mary Barker (2017-01-11) initial version
- Joseph Brown (2017-01-11) initial version

EXAMPLES:

Create a quantum shuffle algebra instance::

    sage: R.<q> = LaurentPolynomialRing(QQ)
    sage: A = QuantumShuffleAlgebra(q, names=build_alphabet(name='lower')); A
    Quantum Shuffle Algebra on 26 generators ['a', ..., 'z'] over
     Univariate Laurent Polynomial Ring in q over Rational Field
    sage: A.is_commutative()
    False

Multiplication::

    sage: A.product_on_basis('abc','de')
    B[abcde] + q*B[abdce] + q*B[abdec] + q*B[adbce] + q*B[adbec]
     + q*B[adebc] + q*B[dabce] + q*B[dabec] + q*B[daebc] + q*B[deabc]

::

    sage: x = A.term('aa') * A.term('bb'); x
    B[aabb] + q*B[abab] + q^2*B[abba] + q^2*B[baab] + q^3*B[baba] + q^4*B[bbaa]
    sage: c = A.term('c'); c
    B[c]
    sage: x * c
    B[aabbc] + q*B[aabcb] + q^2*B[aacbb] + q*B[ababc] + q^2*B[abacb]
     + q^2*B[abbac] + q^2*B[abbca] + q^2*B[abcab] + q^3*B[abcba]
     + q^2*B[acabb] + q^3*B[acbab] + q^4*B[acbba] + q^2*B[baabc]
     + q^3*B[baacb] + q^3*B[babac] + q^3*B[babca] + q^3*B[bacab]
     + q^4*B[bacba] + q^4*B[bbaac] + q^4*B[bbaca] + q^4*B[bbcaa]
     + q^3*B[bcaab] + q^4*B[bcaba] + q^5*B[bcbaa] + q^2*B[caabb]
     + q^3*B[cabab] + q^4*B[cabba] + q^4*B[cbaab] + q^5*B[cbaba]
     + q^6*B[cbbaa]
"""

#*****************************************************************************
#  Copyright (C) 2017 Mary Barker <marybarker103@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.alphabet import Alphabet
from sage.combinat.words.words import Words
from sage.combinat.words.word import Word
from sage.categories.algebras import Algebras

class QuantumShuffleAlgebra(CombinatorialFreeModule):
    r"""
    The quantum shuffle algebra.

    INPUT: 

    - ``q`` -- the quantum parameter
    - ``names`` -- generator names for basis (string or alphabet)
    - ``cartan`` -- (optional) defining Cartan matrix; if not specified,
      then it is of type `A_n`, where `n` are the number of generators
    - ``R`` -- (default: ``q.parent()``) the base ring containing ``q``

    EXAMPLES:: 

        sage: R.<q> = LaurentPolynomialRing(QQ)
        sage: QS = QuantumShuffleAlgebra(q, 'abc', CartanMatrix(['A', 3])); QS
        Quantum Shuffle Algebra on 3 generators ['a', 'b', 'c'] over
         Univariate Laurent Polynomial Ring in q over Rational Field

        sage: mul(QS.gens())
        B[abc] + q*B[acb] + q*B[bac] + q*B[bca] + q*B[cab] + q^2*B[cba]

        sage: QS = QuantumShuffleAlgebra(q, names='abcd')
        sage: QS.cartan_matrix()
        [ 2 -1  0  0]
        [-1  2 -1  0]
        [ 0 -1  2 -1]
        [ 0  0 -1  2]
        
        sage: x, y = QS.term('ab'), QS.term('cc')
        sage: x
        B[ab]
        sage: y
        B[cc]

        sage: x * y
        B[abcc] + q*B[acbc] + q^2*B[accb] + q*B[cabc]
         + q^2*B[cacb] + q^2*B[ccab]

        sage: a,b,c,d = QS.algebra_generators()
        sage: a^2 * b
        (q^-2+1)*B[aab] + (q^-1+q)*B[aba] + (1+q^2)*B[baa]

        sage: term1, term2 = a^2, b * c
        sage: term1 * term2
        (q^-2+1)*B[aabc] + (q^-1+q)*B[aacb]
         + (q^-1+q)*B[abac] + (q^-1+q)*B[abca]
         + (q^-1+q)*B[acab] + (1+q^2)*B[acba]
         + (1+q^2)*B[baac] + (1+q^2)*B[baca]
         + (1+q^2)*B[bcaa] + (q^-1+q)*B[caab]
         + (1+q^2)*B[caba] + (q+q^3)*B[cbaa]

    REFERENCES: 

    - [Lec2004]_
    """
    @staticmethod
    def __classcall_private__(cls, q, names, cartan_matrix=None, R=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: F1 = QuantumShuffleAlgebra(q, 'xyz')
            sage: F2 = QuantumShuffleAlgebra(q, ['x','y','z'])
            sage: F3 = QuantumShuffleAlgebra(q, Alphabet('xyz'))
            sage: F1 is F2 and F1 is F3
            True
        """
        #R = LaurentPolynomialRing(QQ, 'q')
        names = Alphabet(names)

        if cartan_matrix is None:
            cartan_matrix = CartanMatrix(['A', names.cardinality()])
        else:
            cartan_matrix = CartanMatrix(cartan_matrix)
            if names.cardinality() != cartan_matrix.nrows():
                raise ValueError("Cartan matrix is not the same size"
                                 " as the associated generator set")
        if R is None:
            R = q.parent()
        else:
            q = R(q)
        return super(QuantumShuffleAlgebra, cls).__classcall__(cls, q, names,
                                                               cartan_matrix, R)

    def __init__(self, q, names, cartan_matrix, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: QS = QuantumShuffleAlgebra(q, 'ab')
            sage: TestSuite(QS).run()
        """
        self._cartan_matrix = cartan_matrix
        self._alphabet = names
        self._idx = {a: i for i,a in enumerate(self._alphabet)}
        self.__ngens = self._alphabet.cardinality()
        self._q = q

        CombinatorialFreeModule.__init__(self, R, Words(names, infinite=False),
                                         latex_prefix="", string_quotes=False,
                                         category=Algebras(R).WithBasis())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: F = QuantumShuffleAlgebra(q, 'ab')
            sage: F
            Quantum Shuffle Algebra on 2 generators ['a', 'b'] over
             Univariate Laurent Polynomial Ring in q over Rational Field

            sage: QuantumShuffleAlgebra(q, names='x')
            Quantum Shuffle Algebra on one generator ['x'] over
             Univariate Laurent Polynomial Ring in q over Rational Field
        """
        if self.__ngens == 1:
            gen = "one generator"
        else:
            gen = "{} generators".format(self.__ngens)
        return "Quantum Shuffle Algebra on {} {} over {}".format(
            gen, self._alphabet.list(), self.base_ring())

    def _repr_term(self, w):
        """
        Return a string representation of the term indexed by ``w``.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: QS = QuantumShuffleAlgebra(q, 'abc')
            sage: QS._repr_term(QS._indices('abccabbac'))
            'B[abccabbac]'
        """
        return super(QuantumShuffleAlgebra, self)._repr_term(str(w))

    def cartan_matrix(self):
        """
        Return the defining Cartan matrix of ``self``.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: QS = QuantumShuffleAlgebra(q, 'abc', CartanMatrix(['B', 3]))
            sage: QS.cartan_matrix()
            [ 2 -1  0]
            [-1  2 -1]
            [ 0 -2  2]
            sage: QS = QuantumShuffleAlgebra(q, names='abcd')
            sage: QS.cartan_matrix()
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [ 0  0 -1  2]
        """
        return self._cartan_matrix

    def is_commutative(self):
        r"""
        Return if the quantum shuffle algebra is commutative.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: QS = QuantumShuffleAlgebra(q, names='xyz')
            sage: QS.is_commutative()
            False
            sage: QS = QuantumShuffleAlgebra(q, names='uvwxyz')
            sage: QS.is_commutative()
            False
            sage: QS = QuantumShuffleAlgebra(q, names='x')
            sage: QS.is_commutative()
            True
        """
        return len(self._alphabet) < 2

    def variable_names(self):
        r"""
        Return the names of the variables.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: QS = QuantumShuffleAlgebra(q, 'xy')
            sage: QS.variable_names()
            {'x', 'y'}
        """
        return self._alphabet

    @cached_method
    def one_basis(self):
        r"""
        Return the empty word, which index of `1` of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: A = QuantumShuffleAlgebra(q, 'ab')
            sage: A.one_basis()
            word:
            sage: A.one()
            B[]
        """
        return self.basis().keys()([])

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: A = QuantumShuffleAlgebra(q, 'fgh')
            sage: A.algebra_generators()
            Family (B[f], B[g], B[h])

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: A = QuantumShuffleAlgebra(q, ['x1','x2'])
            sage: A.algebra_generators()
            Family (B[x1], B[x2])
        """
        I = self.basis().keys()
        return Family( [self.monomial(I([a])) for a in self._alphabet] )

    @cached_method
    def gens(self):
        """
        Return a tuple of the generators of ``self``.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: A = QuantumShuffleAlgebra(q, 'fgh')
            sage: A.gens()
            (B[f], B[g], B[h])
        """
        return tuple(self.algebra_generators())

    def product_on_terms(self, t1, t2):
        r"""
        Return the quantum shuffle product of two terms.

        For example The product q^-2'ab' * 'a' will be:
        q^-2 ( 'aba' + q'aab' + q*q^-2'aab')
        which is
        (q^-3 + q^-1)'aab' + q^-2'aba'
        in simplified form.

        INPUT: 

        - ``t1`` -- finite sum of laurent polynomial multiples of base elements in the quantum shuffle algebra
        - ``t2`` -- finite sum of laurent polynomial multiples of base elements in the quantum shuffle algebra

        EXAMPLES:: 

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: QS = QuantumShuffleAlgebra(q, names='xyz')

            sage: a = QS.term('xy'); b = QS.term('yy')
            sage: a * b
            (q^-4+q^-2+1)*B[xyyy] + (q^-3+q^-1)*B[yxyy] + (q^-2)*B[yyxy]

            sage: a, b, c = QS.algebra_generators()
            sage: a * b * c
            B[xyz] + q*B[xzy] + q*B[yxz]
             + q*B[yzx] + q*B[zxy] + q^2*B[zyx]

            sage: aa = a * a
            sage: aa * a
            (q^-6+2*q^-4+2*q^-2+1)*B[xxx]
            sage: ab = a * b
            sage: aa * ab
            (q^-6+2*q^-4+2*q^-2+1)*B[xxxy]
             + (q^-5+2*q^-3+2*q^-1+q)*B[xxyx]
             + (q^-4+2*q^-2+2+q^2)*B[xyxx]
             + (q^-3+2*q^-1+2*q+q^3)*B[yxxx]
        """
        mylist = []
        for term1 in t1.terms():
            coef1 = term1.leading_coefficient()
            w1 = str(term1.leading_monomial().support()[0])
            name1 = [[x, 1] for x in w1]
            for term2 in t2.terms():
                coef2 = term2.leading_coefficient()
                w2 = str(term2.leading_monomial().support()[0])
                name2 = [[x, 2] for x in w2]

                name3 = interleave(name2, name1)
                for name in name3: 
                    w = ''.join([y[0] for y in name])
                    nameval = self._calc_laurent(name)
                    mylist.append([w, coef1 * coef2 * (self._q**nameval)])
        return sum(u[1] * self.basis()[u[0]] for u in mylist)

    def product_on_basis(self, w1, w2):
        r"""
        Return the Quantum shuffle product of basis elements indexed
        by words ``w1`` and ``w2`` in ``self``.

        INPUT: 

        - ``w1`` -- a finite word in the defining alphabet
        - ``w2`` -- a finite word in the defining alphabet

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: QS = QuantumShuffleAlgebra(q, names='xyz')
            sage: x, y, z = QS.algebra_generators()
            sage: x * y * z
            B[xyz] + q*B[xzy] + q*B[yxz] + q*B[yzx] + q*B[zxy] + q^2*B[zyx]

            sage: QS.product_on_basis('x', 'y')
            B[xy] + q*B[yx]
        """
        I = self.basis().keys()
        mylist = []
        name1 = [[x, 1] for x in str(w1)]
        name2 = [[x, 2] for x in str(w2)]
        name3 = interleave(name2, name1)
        for name in name3:
            w = ''.join(y[0] for y in name)
            nameval = self._calc_laurent(name)
            mylist.append([I(w) , self._q**nameval])
        return sum(u[1] * self.basis()[u[0]] for u in mylist)

    def _calc_laurent(self, base):
        r"""
        Return the power for Laurent polynomial associated to a product of two 
        words (stored in base) given underlying Cartan matrix. 
    
        INPUT: 

        - ``base`` -- list of lists; in each list is a series of tuples of the form ('a', idx) where idx is an integer referring to the word it came from and 'a' is a character from that word

        EXAMPLES::
    
            sage: from sage.algebras.quantum_shuffle_algebra import interleave
            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: QS = QuantumShuffleAlgebra(q, 'ab')
            sage: a = QS.term('a')
            sage: b = QS.term('ab')
            sage: alist = [[x, 1] for x in a.leading_monomial().support()[0]]
            sage: blist = [[x, 2] for x in b.leading_monomial().support()[0]]
            sage: shuffled = interleave(alist, blist)
            sage: QS._calc_laurent(shuffled[0])
            -1
            sage: QS._calc_laurent(shuffled[1])
            -2
            sage: QS._calc_laurent(shuffled[2])
            0
        """
        power = 0
        for i, char in enumerate(base):
            if char[1] == 2:
                for j, passing in enumerate(base[i:]):
                    if passing[1] == 1:
                        power -= self._cartan_matrix[self._idx[char[0]],
                                                     self._idx[passing[0]]]
        return power

def interleave(str1, str2, min_idx=0):
    r"""
    Return the shuffle product of two words str1 and str2 together 
    using the normal shuffle method see, e.g.
    :mod:`sage.combinat.words.shuffle_product` and
    :meth:`sage.combinat.words.finite_word.FiniteWord_class.shuffle()`,
    but recording the base string each character
    in the final result comes from so that the powers for the
    Laurent polynomial coefficients can be generated.

    INPUT:

    - ``str1`` -- a list with each element of the form ['a', 1] where 'a' is a character from the first string in the shuffle product

    - ``str2`` -- a list with each element of the form ['a', 2] where 'a' is a character from the second string in the shuffle product

    EXAMPLES::

        sage: from sage.algebras.quantum_shuffle_algebra import interleave
        sage: word1 = 'abc'
        sage: word2 = 'def'
        sage: list1 = [[x, 1] for x in word1]
        sage: list2 = [[x, 2] for x in word2]
        sage: interleave(list1, list2)
        [[['d', 2], ['e', 2], ['f', 2], ['a', 1], ['b', 1], ['c', 1]],
         [['d', 2], ['e', 2], ['a', 1], ['f', 2], ['b', 1], ['c', 1]],
         [['d', 2], ['e', 2], ['a', 1], ['b', 1], ['f', 2], ['c', 1]],
         [['d', 2], ['e', 2], ['a', 1], ['b', 1], ['c', 1], ['f', 2]],
         [['d', 2], ['a', 1], ['e', 2], ['f', 2], ['b', 1], ['c', 1]],
         [['d', 2], ['a', 1], ['e', 2], ['b', 1], ['f', 2], ['c', 1]],
         [['d', 2], ['a', 1], ['e', 2], ['b', 1], ['c', 1], ['f', 2]],
         [['d', 2], ['a', 1], ['b', 1], ['e', 2], ['f', 2], ['c', 1]],
         [['d', 2], ['a', 1], ['b', 1], ['e', 2], ['c', 1], ['f', 2]],
         [['d', 2], ['a', 1], ['b', 1], ['c', 1], ['e', 2], ['f', 2]],
         [['a', 1], ['d', 2], ['e', 2], ['f', 2], ['b', 1], ['c', 1]],
         [['a', 1], ['d', 2], ['e', 2], ['b', 1], ['f', 2], ['c', 1]],
         [['a', 1], ['d', 2], ['e', 2], ['b', 1], ['c', 1], ['f', 2]],
         [['a', 1], ['d', 2], ['b', 1], ['e', 2], ['f', 2], ['c', 1]],
         [['a', 1], ['d', 2], ['b', 1], ['e', 2], ['c', 1], ['f', 2]],
         [['a', 1], ['d', 2], ['b', 1], ['c', 1], ['e', 2], ['f', 2]],
         [['a', 1], ['b', 1], ['d', 2], ['e', 2], ['f', 2], ['c', 1]],
         [['a', 1], ['b', 1], ['d', 2], ['e', 2], ['c', 1], ['f', 2]],
         [['a', 1], ['b', 1], ['d', 2], ['c', 1], ['e', 2], ['f', 2]],
         [['a', 1], ['b', 1], ['c', 1], ['d', 2], ['e', 2], ['f', 2]]]
    """
    mylist = []
    n1 = len(str1)
    n2 = len(str2)

    if n2 < 2:
        if n2 > 0: 
            for i in range(min_idx, n1+1):
                mylist.append(str1[0:i] + str2 + str1[i:n1])
        else: 
            mylist = ([str1])
        return mylist

    else:
        minvec = range(min_idx, n1+1)
        mychar = str2[0]
        newlist = [str1[0:i] + [mychar] + str1[i:n1] for i in minvec]

        for i, st in enumerate(newlist):
            ret_val = interleave(st, str2[1:n2], minvec[i]+1)
            mylist.extend(ret_val)
        return mylist

