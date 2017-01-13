r"""
Quantum shuffle algebra on some generators over a base ring. 

Quantum shuffle algebras are noncommutative associative algebras 
with bases generated from words. 

This class is built using parent class ShuffleAlgebra() but uses a 
customized shuffle product for calculating the appropriate Laurent 
polynomial coefficients from the Cartan matrix. 

The effective difference between a quantum shuffle algebra (QSA)
and a normal shuffle algebra (SA) is that the characters in the 
'words' in a QSA are associated to a Cartan matrix that specifies 
a connectivity between them. e.g. if we have the generators 'a','b' 
and associated matrix 
[ 2 -1]
[-1  2]
then the quantum shuffle product 'a'*'b' is given by 
'ab' + q'ba'
powers of the Laurent monomial q are introduced when characters from 
the second word move past characters from the first word. These powers 
are the negative of the corresponding entry in the Cartan matrix

AUTHORS: 

- Mary Barker (2017-01-11) initial version

- Joseph Brown (2017-01-11) initial version


EXAMPLES:

Create a quantum shuffle algebra instance::

    sage: A = QuantumShuffleAlgebra(names=build_alphabet(name='lower'));A
    Quantum Shuffle Algebra on 26 generators ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'] over Univariate Laurent Polynomial Ring in q over Rational Field
    sage: A.is_commutative()
    False

Multiplication::

    sage: A.product_on_basis('abc','de')
    B['abcde'] + q*B['abdce'] + q*B['abdec'] + q*B['adbce'] + q*B['adbec'] + q*B['adebc'] + q*B['dabce'] + q*B['dabec'] + q*B['daebc'] + q*B['deabc']

::

    sage: a = A.term('aa')*A.term('bb');a
    B['aabb'] + q*B['abab'] + q^2*B['abba'] + q^2*B['baab'] + q^3*B['baba'] + q^4*B['bbaa']
    sage: b = A.term('c'); b
    B['c']
    sage: a*b
    B['aabbc'] + q*B['aabcb'] + q^2*B['aacbb'] + q*B['ababc'] + q^2*B['abacb'] + q^2*B['abbac'] + q^2*B['abbca'] + q^2*B['abcab'] + q^3*B['abcba'] + q^2*B['acabb'] + q^3*B['acbab'] + q^4*B['acbba'] + q^2*B['baabc'] + q^3*B['baacb'] + q^3*B['babac'] + q^3*B['babca'] + q^3*B['bacab'] + q^4*B['bacba'] + q^4*B['bbaac'] + q^4*B['bbaca'] + q^4*B['bbcaa'] + q^3*B['bcaab'] + q^4*B['bcaba'] + q^5*B['bcbaa'] + q^2*B['caabb'] + q^3*B['cabab'] + q^4*B['cabba'] + q^4*B['cbaab'] + q^5*B['cbaba'] + q^6*B['cbbaa']

REFERENCES: 

- [Lec2004]_

"""

#*****************************************************************************
#  Copyright (C) 2017 Mary Barker <marybarker103@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.algebras.shuffle_algebra import ShuffleAlgebra
from sage.rings.all import QQ
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.combinat.words.alphabet import Alphabet
from sage.combinat.words.words import Words
from sage.combinat.words.word import Word

class QuantumShuffleAlgebra(ShuffleAlgebra):
    r"""
    INPUT: 

    - ``names`` -- generator names for basis (string or alphabet)

    - ``cartan`` -- Cartan matrix for underlying connectivity structure

    EXAMPLES:: 

        sage: QS = QuantumShuffleAlgebra('abc', CartanMatrix(['A', 3])); QS
        Quantum Shuffle Algebra on 3 generators ['a', 'b', 'c'] over Univariate Laurent Polynomial Ring in q over Rational Field

        sage: mul(QS.gens())
        B['abc'] + q*B['acb'] + q*B['bac'] + q*B['bca'] + q*B['cab'] + q^2*B['cba']

        sage: QS = QuantumShuffleAlgebra(names='abcd'); QS._cartan
        [ 2 -1  0  0]
        [-1  2 -1  0]
        [ 0 -1  2 -1]
        [ 0  0 -1  2]
        
        sage: a, b = QS.term('ab'), QS.term('cc'); a; b
        B['ab']
        B['cc']

        sage: a*b
        B['abcc'] + q*B['acbc'] + q^2*B['accb'] + q*B['cabc'] + q^2*B['cacb'] + q^2*B['ccab']

        sage: (a,b,c,d) = QS.algebra_generators()
        sage: a^2*b
        (q^-2+1)*B['aab'] + (q^-1+q)*B['aba'] + (1+q^2)*B['baa']

        sage: term1, term2 = a^2, b*c
        sage: term1*term2
        (q^-2+1)*B['aabc'] + (q^-1+q)*B['aacb'] + 
        (q^-1+q)*B['abac'] + (q^-1+q)*B['abca'] + 
        (q^-1+q)*B['acab'] + (1+q^2)*B['acba'] + 
        (1+q^2)*B['baac'] + (1+q^2)*B['baca'] + 
        (1+q^2)*B['bcaa'] + (q^-1+q)*B['caab'] + 
        (1+q^2)*B['caba'] + (q+q^3)*B['cbaa']

    """
    def __init__(self, names='ab', cartan=0):
        """
        Initialize quantum shuffle algebra instance. 

        INPUT: 

        - 'names' -- (default: 2) a string or alphabet

        - 'cartan' -- (default: 0) a CartanMatrix() instance

        EXAMPLES::

            sage: QS = QuantumShuffleAlgebra(); QS
            Quantum Shuffle Algebra on 2 generators ['a', 'b'] over Univariate Laurent Polynomial Ring in q over Rational Field
        """
        R = LaurentPolynomialRing(QQ, 'q')

        names = Alphabet(names)

        if cartan != 0:
            self._cartan = cartan
            if names.cardinality() != cartan.nrows():
                print 'Error: Cartan Matrix Is not the same size as the associated generator set'
                return
        else:
            self._cartan = CartanMatrix(['A', names.cardinality()])

        ShuffleAlgebra.__init__(self, R, names)

        self._idx = dict(zip(self.variable_names(),range(names.cardinality())))
        self._alphabet = Alphabet(names)
        self.__ngens = self._alphabet.cardinality()


    def _repr_(self):
        r"""
        Return text representation of this quantum shuffle algebra.

        EXAMPLES::

            sage: F = QuantumShuffleAlgebra()
            sage: F  # indirect doctest
            Quantum Shuffle Algebra on 2 generators ['a', 'b'] over Univariate Laurent Polynomial Ring in q over Rational Field

            sage: QuantumShuffleAlgebra(names='x')
            Quantum Shuffle Algebra on one generator ['x'] over Univariate Laurent Polynomial Ring in q over Rational Field
        """
        if self.__ngens == 1:
            gen = "one generator"
        else:
            gen = "%s generators" %self.__ngens
        return "Quantum Shuffle Algebra on "+ gen +" %s over %s"%(
            self._alphabet.list(), self.base_ring())


    def is_commutative(self):
        r"""
        Return ``False`` as the quantum shuffle algebra is NOT commutative.

        EXAMPLES::

            sage: QS = QuantumShuffleAlgebra(names='xyz')
            sage: QS.is_commutative()
            False
            sage: QS = QuantumShuffleAlgebra(names='uvwxyz')
            sage: QS.is_commutative()
            False
        """
        return False


    def product_on_terms(self, t1, t2):
        r"""
        Returns the quantum shuffle product of two terms. 
        e.g. the product q^-2'ab' * 'a' will be: 
        q^-2 ( 'aba' + q'aab' + q*q^-2'aab') 
        which is  
        (q^-3 + q^-1)'aab' + q^-2'aba'
        in simplified form. 

        INPUT: 

        - 't1' -- finite sum of laurent polynomial multiples of base elements in the quantum shuffle algebra

        - 't2' -- finite sum of laurent polynomial multiples of base elements in the quantum shuffle algebra

        EXAMPLES:: 

            sage: QS = QuantumShuffleAlgebra(names='xyz')

            sage: a = QS.term('xy'); b = QS.term('yy')
            sage: a; b
            B['xy']
            B['yy']
            sage: a*b
            (q^-4+q^-2+1)*B['xyyy'] + (q^-3+q^-1)*B['yxyy'] + (q^-2)*B['yyxy']

            sage: (a, b, c) = QS.algebra_generators(); a; b; c
            B[word: x]
            B[word: y]
            B[word: z]

            sage: a*b*c
            B['xyz'] + q*B['xzy'] + q*B['yxz'] + q*B['yzx'] + q*B['zxy'] + q^2*B['zyx']

            sage: aa = a * a
            sage: aa * a
            (q^-6+2*q^-4+2*q^-2+1)*B['xxx']
            sage: bb = a*b
            sage: aa*bb
            (q^-6+2*q^-4+2*q^-2+1)*B['xxxy'] + 
            (q^-5+2*q^-3+2*q^-1+q)*B['xxyx'] + 
            (q^-4+2*q^-2+2+q^2)*B['xyxx'] + 
            (q^-3+2*q^-1+2*q+q^3)*B['yxxx']
        """
        mylist = list()
        for term1 in t1.terms():
            coef1 = term1.leading_coefficient()
            w1 = str(term1.leading_monomial().support()[0])
            name1 = [[x, 1] for x in str(w1)]
            for term2 in t2.terms():
                coef2 = term2.leading_coefficient()
                w2 = str(term2.leading_monomial().support()[0])
                name2 = [[x, 2] for x in str(w2)]

                name3 = interleave(name2, name1)
                for name in name3: 
                    word = ''.join([y[0] for y in name])
                    nameval = self.calclaurent(name)# self._cartan, self._idx)
                    mylist.append(list((word, coef1 * coef2 * (base_multiple**nameval))))
        return sum(u[1] * self.basis()[u[0]] for u in mylist)


    def product_on_basis(self, w1, w2):
        r"""
        Return the Quantum Shuffle product of two basis 'words' together. 

        INPUT: 

        - 'w1' -- basis element i.e. a finite string composed of generators for the algebra

        - 'w2' -- basis element i.e. a finite string composed of generators for the algebra

        EXAMPLES::

            sage: QS = QuantumShuffleAlgebra(names='xyz')
            sage: (a, b, c) = QS.algebra_generators(); a; b; c
            B[word: x]
            B[word: y]
            B[word: z]

            sage: a*b*c
            B['xyz'] + q*B['xzy'] + q*B['yxz'] + q*B['yzx'] + q*B['zxy'] + q^2*B['zyx']

            sage: QS.product_on_basis('x','y')
            B['xy'] + q*B['yx']
        """
        base_multiple = self.base_ring().gen()
        mylist = list()
        name1 = [[x, 1] for x in str(w1)]
        name2 = [[x, 2] for x in str(w2)]
        name3 = interleave(name2, name1)
        for name in name3:
            word = ''.join([y[0] for y in name])
            nameval = self.calclaurent(name)#, self._cartan, self._idx)
            mylist.append(list((word, base_multiple**nameval)))
        return sum(u[1] * self.basis()[u[0]] for u in mylist)


    def calclaurent(self, base):
        r"""
        Return the power for Laurent polynomial associated to a product of two 
        words (stored in base) given underlying Cartan matrix. 
    
        INPUT: 

        - 'base' -- list of lists. in each list is a series of tuples of the form ('a', idx) where idx is an integer referring to the word it came from and 'a' is a character from that word

        EXAMPLES::
    
            sage: from sage.algebras.quantum_shuffle_algebra import interleave
            sage: QS = QuantumShuffleAlgebra()
            sage: a = QS.term('a')
            sage: b = QS.term('ab')
            sage: alist = [[x, 1] for x in a.leading_monomial().support()[0]]
            sage: blist = [[x, 2] for x in b.leading_monomial().support()[0]]
            sage: shuffled = interleave(alist, blist)
            sage: QS.calclaurent(shuffled[0])
            -1
            sage: QS.calclaurent(shuffled[1])
            -2
            sage: QS.calclaurent(shuffled[2])
            0
        """
        power = 0
        for i, char in enumerate(base):
            if char[1] == 2:
                for j, passing in enumerate(base[i:]):
                    if passing[1] == 1:
                        power -= self._cartan[self._idx[char[0]], self._idx[passing[0]]]
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

        - 'str1' -- a list with each element of the form ['a', 1] where 'a' is a character from the first string in the shuffle product.

        - 'str2' -- a list with each element of the form ['a', 2] where 'a' is a character from the second string in the shuffle product.

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

    if(len(str2) < 2):
        n1 = len(str1)
        n2 = len(str2)
        if n2 > 0: 
            for i in range(min_idx, n1+1):
                mylist.append(str1[0:i] + str2 + str1[i:n1])
        else: 
            mylist = ([str1])
        return mylist

    else:
        mylist = []
        n1 = len(str1)
        n2 = len(str2)
        minvec = range(min_idx, n1+1)
        mychar = str2[0]
        newlist = [str1[0:i] + [mychar] + str1[i:n1] for i in minvec]

        for i, st in enumerate(newlist):
            ret_val = interleave(st, str2[1:n2], minvec[i]+1)
            mylist.extend(ret_val)
        return mylist


