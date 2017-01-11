r"""
The quantum shuffle algebra on some generators over a base ring. 

quantum shuffle algebras are noncommutative associative algebras 
with bases generated from words. 

This class is built using parent class ShuffleAlgebra() but uses a 
customized shuffle product for calculating the appropriate laurent 
polynomial coefficients from the cartan matrix. 

AUTHORS: 

- Mary Barker

- Joseph Brown


REFERENCES: 

Bernard Leclerc, Dual canonical bases, quantum shuffles and q-characters, Mathematische Zeitschrift 246
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

    - ``cartan`` -- cartan matrix for underlying connectivitystructure

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
        initialize quantum shuffle algebra instance

        EXAMPLES: 
            sage: QS = QuantumShuffleAlgebra(); QS
            Quantum Shuffle Algebra on 2 generators ['a', 'b'] over Univariate Laurent Polynomial Ring in q over Rational Field
        """
        R = LaurentPolynomialRing(QQ, 'q')

        names = Alphabet(names)

        if cartan != 0:
            self._cartan = cartan
            if names.cardinality() != cartan.nrows():
                names = Alphabet([str(i) for i in range(0, cartan.nrows())])
        else:
            self._cartan = CartanMatrix(['A', names.cardinality()])

        ShuffleAlgebra.__init__(self, R, names)

        self._idx = dict(zip(self.variable_names(),range(0, names.cardinality())))
        self._alphabet = Alphabet(names)
        self.__ngens = self._alphabet.cardinality()


    def _repr_(self):
        r"""
        Text representation of this quantum shuffle algebra.

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
        Multiply two terms together

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
            w1 = str(term1.leading_monomial().support[0])
            name1 = [[x, 1] for x in str(w1)]
            for term2 in t2.terms():
                coef2 = term2.leading_coefficient()
                w2 = str(term2.leading_monomial().support[0])
                name2 = [[x, 2] for x in str(w2)]

                name3 = interleave(name2, name1)
                for name in name3: 
                    word = ''.join([y[0] for y in name])
                    nameval = calclaurent(name, self._cartan, self._idx)
                    mylist.append(list((word, coef1 * coef2 * (base_multiple**nameval))))
        return sum(u[1] * self.basis()[u[0]] for u in mylist)


    def product_on_basis(self, w1, w2):
        r"""
        Multiply two basis 'words' elements together 

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
            nameval = calclaurent(name, self._cartan, self._idx)
            mylist.append(list((word, base_multiple**nameval)))
        return sum(u[1] * self.basis()[u[0]] for u in mylist)


def calclaurent(base, cartan, idx):
    r"""
    Calculate the power for laurent polynomial given underlying cartan matrix
    """
    power = 0
    for i, char in enumerate(base):
        if char[1] == 2:
            for j, passing in enumerate(base[i:]):
                if passing[1] == 1:
                    power -= cartan[idx[char[0]], idx[passing[0]]]
    return power


def interleave(str1, str2, min_idx=0):
    r"""
    shuffle the two words str1 and str2 together using the normal 
    shuffle method (see, e.g. 
    :mod:`~sage.combinat.words.shuffle_product` and
    :meth:`~sage.combinat.words.finite_word.FiniteWord_class.shuffle()`.
    )
    but recording the base string each character 
    in the final result comes from so that the powers for the 
    laurent polynomial coefficients can be generated
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

