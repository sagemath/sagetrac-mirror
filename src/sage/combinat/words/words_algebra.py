from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.category_types import Category_over_base_ring
from copy import copy

class WordsAlgebra(CombinatorialFreeModule):

    class BaseOfWordsAlgebra(Category_over_base_ring):

        def super_categories(self):
            R = self.base().base_ring()
            return [AlgebrasWithBasis(R)]

        class ParentMethods:

            def product_on_basis(self, left, right):
                return self( left + right )

            def one_basis(self):
                return self.monomial(self._basis_keys([]))

        class ElementMethods:
            def e(self, i):
                """
                Implementation of the e operator.

                EXAMPLES::

                    sage: W = Words(4)
                    sage: A = W.algebra(QQ)
                    sage: w = A(W([4,1,2,2,2,1,3]))
                    sage: w.e(1)
                    B[word: 4112213]
                    sage: w.e(3)
                    B[word: ]
                    sage: A(W([2,1,1,1,2,2,1,1,1,1,2])).e(1)
                    B[word: 21112211111]

                """
                W = self.parent()
                terms = []
                for t in self.monomials():
                    r = t._last_letter_subword(i)
                    w = t.support_of_term()
                    if r == (-1):
                        terms.append((W.one_basis(), self.coefficient(w)))
                    else:
                        res = list(w)
                        res[r] = i
                        terms.append((W.monomial(W._basis_keys(res)), self.coefficient(w)))
                return W.linear_combination(terms)

            def f(self, i):
                """
                Implementation of the f operator.

                TESTS::

                    sage: W = Words(4)
                    sage: A = W.algebra(QQ)
                    sage: A(W([2,1,1,1,2,2,1,1,1,1,2])).f(1)
                    B[word: 21112211122]
                    sage: A(W([2,1,1,1,2,2,1,1,1,1,2])).f(2)
                    B[word: 21112211132]
                """
                W = self.parent()
                terms = []
                for t in self.monomials():
                    r = t._last_letter_subword(i)
                    w = t.support_of_term()
                    if r == 0:
                        terms.append((W.one_basis(), self.coefficient(w)))
                    else:
                        res = list(w)
                        res[r-1] = i+1
                        terms.append((W.monomial(W._basis_keys(res)), self.coefficient(w)))
                return W.linear_combination(terms)

            def sigma(self, i):
                """
                Implementation of the sigma operator.

                TESTS::

                    sage: W = Words(4)
                    sage: A = W.algebra(QQ)
                    sage: A(W([2,1,1,1,2,2,1,1,1,1,2])).f(1)
                    B[word: 21112211122]

                """
                W = self.parent()
                terms = []
                for t in self.monomials():
                    w = t.support_of_term()
                    sub = t._extract_subword(i)

                    n = len(sub)
                    j = 0
                    res = list(w)
                    #count number of a_i
                    r = len([s for s in sub if res[s] == i])
                    s = n - r #number of a_{i+1}
                    for k in range(s):
                        res[sub[k]] = i
                    for k in range(s, r+s):
                        res[sub[k]] = i+1
                    terms.append((W.monomial(W._basis_keys(res)), self.coefficient(w)))
                return W.linear_combination(terms)


            def zeta(self):
                """
                Implementation of the zeta operator.

                TESTS::

                    sage: W = Words(4)
                    sage: A = W.algebra(QQ)
                    sage: A(W([3,4,2,2,1,1,1])).zeta()
                    B[word: 4221113]

                And this operator commute with sigma.

                TESTS::

                    sage: W = Words(4)
                    sage: A = W.algebra(QQ)
                    sage: A(W([3,4,2,2,1,1,1])).zeta().sigma(2) == A(W([3,4,2,2,1,1,1])).sigma(2).zeta()
                    True

                """
                W = self.parent()
                terms = []
                for t in self.monomials():
                    w = t.support_of_term()
                    res = list(w)
                    res = res[1:] + [res[0]]
                    terms.append((W.monomial(W._basis_keys(res)), self.coefficient(w)))
                return W.linear_combination(terms)


            def _extract_subword(self, i):
                """
                Return a list of index of ``self`` restricted on the two letter
                alaphabet i, i+1 where we recursively removed all a_i{i+1} a_i.
                This method has been implemented in order to simplify the code of
                the operators e and f.

                INPUT:

                    An element of the basis
                """
                w = self.support_of_term()
                res = range(len(w))
                #First we take the subword with letters i and i+1
                res = [j for j in range(len(w)) if w[j] in {i, i+1}]
                #Then we recursively found the right subword by deleting
                #consecutives a_{i+1} a_{i}.
                work = True
                while work:
                    work = False
                    to_del = []
                    j = 0
                    #find the indexes...
                    while j < len(res)-1:
                        if (w[res[j]] == i+1) and (w[res[j+1]] == i):
                            to_del.append(j)
                            j = j+1
                            work = True
                        j = j+1
                    #...then remove them in reverse order to preserve indexing
                    for k in reversed(to_del):
                        res.remove(res[k+1])
                        res.remove(res[k])
                return res

            def _last_letter_subword(self, i):
                w = self.support_of_term()
                subword = self._extract_subword(i)
                n = len(subword)
                #Find the last element `i`
                found = False
                j = 0
                while j < n and w[subword[j]] <> i+1:
                    j = j+1
                if j == n:
                    return (-1)
                else:
                    return subword[j]

    def __init__(self, W, base_ring):
        CombinatorialFreeModule.__init__(self, base_ring,
                basis_keys=W,
                category=WordsAlgebra.BaseOfWordsAlgebra(base_ring))
