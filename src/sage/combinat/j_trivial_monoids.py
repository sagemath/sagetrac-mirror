from sage.categories.j_trivial_monoids import JTrivialMonoids
from sage.matrix.constructor import Matrix
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.rings.finite_rings.constructor import GF
from sage.sets.family import Family
from itertools import product


"""
Implement the monoid of Unitriangular Boolean matrices, which is J-trivial
"""
class UnitriangularBooleanMatrices(Parent, UniqueRepresentation):

    def __init__(self, n):
        self.n = n
        Parent.__init__(self, category=JTrivialMonoids().Finite())

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.combinat.j_trivial_monoids import UnitriangularBooleanMatrices
            sage: UnitriangularBooleanMatrices(4)
            Monoid of Unitriangular Boolean Matrices 4 x 4
        """
        return "Monoid of Unitriangular Boolean Matrices %s x %s"%(self.n, self.n)

    @cached_method
    def one(self):
        """
        Return the identity element of the monoid ``self``
        """
        M = Matrix(GF(2), self.n, self.n, lambda i, j: i == j)
        M.set_immutable()
        return self(M)

    def product(self, right, left):
        """
        Return the product of the elements ``left`` and ``right`` in ``self``

        EXAMPLES::

            sage: from sage.combinat.j_trivial_monoids import UnitriangularBooleanMatrices
            sage: U4 = UnitriangularBooleanMatrices(4)
            sage: a = U4(Matrix(GF(2), 4, 4, [[1,0,1,1],[0,1,1,0],[0,0,1,1],[0,0,0,1]]))
            sage: b = U4(Matrix(GF(2), 4, 4, [[1,1,0,1],[0,1,0,0],[0,0,1,1],[0,0,0,1]]))
            sage: a * b
            [1 1 1 1]
            [0 1 1 1]
            [0 0 1 1]
            [0 0 0 1]
            sage: b * a
            [1 1 1 1]
            [0 1 1 0]
            [0 0 1 1]
            [0 0 0 1]
        """
        M = Matrix(GF(2), self.n, self.n,
                lambda i, j: right.value[i][j] or
                left.value[i][j] or
                any(right.value[i][k] and left.value[k][j]
                    for k in range(self.n)))
        M.set_immutable()
        return self(M)

    def semigroup_generators(self):
        """
        Return the generating family of the monoid ``self``, namely all
        unitriangular boolean matrices with n+1 1s

        TODO:: Make immutable matrices in an easier way... This method should
        be only one line.

        EXAMPLES::

            sage: from sage.combinat.j_trivial_monoids import UnitriangularBooleanMatrices
            sage: U4 = UnitriangularBooleanMatrices(4)
            sage: len(U4.semigroup_generators())
            16
        """
        res = [Matrix(GF(2), self.n, self.n,
            lambda i, j: (i == j) or (i <= j and (i == k and j == l))) for k,l in
            product(range(self.n), range(self.n))]
        for M in res:
            M.set_immutable()

        return Family(self(M) for M in res)

    def cardinality(self):
        """
        Return the cardinality of ``self``

        EXAMPLES::

            sage: from sage.combinat.j_trivial_monoids import UnitriangularBooleanMatrices
            sage: UnitriangularBooleanMatrices(6).cardinality()
            32768
        """
        return 2**((self.n**2-self.n)//2)

    def _an_element_(self):
        """
        Return an element of ``self``
        """
        return self.one()

    def random_element(self):
        """
        Return a random element of ``self``.

        EXAMPLES::

            sage: from sage.combinat.j_trivial_monoids import UnitriangularBooleanMatrices
            sage: U = UnitriangularBooleanMatrices(5)
            sage: U.an_element() # random
            [1 1 0 0 0]
            [0 1 0 1 1]
            [0 0 1 1 1]
            [0 0 0 1 0]
            [0 0 0 0 1]
            sage: U.an_element() # random
            [1 0 0 1 0]
            [0 1 1 1 0]
            [0 0 1 1 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
        """
        from sage.misc.prandom import random
        M = Matrix(GF(2), self.n, self.n,
                lambda i, j: (i <= j) and ((i == j) or (random()<.5)))
        M.set_immutable()
        return self(M)

    class Element(ElementWrapper):
        wrapped_class = Matrix
