from __future__ import print_function, absolute_import

import sage.modules.glattice
from sage.modules.free_module_morphism import FreeModuleMorphism
from sage.modules.free_module_element import FreeModuleElement
from sage.structure.element import is_Matrix
from inspect import isfunction
from . import matrix_morphism
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace


class GLatticeMorphism_left(sage.categories.morphism.Morphism):
    def __init__(self, mat, domain, codomain):
        D = domain
        C = codomain
        M = mat
        for e in D.basis():
            for g in D.group().gens():
                if M*((D._action(g)(e))) != C._action(g)(M*e):
                    raise TypeError("The morphism does not preserve the action of the group")
        self._free_morphism = D.hom(M.transpose(),C)
        self._domain = D
        self._codomain = C
        self._matrix = M



    def _repr_(self):
        if self.domain() == self.codomain():
            r = "Lattice endomorphism defined by the left action of the matrix\n{!r}\nDomain: {}\n"
            return r.format(self.matrix(), self.domain())
        else:
            r = "Lattice morphism defined by the left action of the matrix\n{!r}\nDomain: {}\nCodomain: {}"
            return r.format(self.matrix(), self.domain(), self.codomain())

    def __call__(self, elt):
        """
        Compute the image of either an element of a lattice, or a sublattice

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: L = GLattice(G,3)
            sage: m = matrix(3, [2,4,2,1,2,1,3,6,3])
            sage: h = L.hom(m); h
            Lattice morphism defined by the matrix
            [2 4 2]
            [1 2 1]
            [3 6 3]
            Domain: Ambient lattice of rank 3 with the trivial action of a group of order 6
            Codomain: Ambient lattice of rank 3 with the trivial action of a group of order 6
            sage: SL = L.sublattice([sum(b for b in L.basis())])
            sage: h(SL)
            Sublattice of degree 3 and rank 1 with the trivial action of a group of order 6 and echelon basis matrix
            [ 8  4 12]
            sage: [b] = SL.basis()
            sage: m*b
            (8, 4, 12)
        """
        from sage.modules.glattice import SubLattice
        if isinstance(elt, SubLattice):
            B = elt.basis()
            m = self._matrix
            iB = [m*b for b in B]
            return self.codomain().sublattice(iB)
        else:
            return self._matrix*elt

    def __add__(self, other):
        return self.sum(other)

    def sum(self, other):
        if isinstance(other, list):
            if len(other) == 0:
                return self
            else:
                mor = self.sum(other[0])
                return mor.sum(other[1:])
        else:
            A = self.domain()
            B = other.domain()
            TA = self.codomain()
            TB = other.codomain()
            if TA == TB : 
                from sage.matrix.special import block_matrix
                C = A.direct_sum(B)
                M = block_matrix(1,2,[self.matrix(),other.matrix()])
                return C.left_morphism(M,TA)
            if TA != TB :
                raise ValueError("The two morphisms must share the same codomain. You might want to use direct_sum")

    def direct_sum(self, other):
        A = self.domain()
        B = other.domain()
        TA = self.codomain()
        TB = other.codomain()
        C = A.direct_sum(B)
        TC = TA.direct_sum(TB)
        from sage.matrix.special import block_diagonal_matrix
        M = block_diagonal_matrix([self.matrix(),other.matrix()])
        return C.left_morphism(M,TC)

    def group_action(self, element, side="Left"):
        D = self.domain()
        if side == "Left":
            M = D.action_matrix(element).inverse()
        elif side == "Right":
            M = D.action_matrix(element)
        else: 
            raise ValueError("The side of the action must be either Left or Right")
        A = self.matrix()*M
        return D.left_morphism(A, self.codomain())


    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain

    def matrix(self):
        return self._matrix

    def free_module_morphism(self):
        return self._free_morphism

    def kernel(self):
        K = self.matrix().transpose().kernel()
        return self.domain().sublattice(K.basis())

    def image(self):
        I = self.matrix().transpose().image()
        return self.codomain().sublattice(I.basis())

    def dual(self):
        return self.codomain().colattice().left_morphism(self.matrix().transpose(), self.domain().colattice())

    def is_identity(self):
        return self.free_module_morphism().is_identity()    

    def characteristic_polynomial(self):
        return self.matrix().characteristic_polynomial()

    def det(self):
        return self.matrix().det()

    def inverse_image(self, M):
        I = self.free_module_morphism().inverse_image(M)
        return self.domain().sublattice(I.basis())

    def is_bijective(self):
        return self.free_module_morphism().is_bijective()

    def is_endomorphism(self):
        return self.domain() == self.codomain()

    def is_injective(self):
        return self.free_module_morphism().is_injective()

    def is_surjective(self):
        return self.free_module_morphism().is_surjective()

    def is_idempotent(self):
        return self.free_module_morphism().is_idempotent()

    def minpoly(self):
        return self.free_module_morphism().minpoly() 

    def trace(self):
        return self.matrix().trace()

    def pre_compose(self, morphism):
        M = self.matrix()*morphism.matrix()
        return morphism.domain().left_morphism(M, self.codomain())

    def group(self):
        return self.domain().group()

    def cokernel(self):
        I = self.image()
        if self.is_endomorphism():
            return GLattice(self.group(),0)
        return self.codomain().quotient_lattice(I)
