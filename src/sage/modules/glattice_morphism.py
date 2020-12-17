# -*- coding: utf-8 -*-
r"""
Classes of G-Lattices morphisms
===============================

G-Lattices morphismsare defined specifying a matrix acting on the left.


To specify a lattice morphism, you need to specify:

- ``m`` -- a matrix acting on the left, equivariant under the action of the group.

- ``C`` -- the codomain of the morphism to test equivariance. If the matrix is square 
    and the codomain is left blank, it is assumed that the morphism is an endomorphism.


Attributes of a lattice morphism
---------------------------------

- ``morphism._domain`` -- the domain of the morphism.

- ``morphism._codomain`` -- the codomain of the morphism.

- ``morphism._matrix`` -- the matrix corresponding to the morphism.

- ``morphism._free_morphism`` -- the morphism as morphism of free modules acting on the right.

Methods of a lattice morphism
-----------------------------

- :meth:`GLatticeMorphism_left.group` -- the group acting on the domain and codomain.

- :meth:`GLatticeMorphism_left.sum` -- compute the sum of two morphisms.

- :meth:`GLatticeMorphism_left.group_action` -- computes the action of a group element on the morphism.

- :meth:`GLatticeMorphism_left.domain` -- return the domain of the morphism.

- :meth:`GLatticeMorphism_left.codomain` -- the codomain of the morphism.

- :meth:`GLatticeMorphism_left.matrix` -- the matrix corresponding to the morphism.

- :meth:`GLatticeMorphism_left.free_module_morphism` -- the morphism as homomorphism of free modules.

- :meth:`GLatticeMorphism_left.kernel` -- kernel of the morphism.

- :meth:`GLatticeMorphism_left.image` -- image of the morphism.

- :meth:`GLatticeMorphism_left.dual` -- dual of the morphism, form the dual of the codomain, to the dual of the domain.

- :meth:`GLatticeMorphism_left.is_identity` -- determine if the morphism is the identity morphism.

- :meth:`GLatticeMorphism_left.characteristic_polynomial` -- characteristic polynomial of the morphism.

- :meth:`GLatticeMorphism_left.det` -- determinant of the morphism

- :meth:`GLatticeMorphism_left.inverse_image` -- the preimage of a sublattice.

- :meth:`GLatticeMorphism_left.is_bijective` -- determine if the morphism is a bijection.

- :meth:`GLatticeMorphism_left.is_endomorphism` -- determine if the morphism is an endomorphism.

- :meth:`GLatticeMorphism_left.is_injective` -- determine if the morphism is injective.

- :meth:`GLatticeMorphism_left.is_surjective` -- determine if the morphism is surjective.

- :meth:`GLatticeMorphism_left.is_idempotent` -- determine if the morphism is idempotent.

- :meth:`GLatticeMorphism_left.minpoly` -- compute the minimal polynomial of the morphism.

- :meth:`GLatticeMorphism_left.trace` -- compute the trace of the morphism.

- :meth:`GLatticeMorphism_left.pre_compose` -- compute the precomposition of the morphism with another one.

- :meth:`GLatticeMorphism_left.cokernel` -- compute the cokernel of the morphism.


EXAMPLES::

    sage: L = GLattice([2])
    sage: h1 = L.identity_morphism()
    sage: h2 = L.norm()
    sage: h1 + h2
    Lattice endomorphism defined by the left action of the matrix
    [2 1]
    [1 2]
    Domain: Ambient lattice of rank 2 with a faithful action by a group of order 2
    sage: h1 + [h2, h2]
    Lattice endomorphism defined by the left action of the matrix
    [3 2]
    [2 3]
    Domain: Ambient lattice of rank 2 with a faithful action by a group of order 2
    sage: L.zero_sum_sublattice().injection_morphism()
    Lattice morphism defined by the left action of the matrix
    [ 1]
    [-1]
    Domain: Ambient lattice of rank 1 with a faithful action by a group of order 2
    Codomain: Ambient lattice of rank 2 with a faithful action by a group of order 2
    sage: L.diagonal_embedding()
    Lattice morphism defined by the left action of the matrix
    [1 0]
    [0 1]
    [---]
    [1 0]
    [0 1]
    Domain: Ambient lattice of rank 2 with a faithful action by a group of order 2
    Codomain: Ambient lattice of rank 4 with a faithful action by a group of order 2
    sage: L.surjection_from_square()
    Lattice morphism defined by the left action of the matrix
    [1 0|1 0]
    [0 1|0 1]
    Domain: Ambient lattice of rank 4 with a faithful action by a group of order 2
    Codomain: Ambient lattice of rank 2 with a faithful action by a group of order 2
  
    
::

    sage: L = GLattice(SymmetricGroup(3))
    sage: SL = L.zero_sum_sublattice()
    sage: [Q, h] = L.quotient_lattice(SL, True, True)
    sage: h
    Lattice morphism defined by the left action of the matrix
    [1 1 1]
    Domain: Ambient lattice of rank 3 with a faithful action by a group of order 6
    Codomain: Ambient lattice of rank 1 with the trivial action of a group of order 6
    sage: h.is_surjective()
    True
    sage: SL2 = L.sublattice([sum(L.basis())])
    sage: [Q2, h2] = L.quotient_lattice(SL2, True, True); h2
    Lattice morphism defined by the left action of the matrix
    [-1  0  1]
    [-1  1  0]
    Domain: Ambient lattice of rank 3 with a faithful action by a group of order 6
    Codomain: Ambient lattice of rank 2 with a faithful action by a group of order 6




If one enters the matrix manually, one must make sure that the morphism preserves the action of the group.

::

    sage: L = GLattice([4])
    sage: L.left_morphism(matrix(4, range(16)))
    Traceback (most recent call last):
    ...
    TypeError: The morphism does not preserve the action of the group
    sage: L.left_morphism(identity_matrix(4))
    Lattice endomorphism defined by the left action of the matrix
    [1 0 0 0]
    [0 1 0 0]
    [0 0 1 0]
    [0 0 0 1]
    Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
    sage: L2 = GLattice([4], 1)
    sage: L.left_morphism(matrix(1, [1, 1, 1, 1]), L2)
    Lattice morphism defined by the left action of the matrix
    [1 1 1 1]
    Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
    Codomain: Ambient lattice of rank 1 with the trivial action of a group of order 4
"""


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
        """
        Initialize the morphism object and checks that the map is equivariant.

        EXAMPLES::

            sage: L1 = GLattice([2,2])
            sage: L2 = GLattice([2,2], 4)
            sage: L3 = GLattice([2],4)
            sage: m = identity_matrix(4)
            sage: L1.left_morphism(m)
            Lattice endomorphism defined by the left action of the matrix
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4

            sage: L1.left_morphism(m, L2)
            Traceback (most recent call last):
            ...
            TypeError: The morphism does not preserve the action of the group
            sage: L2.left_morphism(m, L3)
            Traceback (most recent call last):
            ...
            ValueError: The two lattices must be acted on by the same group
        """
        D = domain
        C = codomain
        if not (D.group() == C.group()):
            raise ValueError("The two lattices must be acted on by the same group")
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
        """
        The print representation of a lattice morphism.

        EXAMPLES::

            sage: L = GLattice([3])
            sage: h = L.left_morphism(identity_matrix(3)); h
            Lattice endomorphism defined by the left action of the matrix
            [1 0 0]
            [0 1 0]
            [0 0 1]
            Domain: Ambient lattice of rank 3 with a faithful action by a group of order 3
            sage: q = L1.quotient_lattice(L1.zero_sum_sublattice(), True, True); q[1]
            Lattice morphism defined by the left action of the matrix
            [1 1 1]
            Domain: Ambient lattice of rank 3 with a faithful action by a group of order 3
            Codomain: Ambient lattice of rank 1 with the trivial action of a group of order 3
        """
        if self.domain() == self.codomain():
            r = "Lattice endomorphism defined by the left action of the matrix\n{!r}\nDomain: {}"
            return r.format(self.matrix(), self.domain())
        else:
            r = "Lattice morphism defined by the left action of the matrix\n{!r}\nDomain: {}\nCodomain: {}"
            return r.format(self.matrix(), self.domain(), self.codomain())

    def __call__(self, elt):
        """
        Compute the image of either an element of a lattice, or a sublattice.
       

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: L = GLattice(G,3)
            sage: m = matrix(3, [2,4,2,1,2,1,3,6,3])
            sage: h = L.left_morphism(m); h
            Lattice endomorphism defined by the left action of the matrix
            [2 4 2]
            [1 2 1]
            [3 6 3]
            Domain: Ambient lattice of rank 3 with the trivial action of a group of order 6
            sage: SL = L.sublattice([sum(L.basis())])
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

    def __add__(self, other,):
        """
        Adds morphisms with same domain and codomain. The second argument can be a list of morphisms. 
        
        INPUT:

        - ``other`` -- Lattice or list of lattice we want to take the sum with.

        EXAMPLES::

            sage: L = GLattice([2, 2])
            sage: m = identity_matrix(4)
            sage: h = L.left_morphism(m)
            sage: h2 = h + h; h2
            Lattice endomorphism defined by the left action of the matrix
            [2 0 0 0]
            [0 2 0 0]
            [0 0 2 0]
            [0 0 0 2]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: h2.matrix() == L.surjection_from_square().pre_compose(L.diagonal_embedding()).matrix()
            True
            sage: h2 = h + [h, h]; h2
            Lattice endomorphism defined by the left action of the matrix
            [3 0 0 0]
            [0 3 0 0]
            [0 0 3 0]
            [0 0 0 3]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
       
        ::

            sage: G = SymmetricGroup(3)
            sage: L = GLattice([]).norm_one_restriction_of_scalars(G)
            sage: L
            Ambient lattice of rank 5 with a faithful action by a group of order 6
            sage: [P, h] = L.permutation_cover(); h
            Lattice morphism defined by the left action of the matrix
            [ 1| 0| 0| 0|-1| 0]
            [ 0| 0| 0| 1|-1| 0]
            [ 0| 0| 0| 0|-1| 1]
            [ 0| 0| 1| 0|-1| 0]
            [ 0| 1| 0| 0|-1| 0]
            Domain: Ambient lattice of rank 6 with a faithful action by a group of order 6
            Codomain: Ambient lattice of rank 5 with a faithful action by a group of order 6
            sage: h + h
            Lattice morphism defined by the left action of the matrix
            [ 2  0  0  0 -2  0]
            [ 0  0  0  2 -2  0]
            [ 0  0  0  0 -2  2]
            [ 0  0  2  0 -2  0]
            [ 0  2  0  0 -2  0]
            Domain: Ambient lattice of rank 6 with a faithful action by a group of order 6
            Codomain: Ambient lattice of rank 5 with a faithful action by a group of order 6
     
        ::

            sage: L = GLattice([3], 2)
            sage: m1 = identity_matrix(2)
            sage: m2 = matrix(2, [1, 2, 3, 4])
            sage: h1 = L.left_morphism(m1); h1
            Lattice endomorphism defined by the left action of the matrix
            [1 0]
            [0 1]
            Domain: Ambient lattice of rank 2 with the trivial action of a group of order 3
            sage: h2 = L.left_morphism(m2); h2
            Lattice endomorphism defined by the left action of the matrix
            [1 2]
            [3 4]
            Domain: Ambient lattice of rank 2 with the trivial action of a group of order 3
            sage: h1 + h2
            Lattice endomorphism defined by the left action of the matrix
            [2 2]
            [3 5]
            Domain: Ambient lattice of rank 2 with the trivial action of a group of order 3
        """
        return self.sum(other, "inner", "inner")

    def sum(self, other, domainsum = "outer", codomainsum = "outer"):
        """
        Adds morphisms. The second argument can be a list of morphisms. 
        
        INPUT:

        - ``other`` -- Lattice or list of lattice we want to take the sum with.

        - ``domainsum`` -- String (default ``outer``), declares it we want to take the inner or 
          outer sum for the domain. If ``inner``, the domain lattices must match. If the argument
          is ``outer`` then the domain of the sum will be a direct sum of lattices.

        - ``codomainsum`` -- String (default ``outer``), declares it we want to take the inner or 
          outer sum for the codomain. If ``inner``, the codomain lattices must match. If the argument
          is ``outer`` then the codomain of the sum will be a direct sum of lattices.



        EXAMPLES::

            sage: L = GLattice([2, 2])
            sage: m = identity_matrix(4)
            sage: h = L.left_morphism(m)
            sage: h.sum(h, "inner", "inner")
            Lattice endomorphism defined by the left action of the matrix
            [2 0 0 0]
            [0 2 0 0]
            [0 0 2 0]
            [0 0 0 2]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: h.sum(h)
            Lattice endomorphism defined by the left action of the matrix
            [1 0 0 0|0 0 0 0]
            [0 1 0 0|0 0 0 0]
            [0 0 1 0|0 0 0 0]
            [0 0 0 1|0 0 0 0]
            [-------+-------]
            [0 0 0 0|1 0 0 0]
            [0 0 0 0|0 1 0 0]
            [0 0 0 0|0 0 1 0]
            [0 0 0 0|0 0 0 1]
            Domain: Ambient lattice of rank 8 with a faithful action by a group of order 4
            sage: h.sum(h, "inner")
            Lattice morphism defined by the left action of the matrix
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            [-------]
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            Codomain: Ambient lattice of rank 8 with a faithful action by a group of order 4
            sage: h.sum(h, "outer", "inner")
            Lattice morphism defined by the left action of the matrix
            [1 0 0 0|1 0 0 0]
            [0 1 0 0|0 1 0 0]
            [0 0 1 0|0 0 1 0]
            [0 0 0 1|0 0 0 1]
            Domain: Ambient lattice of rank 8 with a faithful action by a group of order 4
            Codomain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: h.sum([h, h], "inner", "inner")
            Lattice endomorphism defined by the left action of the matrix
            [3 0 0 0]
            [0 3 0 0]
            [0 0 3 0]
            [0 0 0 3]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: h.sum([h, h], "inner", "outer")
            Lattice morphism defined by the left action of the matrix
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            [-------]
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            Codomain: Ambient lattice of rank 12 with a faithful action by a group of order 4
            sage: h.sum([h, h], "outer", "inner")
            Lattice morphism defined by the left action of the matrix
            [1 0 0 0 1 0 0 0|1 0 0 0]
            [0 1 0 0 0 1 0 0|0 1 0 0]
            [0 0 1 0 0 0 1 0|0 0 1 0]
            [0 0 0 1 0 0 0 1|0 0 0 1]
            Domain: Ambient lattice of rank 12 with a faithful action by a group of order 4
            Codomain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: h.sum([h, h], "outer", "outer")
            Lattice endomorphism defined by the left action of the matrix
            [1 0 0 0 0 0 0 0|0 0 0 0]
            [0 1 0 0 0 0 0 0|0 0 0 0]
            [0 0 1 0 0 0 0 0|0 0 0 0]
            [0 0 0 1 0 0 0 0|0 0 0 0]
            [0 0 0 0 1 0 0 0|0 0 0 0]
            [0 0 0 0 0 1 0 0|0 0 0 0]
            [0 0 0 0 0 0 1 0|0 0 0 0]
            [0 0 0 0 0 0 0 1|0 0 0 0]
            [---------------+-------]
            [0 0 0 0 0 0 0 0|1 0 0 0]
            [0 0 0 0 0 0 0 0|0 1 0 0]
            [0 0 0 0 0 0 0 0|0 0 1 0]
            [0 0 0 0 0 0 0 0|0 0 0 1]
            Domain: Ambient lattice of rank 12 with a faithful action by a group of order 4
        """
        if isinstance(other, list):
            if len(other) == 0:
                return self
            else:
                mor = self.sum(other[0], domainsum, codomainsum)
                return mor.sum(other[1:], domainsum, codomainsum)
        else:
            A = self.domain()
            B = other.domain()
            TA = self.codomain()
            TB = other.codomain()
            C = A.direct_sum(B)
            TC = TA.direct_sum(TB)
            from sage.matrix.special import block_diagonal_matrix
            M = block_diagonal_matrix([self.matrix(),other.matrix()])
            mor = C.left_morphism(M, TC)
            if domainsum == "inner":
                if not(A == B):
                    raise ValueError("The domains must match.")
                mor = mor.pre_compose(A.diagonal_embedding())
            if codomainsum == "inner":
                if not(TA == TB):
                    raise ValueError("The codomains of the morphisms must match")
                surj = TA.surjection_from_square()
                mor = surj.pre_compose(mor)
            return mor


    def group_action(self, element, side="Left"):
        """ 
        Compute the action of the group on the homomorphism. If the element is not central, the 
        result might not be a lattice morphism.

        INPUT:

        - ``element`` -- group element acting on the morphism.

        - ``side`` -- string (default ``Left``). Can be assigned either ``Right`` or ``Left``. 
          Side on which the element should be acting.

        EXAMPLES::

            sage: L = GLattice([2, 2])
            sage: h = L.left_morphism(identity_matrix(4))
            sage: h
            Lattice endomorphism defined by the left action of the matrix
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: g = L.group()[2]
            sage: g
            (1,2)
            sage: h.group_action(g)
            Lattice endomorphism defined by the left action of the matrix
            [0 1 0 0]
            [1 0 0 0]
            [0 0 1 0]
            [0 0 0 1]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4

        ::

            sage: G =QuaternionGroup()
            sage: L2 = GLattice(1).norm_one_restriction_of_scalars(G)
            sage: L2
            Ambient lattice of rank 7 with a faithful action by a group of order 8
            sage: [P, h] = L2.permutation_cover(); h
            Lattice morphism defined by the left action of the matrix
            [ 1| 0| 0| 0| 0| 0| 0|-1]
            [ 0| 0| 0| 1| 0| 0| 0|-1]
            [ 0| 1| 0| 0| 0| 0| 0|-1]
            [ 0| 0| 1| 0| 0| 0| 0|-1]
            [ 0| 0| 0| 0| 1| 0| 0|-1]
            [ 0| 0| 0| 0| 0| 0| 1|-1]
            [ 0| 0| 0| 0| 0| 1| 0|-1]
            Domain: Ambient lattice of rank 8 with a faithful action by a group of order 8
            Codomain: Ambient lattice of rank 7 with a faithful action by a group of order 8
            sage: g = G.center()[1]
            sage: h.group_action(g)
            Lattice morphism defined by the left action of the matrix
            [ 0  1  0  0  0  0 -1  0]
            [ 0  0  1  0  0  0 -1  0]
            [ 1  0  0  0  0  0 -1  0]
            [ 0  0  0  1  0  0 -1  0]
            [ 0  0  0  0  0  1 -1  0]
            [ 0  0  0  0  0  0 -1  1]
            [ 0  0  0  0  1  0 -1  0]
            Domain: Ambient lattice of rank 8 with a faithful action by a group of order 8
            Codomain: Ambient lattice of rank 7 with a faithful action by a group of order 8

        ::

            sage: G = GLattice([4]).group()
            sage: L = GLattice(1).induced_lattice(G)
            sage: h = L.left_morphism(identity_matrix(4))
            sage: g = G[1]
            sage: h.group_action(g)
            Lattice endomorphism defined by the left action of the matrix
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            [1 0 0 0]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: h.group_action(g, "Right")
            Lattice endomorphism defined by the left action of the matrix
            [0 0 0 1]
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
        """
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
        """
        Return the domain of a morphism

        EXAMPLES::

            sage: G = GLattice([4]).group()
            sage: GM = GLattice(1)
            sage: IL = GM.norm_one_restriction_of_scalars(G)
            sage: IL
            Ambient lattice of rank 3 with a faithful action by a group of order 4
            sage: [P, h] = IL.permutation_cover()
            sage: h.domain()
            Ambient lattice of rank 4 with a faithful action by a group of order 4

        ::

            sage: L = GLattice([4, 2])
            sage: h = L.left_morphism(identity_matrix(6))
            sage: h.domain() == h.codomain()
            True
        """
        return self._domain

    def codomain(self):
        """ 
        Return the codomain of a morphism.

        EXAMPLES::

            sage: L = GLattice([2])
            sage: h = L.diagonal_embedding()
            sage: h.codomain()
            Ambient lattice of rank 4 with a faithful action by a group of order 2
        """
        return self._codomain

    def matrix(self):
        """
        The matrix corresponding to the morphism.

        EXAMPLES::

            sage: L = GLattice([2])
            sage: h = L.surjection_from_square()
            sage: h.matrix()
            [1 0|1 0]
            [0 1|0 1]
            sage: hh = h + h; hh.matrix()
            [2 0 2 0]
            [0 2 0 2]
        """
        return self._matrix

    def free_module_morphism(self):
        """
        Return the corresponding free module morphisms. Those morphisms are acting on the right.

        EXAMPLES::

            sage: L = GLattice([2])
            sage: h = L.surjection_from_square(); h
            Lattice morphism defined by the left action of the matrix
            [1 0|1 0]
            [0 1|0 1]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 2
            Codomain: Ambient lattice of rank 2 with a faithful action by a group of order 2
            sage: h.free_module_morphism()
            Free module morphism defined by the matrix
            [1 0]
            [0 1]
            [1 0]
            [0 1]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 2
            Codomain: Ambient lattice of rank 2 with a faithful action by a group of order 2
        """
        return self._free_morphism

    def kernel(self):
        """
        Returns the Kernel of a morphism

        EXAMPLES::

            sage: L = GLattice([2])
            sage: h = L.surjection_from_square(); h
            Lattice morphism defined by the left action of the matrix
            [1 0|1 0]
            [0 1|0 1]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 2
            Codomain: Ambient lattice of rank 2 with a faithful action by a group of order 2
            sage: h.kernel()
            Sublattice of degree 4 and rank 2 with a faithful action by a group of order 2 and echelon basis matrix
            [ 1  0 -1  0]
            [ 0  1  0 -1]
            sage: L.norm().kernel()
            Sublattice of degree 2 and rank 1 with a faithful action by a group of order 2 and echelon basis matrix
            [ 1 -1]
        """
        K = self.matrix().transpose().kernel()
        return self.domain().sublattice(K.basis())

    def image(self):
        """
        Return the image of a morphism.

        EXAMPLES::

            sage: L = GLattice([2])
            sage: h = L.diagonal_embedding()
            sage: h.image()
            Sublattice of degree 4 and rank 2 with a faithful action by a group of order 2 and echelon basis matrix
            [1 0 1 0]
            [0 1 0 1]
        """
        I = self.matrix().transpose().image()
        return self.codomain().sublattice(I.basis())

    def dual(self):
        """
        Return the dual of a morphism.
        
        EXAMPLES::

            sage: L = GLattice([2])
            sage: h = L.diagonal_embedding()
            sage: h.image()
            Sublattice of degree 4 and rank 2 with a faithful action by a group of order 2 and echelon basis matrix
            [1 0 1 0]
            [0 1 0 1]
            sage: h.dual()
            Lattice morphism defined by the left action of the matrix
            [1 0|1 0]
            [0 1|0 1]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 2
            Codomain: Ambient lattice of rank 2 with a faithful action by a group of order 2
            sage: L.dual().surjection_from_square()
            Lattice morphism defined by the left action of the matrix
            [1 0|1 0]
            [0 1|0 1]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 2
            Codomain: Ambient lattice of rank 2 with a faithful action by a group of order 2
        """
        return self.codomain().colattice().left_morphism(self.matrix().transpose(), self.domain().colattice())

    def is_identity(self):
        """
        Determine if a morphism is the identity morphism.

        EXAMPLES::
        sage: L = GLattice([2], 2)
        sage: h = L.left_morphism(identity_matrix(2))
        sage: h2 = L.left_morphism(matrix(2, [1, 2, 3, 4]))
        sage: [P, h3] = L.permutation_cover()
        sage: h.is_identity()
        True
        sage: h2.is_identity()
        False
        sage: h3.is_identity()
        True
        """
        return self.free_module_morphism().is_identity()    

    def characteristic_polynomial(self):
        """
        Return the characteristic polynomial of a morphism.

        EXAMPLES::

            sage: L = GLattice([4])
            sage: h = L.identity_morphism(); h
            Lattice endomorphism defined by the left action of the matrix
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: h2 = h.group_action(L.group()[1]); h2
            Lattice endomorphism defined by the left action of the matrix
            [0 0 0 1]
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: h.characteristic_polynomial()
            x^4 - 4*x^3 + 6*x^2 - 4*x + 1
            sage: h.characteristic_polynomial().factor()
            (x - 1)^4
            sage: h2.characteristic_polynomial()
            x^4 - 1
            sage: h2.characteristic_polynomial().factor()
            (x - 1) * (x + 1) * (x^2 + 1)
        """
        return self.matrix().characteristic_polynomial()

    def det(self):
        """
        Return the determinant of a morphism.

        EXAMPLES::

            sage: L = GLattice([4])
            sage: h = L.identity_morphism()
            sage: h2 = h.group_action(L.group()[1])
            sage: h.det()
            1
            sage: h2.det()
            -1
            sage: (h+h).det()
            16
        """
        return self.matrix().det()

    def inverse_image(self, M):
        """
        Compute the invese image of a sublattice of the codomain.

        EXAMPLES::

            sage: L = GLattice([4])
            sage: h = L.identity_morphism()
            sage: SL = L.zero_sum_sublattice()
            sage: h.inverse_image(SL)
            Sublattice of degree 4 and rank 3 with a faithful action by a group of order 4 and echelon basis matrix
            [ 1  0  0 -1]
            [ 0  1  0 -1]
            [ 0  0  1 -1]
            sage: SL2 = L.fixed_sublattice()
            sage: h.inverse_image(SL2)
            Sublattice of degree 4 and rank 1 with a faithful action by a group of order 4 and echelon basis matrix
            [1 1 1 1]


        ::

            sage: L = GLattice([2, 2])
            sage: [a, b, c, d] = L.basis()
            sage: SL = L.sublattice([a+b, c+d])
            sage: n = L.norm(); n
            Lattice endomorphism defined by the left action of the matrix
            [2 2 0 0]
            [2 2 0 0]
            [0 0 2 2]
            [0 0 2 2]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: n.inverse_image(SL)
            Sublattice of degree 4 and rank 4 with a faithful action by a group of order 4 and echelon basis matrix
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            sage: n.inverse_image(L.sublattice([a+b+c+d]))
            Sublattice of degree 4 and rank 3 with a faithful action by a group of order 4 and echelon basis matrix
            [ 1  0  0  1]
            [ 0  1  0  1]
            [ 0  0  1 -1]
        """
        I = self.free_module_morphism().inverse_image(M)
        return self.domain().sublattice(I.basis())

    def is_bijective(self):
        """
        Determine if a morphism is bijective.

        EXAMPLES::

            sage: L = GLattice([2, 2])
            sage: h = L.identity_morphism()
            sage: h.is_bijective()
            True
            sage: L.norm().is_bijective()
            False
            sage: L.norm().kernel()
            Sublattice of degree 4 and rank 2 with a faithful action by a group of order 4 and echelon basis matrix
            [ 1 -1  0  0]
            [ 0  0  1 -1]
        """
        return self.free_module_morphism().is_bijective()

    def is_endomorphism(self):
        """
        Return whether a morphism is an endomorphism.

        EXAMPLES::

            sage: L = GLattice([2, 2])
            sage: L.identity_morphism().is_endomorphism()
            True
            sage: L.norm().is_endomorphism()
            True
            sage: L.diagonal_embedding().is_endomorphism()
            False
        """
        return self.domain() == self.codomain()

    def is_injective(self):
        """
        Determine whether a map is injective.

        EXAMPLES::

            sage: L = GLattice([2, 2])
            sage: L.diagonal_embedding().is_injective()
            True
            sage: L.diagonal_embedding().is_bijective()
            False
            sage: L.surjection_from_square().is_injective()
            False
            sage: L.norm().is_injective()
            False
        """
        return self.free_module_morphism().is_injective()

    def is_surjective(self):
        """
        Determine whether the map is surjective.

        EXAMPLES::

            sage: L = GLattice([2, 2])
            sage: L.surjection_from_square().is_surjective()
            True
            sage: L.surjection_from_square().is_bijective()
            False
            sage: L.diagonal_embedding().is_surjective()
            False
            sage: L.norm().is_surjective()
            False
        """
        return self.free_module_morphism().is_surjective()

    def is_idempotent(self):
        """
        Determine whether the map is idempotent.

        EXAMPLES::

            sage: L = GLattice([2, 2])
            sage: h = L.left_morphism(matrix(4, [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
            sage: h.is_identity()
            False
            sage: h.is_idempotent()
            True
        """
        return self.free_module_morphism().is_idempotent()

    def minpoly(self):
        """
        Determine the minimal polynomial of the morphism.

        EXAMPLES::

            sage: L = GLattice([2, 2])
            sage: h = L.left_morphism(matrix(4, [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
            sage: h.is_identity()
            False
            sage: h.is_idempotent()
            True
            sage: h.characteristic_polynomial()
            x^4 - 2*x^3 + x^2
            sage: h.minpoly()
            x^2 - x

        ::

            sage: G = QuaternionGroup()
            sage: L = GLattice(1).norm_one_restriction_of_scalars(G)
            sage: Z = G.center()
            sage: h = L.identity_morphism()
            sage: h2 = h.group_action(Z[1])
            sage: h2.minpoly()
            x^2 - 1
            sage: h2.characteristic_polynomial()
            x^7 + x^6 - 3*x^5 - 3*x^4 + 3*x^3 + 3*x^2 - x - 1
        """
        return self.free_module_morphism().minpoly() 

    def trace(self):
        """
        Compute the trace of a morphism.

        EXAMPLES::

            sage: L = GLattice([2, 2])
            sage: L.identity_morphism().trace()
            4
            sage: L.norm().trace()
            8
        """
        return self.matrix().trace()

    def pre_compose(self, morphism):
        """
        Precompose the morphism with another morphism.

        EXAMPLES::

            sage: L = GLattice([2])
            sage: h1 = L.diagonal_embedding()
            sage: h2 = L.surjection_from_square()
            sage: h3 = L.norm()
            sage: h2.pre_compose(h1)
            Lattice endomorphism defined by the left action of the matrix
            [2 0]
            [0 2]
            Domain: Ambient lattice of rank 2 with a faithful action by a group of order 2
            sage: h2.pre_compose(h1.pre_compose(h3))
            Lattice endomorphism defined by the left action of the matrix
            [2 2]
            [2 2]
            Domain: Ambient lattice of rank 2 with a faithful action by a group of order 2

        The following example builds an injection of a sublattice, and the corresponding
        quotient map. It checks that the composite is trivial, and that it forms an exact sequence.
        ::

            sage: L = GLattice([2, 2])
            sage: [a, b, c, d] = L.basis()
            sage: SL = L.sublattice([a, b])
            sage: i = SL.injection_morphism()
            sage: [Q, s] = L.quotient_lattice(SL, True, True)
            sage: s.pre_compose(i)
            Lattice endomorphism defined by the left action of the matrix
            [0 0]
            [0 0]
            Domain: Ambient lattice of rank 2 with an action by a group of order 4
            sage: s.kernel() == i.image()
            True

        """
        M = self.matrix()*morphism.matrix()
        return morphism.domain().left_morphism(M, self.codomain())

    def group(self):
        """
        Return the group acting on both domain and codomain.

        EXAMPLES::

            sage: L = GLattice([2])
            sage: h = L.identity_morphism()
            sage: h.group()
            Permutation Group with generators [(1,2)]
        """
        return self.domain().group()

    def cokernel(self):
        """
        Compute the cokernel of a morphism.

        EXAMPLES::

            sage: L = GLattice([2])
            sage: h = L.diagonal_embedding()
            sage: h.cokernel()
            Ambient lattice of rank 2 with a faithful action by a group of order 2
            sage: L.norm().cokernel()
            Ambient lattice of rank 1 with a faithful action by a group of order 2
            sage: L.identity_morphism().cokernel()
            Ambient lattice of rank 0 with an action by a group of order 2
        """
        I = self.image()
        return self.codomain().quotient_lattice(I)
