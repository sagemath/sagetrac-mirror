
import sage.modules.matrix_morphism as matrix_morphism
import sage.modules.free_module_morphism as free_module_morphism
from . import glattice_homspace
from sage.structure.element import is_Matrix


def is_GLatticeMorphism(x):

    return isinstance(x, GLatticeMorphism)


class GLatticeMorphism(free_module_morphism.FreeModuleMorphism):

    def __init__(self, homspace, A, side="left"):

        if not glattice_homspace.is_GLatticeHomspace(homspace):
            raise TypeError('homspace must be a GLattice space hom space, not {0}'.format(homspace))
        if isinstance(A, matrix_morphism.MatrixMorphism):
            A = A.matrix()
        if not is_Matrix(A):
            msg = 'input must be a matrix representation or another matrix morphism, not {0}'
            raise TypeError(msg.format(A))
        # now have a homspace, and a matrix, check compatibility and equivariance
        D = homspace.domain()
        C = homspace.codomain()
        if not (D.group() == C.group()):
            raise ValueError("The two lattices must be acted on by the same group")
        if side == "left":
            if D.dimension() != A.nrows():
                raise TypeError('domain rank is incompatible with matrix size')
            if C.dimension() != A.ncols():
                raise TypeError('codomain rank is incompatible with matrix size')
            for e in D.basis():
                for g in D.group().gens():
                    if ((D.action_matrix(g)*e))*A != C.action_matrix(g)*(e*A):
                        raise TypeError("The morphism does not preserve the action of the group")

        if side == "right":
            if C.dimension() != A.nrows():
                raise TypeError('Domain rank is incompatible with matrix size')
            if D.dimension() != A.ncols():
                raise TypeError('codomain rank is incompatible with matrix size')
            for e in D.basis():
                for g in D.group().gens():
                    if A*((D.action_matrix(g)*e)) != C.action_matrix(g)*(A*e):
                        raise TypeError("The morphism does not preserve the action of the group")
        # now check 

        A = homspace._matrix_space(side)(A)
        free_module_morphism.FreeModuleMorphism.__init__(self, homspace, A, side)

    def is_invertible(self):

        m = self.matrix()
        if not m.is_square():
            return False
        return m.rank() == m.ncols()


    def _repr_(self):

        m = self.matrix()
        act = ""
        if self.side() == "right":
            act = "as left-multiplication "
        msg = ("Lattice morphism represented {}by the matrix:\n",
               "{!r}\n",
               "Domain: {}\n",
               "Codomain: {}")
        return ''.join(msg).format(act, m, self.domain(), self.codomain())


    def exterior_sum(self, other, domain_ext = True, codomain_ext = True):
        """
        The exterior sum of morphisms. The second argument can be a list of morphisms. 
        
        INPUT:

        - ``other`` -- Lattice or list of lattice we want to take the sum with.

        - ``domain_ext`` -- Boolean (default ``True``). If True, the domain of the sum 
          will be the direct sum of domains. Else morphisms must share a common domain.

        - ``codomain_ext`` -- Boolean (default ``True``). If True, the codomain of the sum 
          will be the direct sum of codomains. Else morphisms must share a common codomain.



        EXAMPLES::

            sage: L = GLattice([2, 2])
            sage: m = identity_matrix(4)
            sage: h = L.left_morphism(m)
            sage: h.exterior_sum(h, False, False)
            Lattice endomorphism defined by the left action of the matrix
            [2 0 0 0]
            [0 2 0 0]
            [0 0 2 0]
            [0 0 0 2]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: h.exterior_sum(h)
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
            sage: h.exterior_sum(h, domain_ext = False)
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
            sage: h.exterior_sum(h, codomain_ext = False)
            Lattice morphism defined by the left action of the matrix
            [1 0 0 0|1 0 0 0]
            [0 1 0 0|0 1 0 0]
            [0 0 1 0|0 0 1 0]
            [0 0 0 1|0 0 0 1]
            Domain: Ambient lattice of rank 8 with a faithful action by a group of order 4
            Codomain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: h.exterior_sum([h, h], False, False)
            Lattice endomorphism defined by the left action of the matrix
            [3 0 0 0]
            [0 3 0 0]
            [0 0 3 0]
            [0 0 0 3]
            Domain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: h.exterior_sum([h, h], domain_ext = False)
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
            sage: h.exterior_sum([h, h], codomain_ext = False)
            Lattice morphism defined by the left action of the matrix
            [1 0 0 0 1 0 0 0|1 0 0 0]
            [0 1 0 0 0 1 0 0|0 1 0 0]
            [0 0 1 0 0 0 1 0|0 0 1 0]
            [0 0 0 1 0 0 0 1|0 0 0 1]
            Domain: Ambient lattice of rank 12 with a faithful action by a group of order 4
            Codomain: Ambient lattice of rank 4 with a faithful action by a group of order 4
            sage: h.exterior_sum([h, h])
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
                mor = self.exterior_sum(other[0], domain_ext, codomain_ext)
                return mor.exterior_sum(other[1:], domain_ext, codomain_ext)
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
            if not(domain_ext):
                if not(A == B):
                    raise ValueError("The domains must match.")
                mor = mor*A.diagonal_embedding()
            if not(codomain_ext):
                if not(TA == TB):
                    raise ValueError("The codomains of the morphisms must match")
                surj = TA.surjection_from_square()
                mor = surj*mor
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





    def kernel(self):
        """
        Return the Kernel of a morphism

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
        if side=="right":
            K = self.matrix().transpose().kernel()
        else:
            K = self.matrix().kernel()
        return self.domain().sublattice(K.basis(), check=False)

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
        if side=="right":
            I = self.matrix().transpose().image()
        else:
            return self.codomain().sublattice(I.basis(), check=False)

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
        return self.codomain().colattice().hom(self.matrix().transpose(), self.domain().colattice(), side=self.side())

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
        from sage.modules.free_module_morphism import FreeModuleMorphism   
        I = FreeModuleMorphism.inverse_image(self, M)
        return self.domain().sublattice(I.basis(), check=False)


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
