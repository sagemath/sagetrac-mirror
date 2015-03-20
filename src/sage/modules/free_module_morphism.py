"""
Morphisms of free modules

AUTHORS:
    - William Stein: initial version

    - Miguel Marco (2010-06-19): added eigenvalues, eigenvectors and minpoly functions


TESTS::

    sage: V = ZZ^2; f = V.hom([V.1,-2*V.0])
    sage: loads(dumps(f))
    Free module morphism defined by the matrix
    [ 0  1]
    [-2  0]
    Domain: Ambient free module of rank 2 over the principal ideal domain ...
    Codomain: Ambient free module of rank 2 over the principal ideal domain ...
    sage: loads(dumps(f)) == f
    True
"""

####################################################################################
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
####################################################################################

# A matrix morphism is a morphism that is defined by multiplication by a
# matrix.  Elements of domain must either have a method "vector()" that
# returns a vector that the defining matrix can hit from the left, or
# be coercible into vector space of appropriate dimension.

import sage.misc.misc as misc
import sage.modules.free_module as free_module
import matrix_morphism
from sage.structure.sequence import Sequence

import free_module_homspace

def is_FreeModuleMorphism(x):
    """
    EXAMPLES::

        sage: V = ZZ^2; f = V.hom([V.1,-2*V.0])
        sage: sage.modules.free_module_morphism.is_FreeModuleMorphism(f)
        True
        sage: sage.modules.free_module_morphism.is_FreeModuleMorphism(0)
        False
    """
    return isinstance(x, FreeModuleMorphism)

class FreeModuleMorphism(matrix_morphism.MatrixMorphism):
    def __init__(self, parent, A):
        """
        INPUT:

            -  ``parent`` - a homspace in a (sub) category of free modules

            -  ``A`` - matrix

        EXAMPLES::

            sage: V = ZZ^3; W = span([[1,2,3],[-1,2,8]], ZZ)
            sage: phi = V.hom(matrix(ZZ,3,[1..9]))
            sage: type(phi)
            <class 'sage.modules.free_module_morphism.FreeModuleMorphism'>
        """
        if not free_module_homspace.is_FreeModuleHomspace(parent):
            raise TypeError("parent (=%s) must be a free module hom space"%parent)
        if isinstance(A, matrix_morphism.MatrixMorphism):
            A = A.matrix()
        A = parent._matrix_space()(A)
        matrix_morphism.MatrixMorphism.__init__(self, parent, A)

    def pushforward(self, x):
        """
        Compute the image of a sub-module of the domain.

        EXAMPLES::

            sage: V = QQ^3; W = span([[1,2,3],[-1,2,5/3]], QQ)
            sage: phi = V.hom(matrix(QQ,3,[1..9]))
            sage: phi.rank()
            2
            sage: phi(V)   #indirect doctest
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]

        We compute the image of a submodule of a ZZ-module embedded in
        a rational vector space::

            sage: V = QQ^3; W = V.span_of_basis([[2,2,3],[-1,2,5/3]], ZZ)
            sage: phi = W.hom([W.0, W.0-W.1]); phi
            Free module morphism defined by the matrix
            [ 1  0]
            [ 1 -1]...
            sage: phi(span([2*W.1],ZZ))
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [  6   0 8/3]
            sage: phi(2*W.1)
            (6, 0, 8/3)
        """
        if free_module.is_FreeModule(x):
            V = self.domain().submodule(x)
            return self.restrict_domain(V).image()
        raise TypeError("`pushforward` is only defined for submodules")

    def _repr_(self):
        r"""
        Return string representation of this morphism of free modules.

        EXAMPLES::

            sage: V = ZZ^3; W = span([[1,2,3],[-1,2,8]], ZZ)
            sage: phi = V.hom(matrix(ZZ,3,[1..9]))
            sage: phi._repr_()
            'Free module morphism defined by the matrix\n[1 2 3]\n[4 5 6]\n[7 8 9]\nDomain: Ambient free module of rank 3 over the principal ideal domain Integer Ring\nCodomain: Ambient free module of rank 3 over the principal ideal domain Integer Ring'

            sage: V = ZZ^6
            sage: W = ZZ^4
            sage: m = matrix(QQ, [[1, 0, 0 ,0], [0]*4, [0]*4, [0]*4, [0]*4, [0]*4])
            sage: phi = V.hom(m, W)
            sage: rho = phi.restrict_codomain(W.span([W.0]))
            sage: rho
            Free module morphism defined by the matrix
            [1]
            [0]
            [0]
            [0]
            [0]
            [0]
            Domain: Ambient free module of rank 6 over the principal ideal domain Integer Ring
            Codomain: Free module of degree 4 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1 0 0 0]

            sage: V = QQ^40
            sage: m = matrix(QQ, 40, 40, 1600)
            sage: phi = V.hom(m, V)
            sage: phi
            Vector space morphism represented by the matrix:
            40 x 40 dense matrix over Rational Field
            Domain: Vector space of dimension 40 over Rational Field
            Codomain: Vector space of dimension 40 over Rational Field
        """
        r = "Free module morphism defined by the matrix\n{0}\nDomain: {1}\nCodomain: {2}"
        return r.format(self.matrix(), self.domain(), self.codomain())

    def change_ring(self, R):
        """
        Change the ring over which this morphism is defined.  This changes the ring of the
        domain, codomain, and underlying matrix.

        EXAMPLES::

            sage: V0 = span([[0,0,1],[0,2,0]],ZZ); V1 = span([[1/2,0],[0,2]],ZZ); W = span([[1,0],[0,6]],ZZ)
            sage: h = V0.hom([-3*V1.0-3*V1.1, -3*V1.0-3*V1.1])
            sage: h.base_ring()
            Integer Ring
            sage: h
            Free module morphism defined by the matrix
            [-3 -3]
            [-3 -3]...
            sage: h.change_ring(QQ).base_ring()
            Rational Field
            sage: f = h.change_ring(QQ); f
            Vector space morphism represented by the matrix:
            [-3 -3]
            [-3 -3]
            Domain: Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [0 1 0]
            [0 0 1]
            Codomain: Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
            sage: f = h.change_ring(GF(7)); f
            Vector space morphism represented by the matrix:
            [4 4]
            [4 4]
            Domain: Vector space of degree 3 and dimension 2 over Finite Field of size 7
            Basis matrix:
            [0 1 0]
            [0 0 1]
            Codomain: Vector space of degree 2 and dimension 2 over Finite Field of size 7
            Basis matrix:
            [1 0]
            [0 1]
        """
        D = self.domain().change_ring(R)
        C = self.codomain().change_ring(R)
        A = self.matrix().change_ring(R)
        return D.hom(A, C)

    def inverse_image(self, V):
        """
        Given a submodule V of the codomain of self, return the
        inverse image of V under self, i.e., the biggest submodule of
        the domain of self that maps into V.

        EXAMPLES:

        We test computing inverse images over a field::

            sage: V = QQ^3; W = span([[1,2,3],[-1,2,5/3]], QQ)
            sage: phi = V.hom(matrix(QQ,3,[1..9]))
            sage: phi.rank()
            2
            sage: I = phi.inverse_image(W); I
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0    0]
            [   0    1 -1/2]
            sage: phi(I.0) in W
            True
            sage: phi(I.1) in W
            True
            sage: W = phi.image()
            sage: phi.inverse_image(W) == V
            True

        We test computing inverse images between two spaces embedded in different
        ambient spaces.::

            sage: V0 = span([[0,0,1],[0,2,0]],ZZ); V1 = span([[1/2,0],[0,2]],ZZ); W = span([[1,0],[0,6]],ZZ)
            sage: h = V0.hom([-3*V1.0-3*V1.1, -3*V1.0-3*V1.1])
            sage: h.inverse_image(W)
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [0 2 1]
            [0 0 2]
            sage: h(h.inverse_image(W)).is_submodule(W)
            True
            sage: h(h.inverse_image(W)).index_in(W)
            +Infinity
            sage: h(h.inverse_image(W))
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [ 3 12]


        We test computing inverse images over the integers::

            sage: V = QQ^3; W = V.span_of_basis([[2,2,3],[-1,2,5/3]], ZZ)
            sage: phi = W.hom([W.0, W.0-W.1])
            sage: Z = W.span([2*W.1]); Z
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [    2    -4 -10/3]
            sage: Y = phi.inverse_image(Z); Y
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [  6   0 8/3]
            sage: phi(Y) == Z
            True
        """
        if self.rank() == 0:
            # Special case -- if this is the 0 map, then the only possibility
            # for the inverse image is that it is the whole domain.
            return self.domain()

        R = self.base_ring()
        A = self.matrix()

        # Replace the module V that we are going to pullback by a
        # submodule that is contained in the image of self, since our
        # plan is to lift all generators of V.
        V = self.image().intersection(V)
        # Write V in terms of the basis for the codomain.
        V = self.codomain().coordinate_module(V)
        B = V.basis_matrix()

        # Compute the kernel, which is contained in the inverse image.
        K = self.kernel()

        if R.is_field():
            # By solving, find lifts of each of the basis elements of V.
            # Each row of C gives a linear combination of the basis for the domain
            # that maps to one of the basis elements V.
            C = A.solve_left(B)

        else:
            if not hasattr(A, 'hermite_form'):
                raise NotImplementedError("base ring (%s) must have hermite_form algorithm in order to compute inverse image"%R)

            # 1. Compute H such that U*A = H = hnf(A) without zero
            # rows. What this "does" is find a basis for the image of
            # A and explicitly represents each element in this basis
            # as the image of some element of the domain (the rows of
            # U give these elements of the domain).
            H, U = A.hermite_form(transformation=True,include_zero_rows=False)

            # 2. Next we find the unique solution to the equation
            #    Y*H = B.  This writes each basis element of V in
            # terms of our image basis found in the previous step.
            Y = H.solve_left(B)

            # 3. Multiply Y by U then takes those same linear combinations
            # from step 2 above and lifts them to coefficients that define
            # linear combinations of the basis for the domain.
            C = Y*U

        # Finally take the linear combinations of the basis for the
        # domain defined by C. Together with the kernel K, this spans
        # the inverse image of V.
        if self.domain().is_ambient():
            L = C.row_module(R)
        else:
            L = (C*self.domain().basis_matrix()).row_module(R)

        return K + L

    def lift(self, x):
        r"""
        Given an element of the image, return an element of the codomain that maps onto it.

        Note that ``lift`` and ``preimage_representative`` are
        equivalent names for this method, with the latter suggesting
        that the return value is a coset representative of the domain
        modulo the kernel of the morphism.

        EXAMPLE::

            sage: X = QQ**2
            sage: V = X.span([[2, 0], [0, 8]], ZZ)
            sage: W = (QQ**1).span([[1/12]], ZZ)
            sage: f = V.hom([W([1/3]), W([1/2])], W)
            sage: f.lift([1/3])
            (8, -16)
            sage: f.lift([1/2])
            (12, -24)
            sage: f.lift([1/6])
            (4, -8)
            sage: f.lift([1/12])
            Traceback (most recent call last):
            ...
            ValueError: element is not in the image
            sage: f.lift([1/24])
            Traceback (most recent call last):
            ...
            TypeError: element [1/24] is not in free module

        This works for vector spaces, too::

            sage: V = VectorSpace(GF(3), 2)
            sage: W = VectorSpace(GF(3), 3)
            sage: f = V.hom([W.1, W.1 - W.0])
            sage: f.lift(W.1)
            (1, 0)
            sage: f.lift(W.2)
            Traceback (most recent call last):
            ...
            ValueError: element is not in the image
            sage: w = W((17, -2, 0))
            sage: f(f.lift(w)) == w
            True

        This example illustrates the use of the ``preimage_representative``
        as an equivalent name for this method.  ::

            sage: V = ZZ^3
            sage: W = ZZ^2
            sage: w = vector(ZZ, [1,2])
            sage: f = V.hom([w, w, w], W)
            sage: f.preimage_representative(vector(ZZ, [10, 20]))
            (0, 0, 10)
        """
        from free_module_element import vector
        x = self.codomain()(x)
        A = self.matrix()
        R = self.base_ring()
        if R.is_field():
            try:
                C = A.solve_left(x)
            except ValueError:
                raise ValueError("element is not in the image")
        else:
            # see inverse_image for similar code but with comments
            if not hasattr(A, 'hermite_form'):
                raise NotImplementedError("base ring (%s) must have hermite_form algorithm in order to compute inverse image"%R)
            H, U = A.hermite_form(transformation=True,include_zero_rows=False)
            Y = H.solve_left(vector(self.codomain().coordinates(x)))
            C = Y*U
        try:
            t = self.domain().linear_combination_of_basis(C)
        except TypeError:
            raise ValueError("element is not in the image")
        assert self(t) == x
        return t

    preimage_representative = lift

    def eigenvalues(self,extend=True):
        r"""
        Returns a list with the eigenvalues of the endomorphism of vector spaces.

        INPUT:

        - ``extend`` -- boolean (default: True) decides if base field
          extensions should be considered or not.

        EXAMPLES:

        We compute the eigenvalues of an endomorphism of `\QQ^3`::

            sage: V=QQ^3
            sage: H=V.endomorphism_ring()([[1,-1,0],[-1,1,1],[0,3,1]])
            sage: H.eigenvalues()
            [3, 1, -1]

        Note the effect of the ``extend`` option::

            sage: V=QQ^2
            sage: H=V.endomorphism_ring()([[0,-1],[1,0]])
            sage: H.eigenvalues()
            [-1*I, 1*I]
            sage: H.eigenvalues(extend=False)
            []
        """
        if self.base_ring().is_field():
            if self.is_endomorphism():
                return self.matrix().eigenvalues(extend=extend)
            else:
                raise TypeError("not an endomorphism")
        else:
            raise NotImplementedError("module must be a vector space")

    def eigenvectors(self,extend=True):
        """
        Computes the subspace of eigenvectors of a given eigenvalue.

        INPUT:

        - ``extend`` -- boolean (default: True) decides if base field
          extensions should be considered or not.

        OUTPUT:

        A sequence of tuples. Each tuple contains an eigenvalue, a sequence
        with a basis of the corresponding subspace of eigenvectors, and the
        algebraic multiplicity of the eigenvalue.

        EXAMPLES::

            sage: V=(QQ^4).subspace([[0,2,1,4],[1,2,5,0],[1,1,1,1]])
            sage: H=(V.Hom(V))(matrix(QQ, [[0,1,0],[-1,0,0],[0,0,3]]))
            sage: H.eigenvectors()
            [(3, [
            (0, 0, 1, -6/7)
            ], 1), (-1*I, [
            (1, 1*I, 0, -0.571428571428572? + 2.428571428571429?*I)
            ], 1), (1*I, [
            (1, -1*I, 0, -0.571428571428572? - 2.428571428571429?*I)
            ], 1)]
            sage: H.eigenvectors(extend=False)
            [(3, [
            (0, 0, 1, -6/7)
            ], 1)]
            sage: H1=(V.Hom(V))(matrix(QQ, [[2,1,0],[0,2,0],[0,0,3]]))
            sage: H1.eigenvectors()
            [(3, [
            (0, 0, 1, -6/7)
            ], 1), (2, [
            (0, 1, 0, 17/7)
            ], 2)]
            sage: H1.eigenvectors(extend=False)
            [(3, [
            (0, 0, 1, -6/7)
            ], 1), (2, [
            (0, 1, 0, 17/7)
            ], 2)]
        """
        if self.base_ring().is_field():
            if self.is_endomorphism():
                seigenvec=self.matrix().eigenvectors_left(extend=extend)
                resu=[]
                for i in seigenvec:
                    V=self.domain().base_extend(i[0].parent())
                    svectors=Sequence(map(lambda j: V(j * V.basis_matrix()),i[1]), cr=True)
                    resu.append((i[0],svectors,i[2]))
                return resu
            else:
                raise TypeError("not an endomorphism")
        else:
            raise NotImplementedError("module must be a vector space")

    def eigenspaces(self,extend=True):
        """
        Compute a list of subspaces formed by eigenvectors of ``self``.

        INPUT:

        - ``extend`` -- (default: ``True``) determines if field
          extensions should be considered

        OUTPUT:

        - a list of pairs ``(eigenvalue, eigenspace)``

        EXAMPLES::

            sage: V = QQ^3
            sage: h = V.hom([[1,0,0],[0,0,1],[0,-1,0]], V)
            sage: h.eigenspaces()
            [(1,
              Vector space of degree 3 and dimension 1 over Rational Field
              Basis matrix:
              [1 0 0]),
             (-1*I,
              Vector space of degree 3 and dimension 1 over Algebraic Field
              Basis matrix:
              [  0   1 1*I]),
             (1*I,
              Vector space of degree 3 and dimension 1 over Algebraic Field
              Basis matrix:
              [   0    1 -1*I])]

            sage: h.eigenspaces(extend=False)
            [(1,
              Vector space of degree 3 and dimension 1 over Rational Field
              Basis matrix:
              [1 0 0])]

            sage: h = V.hom([[2,1,0], [0,2,0], [0,0,-1]], V)
            sage: h.eigenspaces()
            [(-1, Vector space of degree 3 and dimension 1 over Rational Field
              Basis matrix:
              [0 0 1]),
             (2, Vector space of degree 3 and dimension 1 over Rational Field
              Basis matrix:
              [0 1 0])]

            sage: h = V.hom([[2,1,0], [0,2,0], [0,0,2]], V)
            sage: h.eigenspaces()
            [(2, Vector space of degree 3 and dimension 2 over Rational Field
              Basis matrix:
              [0 1 0]
              [0 0 1])]
        """
        ev = self.eigenvectors(extend)
        return [(vec[0], Sequence(vec[1]).universe().subspace(vec[1]))
                for vec in ev]

    def minimal_polynomial(self,var='x'):
        r"""
        Computes the minimal polynomial.

        ``minpoly()`` and ``minimal_polynomial()`` are the same method.

        INPUT:

        - ``var`` - string (default: 'x') a variable name

        OUTPUT:

        polynomial in var - the minimal polynomial of the endomorphism.

        EXAMPLES:

        Compute the minimal polynomial, and check it. ::

            sage: V=GF(7)^3
            sage: H=V.Hom(V)([[0,1,2],[-1,0,3],[2,4,1]])
            sage: H
            Vector space morphism represented by the matrix:
            [0 1 2]
            [6 0 3]
            [2 4 1]
            Domain: Vector space of dimension 3 over Finite Field of size 7
            Codomain: Vector space of dimension 3 over Finite Field of size 7

            sage: H.minpoly()
            x^3 + 6*x^2 + 6*x + 1

            sage: H.minimal_polynomial()
            x^3 + 6*x^2 + 6*x + 1

            sage: H^3 + (H^2)*6 + H*6 + 1
            Vector space morphism represented by the matrix:
            [0 0 0]
            [0 0 0]
            [0 0 0]
            Domain: Vector space of dimension 3 over Finite Field of size 7
            Codomain: Vector space of dimension 3 over Finite Field of size 7
        """
        if self.is_endomorphism():
            return self.matrix().minpoly(var)
        else:
            raise TypeError("not an endomorphism")

    minpoly = minimal_polynomial
