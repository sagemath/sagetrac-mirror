"""
Fusion Rings
"""
# ****************************************************************************
#  Copyright (C) 2019 Daniel Bump <bump at match.stanford.edu>
#                     Nicolas Thiery <nthiery at users.sf.net>
#                     Guillermo Aboumrad <gh_willieab>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.combinat.root_system.weyl_characters import WeylCharacterRing
from functools import reduce
from operator import mul
from sage.combinat.q_analogues import q_int
from sage.matrix.special import diagonal_matrix
from sage.matrix.constructor import matrix
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc import inject_variable
from sage.rings.all import ZZ
from sage.rings.number_field.number_field import CyclotomicField
from sage.misc.cachefunc import cached_method

class FusionRing(WeylCharacterRing):
    r"""
    Return the Fusion Ring (Verlinde Algebra) of level ``k``.

    INPUT:

    - ``ct`` -- the Cartan type of a simple (finite-dimensional) Lie algebra
    - ``k`` -- a nonnegative integer

    This algebra has a basis (sometimes called *primary fields*)
    indexed by the weights of level `\leq k`. These arise as
    the fusion algebras of WZW conformal field theories, or from
    quantum groups at roots of unity. The :class:`FusionRing` class is
    implemented as a variant of the :class:`WeylCharacterRing`.

    REFERENCES:

    - [BK]_ Chapter 3
    - [DFMS1996]_ Chapter 16
    - [EGNO]_ Chapter 8
    - [Feingold2004]_
    - [Fuchs1994]_
    - [Row2006]_
    - [Walton1990]_

    EXAMPLES::

        sage: A22 = FusionRing("A2",2)
        sage: [f1, f2] = A22.fundamental_weights()
        sage: M = [A22(x) for x in [0*f1, 2*f1, 2*f2, f1+f2, f2, f1]]
        sage: [M[3] * x for x in M]
        [A22(1,1),
         A22(0,1),
         A22(1,0),
         A22(0,0) + A22(1,1),
         A22(0,1) + A22(2,0),
         A22(1,0) + A22(0,2)]

    You may assign your own labels to the basis elements. In the next
    example, we create the `SO(5)` fusion ring of level `2`, check the
    weights of the basis elements, then assign new labels to them::

        sage: B22 = FusionRing("B2", 2)
        sage: b = [B22(x) for x in B22.get_order()]; b
        [B22(0,0), B22(1,0), B22(0,1), B22(2,0), B22(1,1), B22(0,2)]
        sage: [x.weight() for x in b]
        [(0, 0), (1, 0), (1/2, 1/2), (2, 0), (3/2, 1/2), (1, 1)]
        sage: B22.fusion_labels(['I0','Y1','X','Z','Xp','Y2'])
        sage: b = [B22(x) for x in B22.get_order()]; b
        [I0, Y1, X, Z, Xp, Y2]
        sage: [(x, x.weight()) for x in b]
        [(I0, (0, 0)),
         (Y1, (1, 0)),
         (X, (1/2, 1/2)),
         (Z, (2, 0)),
         (Xp, (3/2, 1/2)),
         (Y2, (1, 1))]
        sage: X*Y1
        X + Xp
        sage: Z*Z
        I0

    A fixed order of the basis keys is avalailable with :meth:`get_order`.
    This is the order used by methods such as :meth:`s_matrix`.
    You may use :meth:`set_order` to reorder the basis::

        sage: B22.set_order([x.weight() for x in [I0,Y1,Y2,X,Xp,Z]])
        sage: [B22(x) for x in B22.get_order()]
        [I0, Y1, Y2, X, Xp, Z]

    To reset the labels and the order to their defaults,
    you may run :meth:`fusion_labels` with no parameter::

        sage: B22.fusion_labels()
        sage: [B22(x) for x in B22.get_order()]
        [B22(0,0), B22(1,0), B22(0,1), B22(2,0), B22(1,1), B22(0,2)]

    The fusion ring has a number of methods that reflect its role
    as the Grothendieck ring of a modular tensor category. These
    include a twist method :meth:twist for its elements related
    to the ribbon element, and the S-matrix :meth:`s_ij`.

    Let us check the Verlinde formula. This famous identity states
    that

        .. MATH::

            N_{ijk} = \sum_l \frac{S(i,\ell)\,S(j,\ell)\,S(k,\ell)}{S(I,\ell)}

    where ``I`` is the unit object. We may define a function that corresponds
    to the right-hand side::

        sage: def V(i,j,k):
        ....:     R = i.parent()
        ....:     return sum(R.s_ij(i,l)*R.s_ij(j,l)*R.s_ij(k,l)/R.s_ij(R.one(),l) for l in R.basis())

    This does not produce ``self.N_ijk(i,j,k)`` exactly, because our S-matrix
    is omitting a normalization factor. The following code to check the
    Verlinde formula takes this into account::

       sage: def test_verlinde(R):
       ....:     b0 = R.one()
       ....:     c = V(b0, b0, b0)
       ....:     return all(V(i,j,k)==c*R.N_ijk(i,j,k) for i in R.basis() for j in R.basis() for k in R.basis())

    Every fusion ring should pass this test::

        sage: test_verlinde(FusionRing("A2",1))
        True
        sage: test_verlinde(FusionRing("B4",2))
        True

    """
    @staticmethod
    def __classcall__(cls, ct, k, base_ring=ZZ, prefix=None, style="coroots", conjugate=False):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: F1 = FusionRing('B3', 2)
            sage: F2 = FusionRing(CartanType('B3'), QQ(2), ZZ)
            sage: F3 = FusionRing(CartanType('B3'), int(2), style="coroots")
            sage: F1 is F2 and F2 is F3
            True

            sage: A23 = FusionRing('A2', 3)
            sage: TestSuite(A23).run()

            sage: B22 = FusionRing('B2', 2)
            sage: TestSuite(B22).run()

            sage: C31 = FusionRing('C3', 1)
            sage: TestSuite(C31).run()

            sage: D41 = FusionRing('D4', 1)
            sage: TestSuite(D41).run()
        """
        return super(FusionRing, cls).__classcall__(cls, ct, base_ring=base_ring,
                                                    prefix=prefix, style=style, k=k, conjugate=conjugate)

    def field(self):
        """
        This returns a cyclotomic field large enough to
        contain the `2l`-th roots of unity, as well as
        all the S-matrix entries.

        EXAMPLES::

            sage: FusionRing("A2",2).field()
            Cyclotomic Field of order 60 and degree 16
            sage: FusionRing("B2",2).field()
            Cyclotomic Field of order 40 and degree 16

        """

        return CyclotomicField(4*self._fg*self._l)
                
    def get_order(self):
        """
        This returns the weights of the basis vectors in a fixed order.
        You may change the order of the basis using :meth:`set_order`

        EXAMPLES::

            sage: A14=FusionRing("A1",4)
            sage: w = A14.get_order(); w
            [(0, 0), (1/2, -1/2), (1, -1), (3/2, -3/2), (2, -2)]
            sage: A14.set_order([w[k] for k in [0,4,1,3,2]])
            sage: [A14(x) for x in A14.get_order()]
            [A14(0), A14(4), A14(1), A14(3), A14(2)]

        This duplicates :meth:`get_order` from :mod:`combinat.free_module`.
        However unlike the :mod:`combinat.free_module` method with the same
        name this :meth:`get_order` is not cached. Caching of :meth:`get_order` causes
        inconsistent results after calling :meth:`set_order`.
        """
        if self._order is None:
            self.set_order(self.basis().keys().list())
        return self._order

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: D41 = FusionRing('D4', 1)
            sage: D41.some_elements()
            [D41(1,0,0,0), D41(0,0,1,0), D41(0,0,0,1)]
        """
        return [self.monomial(x) for x in self.fundamental_weights()
                if self.level(x) <= self._k]

    def fusion_k(self):
        r"""
        Return the level of the FusionRing.

        EXAMPLES::

            sage: B22=FusionRing('B2',2)
            sage: B22.fusion_k()
            2
        """
        return self._k

    def fusion_l(self):
        r"""
        This Returns the product `m_g(k + h^\vee)`, where `m_g` denotes the
        square of the ratio of the lengths of long to short roots of 
        the underlying Lie algebra, `k` denotes the level of the FusionRing,
        and `h^\vee` denotes the dual Coxeter number of the underlying Lie
        algebra.

        This value is used to define the associated root `2\ell`-th
        of unity `q=e^{i\pi/\ell}`.

        EXAMPLES::

            sage: B22=FusionRing('B2',2)
            sage: B22.fusion_l()
            10
            sage: D52=FusionRing('D5',2)
            sage: D52.fusion_l()
            10
        """
        return self._l

    def twists_matrix(self):
        r"""
        Return a diagonal matrix describing the twist corresponding to
        each simple object in the ``FusionRing``.

        EXAMPLES::

            sage: B22=FusionRing('B2',2)
            sage: B22.twists_matrix()
            [  0   0   0   0   0   0]
            [  0 4/5   0   0   0   0]
            [  0   0 1/2   0   0   0]
            [  0   0   0   2   0   0]
            [  0   0   0   0 3/2   0]
            [  0   0   0   0   0 6/5]
        """
        return diagonal_matrix(self.basis()[x].twist() for x in self.get_order())

    def q_dims(self):
        r"""
        Return a list of quantum dimensions of the simple objects.

        EXAMPLES::

            sage: F41=FusionRing('F4',1,conjugate=True)
            sage: F41.q_dims()
            [1, -zeta80^24 + zeta80^16 + 1]
            sage: B22=FusionRing("B2",2)
            sage: B22.q_dims()
            [1, 2, -2*zeta40^12 + 2*zeta40^8 + 1, 1, -2*zeta40^12 + 2*zeta40^8 + 1, 2]
        """
        b = self.basis()
        return [b[x].q_dimension() for x in self.get_order()]

    @cached_method
    def N_ijk(self, elt_i, elt_j, elt_k):
        """
        INPUT:

        - ``elt_i``, ``elt_j``, ``elt_k`` -- elements of the fusion basis

        Returns the symmetric fusion coefficient `N_{ijk}`. This
        is the same as `N_{ij}^{k\\ast}`, where $N_{ij}^k$ are the
        structure coefficients of the ring (see :meth:`Nk_ij`),
        and `k\\ast`` denotes the dual element. The coefficient `N_{ijk}`
        is unchanged under permutations of the three basis vectors.

        EXAMPLES::

            sage: G23=FusionRing("G2",3)
            sage: G23.fusion_labels("g")
            sage: b = G23.basis().list(); b
            [g0, g1, g2, g3, g4, g5]
            sage: [(x,y,z) for x in b for y in b for z in b if G23.N_ijk(x,y,z)>1]
            [(g3, g3, g3), (g3, g3, g4), (g3, g4, g3), (g4, g3, g3)]
            sage: all([G23.N_ijk(x,y,z)==G23.N_ijk(y,z,x) for x in b for y in b for z in b])
            True
            sage: all([G23.N_ijk(x,y,z)==G23.N_ijk(y,x,z) for x in b for y in b for z in b])
            True

        """
        return (elt_i*elt_j).monomial_coefficients().get(elt_k.dual().weight(),0)

    @cached_method
    def Nk_ij(self, elt_i, elt_j, elt_k):
        """
        Returns the fusion coefficient `N^k_{ij}`. These are
        the structure coefficients of the fusion ring, so
        
        .. MATH::

            i*j = \sum_{k} N_{ij}^k k

        as the example shows.

        EXAMPLES::

            sage: A22=FusionRing("A2",2)
            sage: b = A22.basis().list()
            sage: all(x*y == sum(A22.Nk_ij(x,y,k)*k for k in b) for x in b for y in b)
            True

        """
        return (elt_i*elt_j).monomial_coefficients().get(elt_k.weight(),0)

    @cached_method
    def s_ij(self, elt_i, elt_j):
        """
        Return the element of the S-matrix of this FusionRing corresponding to 
        the given elements. This is computed using the formula

        .. MATH::

            s_{i,j} = \\frac{1}{\\theta_i\\theta_j}\\sum_k N_{ik}^j d_k\\theta_k

        where `\\theta_k` is the twist and `d_k` is the quantum
        dimension. See [Row2006]_ equation (2.2) or [EGNO]_ Proposition 8.13.8.

        INPUT:

        - ``elt_i``, ``elt_j`` -- elements of the fusion basis

        EXAMPLES:: 

            sage: G21=FusionRing("G2",1)
            sage: b=G21.basis()
            sage: [G21.s_ij(x,y) for x in b for y in b]
            [1, -zeta60^14 + zeta60^6 + zeta60^4, -zeta60^14 + zeta60^6 + zeta60^4, -1]

        """
        l = self.fusion_l()*self._fg
        K = self.field()
        q = K.gen()
        ijtwist = -2*l*(elt_i.twist() + elt_j.twist())
        mc = (elt_i.dual()*elt_j).monomial_coefficients()
        return sum(self(k).q_dimension()*self.Nk_ij(elt_i,self(k),elt_j)*q**(2*l*self(k).twist()+ijtwist) for k in mc)

    def s_matrix(self):
        r"""
        Return the S-matrix of this FusionRing. This is the matrix denoted
        `\widetilde{s}` in [BK]_. It is not normalized to be unitary.

        EXAMPLES:: 

            sage: D91=FusionRing('D9',1)
            sage: D91.s_matrix()
            [         1          1          1          1]
            [         1          1         -1         -1]
            [         1         -1 -zeta68^17  zeta68^17]
            [         1         -1  zeta68^17 -zeta68^17]

            sage: D41=FusionRing('D4',1)
            sage: D41.s_matrix()
            [ 1  1  1  1]
            [ 1  1 -1 -1]
            [ 1 -1  1 -1]
            [ 1 -1 -1  1]
        """
        b = self.basis()
        return matrix([[self.s_ij(b[x],b[y]) for x in self.get_order()] for y in self.get_order()])

    def fusion_labels(self, labels=None):
        """
        Set the labels of the basis.

        INPUT:

        - ``labels`` -- (default: ``None``) a list of strings or string

        If ``labels`` is a list, the length of the list must equal the
        number of basis elements. These become the names of
        the basis elements. 

        If ``labels`` is a string, this is treated as a prefix and a
        list of names is generated.

        If ``labels`` is ``None``, then this resets the labels to the default.

        EXAMPLES::

            sage: A13 = FusionRing("A1", 3)
            sage: A13.fusion_labels("x")
            sage: fb = list(A13.basis()); fb
            [x0, x1, x2, x3]
            sage: Matrix([[x*y for y in A13.basis()] for x in A13.basis()])
            [     x0      x1      x2      x3]
            [     x1 x0 + x2 x1 + x3      x2]
            [     x2 x1 + x3 x0 + x2      x1]
            [     x3      x2      x1      x0]

        We reset the labels to the default::

            sage: A13.fusion_labels()
            sage: fb
            [A13(0), A13(1), A13(2), A13(3)]
        """
        if labels is None:
            self._order = None
            self._fusion_labels = None
            return
        elif type(labels) == str:
            labels = [labels+"%d"%k for k in range(len(self.basis()))]
        elif len(labels) != len(self.basis()):
            raise ValueError('invalid data')
        d = {}
        fb = list(self.get_order())
        for j, b in enumerate(fb):
            t = tuple([b.inner_product(x) for x in self.simple_coroots()])
            d[t] = labels[j]
            inject_variable(labels[j], self(b))
        self._fusion_labels = d

    class Element(WeylCharacterRing.Element):
        """
        A class for FusionRing elements.
        """
        def is_simple_object(self):
            """
            Determine whether element is a simple object of the FusionRing.

            EXAMPLES::

                sage: A22=FusionRing("A2",2)
                sage: x = A22(1,0); x
                A22(1,0)
                sage: x.is_simple_object()
                True
                sage: x^2
                A22(0,1) + A22(2,0)
                sage: (x^2).is_simple_object()
                False
            """
            return self.parent()._k is not None and len(self.monomial_coefficients())==1

        def weight(self):
            """
            This method is only available for basis elements. Returns the
            parametrizing dominant weight in the level `k` alcove.

            EXAMPLES::

                sage: A21 = FusionRing("A2",1)
                sage: sorted([x.weight() for x in A21.basis()])
                [(0, 0, 0), (1/3, 1/3, -2/3), (2/3, -1/3, -1/3)]
            """
            if len(self.monomial_coefficients()) != 1:
                raise ValueError("fusion weight is valid for basis elements only")
            return self.leading_support()

        def twist(self):
            r"""
            Compute the object's twist. This returns a rational number `h` such that 
            `\theta = e^{i \pi h}` is the twist of ``self``. This method is
            only available for simple objects.

            We compute the twists following p.7 of [Row2006]_, noting that the bilinear form
            is normalized so that `\langle\alpha, \alpha\rangle = 2` for short roots.

            EXAMPLES::

                sage: G21=FusionRing('G2',1)
                sage: G21.basis()
                Finite family {(0, 0, 0): G21(0,0), (1, 0, -1): G21(1,0)}
                sage: G21(1,0).twist()
                4/5
                sage: F41=FusionRing('F4',1,conjugate=True)
                sage: F41.basis()
                Finite family {(0, 0, 0, 0): F41(0,0,0,0), (1, 0, 0, 0): F41(0,0,0,1)}
                sage: F41(0,0,0,1).twist()
                4/5
            """
            if not self.is_simple_object():
                raise ValueError("quantum twist is only available for simple objects of a FusionRing")
            rho = sum(self.parent().positive_roots())/2
            lam = self.weight()
            inner = lam.inner_product(lam+2*rho)
            twist = self.parent()._conj*self.parent()._nf*inner/self.parent().fusion_l()
            #Reduce to canonical form 
            while twist > 2:
                twist -= 2
            while twist < 0:
                twist += 2
            return twist

        @cached_method
        def q_dimension(self):
            r"""
            This returns the quantum dimension as an element of the cyclotomic
            field of the `2l`-th roots of unity, where `l = m(k+h^\vee)`
            with `m=1,2,3` depending on whether type is simply, doubly or
            triply laced, `k` is the level and `h^\vee` is the dual Coxeter number.

            EXAMPLES::

                sage: B22=FusionRing("B2",2)
                sage: [(b.q_dimension())^2 for b in B22.basis()]
                [1, 4, 5, 1, 5, 4]
            """
            if not self.is_simple_object():
                raise ValueError("quantum dimension is only available for simple objects of a FusionRing")
            lam = self.weight()
            space = self.parent().space()
            rho = space.rho()
            num = reduce(mul, [q_int(self.parent()._nf*alpha.inner_product(lam+rho)) for alpha in space.positive_roots()], 1)
            den = reduce(mul, [q_int(self.parent()._nf*alpha.inner_product(rho)) for alpha in space.positive_roots()], 1)
            expr = num/den
            pr = expr.parent().ring()
            q = pr.gen()**2
            expr = pr(expr)
            expr = expr.substitute(q=q**2)/q**(expr.degree())
            zet = self.parent().field().gen()**(self.parent()._fg)
            return expr.substitute(q=zet)
