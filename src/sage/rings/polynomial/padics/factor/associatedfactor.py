r"""
Irreduble factors of associate polynomials needed for OM computations

AUTHORS:

- Brian Sinclair and Sebastian Pauli (2012-02-22): initial version

"""
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.constructor import GF
from sage.matrix.constructor import Matrix
from sage.rings.polynomial.padics.factor.frameelt import FrameElt
from sage.rings.infinity import infinity

class AssociatedFactor:
    r"""
    An irreducible factor of the associated polynomials of higher order
    newton polygon segments needed for OM computation.

    For each distinct irreducible factor of the associated polynomial,
    the tree of OM representations branches, thus producing distinct factors
    of the original polynomial.

    If ``rho`` is not linear, then we have found inertia. Future associated
    polynomials wll need to be produced over an extension over our ground
    field by ``rho``.  This can produce a tower of finite field extensions
    to be worked in.

    INPUT:

    - ``segment`` -- The segment whose associated polynomial self is a factor of.

    - ``rho`` -- The irreducible finite field polynomial factor of the
      associated polynomial.

    - ``rhoexp`` -- The multiplicity of the factor.

    """
    def __init__(self,segment,rho,rhoexp):
        """
        Initialization.

        See ``AssociatedFactor`` for full documentation.

        TESTS::

            sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: t = OM_tree(x^4+20*x^3+44*x^2+80*x+1040)
            sage: TestSuite(t[0].prev).run()
        """
        self.segment = segment
        self.rho = rho
        self.rhoexp = rhoexp
        self.Fplus = self.rho.degree()

        if self.segment.frame.is_first():
            # In the first frame, so FFbase is the residue class field of O
            self.FFbase = self.segment.frame.R
        else:
            # Not the first frame
            self.FFbase = self.segment.frame.prev.FF

        if self.Fplus == 1:
            self.FF = self.FFbase
            self.FFz = PolynomialRing(self.FF,'z'+str(self.segment.frame.depth))
            # rho is linear delta is the root of rho
            self.delta = self.rho.roots()[0][0]
        else:
            self.FF = GF(self.FFbase.order()**self.Fplus,'a'+str(self.segment.frame.depth))
            self.FFz = PolynomialRing(self.FF,'z'+str(self.segment.frame.depth))
            self.FFbase_gamma = (self.FFz(self.FFbase.modulus())).roots()[0][0]
            FFrho = self.FFz([self.FFbase_elt_to_FF(a) for a in list(rho)])
            self.gamma = FFrho.roots()[0][0]
            basis = [(self.gamma**j*self.FFbase_gamma**i).polynomial() for j in range(0,self.Fplus) for i in range(0,self.FFbase.degree())]
            self.basis_trans_mat = Matrix([self.FF(b)._vector_() for b in basis])

    def __cmp__(self, other):
        """
        Comparison.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: t = OM_tree(x^4+20*x^3+44*x^2+80*x+1040)
            sage: t[0].prev == t[0].polygon[0].factors[0]
            False
        """
        c = cmp(type(self), type(other))
        if c: return c
        return cmp((self.segment, self.rho, self.rhoexp), (other.segment, other.rho, other.rhoexp))

    def FF_elt_to_FFbase_vector(self,a):
        """
        Represents an element in our current extended residue field as a
        vector over its ground residue field.

        INPUT:

        - ``a`` -- Element of our extended residue field

        OUTPUT:

        - A list representing a vector of ``a`` over the ground field of
          the latest extension.

        EXAMPLES::

        First we set up AssociatedFactors building a tower of extensions::

            sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: t = OM_tree(x^4+20*x^3+44*x^2+80*x+1040)
            sage: t[0].prev
            AssociatedFactor of rho z^2 + z + 1
            sage: t[0].polygon[0].factors[0]
            AssociatedFactor of rho z0^2 + a0*z0 + 1

        Then we take elements in the different finite fields and represent
        them as vectors over their base residue field::

            sage: K.<a0> = t[0].prev.FF;K
            Finite Field in a0 of size 2^2
            sage: t[0].prev.FF_elt_to_FFbase_vector(a0+1)
            [1, 1]
            sage: L.<a1> = t[0].polygon[0].factors[0].FF;L
            Finite Field in a1 of size 2^4
            sage: t[0].polygon[0].factors[0].FF_elt_to_FFbase_vector(a1)
            [1, a0 + 1]

        """
        if self.segment.frame.is_first() and self.Fplus == 1:
            return a
        elif self.Fplus == 1:
            return self.segment.frame.prev.FF_elt_to_FFbase_vector(a)
        else:
            basedeg = self.FFbase.degree()
            avec = self.FF(a)._vector_()
            svector = self.basis_trans_mat.solve_left(Matrix(self.FF.prime_subfield(),avec))
            s_list = svector.list()
            s_split = [ s_list[i*basedeg:(i+1)*basedeg] for i in range(0,self.Fplus)]
            s = [sum([ss[i]*self.FFbase.gen()**i for i in range(0,len(ss))]) for ss in s_split]
            return s

    def FFbase_elt_to_FF(self,b):
        """
        Lifts an element up from the previous residue field to the current
        extended residue field.

        INPUT:

        - ``b`` -- Element in the previous residue field.

        OUTPUT:

        - An element in the current extended residue field.

        EXAMPLES::

        First we set up AssociatedFactors building a tower of extensions::

            sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: t = OM_tree(x^4+20*x^3+44*x^2+80*x+1040)

        Then we take elements in the different finite fields and lift them
        to the next residue field upward in the extension tower::

            sage: K.<a0> = t[0].prev.FF;K
            Finite Field in a0 of size 2^2
            sage: L.<a1> = t[0].polygon[0].factors[0].FF;L
            Finite Field in a1 of size 2^4
            sage: t[0].prev.FFbase_elt_to_FF(1)
            1
            sage: t[0].polygon[0].factors[0].FFbase_elt_to_FF(a0+1)
            a1^2 + a1 + 1

        """
        if self.segment.frame.is_first():
            return b
        elif self.Fplus == 1:
            return self.segment.frame.prev.FFbase_elt_to_FF(b)
        elif self.segment.frame.F == 1:
            return b * self.FFbase_gamma
        else:
            bvec = b._vector_()
            return sum([ bvec[i]*self.FFbase_gamma**i for i in range(len(bvec))])

    def __repr__(self):
        """
        Representation of self.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: t = OM_tree(x^4+20*x^3+44*x^2+80*x+1040)
            sage: t[0].prev.__repr__()
            'AssociatedFactor of rho z^2 + z + 1'
            sage: t[0].polygon[0].factors[0].__repr__()
            'AssociatedFactor of rho z0^2 + a0*z0 + 1'

        """
        return "AssociatedFactor of rho "+repr(self.rho)

    def lift(self,delta):
        """
        FrameElt representation of a lift of residue field element ``delta``.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: t = OM_tree(x^4+20*x^3+44*x^2+80*x+1040)
            sage: K.<a0> = t[0].prev.FF;K
            Finite Field in a0 of size 2^2
            sage: t[0].polygon[0].factors[0].lift(a0+1)
            [[1*2^0]phi1^0, [1*2^-1]phi1^1]

        """
        if self.segment.frame.F == 1:
            return FrameElt(self.segment.frame,self.segment.frame.Ox(delta))
        elif self.segment.frame.prev.Fplus == 1:
            return FrameElt(self.segment.frame,self.segment.frame.prev.lift(delta),this_exp=0)
        else:
            dvec = self.segment.frame.prev.FF_elt_to_FFbase_vector(delta)
            return sum([self.segment.frame.prev.gamma_frameelt**i*FrameElt(self.segment.frame,self.segment.frame.prev.lift(dvec[i]),this_exp=0) for i in range(len(dvec)) if dvec[i] != 0])

    def next_frame(self,length=infinity):
        """
        Produce the child Frame in the tree of OM representations with the
        partitioning from self.

        This method generates a new Frame with the ``self`` as previous and
        seeds it with a new approximation with strictly greater valuation
        than the current one.

        INPUT:

        - ``length`` -- Integer or infinity, default infinity; The length of
          the segment generating this factor.  This is used to reduce the
          total number of quotient with remainder operations needed in the
          resulting Frame.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: Phi = ZpFM(2,20,'terse')['x'](x^32+16)
            sage: f = Frame(Phi)
            sage: f.seed(Phi.parent().gen());f
            Frame with phi (1 + O(2^20))*x + (0 + O(2^20))
            sage: f = f.polygon[0].factors[0].next_frame();f
            Frame with phi (1 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (0 + O(2^20))*x^5 + (0 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (0 + O(2^20))*x^2 + (0 + O(2^20))*x + (1048574 + O(2^20))
            sage: f = f.polygon[0].factors[0].next_frame();f
            Frame with phi (1 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (0 + O(2^20))*x^5 + (0 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (1048574 + O(2^20))*x^2 + (0 + O(2^20))*x + (1048574 + O(2^20))
            sage: f = f.polygon[0].factors[0].next_frame();f
            Frame with phi (1 + O(2^20))*x^16 + (0 + O(2^20))*x^15 + (0 + O(2^20))*x^14 + (0 + O(2^20))*x^13 + (0 + O(2^20))*x^12 + (0 + O(2^20))*x^11 + (1048572 + O(2^20))*x^10 + (0 + O(2^20))*x^9 + (1048572 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (1048572 + O(2^20))*x^5 + (4 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (8 + O(2^20))*x^2 + (0 + O(2^20))*x + (4 + O(2^20))

        """
        from frame import Frame
        if self.segment.slope == infinity:
            next = Frame(self.segment.frame.Phi,self,self.segment.frame.iteration)
            self.next = next
            next.seed(self.segment.frame.phi,length=length)
            return next
        if self.Fplus == 1 and self.segment.Eplus == 1:
            next = Frame(self.segment.frame.Phi,self.segment.frame.prev,self.segment.frame.iteration)
        else:
            next = Frame(self.segment.frame.Phi,self,self.segment.frame.iteration)
        self.next = next
        self.gamma_frameelt = FrameElt(next,self.segment.psi**-1,self.segment.Eplus)
        if self.Fplus == 1 and self.segment.frame.F == 1:
            next_phi = self.segment.frame.phi**self.segment.Eplus-(self.segment.psi.polynomial()*self.segment.frame.Ox(self.delta))
            self.reduce_elt = FrameElt(next,self.segment.psi*self.lift(self.delta),0)
            next.seed(next_phi,length=length)
        elif self.Fplus == 1 and self.segment.Eplus == 1:
            delta_elt = self.lift(self.delta)
            next_phi_tail = self.segment.psi*delta_elt.reduce()
            next_phi = self.segment.frame.phi-next_phi_tail.polynomial()
            self.reduce_elt = FrameElt(next,next_phi_tail,0)
            next.seed(next_phi,length=length)
        else:
            lifted_rho_coeffs = [self.lift(r) for r in list(self.rho)]
            lifted_rho_coeffs_with_psi = [FrameElt(next,(self.segment.psi**(self.Fplus-i)*lifted_rho_coeffs[i]).reduce(),0) for i in range(len(lifted_rho_coeffs))]
            phi_elt = FrameElt(next,self.segment.frame.Ox(1),1)
            next_phi_tail = sum([phi_elt**(self.segment.Eplus*i)*lifted_rho_coeffs_with_psi[i] for i in range(len(lifted_rho_coeffs_with_psi)-1)])
            next_phi = (phi_elt**(self.segment.Eplus*self.Fplus)+next_phi_tail).polynomial()
            self.reduce_elt = FrameElt(next)+(-next_phi_tail) # that is -next_phi_tail
            next.seed(next_phi,length=length)
        return next
