r"""
Implements representations of a Iwahori-Hecke algebra as
:class:`CombinatorialFreeModules`, including cell representations and in types A, B and,
more generally, `G(r,1,n)`, the Murphy and seminormal forms of the Specht
modules.

<Paragraph description>

AUTHORS:

- Andrew Mathas (2014-12-03): initial version

EXAMPLES::

<Lots and lots of examples>
"""

#*****************************************************************************
#    Copyright (C) 2014 Andrew Mathas <andrew.mathas@sydney.edu.au>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version. See http://www.gnu.org/licenses/
#    for more details.
#*****************************************************************************


from sage.algebras.algebra import is_Algebra
from sage.categories.action import Action
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.cachefunc import cached_method
from sage.rings.fraction_field import FractionField_generic
from sage.rings.integer_ring import IntegerRing
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.polynomial.polynomial_ring import polygen
import operator

class IwahoriHeckeAlgebraRepresentation(CombinatorialFreeModule):
    r"""
    Implement a particular seminormal representation of the Iwahori-Hecke algebras of type A that 
    restricts well to the alternating group. This module has a basis $\{f_t\}$ indexed by
    the set of standard tableaux of the corresponding shape. If $1\le r\le n$
    then the action of `T_r` on the basis element `f_t` is given by

    .. MATH::

       f_t T_r = \frac{-q^{-1}}{[d]}f_t + \alpha_d f_s

    where $d=c_r(t)-c_{r+1}(t)$  is the axial distance from $r$ to $r+1$ in `t` and

    .. MATH::

       \alpha_d = \begin{cases}
                     \frac{\sqrt{-1}\sqrt{[d+1]}\sqrt{[d-1]}}{[d]},&\text{if }d>0,\\
                    -\frac{\sqrt{-1}\sqrt{[-d+1]}\sqrt{[-d-1]}}{[d]},&\text{if }d<0.
       \end{cases}

    The whole point of choosing this particular basis is that if the module
    is indexed by a self-conjugate partition then it splits when it is
    restricted to the alternating Hecke algebra.

    Finally, it should be mentioned that the particular Iwahori-Hecke algebra
    being used here is the one with quadratic relations

    .. MATH::

        (T_r - q)(T_r + q^{-1}) = 0, \qquad\text{ for }1\le r<n.

    EXAMPLES::

        sage: rep=SeminormalRepresentation([3,2,1]); rep
        SeminormalRepresentation([3,2,1])
        sage: rep([[1,2,3],[4,5],[6]]).T(1)
        q*f(1,2,3/4,5/6)
        sage: rep([[1,2,3],[4,5],[6]]).T(1,2)
        q^2*f(1,2,3/4,5/6)
        sage: rep([[1,2,3],[4,5],[6]]).T(1,2,3)
        (q^2/(-q-q^3-q^5))*f(1,2,3/4,5/6) + (q^2*r2*r4*I/(1+q^2+q^4))*f(1,2,4/3,5/6)
        sage: rep.character(1)
        (-8*q^7 + 8*q^9)/q^8
        sage: rep.an_element()
        3*f(1,2,6/3,5/4) + 2*f(1,3,6/2,5/4) + 2*f(1,4,6/2,5/3)
        sage: rep([[1,3,5],[2,6],[4]])
        f(1,3,5/2,6/4)

        sage: SeminormalRepresentation([3,2,2], prefix='v').an_element()
        3*v(1,2,7/3,5/4,6) + 2*v(1,3,7/2,5/4,6) + 2*v(1,4,7/2,5/3,6)
    """
    def __init__(self, R, basis_keys, prefix, q1,q2, left_module=True):
        pass

    def hecke_algebra_action(self, H):
        r"""
        Return the action of the Hecke algebra `H` upon `self`.
        """
        return IwahoriHeckeAlgebraAction(H,self)

    def _get_action_(self, H, op, self_on_left):
        print 'H=%s, op=%s, self_on_left=%s'%(H,op,self_on_left)
        if is_Algebra(H) and op==operator.mul and self_on_left:
            return self.hecke_algebra_action(H)
        return None

    def character(self, *arg):
        r"""
        Return the trace of `T_w` acting on the representation, where `w` is the
        permutation corresponding to the (reduced) expression given by `*arg`.

        EXAMPLES::

                sage: SeminormalRepresentation([3,2,1]).character()
                16
                sage: SeminormalRepresentation([3,2,1]).character(1)
                (-8*q^7 + 8*q^9)/q^8
        """
        trace=0
        for t in self._basis_keys:
            trace+=self(t).T(*arg).coefficient(t)
        return trace

    @cached_method
    def _Tr_matrix(self,r):
        r"""
        Return the matrix that describes the action of `T_r` on the seminormal representation `self`.
        """
        mat={} # create a dictionary of the non-zero entries in the matrix for passing to matrix()
        for s in self._basis_keys:
            for (tab,coeff) in self(s).T(r):
                mat[(self._basis_keys.rank(s),self._basis_keys.rank(tab))] = coeff
        return matrix(self._base, mat)

    @cached_method
    def T_matrix(self,*arg):
        r"""
        Return the matrix that describes the action of `T_w` on the seminormal representation `self`.

        EXAMPLES::

            sage: SeminormalRepresentation([2,2]).T_matrix(1)
            [1/(-q)      0]
            [     0      q]
            sage: SeminormalRepresentation([2,2]).T_matrix(2)
            [    q^3/(1 + q^2) (-r3*I)/(1 + q^2)]
            [   r3*I/(1 + q^2)      1/(-q - q^3)]
        """
        return prod(self._Tr_matrix(r) for r in arg)

    def _T_on_basis(self, b, r):
        r"""
        Returns the result of `T_r` acting on `self(b)`.
        """
        raise NotImplementedError('the action on the basis has not been implemented yet!')

    def action(self, a):
        return sum(c*self.T_matrix( *w.reduced_word() ) for (w,c) in a.to_T_basis())

    def check_relations(self, verbose=False):
        r"""
        Return `True` if all of the relations are satisfied on `self` and
        `False` otherwise.

        EXAMPLES::

            sage: SeminormalRepresentation([2,1]).check_relations()
            True
        """
        alpha=self.q-self.q**-1#  => T_r^2 = \alpha T_r+1
        for t in self._basis_keys:
            if verbose: print 'Checking relations on ({term})'.format(term=self(t))
            for r in range(1,self.size()):
                # check quadratic relations on self(t)
                assert self(t).T(r,r)==alpha*self(t).T(r)+self(t), (
                        'Quadratic relation failed for T_{r} on {term}'.format(r=r,term=self(t))
                )
            if r<self.size()-1:
                # check length 3 braid relation on self(t): T_rT_{r+1}T_r = T_{r+1}T_rT_{r+1}
                assert self(t).T(r,r+1,r)==self(t).T(r+1,r,r+1), (
                        'Length 3 braid relation failed for T_{r} on ({term})'.format(r=r,term=self(t))
                )
            for s in range(r+2,self.size()):
                # check length 2 braid relation on self(t): T_rT_s=T_sT_r if |r-s|>1
                assert self(t).T(r,s)==self(t).T(s,r), (
                        'Length 2 braid relation failed for T_{r} on ({term})'.format(r=r,term=self(t))
                )
        return True


    class Element(CombinatorialFreeModule.Element):

        def _Tr(self,r):
            r"""
            This is a helper function for computing the action of `T_r` upon the
            element `self` of the seminormal representation.

            EXAMPLES::

                sage: SeminormalRepresentation([2,2]).an_element().T(1)
                2*q*f(1,2/3,4) + (2/(-q))*f(1,3/2,4)
            """

            return self.parent().sum(coeff*self.parent()._T_on_basis(b,r) for (b,coeff) in self)

        def T(self,*args):
            r"""
            Compute the action of `T_w` on self, where `w` is the permutation
            corresponding to the (reduced) word `*args`.

            EXAMPLES::

                sage: SeminormalRepresentation([2,2]).an_element().T(1)
                2*q*f(1,2/3,4) + (2/(-q))*f(1,3/2,4)
            """
            v=self
            for r in args:
                v=v._Tr(r)
            return v

        def action(self,a):
            return self.parent().sum(c*self.T( *w.reduced_word() ) for (w,c) in a.to_T_basis)



class IwahoriHeckeAlgebraAction(Action):
    r"""
    Implements the action of a Hecke algebra, or symmetric group, on a Specht
    module.
    """
    def __init__(self, H, rep, is_left=1):
        super(IwahoriHeckeAlgebraAction,self).__init__(rep, H, is_left, operator.mul)

    def _call_(self, v, a):
        r"""
        Return the action of the algebra element `a` on the seminormal basis vector `v`.

        EXAMPLES::
        """
        return v.action(a)


###############################################################################
##       Kazhdan-Lusztig cell representations of Iwahori-Hecke algebras      ##
###############################################################################

class LeftCellRepresentationOfHeckeAlgebra(IwahoriHeckeAlgebraRepresentation):
    pass

class RightCellRepresentationOfHeckeAlgebra(IwahoriHeckeAlgebraRepresentation)
    pass

###############################################################################
##     Seminormal representations of (cyclotomic) Iwahori-Hecke algebras     ##
###############################################################################

class SeminormalRepresentation_generic(IwahoriHeckeAlgebraRepresentation):

    @staticmetho
    def __classcall__(cls, shape, prefix='f'):
        """
        Magic to allow class to accept a list (which is not hashable) instead
        of a composition (which apparently is).

        EXAMPLES::

            sage: SeminormalRepresentation([3,2,1])    # indirect doctest
            SeminormalRepresentation([3,2,1])
        """
        return super(SeminormalRepresentation, cls).__classcall__(cls, Partition(shape),prefix)

    def __init__(self, shape, prefix):
        r"""
        Initialisation of a seminormal representation of the symmetric group
        indexed by the partition `shape` and where the basis elements are
        labelled by `prefix`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2,1])    # indirect doctest
            SeminormalRepresentation([3,2,1])
        """
        # We need to build the base ring that splits our module. It needs to
        # contain the inverses and square roots of the quantum integers [k]_q,
        # for 0\le k\le n. To achieve this we start with the ring of Laurent
        # polynomials and successively extend this ring by throwing in the
        # square roots that we need. At the end we throw in a square of -1. We
        # don't use the ComplexField() because this seems to drag along with it
        # unwanted precision when printing.
        self._shape=shape
        self._prefix=prefix
        Fq=FractionField(PolynomialRing(IntegerRing(),'q')) 
        self.q=Fq.gen()
        roots=['r%d'%(d+1) for d in range(0,shape.size())] # names for the square roots
        roots[0]='I'
        poly_ring=PolynomialRing(Fq,roots)
        relations=[poly_ring.gen(0)**2+1]
        self.__quantum_integer={0:0, 1:1}       # dictionary of positive quantum integers
        self.__quantum_integer_root={0:0, 1:1}  # dictionary of square roots of positive quantum integers
        qint=1                                  # a quantum integer in waiting
        for k in range(1,shape.size()):
            qint+=self.q**(2*k)                 # now equal to 1+q^2+...+q^2k = [k+1]_{q^2}
            self.__quantum_integer[k+1]=qint    # we remember these because we need them later
            relations.append(poly_ring.gen(k)**2-qint)
        F=poly_ring.quotient(poly_ring.ideal(relations),roots)
        for k in range(1,shape.size()):
            self.__quantum_integer_root[k+1]=F.gen(k)
        self.I=F.gen(0)                         # I=sqrt(-1). We don't use CC as we want issues with precision
        CombinatorialFreeModule.__init__(self,F,shape.standard_tableaux(), prefix='f')

    def shape(self):
        r"""
        Return the partition that indexes the seminormal representation `self`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2]).shape()
            [3, 2]
        """
        return self._shape

    def size(self):
        r"""
        Return the size of the partition that indexes the seminormal representation `self`.

        EXAMPLES::

        sage: SeminormalRepresentation([3,2]).size()
        5
        """
        return self._shape.size()

    def _quantum_integer(self,d):
        r"""
        Return the quantum integer `[d]_{q^2}=1+q^2+\dots+q^{2d-2}`, where `d\ge0`.

        EXAMPLES::

        sage: SeminormalRepresentation([3,2])._quantum_integer(3)
        1 + q^2 + q^4
        sage: SeminormalRepresentation([3,2])._quantum_integer(-3)
        -q^-6 - q^-4 - q^-2
        """
        return self.__quantum_integer[d] if d>=0 else -self.q**(2*d)*self.__quantum_integer[-d]

    def _repr_(self):
        r"""
        Return a string representation of the seminormal representation `self`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2,1])    # indirect doctest
            SeminormalRepresentation([3,2,1])
        """
        return 'SeminormalRepresentation([%s])' % self._shape._repr_compact_high()

    def _repr_term(self,t):
        r"""
        Override the default printing of terms so that compact tableaux are
        used.

        EXAMPLES::

            sage: SeminormalRepresentation([2,1]).an_element()
            2*f(1,2/3) + 2*f(1,3/2)
            sage: SeminormalRepresentation([2,1], prefix='v').an_element()
            2*v(1,2/3) + 2*v(1,3/2)
        """
        return '%s(%s)' % ( self._prefix, t._repr_compact() )

    def _element_constructor_(self,t):
        r"""
        We override the fault element constructor in order to coerce lists into standard tableaux.
        
        EXAMPLES::

            sage: SeminormalRepresentation([3,2])([[1,2,3],[4,5]])
            f(1,2,3/4,5)
        """
        return self.monomial(self._basis_keys(t))

class SeminormalRepresentation(SeminormalRepresentation_generic):

    def __init__(self, shape, prefix='v'):
        r"""
        Initialisation of a seminormal representation of the symmetric group
        indexed by the partition `shape` and where the basis elements are
        labelled by `prefix`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2,1])    # indirect doctest
            SeminormalRepresentation([3,2,1])
        """
        # We need to build the base ring that splits our module. It needs to
        # contain the inverses and square roots of the quantum integers [k]_q,
        # for 0\le k\le n. To achieve this we start with the ring of Laurent
        # polynomials and successively extend this ring by throwing in the
        # square roots that we need. At the end we throw in a square of -1. We
        # don't use the ComplexField() because this seems to drag along with it
        # unwanted precision when printing.
        self._shape=shape
        self._prefix=prefix
        Fq=FractionField(PolynomialRing(IntegerRing(),'q')) 
        self.q=Fq.gen()
        roots=['r%d'%(d+1) for d in range(0,shape.size())] # names for the square roots
        roots[0]='I'
        poly_ring=PolynomialRing(Fq,roots)
        relations=[poly_ring.gen(0)**2+1]
        self.__quantum_integer={0:0, 1:1}       # dictionary of positive quantum integers
        self.__quantum_integer_root={0:0, 1:1}  # dictionary of square roots of positive quantum integers
        qint=1                                  # a quantum integer in waiting
        for k in range(1,shape.size()):
            qint+=self.q**(2*k)                 # now equal to 1+q^2+...+q^2k = [k+1]_{q^2}
            self.__quantum_integer[k+1]=qint    # we remember these because we need them later
            relations.append(poly_ring.gen(k)**2-qint)
        F=poly_ring.quotient(poly_ring.ideal(relations),roots)
        for k in range(1,shape.size()):
            self.__quantum_integer_root[k+1]=F.gen(k)
        self.I=F.gen(0)                         # I=sqrt(-1). We don't use CC as we want issues with precision
        CombinatorialFreeModule.__init__(self,F,shape.standard_tableaux(), prefix='f')

    @cached_method
    def _T_on_basis(self,tab,r):
        r"""
        Compute the action of `T_r` on the basis element `f_tab`, where the action
        of `T_r` is as follows:

        .. MATH::

           f_t T_r = \frac{-q^{-1}}{[d]}f_t + \alpha_d f_s

        where $d=c_r(t)-c_{r+1}(t)$  is the axial distance from $r$ to $r+1$ in `t` and

        .. MATH::

           \alpha_d = \begin{cases}
                         \frac{\sqrt{-1}\sqrt{[d+1]}\sqrt{[d-1]}}{[d]},&\text{if }d>0,\\
                        -\frac{\sqrt{-1}\sqrt{[-d+1]}\sqrt{[-d-1]}}{[d]},&\text{if }d<0.
           \end{cases}

         The whole point of choosing this particular basis is that if the module
         is indexed by a self-conjugate partition then it splits when it is
         restricted to the alternating Hecke algebra.

         EXAMPLES::

             sage: rep=SeminormalRepresentation([3,2,1])
             sage: rep([[1,2,3],[4,5],[6]]).T(1)   # indirect doctest
             q*f(1,2,3/4,5/6)
             sage: rep([[1,2,3],[4,5],[6]]).T(1,2)
             q^2*f(1,2,3/4,5/6)
             sage: rep([[1,2,3],[4,5],[6]]).T(1,2,3)
             (q^2/(-q-q^3-q^5))*f(1,2,3/4,5/6) + (q^2*r2*r4*I/(1+q^2+q^4))*f(1,2,4/3,5/6)
        """
        rtab=tab.symmetric_group_action_on_entries( Permutation((r,r+1)) )
        d=tab.content(r)-tab.content(r+1)
        denom=1/(self.q**(2*d+1)*self._quantum_integer(-d))
        if abs(d)==1:  # r and r+1 in same row or column <=> rtab is not standard
            return self.term(tab,denom)
        else:
            if d>0:
                coeff=1/self._quantum_integer(d)
                coeff*=self.I*self._quantum_integer_root(d+1)*self._quantum_integer_root(d-1)
            else:
                coeff=1/self._quantum_integer(-d)
                coeff*=-self.I*self._quantum_integer_root(-d+1)*self._quantum_integer_root(-d-1)
            return self.sum_of_terms([(tab,denom), (rtab,coeff)],distinct=True)

class SeminormalRepresentation_orthogonal(SeminormalRepresentation_generic):

    def __init__(self, shape, prefix='y'):
        r"""
        Initialisation of a seminormal representation of the symmetric group
        indexed by the partition `shape` and where the basis elements are
        labelled by `prefix`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2,1])    # indirect doctest
            SeminormalRepresentation([3,2,1])
        """
        # We need to build the base ring that splits our module. It needs to
        # contain the inverses and square roots of the quantum integers [k]_q,
        # for 0\le k\le n. To achieve this we start with the ring of Laurent
        # polynomials and successively extend this ring by throwing in the
        # square roots that we need. At the end we throw in a square of -1. We
        # don't use the ComplexField() because this seems to drag along with it
        # unwanted precision when printing.
        self._shape=shape
        self._prefix=prefix
        Fq=FractionField(PolynomialRing(IntegerRing(),'q')) 
        self.q=Fq.gen()
        roots=['r%d'%(d+1) for d in range(0,shape.size())] # names for the square roots
        roots[0]='I'
        poly_ring=PolynomialRing(Fq,roots)
        relations=[poly_ring.gen(0)**2+1]
        self.__quantum_integer={0:0, 1:1}       # dictionary of positive quantum integers
        self.__quantum_integer_root={0:0, 1:1}  # dictionary of square roots of positive quantum integers
        qint=1                                  # a quantum integer in waiting
        for k in range(1,shape.size()):
            qint+=self.q**(2*k)                 # now equal to 1+q^2+...+q^2k = [k+1]_{q^2}
            self.__quantum_integer[k+1]=qint    # we remember these because we need them later
            relations.append(poly_ring.gen(k)**2-qint)
        F=poly_ring.quotient(poly_ring.ideal(relations),roots)
        for k in range(1,shape.size()):
            self.__quantum_integer_root[k+1]=F.gen(k)
        self.I=F.gen(0)                         # I=sqrt(-1). We don't use CC as we want issues with precision
        CombinatorialFreeModule.__init__(self,F,shape.standard_tableaux(), prefix='f')

    def _quantum_integer_root(self,d):
        r"""
        Return the square root of the quantum integer `[d]_{q^2}`. These square
        roots are chosen so that `\sqrt{[-d]}=\sqrt{-1}q^d\sqrt{[d]` if `d<0`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2])._quantum_integer_root(3)
            r3
            sage: SeminormalRepresentation([3,2])._quantum_integer_root(-3)
            ((q^-3)*r3)*I

        """
        return self.__quantum_integer_root[d] if d>=0 else -self.q**(2*d)*self.__quantum_integer_root[-d]

class SeminormalRepresentation_murphy(SeminormalRepresentation_generic):

    def __init__(self, shape, prefix):
        r"""
        Initialisation of a seminormal representation of the symmetric group
        indexed by the partition `shape` and where the basis elements are
        labelled by `prefix`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2,1])    # indirect doctest
            SeminormalRepresentation([3,2,1])
        """
        # We need to build the base ring that splits our module. It needs to
        # contain the inverses and square roots of the quantum integers [k]_q,
        # for 0\le k\le n. To achieve this we start with the ring of Laurent
        # polynomials and successively extend this ring by throwing in the
        # square roots that we need. At the end we throw in a square of -1. We
        # don't use the ComplexField() because this seems to drag along with it
        # unwanted precision when printing.
        self._shape=shape
        self._prefix=prefix
        Fq=FractionField(PolynomialRing(IntegerRing(),'q')) 
        self.q=Fq.gen()
        roots=['r%d'%(d+1) for d in range(0,shape.size())] # names for the square roots
        roots[0]='I'
        poly_ring=PolynomialRing(Fq,roots)
        relations=[poly_ring.gen(0)**2+1]
        self.__quantum_integer={0:0, 1:1}       # dictionary of positive quantum integers
        self.__quantum_integer_root={0:0, 1:1}  # dictionary of square roots of positive quantum integers
        qint=1                                  # a quantum integer in waiting
        for k in range(1,shape.size()):
            qint+=self.q**(2*k)                 # now equal to 1+q^2+...+q^2k = [k+1]_{q^2}
            self.__quantum_integer[k+1]=qint    # we remember these because we need them later
            relations.append(poly_ring.gen(k)**2-qint)
        F=poly_ring.quotient(poly_ring.ideal(relations),roots)
        for k in range(1,shape.size()):
            self.__quantum_integer_root[k+1]=F.gen(k)
        self.I=F.gen(0)                         # I=sqrt(-1). We don't use CC as we want issues with precision
        CombinatorialFreeModule.__init__(self,F,shape.standard_tableaux(), prefix='f')

    @cached_method
    def _T_on_basis(self,tab,r):
        r"""
        Compute the action of `T_r` on the basis element `f_tab`, where the action
        of `T_r` is as follows:

        .. MATH::

           f_t T_r = \frac{-q^{-1}}{[d]}f_t + \alpha_d f_s

        where $d=c_r(t)-c_{r+1}(t)$  is the axial distance from $r$ to $r+1$ in `t` and

        .. MATH::

           \alpha_d = \begin{cases}
                         \frac{\sqrt{-1}\sqrt{[d+1]}\sqrt{[d-1]}}{[d]},&\text{if }d>0,\\
                        -\frac{\sqrt{-1}\sqrt{[-d+1]}\sqrt{[-d-1]}}{[d]},&\text{if }d<0.
           \end{cases}

         The whole point of choosing this particular basis is that if the module
         is indexed by a self-conjugate partition then it splits when it is
         restricted to the alternating Hecke algebra.

         EXAMPLES::

             sage: rep=SeminormalRepresentation([3,2,1])
             sage: rep([[1,2,3],[4,5],[6]]).T(1)   # indirect doctest
             q*f(1,2,3/4,5/6)
             sage: rep([[1,2,3],[4,5],[6]]).T(1,2)
             q^2*f(1,2,3/4,5/6)
             sage: rep([[1,2,3],[4,5],[6]]).T(1,2,3)
             (q^2/(-q-q^3-q^5))*f(1,2,3/4,5/6) + (q^2*r2*r4*I/(1+q^2+q^4))*f(1,2,4/3,5/6)
        """
        rtab=tab.symmetric_group_action_on_entries( Permutation((r,r+1)) )
        d=tab.content(r)-tab.content(r+1)
        denom=1/(self.q**(2*d+1)*self._quantum_integer(-d))
        if abs(d)==1:  # r and r+1 in same row or column <=> rtab is not standard
            return self.term(tab,denom)
        else:
            if d>0:
                coeff=1/self._quantum_integer(d)
                coeff*=self.I*self._quantum_integer_root(d+1)*self._quantum_integer_root(d-1)
            else:
                coeff=1/self._quantum_integer(-d)
                coeff*=-self.I*self._quantum_integer_root(-d+1)*self._quantum_integer_root(-d-1)
            return self.sum_of_terms([(tab,denom), (rtab,coeff)],distinct=True)

class SeminormalRepresentation_alternating(SeminormalRepresentation_generic):

    def __init__(self, shape, prefix='a'):
        r"""
        Initialisation of a seminormal representation of the symmetric group
        indexed by the partition `shape` and where the basis elements are
        labelled by `prefix`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2,1])    # indirect doctest
            SeminormalRepresentation([3,2,1])
        """
        # We need to build the base ring that splits our module. It needs to
        # contain the inverses and square roots of the quantum integers [k]_q,
        # for 0\le k\le n. To achieve this we start with the ring of Laurent
        # polynomials and successively extend this ring by throwing in the
        # square roots that we need. At the end we throw in a square of -1. We
        # don't use the ComplexField() because this seems to drag along with it
        # unwanted precision when printing.
        self._shape=shape
        self._prefix=prefix
        Fq=FractionField(PolynomialRing(IntegerRing(),'q')) 
        self.q=Fq.gen()
        roots=['r%d'%(d+1) for d in range(0,shape.size())] # names for the square roots
        roots[0]='I'
        poly_ring=PolynomialRing(Fq,roots)
        relations=[poly_ring.gen(0)**2+1]
        self.__quantum_integer={0:0, 1:1}       # dictionary of positive quantum integers
        self.__quantum_integer_root={0:0, 1:1}  # dictionary of square roots of positive quantum integers
        qint=1                                  # a quantum integer in waiting
        for k in range(1,shape.size()):
            qint+=self.q**(2*k)                 # now equal to 1+q^2+...+q^2k = [k+1]_{q^2}
            self.__quantum_integer[k+1]=qint    # we remember these because we need them later
            relations.append(poly_ring.gen(k)**2-qint)
        F=poly_ring.quotient(poly_ring.ideal(relations),roots)
        for k in range(1,shape.size()):
            self.__quantum_integer_root[k+1]=F.gen(k)
        self.I=F.gen(0)                         # I=sqrt(-1). We don't use CC as we want issues with precision
        CombinatorialFreeModule.__init__(self,F,shape.standard_tableaux(), prefix='f')

    def _quantum_integer_root(self,d):
        r"""
        Return the square root of the quantum integer `[d]_{q^2}`. These square
        roots are chosen so that `\sqrt{[-d]}=\sqrt{-1}q^d\sqrt{[d]` if `d<0`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2])._quantum_integer_root(3)
            r3
            sage: SeminormalRepresentation([3,2])._quantum_integer_root(-3)
            ((q^-3)*r3)*I

        """
        return self.__quantum_integer_root[d] if d>=0 else -self.q**(2*d)*self.__quantum_integer_root[-d]

    @cached_method
    def _T_on_basis(self,tab,r):
        r"""
        Compute the action of `T_r` on the basis element `f_tab`, where the action
        of `T_r` is as follows:

        .. MATH::

           f_t T_r = \frac{-q^{-1}}{[d]}f_t + \alpha_d f_s

        where $d=c_r(t)-c_{r+1}(t)$  is the axial distance from $r$ to $r+1$ in `t` and

        .. MATH::

           \alpha_d = \begin{cases}
                         \frac{\sqrt{-1}\sqrt{[d+1]}\sqrt{[d-1]}}{[d]},&\text{if }d>0,\\
                        -\frac{\sqrt{-1}\sqrt{[-d+1]}\sqrt{[-d-1]}}{[d]},&\text{if }d<0.
           \end{cases}

         The whole point of choosing this particular basis is that if the module
         is indexed by a self-conjugate partition then it splits when it is
         restricted to the alternating Hecke algebra.

         EXAMPLES::

             sage: rep=SeminormalRepresentation([3,2,1])
             sage: rep([[1,2,3],[4,5],[6]]).T(1)   # indirect doctest
             q*f(1,2,3/4,5/6)
             sage: rep([[1,2,3],[4,5],[6]]).T(1,2)
             q^2*f(1,2,3/4,5/6)
             sage: rep([[1,2,3],[4,5],[6]]).T(1,2,3)
             (q^2/(-q-q^3-q^5))*f(1,2,3/4,5/6) + (q^2*r2*r4*I/(1+q^2+q^4))*f(1,2,4/3,5/6)
        """
        rtab=tab.symmetric_group_action_on_entries( Permutation((r,r+1)) )
        d=tab.content(r)-tab.content(r+1)
        denom=1/(self.q**(2*d+1)*self._quantum_integer(-d))
        if abs(d)==1:  # r and r+1 in same row or column <=> rtab is not standard
            return self.term(tab,denom)
        else:
            if d>0:
                coeff=1/self._quantum_integer(d)
                coeff*=self.I*self._quantum_integer_root(d+1)*self._quantum_integer_root(d-1)
            else:
                coeff=1/self._quantum_integer(-d)
                coeff*=-self.I*self._quantum_integer_root(-d+1)*self._quantum_integer_root(-d-1)
            return self.sum_of_terms([(tab,denom), (rtab,coeff)],distinct=True)

    def tau_character(self, *arg):
        r"""
        Return the trace of `tau*T_w` acting on the representation. Here `tau`
        is the automorphsim of the Specht module that sends the basis element
        `f_t` to `f_{t'}`, where `t'` is the tableau conjugate to `t` and `w` is
        the permutation corresponding to the  (reduced) expression `*arg`.

        In fact, if `w` is a minimal length coset representative then this
        character value is zero unless `w` belongs to the 
        conjugacy class of the symmetric group that splits when restricted to
        the alternating group.

        EXAMPLES

            sage: SeminormalRepresentation([2,2]).tau_character(1)
            0
            sage: SeminormalRepresentation([2,2]).tau_character(1,2)
            (((-1 - 2*q^2 - q^4)*r3)*I)/(-q - 2*q^3 - q^5)
            sage: SeminormalRepresentation([2,2]).tau_character(2,1)
            (((1 + 2*q^2 + q^4)*r3)*I)/(-q - 2*q^3 - q^5)
        """
        if self._shape<> self._shape.conjugate():
            return 0
        trace=0
        for t in self._basis_keys:
            trace+=self(t.conjugate()).T(*arg).coefficient(t)
        return trace

class SpechtModuleWithMurphyBasis(SeminormalRepresentation_generic):

    def __init__(self, shape, prefix='m'):
        r"""
        Initialisation of a seminormal representation of the symmetric group
        indexed by the partition `shape` and where the basis elements are
        labelled by `prefix`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2,1])    # indirect doctest
            SeminormalRepresentation([3,2,1])
        """
        # We need to build the base ring that splits our module. It needs to
        # contain the inverses and square roots of the quantum integers [k]_q,
        # for 0\le k\le n. To achieve this we start with the ring of Laurent
        # polynomials and successively extend this ring by throwing in the
        # square roots that we need. At the end we throw in a square of -1. We
        # don't use the ComplexField() because this seems to drag along with it
        # unwanted precision when printing.
        self._shape=shape
        self._prefix=prefix
        self._generic_seminormalform=SeminormalRepresentation_murphy(shape)
        self._generic_Murphy=lambda t: self._generic_seminormalform().gen().T( *t.reduced_word())
        Fq=FractionField(PolynomialRing(IntegerRing(),'q')) 
        self.q=Fq.gen()
        roots=['r%d'%(d+1) for d in range(0,shape.size())] # names for the square roots
        roots[0]='I'
        poly_ring=PolynomialRing(Fq,roots)
        relations=[poly_ring.gen(0)**2+1]
        self.__quantum_integer={0:0, 1:1}       # dictionary of positive quantum integers
        self.__quantum_integer_root={0:0, 1:1}  # dictionary of square roots of positive quantum integers
        qint=1                                  # a quantum integer in waiting
        for k in range(1,shape.size()):
            qint+=self.q**(2*k)                 # now equal to 1+q^2+...+q^2k = [k+1]_{q^2}
            self.__quantum_integer[k+1]=qint    # we remember these because we need them later
            relations.append(poly_ring.gen(k)**2-qint)
        F=poly_ring.quotient(poly_ring.ideal(relations),roots)
        for k in range(1,shape.size()):
            self.__quantum_integer_root[k+1]=F.gen(k)
        self.I=F.gen(0)                         # I=sqrt(-1). We don't use CC as we want issues with precision
        CombinatorialFreeModule.__init__(self,F,shape.standard_tableaux(), prefix='m')

    @cached_method
    def _T_on_basis(self,tab,r):
        r"""
        Compute the action of `T_r` on the basis element `f_tab`, where the action
        of `T_r` is as follows:

        .. MATH::

           f_t T_r = \frac{-q^{-1}}{[d]}f_t + \alpha_d f_s

        where $d=c_r(t)-c_{r+1}(t)$  is the axial distance from $r$ to $r+1$ in `t` and

        .. MATH::

           \alpha_d = \begin{cases}
                         \frac{\sqrt{-1}\sqrt{[d+1]}\sqrt{[d-1]}}{[d]},&\text{if }d>0,\\
                        -\frac{\sqrt{-1}\sqrt{[-d+1]}\sqrt{[-d-1]}}{[d]},&\text{if }d<0.
           \end{cases}

         The whole point of choosing this particular basis is that if the module
         is indexed by a self-conjugate partition then it splits when it is
         restricted to the alternating Hecke algebra.

         EXAMPLES::

             sage: rep=SeminormalRepresentation([3,2,1])
             sage: rep([[1,2,3],[4,5],[6]]).T(1)   # indirect doctest
             q*f(1,2,3/4,5/6)
             sage: rep([[1,2,3],[4,5],[6]]).T(1,2)
             q^2*f(1,2,3/4,5/6)
             sage: rep([[1,2,3],[4,5],[6]]).T(1,2,3)
             (q^2/(-q-q^3-q^5))*f(1,2,3/4,5/6) + (q^2*r2*r4*I/(1+q^2+q^4))*f(1,2,4/3,5/6)
        """
        rtab=tab.symmetric_group_action_on_entries( Permutation((r,r+1)) )
        d=tab.content(r)-tab.content(r+1)
        if rtab.is_standard():
            if d>0: 
                # tab dominates rtab => m_tab T_r = m_{rtab}
                return self.monomial(rtab)
            else:
                # rtab dominates tab => m_tab T_r = (q_1+q_2)m_{rtab} + q_1q_2 m_tab
                return self.sum_of_terms([(tab,self.parent()._q_prod), (rtab,rtab,self.parent()._qsum)],distinct=True)
        else:
           # r and r+1 are in the same column so we compute via the generic
           # seminormal form using the fact the transition matrix between the
           # Murphy basis the SeminormalRepresentation_murphy basis is unitriangular
           specht=self.parent()
           mtabr=specht.zero()
           generic_prod=specht._generic_Murphy(t)
           while generic_prod != self._generic_seminormalform().zero():
                (s,c) = generic_prod.trailing_item(index_cmp)
                generic_prod -= c * specht._generic_Murphy(s)
                mtabr += specht._specialise_coeff(c) * self.parent().monomial(s)

           return mtabr


