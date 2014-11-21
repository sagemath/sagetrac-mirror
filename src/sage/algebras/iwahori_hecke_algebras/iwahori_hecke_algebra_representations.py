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
#    Copyright (C) 2014 Andrew Mathas <andrew mathas at sydney edu au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at: http://www.gnu.org/licenses/.
#
#*****************************************************************************


from sage.algebras.algebra import is_Algebra
from sage.categories.action import Action
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.partition_tuple import PartitionTuple
from sage.combinat.permutation import Permutation
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.rings.fraction_field import FractionField
from sage.rings.integer_ring import IntegerRing
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import RationalField


import operator

class IwahoriHeckeAlgebraRepresentations(object):
    r"""
    This is a shortcut that returns the possible representations of this
    algebra.

    EXAMPLES::

        sage: H=IwahoriHeckeAlgebra("A4",1)
        sage: H.representations.<tab>             # not tested - tab-completion
        H.representations.LeftCellRepresentation       H.representations.SeminormalForm_Murphy
        H.representations.RightCellRepresentation      H.representations.SeminormalForm_Orthogonal
        H.representations.SeminormalForm               H.representations.SpechtModuleWithMurphyBasis
        sage: H.representations.SeminormalForm([3,2])

    """
    pass

class HeckeAlgebraRepresentation(CombinatorialFreeModule):
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
    def __init__(self, base_ring, basis_keys, prefix, q1, q2, cartan_type=None, module_action="left"):

        # check that q1 and q2 belong to the base ring
        if not (q1 in base_ring and q2 in base_ring):
            raise ValueError('q1={}, q2={} must be in the same ring as charge={}'.format(q1,q2,self._charge))

        # Record the Hecke algebra parameters
        self._q1=q1
        self._q2=q2

        # This gives a speed-up as it avoids the need to constantly add and
        # multiply the parameters when applying the quadratic relation: 
        #   T^2 = (q1+q2)T - q1*q2
        self._q_sum  = base_ring(q1+q2)
        self._q_prod = base_ring(-q1*q2)

        # in order to allow the corresponding Hecke algebra to act directly on
        # this representation we need to know when something belongs to this algebra
        try:
            from sage.combinat.root_system.cartan_type import CartanType
            cartan_type=CartanType(cartan_type)
            from sage.algebras.iwahori_hecke_algebras.iwahori_hecke_algebra import IwahoriHeckeAlgebra
            self._algebra=IwahoriHeckeAlgebra(cartan_type, base_ring(q1), base_ring(q2))
            self._cartan_type=cartan_type
        except TypeError:
            # if this is not a Hecke algebra rep then we still need self._algebra
            # to be an iterable object so that `a in self._algebra` is valid
            # syntax (and returns false)
            self._algebra=[]

        if module_action=="left":
            self.Element.T=self.Element._T_left_action
        elif module_action=="right":
            self.Element.T=self.Element._T_right_action
        elif module_action=="bimodule": # TODO: what do we need to do to allow this?
            self.Element.T_left=self.Element._T_left_action
            self.Element.T_right=self.Element._T_right_action
        else:
            raise ValueError('module_action must be left, right or bimodule!')
 
        print 'Prefix = {}, basis keys = {}'.format(prefix, basis_keys)
        super(HeckeAlgebraRepresentation,self).__init__(
                base_ring, basis_keys, prefix=prefix
        )

    ## compute the action of the Iwahori-Hecke algebra generators ---------------------------

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

    @abstract_method
    def _T_on_basis(self, b, r):
        r"""
        Returns the result of `T_r` acting on `b`.
        """
        pass

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

    ## Hecke algebra action ---------------------------

    def hecke_algebra_action(self, H):
        r"""
        Return the action of the Hecke algebra `H` upon `self`.
        """
        return IwahoriHeckeAlgebraAction(H,self)

    def _get_action_(self, H, op, self_on_left):
        r"""

        TODO: bimodule action!

        """
        print 'H=%s, op=%s, self_on_left=%s'%(H,op,self_on_left)
        if is_Algebra(H) and op==operator.mul and self_on_left:
            return self.hecke_algebra_action(H)
        return None

    def action(self,a):
        if not a in self._algebra:
            raise ValueError('%{a} does not act on {rep}'.format(a=a,rep=self))

        return self.parent().sum(c*self.T( *w.reduced_word() ) for (w,c) in a.to_T_basis())

    ## check that the algebra's relations are satisfied on the representation ---------------------------

    def check_relations(self, verbose=False):
        r"""
        Return `True` if all of the relations are satisfied on `self` and
        `False` otherwise.

        EXAMPLES::

            sage: SeminormalRepresentation([2,1]).check_relations()
            True
        """
        for t in self._basis_keys:
            if verbose: print 'Checking relations on ({term})'.format(term=self(t))
            for r in range(1,self.size()):
                # check quadratic relations on self(t)
                assert self(t).T(r,r)==self._q_sum*self(t).T(r)+self._q_prod*self(t), (
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
        r"""
        Element methods.
        """

        def _Tr(self,r):
            r"""
            This is a helper function for computing the action of `T_r` upon the
            element `self` of the seminormal representation.

            EXAMPLES::

                sage: SeminormalRepresentation([2,2]).an_element().T(1)
                2*q*f(1,2/3,4) + (2/(-q))*f(1,3/2,4)
            """

            return self.parent().sum(coeff*self.parent()._T_on_basis(b,r) for (b,coeff) in self)

        def _T_left_action(self,*args):
            r"""
            Compute the action of `T_w` on self, where `w` is the permutation
            corresponding to the (reduced) word `*args`.

            EXAMPLES::

                sage: SeminormalRepresentation([2,2]).an_element().T(1)
                2*q*f(1,2/3,4) + (2/(-q))*f(1,3/2,4)
            """
            v=self
            for r in args[::-1]:
                v=v._Tr(r)
            return v

        def _T_right_action(self,*args):
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


##################################################################
##  Kazhdan-Lusztig cell representations of Iwahori-Hecke algebras
##################################################################

class LeftCellRepresentationOfHeckeAlgebra(HeckeAlgebraRepresentation):
    pass

class RightCellRepresentationOfHeckeAlgebra(HeckeAlgebraRepresentation):
    pass

############################################################
## Seminormal representations of (cyclotomic) Hecke algebras
############################################################

# class degenerate_or_nondegenerate(object):
#     '''
#     A custom decorator that, on the first use, replaces a method function with the
#     appropriate version depending on the degeneracy and the level. This is a
#     little like :meth:`lazy_attribute` except that it defines a method rather
#     than an attribute.
# 
#      EXAMPLES::
#
#          sage: class A(SeminormalRepresentation_generic):
#          ....:     @degenerate_or_nondegenerate
#          ....:     def _some_method(self, *args): pass
#          ....:     def _degenerate_some_method(self, *args): pass
#          ....:     def _nondegenerate_some_method(self, *args): pass
#          ....:     def _degenerate_some_method_level_one(self, *args): pass
#          ....:     def _nondegenerate_some_method_level_one(self, *args): pass
#
#     When the `_some_method` method is first called this method well be set to
#     the appropriate one of the four methods below and the value computed and
#     returned. All subsequent calls now use the appropriate method.
#
#     In the end I decided to do this a different way, but I'm leaving the code
#     here as this approach might be useful elsewhere.
#     '''
# 
#     def __init__(self,func):
#         self.func =func
#         self.func_name = func.__name__
# 
#     def __get__(self,seminormal,cls):
#         if seminormal is None:
#             return None
#         function_name='_{}degenerate{}{}'.format('' if seminormal._degenerate else 'non',
#              self.func_name, '_level_one' if seminormal._shape.level()==1 else '')
#         try:
#             real_func=getattr(seminormal, function_name)
#         except:
#             raise AttributeError('{} has no attribute {}'.format(seminormal, function_name))
#         setattr(seminormal,self.func_name,real_func)
#         return real_func

class SeminormalRepresentation_generic(HeckeAlgebraRepresentation):
    r"""
    This is a private class that implements the broad skeleton of the seminormal
    representations of the (cyclotomic) Hecke algebras of type A. The derived 
    classes of this class need to::

    - define the `base_ring` for the representation and set q1 and q1
    - define the `charge`, which determines the eigenvalues of $L_1$.
    - define a `prefix` for printing the basis elements
    - implement the four methods::
        - _rho_and_beta_degenerate
        - _rho_and_beta_degenerate_level_one
        - _rho_and_beta_nondegenerate
        - _rho_and_beta_nondegenerate_level_one

    """
    @staticmethod
    def __classcall__(cls, shape, q1=None, q2=-1, charge=None, degenerate=False, prefix='?', **kwargs):
        """
        Magic to allow class to accept a list (which is not hashable) instead
        of a :class:`PartitionTuple`, which apparently is.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2,1])    # indirect doctest
            SeminormalRepresentation([3,2,1])
        """
        try:
            shape=PartitionTuple(shape)
            if charge is not None and len(charge)!=shape.level():
                raise ValueError('the shape and charge have different lengths')
        except (TypeError, ValueError):
             raise ValueError('the shape must be a partition or partition tuple')

        return super(SeminormalRepresentation_generic, cls).__classcall__(cls, shape=PartitionTuple(shape),
                           q1=q1, q2=q2, charge=charge, degenerate=degenerate, **kwargs)

    def __init__(self, shape, q1, q2, prefix, degenerate, base_ring, **kwargs):
        r"""
        Initialisation of a seminormal representation of the symmetric group
        indexed by the partition `shape` and with the basis elements labelled by
        `prefix`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2,1])    # indirect doctest
            SeminormalRepresentation([3,2,1])
        """
        self._shape=shape            # shape of the partition defining the representation
        self._degenerate=degenerate  # true/false depending on whether we use a
                                     # degenerate or non-degenerate seminormal form.

        # Iwahori-Hecke algebras are currently only implemented in type A (level 1)
        # and in type B (level 2) in the equal parameter case
        level=self._shape.level()
        from sage.algebras.iwahori_hecke_algebras.iwahori_hecke_algebra import IwahoriHeckeAlgebra
        if level==1:
            cartan_type=['A', shape.size()]
        elif level==2 and kwargs['q1']==charge[0] and kwargs['q2']==charge[1]:
            cartan_type=['B', shape.size()]
        else:
            cartan_type=['G', shape.level(), 1, shape.size()]

        # set the defaults for the methods that depend on degeneracy and level
        for func in ['rho_and_beta','tableau_content']:
            real_func='_{}_{}degenerate{}'.format(func, '' if degenerate else 'non',
                                            '_level_one' if shape.level()==1 else '')
            try:
                setattr(self,'_'+func,getattr(self,real_func))
            except:
                raise AttributeError('{} has no attribute {}'.format(self, real_func))

        super(SeminormalRepresentation_generic,self).__init__(
                    basis_keys=shape.standard_tableaux(), # indexing set for basis
                    base_ring=base_ring,                  # coefficient field for representation
                    q1=q1, q2=q2,                         # parameters
                    prefix=prefix,                        # prefix for printing basis
                    cartan_type=cartan_type,              # Cartan type of representation
                    **kwargs                              # anything used upstream
        )

    def check_relations(self, verbose):
        r"""
        Return `True` if the action of the algebra of the representation `self`
        repects all of the algebras relations and `False` otherwise.

        """
        # we only need to check the relations involving L_1,...,L_n are the
        # remaining relations are verified by HeckeAlgebraRepresentation.
        for t in self._basis_keys:
            if verbose: print 'Checking cyclotomic relations on ({term})'.format(term=self(t))
            for r in range(1,self.size()+1):
                # check quadratic relations on self(t)
                assert self(t).L(r)==self._tableau_content(t,r)*self(t), (
                        'L_{r} is not acting with the expected eigenvalue on {term}'.format(r=r,term=self(t))
                )
        return super(SeminormalRepresentation_generic,self).check_relations(verbose)

    def shape(self):
        r"""
        Return the partition that indexes the seminormal representation `self`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2]).shape()
            [3, 2]
        """
        return self._shape

    @lazy_attribute
    def size(self):
        r"""
        Return the size of the partition that indexes the seminormal representation `self`.

        EXAMPLES::

        sage: SeminormalRepresentation([3,2]).size()
        5
        """
        return self._shape.size()

    def charge(self):
        return self._charge


    @lazy_attribute
    def generator(self):
        r"""
        Returna generator, as a module for the Hecke algebra, of the Specht module `self`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2]).generator()
        """
        return self.monomial(self._shape.initial_tableau())

    @lazy_attribute
    def is_degenerate(self):
        r"""
        Return `True` if `self` is a degenerate seminormal representation and
        `False` otherwise.
        """
        return self._degenerate

    def _repr_(self):
        r"""
        Return a string representation of the seminormal representation `self`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2,1])    # indirect doctest
            SeminormalRepresentation([3,2,1])
        """
        if self._degenerate:
            return 'SeminormalRepresentation([{}], degenerate=True)'.format(self._shape._repr_compact_high())
        else:
            return 'SeminormalRepresentation([{}])'.format(self._shape._repr_compact_high())

    def _repr_term(self,tab):
        r"""
        Override the default printing of terms so that compact tableaux are
        used.

        EXAMPLES::

            sage: SeminormalRepresentation([2,1]).an_element()  # indirect doctest
            2*f(1,2/3) + 2*f(1,3/2)
            sage: SeminormalRepresentation([2,1], prefix='v').an_element()
            2*v(1,2/3) + 2*v(1,3/2)
        """
        return '{prefix}({tab})'.format(prefix=self._prefix, tab=tab._repr_compact() )

    def _element_constructor_(self,t):
        r"""
        We override the fault element constructor in order to coerce lists into standard tableaux.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2])([[1,2,3],[4,5]])
            f(1,2,3/4,5)
        """
        return self.monomial(self._basis_keys(t))

    def _tableau_content_degenerate(self, t, r):
        r"""
        Return the content of the integer `r` in the tableau `t` in the degenerate case.
        """
        c=t.content(r)
        return self._quantum_integer(self._charge[c[0]]+c[1])

    def _tableau_content_degenerate_level_one(self, t, r):
        r"""
        Return the content of the integer `r` in the tableau `t` in the degenerate case.
        """
        return self._quantum_integer(t.content(r))

    def _tableau_content_nondegenerate(self, t, r):
        r"""
        Return the content of the integer `r` in the tableau `t` in the degenerate case.
        """
        c=t.content(r)
        return self._charge[c[0]]*self._q_prod**c[1]

    def _tableau_content_nondegenerate_level_one(self, t, r):
        r"""
        Return the content of the integer `r` in the tableau `t` in the degenerate case.
        """
        return self._q_prod**t.content(r)

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
        c = tab.content(r)     # content of r
        d = tab.content(r+1)   # content of r+1
        rho, beta = self._rho_and_beta(c, d)
        if rtab.is_standard():  # r and r+1 in different row or column <=> rtab is standard
            return self.sum_of_terms([(tab,rho), (rtab,beta)],distinct=True)
        else:                   # r and r+1 in same different row or column <=> rtab is not standard
            return self.term(tab,rho)

    class Element(HeckeAlgebraRepresentation.Element):

        def L(self, *arg):
            r"""
            Return the action of `L(r_1, r_2,...)` on the module element `self`.
            """
            return self.parent().sum_of_terms(
                       [(tab,c*prod([self._tableau_content(tab,r) for r in arg])) for tab,c in self],
                       distinct=True)


class SeminormalRepresentation(SeminormalRepresentation_generic):
    def __init__(self, shape, q1, q2, charge, degenerate, **kwargs):
        r"""
        Initialisation of a seminormal representation of the symmetric group
        indexed by the partition `shape` and where the basis elements are
        labelled by `prefix`.

        EXAMPLES::

            sage: SeminormalRepresentation([3,2,1])    # indirect doctest
            SeminormalRepresentation([3,2,1])
        """
        if q1 is None:
            if q2 is None:
                # the generic relation (T_r-u)(T_r-v)=0
                F=FractionField(PolynomialRing(IntegerRing(),'u,v'))
                q1,q2=F.gens()
            else:
                # the default parameters for the relation (T_r-q)(T_r+1)=0
                F=FractionField(PolynomialRing(IntegerRing(),'q'))
                q1=F.gen()
        else:
            F=FractionField(q1.parent())
            if q2 not in F and q1 in q2.parent():
                F=FractionField(q2.parent())
            else:
                raise ValueError('q1={} and q2={} do not belong to a common ring'.format(q1,q2))

        if charge is None:
            if shape.level()==1:
                self._charge=[1]
                base_ring=F
            else:
                base_ring=PolynomialRing(F,['Q%s'%(r+1) for r in range(shape.level())])
                self._charge=base_ring.gens()
        else:
            self._charge=charge
            base_ring=self._charge[-1].parent()

        super(SeminormalRepresentation,self).__init__(
                shape=shape,                       # for creating the indexing set of the basis
                q1=q1, q2=q2, base_ring=base_ring, # set the representations parameters upstream
                prefix='v',                        # seminormal basis prints as v(...)
                degenerate=degenerate,
                **kwargs
        )

        self.__quantum_integers={0:0} # quick look up table for quantum integers
        if degenerate:
            self.__quantum_integers={0:0} # quick look up table for quantum integers
            for k in range(shape.size()):
                self.__quantum_integers[k+1]=self._q2+self._q2/self._q1*self.__quantum_integers[k]
                self.__quantum_integers[-k-1]=self._q1+self._q1/self._q2*self.__quantum_integers[-k]
        else:
            self.__quantum_integers={0:1, 1: self._q2, -1: self._q1} # quick look up table for quantum integers
            for k in range(1,shape.size()):
                self.__quantum_integers[k+1]=self._q2/self._q1*self.__quantum_integers[k]
                self.__quantum_integers[-k-1]=self._q1/self._q2*self.__quantum_integers[-k]

    def _rho_and_beta_degenerate_level_one(self, c, d):
        return -1/self.__quantum_integer(c-d), self.__quantum_integer(c-d+1)/self.__quantum_integer(c-d)

    def _rho_and_beta_nondegenerate_level_one(self, c, d):
        rho  = (self._q1-self._q2)*self._q1**c/(self._q1**c-self._q2**d)
        beta = ( self._q1**(c+1)-self._q2**d)/(self._q1**c-self._q2**d)
        return rho, beta

    def _rho_and_beta_degenerate(self, c, d):
        return -1/self.__quantum_integer(c-d), self.__quantum_integer(c-d+1)/self.__quantum_integer(c-d)

    def _rho_and_beta_nondegenerate(self, c, d):
        rho  = (self._q1-self._q2)*self._q1**c/(self._q1**c-self._q2**d)
        beta = ( self._q1**(c+1)-self._q2**d)/(self._q1**c-self._q2**d)
        return rho, beta


class SeminormalRepresentation_Orthogonal(SeminormalRepresentation_generic):
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

    def _rho_and_beta_degenerate_level_one(self, c, d):
        return -1/self.__quantum_integer(c-d), self.__quantum_integer_root(c-d+1)/self.__quantum_integer(c-d)

    def _rho_and_beta_nondegenerate_level_one(self, c, d):
        rho  = (self._q1-self._q2)*self._q1**c/(self.__q1**c-self._q2**d)
        beta = ( self._q1**(c+1)-self._q2**d)/(self.__q1**c-self._q2**d)
        return rho, beta

    def _rho_and_beta_degenerate(self, c, d):
        return -1/self.__quantum_integer(c-d), self.__quantum_integer(c-d+1)/self.__quantum_integer(c-d)

    def _rho_and_beta_nondegenerate(self, c, d):
        rho  = (self._q1-self._q2)*self._q1**c/(self.__q1**c-self._q2**d)
        beta = ( self._q1**(c+1)-self._q2**d)/(self.__q1**c-self._q2**d)
        return rho, beta

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

class SeminormalRepresentation_Murphy(SeminormalRepresentation_generic):
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

    def _rho_and_beta_degenerate_level_one(self, c, d):
        if c>d:
            return -1/self.__quantum_integer(c-d), 1
        else:
            beta= self.__quantum_integer_root(c-d+1)*self.__quantum_integer_root(c-d+1)/self.__quantum_integer(c-d)**2
            return -1/self.__quantum_integer(c-d), beta

    def _rho_and_beta_nondegenerate_level_one(self, c, d):
        rho  = (self._q1-self._q2)*self._q1**c/(self.__q1**c-self._q2**d)
        beta = ( self._q1**(c+1)-self._q2**d)/(self.__q1**c-self._q2**d)
        return rho, beta

    def _rho_and_beta_degenerate(self, c, d):
        return -1/self.__quantum_integer(c-d), self.__quantum_integer(c-d+1)/self.__quantum_integer(c-d)

    def _rho_and_beta_nondegenerate(self, c, d):
        rho  = (self._q1-self._q2)*self._q1**c/(self.__q1**c-self._q2**d)
        beta = ( self._q1**(c+1)-self._q2**d)/(self.__q1**c-self._q2**d)
        return rho, beta

    @cached_method
    def _murphy(self,tab):
        r"""
        Return the Murphy basis element indexed by the tableau `tab` as a linear
        combination of seminormal basis elements.
        """
        initial_tableau=self._shape.initial_tableau()
        if tab==initial_tableau: 
            return self.monomial(initial_tableau)
        else:
            r=tab.reduced_row_word()[-1]
            s=tab.symmetric_group_action_on_entries( Permutation((r,r+1)) )
            return self._murphy(s).T(r)

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
    def _seminormal(self, tab):
        r"""
        Return the seminormal basis element indexed by the tableau `tab` as a
        linearcombination of Murphy basis elements.

        EXAMPES::

        TODO: problematic because, in general, the coefficients will not belong
        to the current base-ring.
        """
        pass

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

class SeminormalRepresentation_Alternating(SeminormalRepresentation_generic):
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

    def _rho_and_beta_degenerate_level_one(self, c, d):
        rho=1/(self.q**(2*d+1)*self._quantum_integer(-d))
        if d>0:
            beta=1/self._quantum_integer(d)
            beta*=self.I*self._quantum_integer_root(d+1)*self._quantum_integer_root(d-1)
        else:
            beta=1/self._quantum_integer(-d)
            beta*=-self.I*self._quantum_integer_root(-d+1)*self._quantum_integer_root(-d-1)
        return rho, beta

    def _rho_and_beta_degenerate(self, c, d):
        pass

    def _rho_and_beta_nondegenerate_level_one(self, c, d):
        pass

    def _rho_and_beta_nondegenerate(self, c, d):
        pass

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
