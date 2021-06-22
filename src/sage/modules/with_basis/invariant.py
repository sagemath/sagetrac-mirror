r"""
Invariant modules
"""

# ****************************************************************************
#       Copyright (C) 2021 Trevor K. Karn <karnx018 at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

## TODO: COMB THORUGH TO MAKE SURE ALL STUFF IS HAPPENING IN REPN NOT IN MODULE

from sage.modules.with_basis.subquotient import SubmoduleWithBasis
from sage.categories.finitely_generated_semigroups import FinitelyGeneratedSemigroups
from sage.categories.finite_dimensional_modules_with_basis import FiniteDimensionalModulesWithBasis
from sage.categories.groups import Groups
from sage.sets.family import Family

class FiniteDimensionalInvariantModule(SubmoduleWithBasis):
    r"""
    Construct the `S`-invariant submodule of `M`. When a semigroup `S` acts on a module
    `M`, the invariant module is the collection of elements `m` in `M` such that
    `s \cdot m = m` for all `s \in S.

    ..MATH::

        M^S = \{m \in M : s\cdot m = m,\, \forall s \in S \}

    INPUTS:

    - ``R`` -- an instance of a ``Representation`` of a semigroup `S`
               acting on the module `M`.

    OUTPUTS:

    - ``MS`` -- the invariant algebra of the semigroup action of `S` on `M`, or
                equivalently, the isotypic component of the representation of
                `S` carried by `M` corresponding to the trivial character.

    EXAMPLES::

        sage: G = CyclicPermutationGroup(3)
        sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
        sage: from sage.modules.with_basis.representation import Representation
        sage: action = lambda g, m: M.term(g(m)) #cyclically permute coordinates
        sage: R = Representation(G, M, action)
        sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
        sage: I = FiniteDimensionalInvariantModule(R)
        sage: [I.lift(b) for b in I.basis()]
        [M[1] + M[2] + M[3]]

        sage: G = CyclicPermutationGroup(3)
        sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
        sage: from sage.modules.with_basis.representation import Representation
        sage: action = lambda g, m: M.term(g(m)) #cyclically permute coordinates
        sage: R = Representation(G, M, action, side='right') #same as last but on right
        sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
        sage: g = G.an_element(); g
        (1,2,3)
        sage: r = R.an_element(); r
        2*M[1] + 2*M[2] + 3*M[3]
        sage: R.side()
        'right'
        sage: r*g
        3*M[1] + 2*M[2] + 2*M[3]
        sage: I = FiniteDimensionalInvariantModule(R)
        sage: [I.lift(b) for b in I.basis()]
        [M[1] + M[2] + M[3]]

        sage: G = SymmetricGroup(3)
        sage: R = G.regular_representation(QQ)
        sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
        sage: I = FiniteDimensionalInvariantModule(R)
        sage: [I.lift(b).to_vector() for b in I.basis()]
        [(1, 1, 1, 1, 1, 1)]

    .. NOTE:: 

        The current implementation works when `S` is a finitely-generated semigroup,
        and when `M` is a finite-dimensional free module with a distinguished basis.

    .. TODO::

        Extend when `M` does not have a basis and `S` is a permutation group using:
        - https://arxiv.org/abs/0812.3082
        - https://www.dmtcs.org/pdfpapers/dmAA0123.pdf

    """

    def __init__(self, R, *args, **kwargs):
        """
        TESTS::
            
            sage: G = QQ # a group that is not finitely generated
            sage: M = CombinatorialFreeModule(QQ, [1,2,3])
            sage: on_basis = lambda g,m: M.term(m) # trivial rep'n
            sage: from sage.modules.with_basis.representation import Representation
            sage: R = Representation(G,M,on_basis)
            sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
            sage: I = FiniteDimensionalInvariantModule(R)
            Traceback (most recent call last):
            ...
            ValueError: Rational Field is not finitely generated

            sage: G = SymmetricGroup(3)
            sage: M = PolynomialRing(QQ, ['x1','x2','x3'])
            sage: on_basis = lambda g,m: M.prod(M('x'+str(g(Integer(x[-1])))) for x in m.support())
            sage: f = M('x1')*M('x3')
            sage: g = G.an_element(); g
            (2,3)
            sage: on_basis(g,f)
            x1*x2
            sage: R = Representation(G, M, on_basis)
            sage: I = FiniteDimensionalInvariantModule(R)
            Traceback (most recent call last):
            ...
            ValueError: Multivariate Polynomial Ring in x1, x2, x3
             over Rational Field is not finite-dimensional

            sage: #from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
            sage: #I
            sage: #TestSuite(I).run()
        """
        self._semigroup_representation = R
        self._semigroup = R.semigroup()

        if self._semigroup not in FinitelyGeneratedSemigroups:
            raise ValueError(f'{self._semigroup} is not finitely generated')

        if self._semigroup_representation._module not in (
            FiniteDimensionalModulesWithBasis):
            raise ValueError(
                f'{self._semigroup_representation._module
                } is not finite dimensional')

        # The left/right multiplication is taken care of
        # by self._semigroup_representation, so here
        # we can just pass the left multiplication.
        # This means that the side argument of annihilator_basis
        # (see below) will always be side = 'left'
        if self._semigroup_representation.side() == 'left':
            def _invariant_map(g, x):
                return g*x - x
        elif self._semigroup_representation.side() == 'right':
            def _invariant_map(g, x):
                return x*g - x
        
        self._invariant_map = _invariant_map
        
        category = kwargs.pop('category', R.category().Subobjects())

        # Give the intersection of kernels of the map `s*x-x` to determine when
        # `s*x = x` for all generators `s` of `S`
        basis = self._semigroup_representation.annihilator_basis(
                    self._semigroup.gens(),
                    action = self._invariant_map,
                    side = 'left')

        super().__init__(Family(basis),
                        support_order = self._semigroup_representation._compute_support_order(basis),
                        ambient = self._semigroup_representation,
                        unitriangular = False,#is this right?
                        category = category,
                        *args, **kwargs)

    def _repr_(self):
        """
        EXAMPLES::

            sage: M = CombinatorialFreeModule(QQ,[1,2,3])
            sage: G = CyclicPermutationGroup(3)
            sage: from sage.modules.with_basis.representation import Representation
            sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
            sage: R = Representation(G,M,lambda g,x: M.monomial(g(x)))
            sage: FiniteDimensionalInvariantModule(R)
            (Cyclic group of order 3 as a permutation group)-invariant submodule of
             Free module generated by {1, 2, 3} over Rational Field

        """

        return f"({self._semigroup})-invariant submodule of {self._semigroup_representation._module}"


    def semigroup(self):
        """
        Return the semigroup `S` whose action ``self`` is invariant under.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
            sage: action = lambda g,x: M.term(g(x))
            sage: from sage.modules.with_basis.representation import Representation
            sage: R = Representation(G, M, action)
            sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
            sage: I = FiniteDimensionalInvariantModule(R)
            sage: I.semigroup()
            Symmetric group of order 3! as a permutation group

        """

        return self._semigroup

    def action_on_basis(self,*args):
        if len(args) == 0:
            return self._action_on_basis

        return self._action_on_basis(args)

    # def _test_

    class Element(SubmoduleWithBasis.Element):

        def _mul_(self, other):
            P = self.parent()
            return P.retract(P.lift(self) * P.lift(other))

        # lmul -- self*right
        def _lmul_(self, right):
            """
            EXAMPLES::

                sage:        #### create invariant module
                sage:        #### multiply group element
            """

            if right in self.parent()._semigroup and self.parent()._side == 'right':
                return self

            return super()._lmul_(right)

        # rmul -- left * self
        def _rmul_(self, left):
            """
            EXAMPLES::

                sage: ### create invariant module
                sage: ### multiply group element
            """
            if left in self.parent()._semigroup and self.parent()._side == 'left':
                return self

            return super()._rmul_(left)

        def _acted_upon_(self, scalar, self_on_left = False):
            """
            EXAMPLES::

            """

            if scalar in self.parent()._semigroup and self_on_left == (self.parent()._side == 'right'):

                return self

            return None


# Group action should be lifted on the module 

class FiniteDimensionalTwistedInvariantModule(SubmoduleWithBasis):
    r"""
    Construct the `\chi`-twisted invariant submodule of `M`. When a semigroup `S` acts on a module
    `M`, the `\chi`-twisted invariant submodule of `M` is the isotypic component of the representation
    `M` corresponding to the irreducible character `\chi`.

    .. MATH::

        M^S = \{m \in M : s\cdot m = m,\, \forall s \in S \}

    NOTE: The current implementation works when `S` is a finitely-generated semigroup,
    and when `M` is a finite-dimensional free module with a distinguished basis.
    """

    def __init__(self, R, character = 'trivial'):

        super.__init__(R)

        if character != 'trivial':

            pass


    def projection(self, element):
        """
        Give the projection of element (in `self.module()`) onto self
        """
        pass

    # class Element

    #     _lmul_
    #     _rmul_
    #     _acted_upon_