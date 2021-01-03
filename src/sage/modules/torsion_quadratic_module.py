r"""
Finite `\ZZ`-modules with with bilinear and quadratic forms.

AUTHORS:

- Simon Brandhorst (2017-09): First created
"""

# ****************************************************************************
#       Copyright (C) 2017 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.modules.fg_pid.fgp_module import FGP_Module_class
from sage.modules.fg_pid.fgp_element import FGP_Element
from sage.modules.free_quadratic_module import FreeQuadraticModule
from sage.modules.free_quadratic_module_integer_symmetric import FreeQuadraticModule_integer_symmetric
from sage.arith.misc import gcd
from sage.rings.all import ZZ, Zp, QQ, IntegerModRing
from sage.groups.additive_abelian.qmodnz import QmodnZ
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix
from sage.misc.cachefunc import cached_method
from sage.arith.misc import legendre_symbol
from sage.structure.unique_representation import CachedRepresentation
from sage.rings.finite_rings.integer_mod import mod
from sage.libs.gap.libgap import libgap

def TorsionQuadraticForm(q):
    r"""
    Create a torsion quadratic form module from a rational matrix.

    The resulting quadratic form takes values in `\QQ / \ZZ`
    or `\QQ / 2 \ZZ` (depending on ``q``).
    If it takes values modulo `2`, then it is non-degenerate.
    In any case the bilinear form is non-degenerate.

    INPUT:

    - ``q`` -- a symmetric rational matrix

    EXAMPLES::

        sage: q1 = Matrix(QQ,2,[1,1/2,1/2,1])
        sage: TorsionQuadraticForm(q1)
        Finite quadratic module over Integer Ring with invariants (2, 2)
        Gram matrix of the quadratic form with values in Q/2Z:
        [  1 1/2]
        [1/2   1]

    In the following example the quadratic form is degenerate.
    But the bilinear form is still non-degenerate::

        sage: q2 = diagonal_matrix(QQ,[1/4,1/3])
        sage: TorsionQuadraticForm(q2)
        Finite quadratic module over Integer Ring with invariants (12,)
        Gram matrix of the quadratic form with values in Q/Z:
        [7/12]

    TESTS::

        sage: TorsionQuadraticForm(matrix.diagonal([3/8,3/8,3/4]))
        Finite quadratic module over Integer Ring with invariants (4, 8, 8)
        Gram matrix of the quadratic form with values in Q/2Z:
        [3/4   0   0]
        [  0 3/8   0]
        [  0   0 3/8]
    """
    q = matrix(QQ, q)
    if q.nrows() != q.ncols():
        raise ValueError("the input must be a square matrix")
    if q != q.transpose():
        raise ValueError("the input must be a symmetric matrix")

    Q, d = q._clear_denom()
    S, U, V = Q.smith_form()
    D = U * q * V
    Q = FreeQuadraticModule(ZZ, q.ncols(), inner_product_matrix=d**2 * q)
    denoms = [D[i, i].denominator() for i in range(D.ncols())]
    rels = Q.span(diagonal_matrix(ZZ, denoms) * U)
    return TorsionQuadraticModule((1/d)*Q, (1/d)*rels, modulus=1)


class TorsionQuadraticModuleElement(FGP_Element):
    r"""
    An element of a torsion quadratic module.

    INPUT:

    - ``parent`` -- parent

    - ``x`` -- element of ``parent.V()``

    - ``check`` -- bool (default: ``True``)

    TESTS::

        sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
        sage: T = TorsionQuadraticModule(ZZ^3, 6*ZZ^3)
        sage: loads(dumps(T)) == T
        True
        sage: t = T.gen(0)
        sage: loads(dumps(t)) == t
        True

        sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
        sage: V = span([[1/2,1,1], [3/2,2,1], [0,0,1]], ZZ)
        sage: b = V.basis()
        sage: W = V.span([2*b[0]+4*b[1], 9*b[0]+12*b[1], 4*b[2]])
        sage: Q = TorsionQuadraticModule(V, W)
        sage: x = Q(b[0] - b[1])
        sage: TestSuite(x).run()
        """

    def _mul_(self, other):
        r"""
        Compute the inner product of two elements.

        OUTPUT:

        - an element of `\QQ / m\ZZ` with `m\ZZ = (V, W)`

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = (1/2)*ZZ^2; W = ZZ^2
            sage: T = TorsionQuadraticModule(V, W)
            sage: g = T.gens()
            sage: x = g[0]
            sage: y = g[0] + g[1]
            sage: x
            (1, 0)
            sage: x*y
            1/4

        The inner product has further aliases::

            sage: x.inner_product(y)
            1/4
            sage: x.b(y)
            1/4
        """
        value_module = self.parent().value_module()
        return value_module( self.lift().inner_product(other.lift()) )

    inner_product = _mul_
    b = _mul_

    def quadratic_product(self):
        r"""
        Compute the quadratic_product of ``self``.

        OUTPUT:

        - an element of `\QQ / n\ZZ` where `n\ZZ = 2(V,W) +
          \ZZ \{ (w,w) | w \in W \}`

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: W = FreeQuadraticModule(ZZ, 2, 2*matrix.identity(2))
            sage: V = (1/2) * W
            sage: T = TorsionQuadraticModule(V,W)
            sage: x = T.gen(0)
            sage: x
            (1, 0)
            sage: x.quadratic_product()
            1/2
            sage: x.quadratic_product().parent()
            Q/2Z
            sage: x*x
            1/2
            sage: (x*x).parent()
            Q/Z
        """
        value_module_qf = self.parent().value_module_qf()
        lift = self.lift()
        return value_module_qf(lift.inner_product(lift))

    q = quadratic_product


class TorsionQuadraticModule(FGP_Module_class, CachedRepresentation):
    r"""
    Finite quotients with a bilinear and a quadratic form.

    Let `V` be a symmetric FreeQuadraticModule and `W \subseteq V` a
    submodule of the same rank as `V`. The quotient `V / W` is a torsion
    quadratic module. It inherits a bilinear form `b` and a quadratic
    form `q`.

    `b: V \times V \to \QQ / m\ZZ`, where  `m\ZZ = (V,W)`
    and `b(x,y) = (x,y) + m\ZZ`

    `q: V \to \QQ / n\ZZ`, where `n\ZZ = 2(V,W) + \ZZ \{ (w,w) | w \in W \}`

    INPUT:

    - ``V`` -- a :class:`FreeModule` with a symmetric inner product matrix

    - ``W`` -- a submodule of ``V`` of the same rank as ``V``

    - ``check`` -- bool (default: ``True``)

    - ``modulus`` -- a rational number dividing `m` (default: `m`);
      the inner product `b` is defined in `\QQ /` ``modulus`` `\ZZ`

    - ``modulus_qf`` -- a rational number dividing `n` (default: `n`);
      the quadratic form `q` is defined in `\QQ /` ``modulus_qf`` `\ZZ`

    EXAMPLES::

        sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
        sage: V = FreeModule(ZZ, 3)
        sage: T = TorsionQuadraticModule(V, 5*V)
        sage: T
        Finite quadratic module over Integer Ring with invariants (5, 5, 5)
        Gram matrix of the quadratic form with values in Q/5Z:
        [1 0 0]
        [0 1 0]
        [0 0 1]
    """
    Element = TorsionQuadraticModuleElement

    @staticmethod
    def __classcall__(cls, V, W, gens=None, modulus=None, modulus_qf=None, check=True):
        r"""
        Return a :class:`TorsionQuadraticModule`.

        This method does the preprocessing for :meth:`sage.structure.CachedRepresentation`.

        TESTS::

            sage: q = matrix([1/2])
            sage: D1 = TorsionQuadraticForm(q)
            sage: D2 = TorsionQuadraticForm(q)
            sage: D1 is D2
            True
        """
        if check:
            if V.rank() != W.rank():
                raise ValueError("modules must be of the same rank")
            if V.base_ring() is not ZZ:
                raise NotImplementedError("only currently implemented over ZZ")
            if V.inner_product_matrix() != V.inner_product_matrix().transpose():
                raise ValueError("the cover must have a symmetric inner product")

            if gens is not None and V.span(gens) + W != V:
                raise ValueError("provided gens do not generate the quotient")

        # compute the modulus - this may be expensive
        if modulus is None or check:
            # The inner product of two elements `b(v1+W,v2+W)`
            # is defined `mod (V,W)`
            num = V.basis_matrix() * V.inner_product_matrix() * W.basis_matrix().T
            max_modulus = gcd(num.list())

        if modulus is None:
            modulus = max_modulus
        elif check and max_modulus / modulus not in V.base_ring():
            raise ValueError("the modulus must divide (V, W)")

        if modulus_qf is None or check:
            # The quadratic_product of an element `q(v+W)` is defined
            # `\mod 2(V,W) + ZZ\{ (w,w) | w in w\}`
            norm = gcd(W.gram_matrix().diagonal())
            max_modulus_qf = gcd(norm, 2 * modulus)

        if modulus_qf is None:
            modulus_qf = max_modulus_qf
        elif check and max_modulus_qf / modulus_qf not in V.base_ring():
            raise ValueError("the modulus_qf must divide (V, W)")
        return super(TorsionQuadraticModule, cls).__classcall__(cls, V, W, gens, modulus, modulus_qf)

    def __init__(self, V, W, gens, modulus, modulus_qf):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: T = TorsionQuadraticModule(ZZ^3, 6*ZZ^3)
            sage: TestSuite(T).run()
        """

        FGP_Module_class.__init__(self, V, W, check=True)
        if gens is not None:
            self._gens_user = tuple(self(g) for g in gens)
        else:
            # this is taken care of in the .gens method
            # we do not want this at initialization
            self._gens_user = None
        self._modulus = modulus
        self._modulus_qf = modulus_qf


    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeModule(ZZ,3)
            sage: T = TorsionQuadraticModule(V, 5*V,modulus=1)
            sage: T
            Finite quadratic module over Integer Ring with invariants (5, 5, 5)
            Gram matrix of the quadratic form with values in Q/Z:
            [0 0 0]
            [0 0 0]
            [0 0 0]
        """
        return ( "Finite quadratic module over %s with invariants %s\n"
                 % (self.base_ring(), self.invariants()) +
                 "Gram matrix of the quadratic form with values in %r:\n%r"
                 % (self.value_module_qf(), self.gram_matrix_quadratic()))

    def _module_constructor(self, V, W, check=False):
        r"""
        Construct a torsion quadratic module ``V / W``.

        INPUT:

        - ``V`` -- an module

        - ``W`` -- an submodule of ``V`` over the same base ring

        - ``check`` -- bool (default: ``False``);

          * if ``False``, then the value modulus is inherited from ``self``
          * if ``True``, it figures it out on its own. But that is expensive

        OUTPUT:

        The quotient ``V / W`` as a :class:`TorsionQuadraticModule`.

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = span([[1/2,1,1], [3/2,2,1], [0,0,1]], ZZ)
            sage: b = V.basis()
            sage: W = V.span([2*b[0]+4*b[1], 9*b[0]+12*b[1], 4*b[2]])
            sage: Q = TorsionQuadraticModule(V, W); Q
            Finite quadratic module over Integer Ring with invariants (4, 12)
            Gram matrix of the quadratic form with values in Q/(1/4)Z:
            [0 0]
            [0 0]
            sage: Q._module_constructor(V,W)
            Finite quadratic module over Integer Ring with invariants (4, 12)
            Gram matrix of the quadratic form with values in Q/(1/4)Z:
            [0 0]
            [0 0]
        """
        if check:
            # figuring out the modulus can be expensive
            return TorsionQuadraticModule(V, W, check=check)
        else:
            return TorsionQuadraticModule(V, W, check=check,
                                          modulus=self._modulus,
                                          modulus_qf=self._modulus_qf)

    def all_submodules(self):
        r"""
        Return a list of all submodules of ``self``.

        .. WARNING::

            This method creates all submodules in memory. The number of submodules
            grows rapidly with the number of generators. For example consider a
            vector space of dimension `n` over a finite field of prime order `p`.
            The number of subspaces is (very) roughly `p^{(n^2-n)/2}`.

        EXAMPLES::

            sage: D = IntegralLattice("D4").discriminant_group()
            sage: D.all_submodules()
            [Finite quadratic module over Integer Ring with invariants ()
              Gram matrix of the quadratic form with values in Q/2Z:
              [],
             Finite quadratic module over Integer Ring with invariants (2,)
              Gram matrix of the quadratic form with values in Q/2Z:
              [1],
             Finite quadratic module over Integer Ring with invariants (2,)
              Gram matrix of the quadratic form with values in Q/2Z:
              [1],
             Finite quadratic module over Integer Ring with invariants (2,)
              Gram matrix of the quadratic form with values in Q/2Z:
              [1],
             Finite quadratic module over Integer Ring with invariants (2, 2)
              Gram matrix of the quadratic form with values in Q/2Z:
              [  1 1/2]
              [1/2   1]]
        """
        from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        invs = self.invariants()
        # knows how to compute all subgroups
        A = AbelianGroupGap(invs)
        S = A.all_subgroups()
        # over ZZ submodules and subgroups are the same thing.
        submodules = []
        for sub in S:
            gen = [A(g).exponents() for g in sub.gens()]
            gen = [self.linear_combination_of_smith_form_gens(g) for g in gen]
            submodules.append(self.submodule(gen))
        return submodules

    def all_totally_isotropic_submodules(self, bilinear=False):
        r"""
        Return a list of all totally isotropic submodules.

        Be aware that the number of such submodules grows quickly.
        The method works by first creating all subgroups and then discarding the
        non - isotropic ones. This can be time and memory consuming.

        EXAMPLES::

            sage: q = Matrix(QQ,2,[2,1,1,2])/2
            sage: T = TorsionQuadraticForm(q)
            sage: T.all_totally_isotropic_submodules()
            [Finite quadratic module over Integer Ring with invariants ()
            Gram matrix of the quadratic form with values in Q/2Z:
            []]
            sage: T.all_totally_isotropic_submodules(bilinear=True)
            [Finite quadratic module over Integer Ring with invariants ()
            Gram matrix of the quadratic form with values in Q/2Z:
            [], Finite quadratic module over Integer Ring with invariants (2,)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1], Finite quadratic module over Integer Ring with invariants (2,)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1], Finite quadratic module over Integer Ring with invariants (2,)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1]]
        """
        isotropic = []
        n = self.value_module_qf().n
        for sub in self.all_submodules():
            if bilinear:
                q = sub.gram_matrix_bilinear()
            else:
                q = sub.gram_matrix_quadratic()
            if q == 0:
                isotropic.append(sub)
        return isotropic

    @cached_method
    def brown_invariant(self):
        r"""
        Return the Brown invariant of this torsion quadratic form.

        Let `(D,q)` be a torsion quadratic module with values in `\QQ / 2 \ZZ`.
        The Brown invariant `Br(D,q) \in \Zmod{8}` is defined by the equation

        .. MATH::

            \exp \left( \frac{2 \pi i }{8} Br(q)\right) =
            \frac{1}{\sqrt{D}} \sum_{x \in D} \exp(i \pi q(x)).

        The Brown invariant is additive with respect to direct sums of
        torsion quadratic modules.

        OUTPUT:

        - an element of `\Zmod{8}`

        EXAMPLES::

            sage: L = IntegralLattice("D4")
            sage: D = L.discriminant_group()
            sage: D.brown_invariant()
            4

        We require the quadratic form to be defined modulo `2 \ZZ`::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ,3,matrix.identity(3))
            sage: T = TorsionQuadraticModule((1/10)*V, V)
            sage: T.brown_invariant()
            Traceback (most recent call last):
            ...
            ValueError: the torsion quadratic form must have values in QQ / 2 ZZ
        """
        if self._modulus_qf != 2:
            raise ValueError("the torsion quadratic form must have values in "
                             "QQ / 2 ZZ")
        from sage.quadratic_forms.genera.normal_form import collect_small_blocks
        brown = IntegerModRing(8).zero()
        for p in self.annihilator().gen().prime_divisors():
            q = self.primary_part(p).normal_form()
            q = q.gram_matrix_quadratic()
            L = collect_small_blocks(q)
            for qi in L:
                brown += _brown_indecomposable(qi, p)
        return brown

    @cached_method
    def gram_matrix_bilinear(self):
        r"""
        Return the Gram matrix with respect to the generators.

        OUTPUT:

        A rational matrix ``G`` with ``G[i,j]`` given by the inner product
        of the `i`-th and `j`-th generator. Its entries are only well
        defined `\mod (V, W)`.

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ, 3, matrix.identity(3)*5)
            sage: T = TorsionQuadraticModule((1/5)*V, V)
            sage: T.gram_matrix_bilinear()
            [1/5   0   0]
            [  0 1/5   0]
            [  0   0 1/5]
        """
        gens = self.gens()
        n = len(gens)
        Q = self.base_ring().fraction_field()
        G = matrix.zero(Q, n)
        for i in range(n):
            for j in range(i + 1):
                G[i, j] = G[j, i] = (gens[i] * gens[j]).lift()
        return G

    @cached_method
    def gram_matrix_quadratic(self):
        r"""
        The Gram matrix of the quadratic form with respect to the generators.

        OUTPUT:

        - a rational matrix ``Gq`` with ``Gq[i,j] = gens[i]*gens[j]``
          and ``G[i,i] = gens[i].q()``

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: D4_gram = Matrix(ZZ, [[2,0,0,-1],[0,2,0,-1],[0,0,2,-1],[-1,-1,-1,2]])
            sage: D4 = FreeQuadraticModule(ZZ, 4, D4_gram)
            sage: D4dual = D4.span(D4_gram.inverse())
            sage: discrForm = TorsionQuadraticModule(D4dual, D4)
            sage: discrForm.gram_matrix_quadratic()
            [  1 1/2]
            [1/2   1]
            sage: discrForm.gram_matrix_bilinear()
            [  0 1/2]
            [1/2   0]
        """
        gens = self.gens()
        n = len(gens)
        Q = self.base_ring().fraction_field()
        G = matrix.zero(Q, n)
        for i in range(n):
            for j in range(i):
                G[i, j] = G[j, i] = (gens[i] * gens[j]).lift()
            G[i, i] = gens[i].q().lift()
        return G

    def gens(self):
        r"""
        Return generators of ``self``.

        There is no assumption on the generators except that they
        generate the module.

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeModule(ZZ, 3)
            sage: T = TorsionQuadraticModule(V, 5*V)
            sage: T.gens()
            ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        """
        if self._gens_user is None:
            return self.smith_form_gens()
        return self._gens_user

    def genus(self, signature_pair):
        r"""
        Return the genus defined by ``self`` and the ``signature_pair``.

        If no such genus exists, raise a ``ValueError``.

        REFERENCES:

        [Nik1977]_ Corollary 1.9.4 and 1.16.3.

        EXAMPLES::

            sage: L = IntegralLattice("D4").direct_sum(IntegralLattice("A2"))
            sage: D = L.discriminant_group()
            sage: genus = D.genus(L.signature_pair())
            sage: genus
            Genus of
            None
            Signature:  (6, 0)
            Genus symbol at 2:    1^4:2^-2
            Genus symbol at 3:     1^-5 3^-1
            sage: genus == L.genus()
            True

        Let `H` be an even unimodular lattice of signature `(9, 1)`.
        Then `L = D_4 + A_2` is primitively embedded in `H`. We compute the discriminant
        form of the orthogonal complement of `L` in `H`::

            sage: DK = D.twist(-1)
            sage: DK
            Finite quadratic module over Integer Ring with invariants (2, 6)
            Gram matrix of the quadratic form with values in Q/2Z:
            [  1 1/2]
            [1/2 1/3]

        We know that  `K` has signature `(5, 1)` and thus we can compute
        the genus of `K` as::

            sage: DK.genus((3,1))
            Genus of
            None
            Signature:  (3, 1)
            Genus symbol at 2:    1^2:2^-2
            Genus symbol at 3:     1^-3 3^1

        We can also compute the genus of an odd lattice
        from its discriminant form::

            sage: L = IntegralLattice(matrix.diagonal(range(1,5)))
            sage: D = L.discriminant_group()
            sage: D.genus((4,0))
            Genus of
            None
            Signature:  (4, 0)
            Genus symbol at 2:    [1^-2 2^1 4^1]_6
            Genus symbol at 3:     1^-3 3^1

        TESTS::

            sage: L.genus() == D.genus((4,0))
            True
            sage: D.genus((1,0))
            Traceback (most recent call last):
            ...
            ValueError: this discriminant form and signature do not define a genus

        A systematic test of lattices of small ranks and determinants::

            sage: from sage.quadratic_forms.genera.genus import genera
            sage: signatures = [(1,0),(1,1),(1,2),(3,0),(0,4)]
            sage: dets = range(1,33)
            sage: genera = flatten([genera(s, d, even=False) for d in dets for s in signatures])    # long time
            sage: all(g == g.discriminant_form().genus(g.signature_pair()) for g in genera)  # long time
            True
            """
        from sage.quadratic_forms.genera.genus import (Genus_Symbol_p_adic_ring,
                                                       GenusSymbol_global_ring,
                                                       p_adic_symbol,
                                                       is_GlobalGenus,
                                                       _blocks)
        from sage.misc.misc_c import prod
        s_plus = ZZ(signature_pair[0])
        s_minus = ZZ(signature_pair[1])
        rank = s_plus + s_minus
        if s_plus < 0 or s_minus < 0:
            raise ValueError("signatures must be non-negative")
        if len(self.invariants()) > rank:
            raise ValueError("this discriminant form and " +
                             "signature do not define a genus")
        disc = self.cardinality()
        determinant = ZZ(-1)**s_minus * disc
        local_symbols = []
        for p in (2 * disc).prime_divisors():
            D = self.primary_part(p)
            if len(D.invariants()) != 0:
                G_p = D.gram_matrix_quadratic().inverse()
                # get rid of denominators without changing the local equivalence class
                G_p *= G_p.denominator()**2
                G_p = G_p.change_ring(ZZ)
                local_symbol = p_adic_symbol(G_p, p, D.invariants()[-1].valuation(p))
            else:
                local_symbol = []

            rk = rank - len(D.invariants())
            if rk > 0:
                if p == 2:
                    det = determinant.prime_to_m_part(2)
                    det *= prod([di[2] for di in local_symbol])
                    det = det % 8
                    local_symbol.append([ZZ(0), rk, det, ZZ(0), ZZ(0)])
                else:
                    det = legendre_symbol(determinant.prime_to_m_part(p), p)
                    det = (det * prod([di[2] for di in local_symbol]))
                    local_symbol.append([ZZ(0), rk, det])
            local_symbol.sort()
            local_symbol = Genus_Symbol_p_adic_ring(p, local_symbol)
            local_symbols.append(local_symbol)

        # This genus has the right discriminant group
        # but it may be empty
        sym2 = local_symbols[0].symbol_tuple_list()

        if sym2[0][0] != 0:
            sym2 = [[ZZ(0), ZZ(0), ZZ(1), ZZ(0), ZZ(0)]] + sym2
        if len(sym2) <= 1 or sym2[1][0] != 1:
            sym2 = sym2[:1] + [[ZZ(1), ZZ(0), ZZ(1), ZZ(0), ZZ(0)]] + sym2[1:]
        if len(sym2) <= 2 or sym2[2][0] != 2:
            sym2 = sym2[:2] + [[ZZ(2), ZZ(0), ZZ(1), ZZ(0), ZZ(0)]] + sym2[2:]

        if self.value_module_qf().n == 1:
            # in this case the blocks of scales 1, 2, 4 are under determined
            # make sure the first 3 symbols are of scales 1, 2, 4
            # i.e. their valuations are 0, 1, 2

            # the form is odd
            block0 = [b for b in _blocks(sym2[0]) if b[3] == 1]

            o = sym2[1][3]
            # no restrictions on determinant and
            # oddity beyond existence
            # but we know if even or odd
            block1 = [b for b in _blocks(sym2[1]) if b[3] == o]

            d = sym2[2][2]
            o = sym2[2][3]
            t = sym2[2][4]
            # if the jordan block of scale 2 is even we know it
            if o == 0:
                block2 = [sym2[2]]
            # if it is odd we know det and oddity mod 4 at least
            else:
                block2 = [b for b in _blocks(sym2[2]) if b[3] == o
                          and (b[2] - d) % 4 == 0
                          and (b[4] - t) % 4 == 0
                          and (b[2] - d) % 8 == (b[4] - t) % 8  # if the oddity is altered by 4 then so is the determinant
                         ]
        elif self.value_module_qf().n == 2:
            # the form is even
            block0 = [b for b in _blocks(sym2[0]) if b[3] == 0]

            # if the jordan block of scale 2 is even we know it
            d = sym2[1][2]
            o = sym2[1][3]
            t = sym2[1][4]
            if o == 0:
                block1 = [sym2[1]]
            else:
                # the block is odd and we know det and oddity mod 4
                block1 = [b for b in _blocks(sym2[1])
                          if b[3] == o
                          and (b[2] - d) % 4 == 0
                          and (b[4] - t) % 4 == 0
                          and (b[2] - d) % 8 == (b[4] - t) % 8 # if the oddity is altered by 4 then so is the determinant
                         ]
            # this is completely determined
            block2 = [sym2[2]]
        else:
            raise ValueError("this is not a discriminant form")

        # figure out which symbol defines a genus and return that
        for b0 in block0:
            for b1 in block1:
                for b2 in block2:
                    sym2[:3] = [b0, b1, b2]
                    local_symbols[0] = Genus_Symbol_p_adic_ring(2, sym2)
                    genus = GenusSymbol_global_ring(signature_pair, local_symbols)
                    if is_GlobalGenus(genus):
                        # make the symbol sparse again.
                        i = 0
                        k = 0
                        while i < 3:
                            if sym2[k][1] == 0:
                                sym2.pop(k)
                            else:
                                k = k + 1
                            i = i + 1
                        local_symbols[0] = Genus_Symbol_p_adic_ring(2, sym2)
                        genus = GenusSymbol_global_ring(signature_pair, local_symbols)
                        return genus
        raise ValueError("this discriminant form and signature do not define a genus")

    @cached_method
    def is_degenerate(self):
        r"""
        Return if the underlying bilinear form is degenerate.

        EXAMLES::

            sage: T = TorsionQuadraticForm(matrix([1/27]))
            sage: D = T.submodule([T.gen(0)*3])
            sage: D.is_degenerate()
            True
        """
        return 1 != self.orthogonal_submodule_to(self.gens()).cardinality()

    def is_genus(self, signature_pair, even=True):
        r"""
        Return ``True`` if there is a lattice with this signature and discriminant form.

        .. TODO::

            implement the same for odd lattices

        INPUT:

        - signature_pair -- a tuple of non negative integers ``(s_plus, s_minus)``
        - even -- bool (default: ``True``)

        EXAMPLES::

            sage: L = IntegralLattice("D4").direct_sum(IntegralLattice(3 * Matrix(ZZ,2,[2,1,1,2])))
            sage: D = L.discriminant_group()
            sage: D.is_genus((6,0))
            True

        Let us see if there is a lattice in the genus defined by the same discriminant form
        but with a different signature::

            sage: D.is_genus((4,2))
            False
            sage: D.is_genus((16,2))
            True
        """
        s_plus = ZZ(signature_pair[0])
        s_minus = ZZ(signature_pair[1])
        if s_plus < 0 or s_minus < 0:
            raise ValueError("signature invariants must be non negative")
        rank = s_plus + s_minus
        signature = s_plus - s_minus
        D = self.cardinality()
        det = (-1)**s_minus * D
        if rank < len(self.invariants()):
            return False
        if even and self._modulus_qf != 2:
            raise ValueError("the discriminant form of an even lattice has"
                             "values modulo 2.")
        if (not even) and not (self._modulus == self._modulus_qf == 1):
            raise ValueError("the discriminant form of an odd lattice has"
                             "values modulo 1.")
        if not even:
            raise NotImplementedError("at the moment sage knows how to do this only for even genera. " +
                                      " Help us to implement this for odd genera.")
        for p in D.prime_divisors():
            # check the determinant conditions
            Q_p = self.primary_part(p)
            gram_p = Q_p.gram_matrix_quadratic()
            length_p = len(Q_p.invariants())
            u = det.prime_to_m_part(p)
            up = gram_p.det().numerator().prime_to_m_part(p)
            if p != 2 and length_p == rank:
                if legendre_symbol(u, p) != legendre_symbol(up, p):
                    return False
            if p == 2:
                if rank % 2 != length_p % 2:
                    return False
                n = (rank - length_p) / 2
                if u % 4 != (-1)**(n % 2) * up % 4:
                    return False
                if rank == length_p:
                    a = QQ(1) / QQ(2)
                    b = QQ(3) / QQ(2)
                    diag = gram_p.diagonal()
                    if not (a in diag or b in diag):
                        if u % 8 != up % 8:
                            return False
        if self.brown_invariant() != signature:
            return False
        return True


    def is_isomorphic_to(self, other):
        r"""
        Return if the underlying quadratic forms are isomorphic.

        INPUT:

        - ``other`` -- a torsion quadratic form
        """
        if self._modulus_qf != 2:
            raise NotImplementedError()
        if self.value_module_qf() != other.value_module_qf():
            return False, None
        if self.invariants() != other.invariants():
            return False, None
        if self.is_degenerate() != other.is_degenerate():
            return False, None
        if len(self.invariants())==0:
            return True, (self.gens(), other.gens())
        qf1 = self.normal_form().gram_matrix_quadratic()
        qf2 = other.normal_form().gram_matrix_quadratic()
        isom = (self.normal_form().gens(), other.normal_form().gens())
        if (qf1 - qf2).denominator() != 1: # the bilinear forms agree
            return False, None
        if not self.is_degenerate():
            # check the diagonal agrees too
            if ZZ(2).divides(self.invariants()[-1]) and qf1 != qf2:
                return False, None
            return True, isom
        if len(set(self.invariants())) == 1 and not ZZ(2).divides(self.invariants()[-1]):
            return True, isom
        # now the normal form is not unique anymore.
        ker1 = self.orthogonal_submodule_to(self.gens())
        ker2 = other.orthogonal_submodule_to(other.gens())
        if ker1.invariants() != ker2.invariants():
            return False, None
        #if 2 == self.invariants()[-1]:
        #    d1 = set(ker1.gram_matrix_quadratic().diagonal())
        #    d2 = set(ker2.gram_matrix_quadratic().diagonal())
        #    return d1 == d2
        from sage.groups.fqf_orthogonal import _isom_fqf
        try:
            isom = _isom_fqf(self,other)
            return True, (self.smith_form_gens(), isom)
        except ValueError:
            return False, None

    def orthogonal_group(self, gens=None, check=False):
        r"""
        Orthogonal group of the associated torsion quadratic form.

        .. WARNING::

            This is can be smaller than the orthogonal group of the bilinear form.

        INPUT:

        - ``gens`` --  a list of generators, for instance square matrices,
                       something that acts on ``self``, or an automorphism
                       of the underlying abelian group
        - ``check`` -- perform additional checks on the generators

        EXAMPLES:

        You can provide generators to obtain a subgroup of the full orthogonal group::

            sage: D = TorsionQuadraticForm(matrix.identity(2)/2)
            sage: f = matrix(2,[0,1,1,0])
            sage: D.orthogonal_group(gens=[f]).order()
            2

        If no generators are given a slow brute force approach is used to calculate the full orthogonal group::

            sage: D = TorsionQuadraticForm(matrix.identity(3)/2)
            sage: OD = D.orthogonal_group()
            sage: OD.order()
            6
            sage: fd = D.hom([D.1,D.0,D.2])
            sage: OD(fd)
            [0 1 0]
            [1 0 0]
            [0 0 1]

        We compute the kernel of the action of the orthogonal group of `L` on the discriminant group.

            sage: L = IntegralLattice('A4')
            sage: O = L.orthogonal_group()
            sage: D = L.discriminant_group()
            sage: Obar = D.orthogonal_group(O.gens())
            sage: O.order()
            240
            sage: Obar.order()
            2
            sage: phi = O.hom([Obar(g) for g in O.gens()])
            sage: phi.kernel().order()
            120
        """
        from sage.groups.fqf_orthogonal import FqfOrthogonalGroup,_isom_fqf
        try:
            from sage.groups.fqf_orthogonal_gens import _compute_gens as _get_gens
        except ImportError:
            from sage.groups.fqf_orthogonal import _isom_fqf as _get_gens

        from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap

        if self.value_module_qf().n != 2:
            raise NotImplementedError("orthogonal groups are implemented/tested only for even forms")
            # slow brute force implementation
            # gens = [matrix(g) for g in _isom_fqf(self)]

        flag = False
        if gens is None:
            try:
                return self._orthogonal_group
            except AttributeError:
                flag = True
                gens = _get_gens(self)
        else:
            # see if there is an action
            try:
                gens = [matrix(x*g for x in self.smith_form_gens()) for g in gens]
            except TypeError:
                pass
        ambient = AbelianGroupGap(self.invariants()).aut()
        # the ambient knows what to do with the generators
        gens = tuple(ambient(g) for g in gens)
        gens = tuple(g for g in gens if not g.is_one())
        Oq =  FqfOrthogonalGroup(ambient, gens, self, check=check)
        if flag:
            self._orthogonal_group = Oq
        return Oq

    O = orthogonal_group


    def orthogonal_submodule_to(self, S):
        r"""
        Return the submodule orthogonal to ``S``.

        INPUT:

        - ``S`` -- a submodule, list, or tuple of generators

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeModule(ZZ, 10)
            sage: T = TorsionQuadraticModule(V, 3*V)
            sage: S = T.submodule(T.gens()[:5])
            sage: O = T.orthogonal_submodule_to(S)
            sage: O
            Finite quadratic module over Integer Ring with invariants (3, 3, 3, 3, 3)
            Gram matrix of the quadratic form with values in Q/3Z:
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            sage: O.V() + S.V() == T.V()
            True
        """
        if not isinstance(S, TorsionQuadraticModule):
            S = self.submodule(S)
        else:
            if not S.is_submodule(self):
                raise ValueError("S must be a submodule of this module")

        G = self.V().inner_product_matrix()
        T = self.V().basis_matrix()
        S = S.V().basis_matrix()
        m = self._modulus

        Y = T * G * S.transpose()
        # Elements of the ambient module which pair integrally with self.V()
        integral = Y.inverse() * T
        # Element of the ambient module which pair in mZZ with self.V()
        orthogonal = m * integral
        orthogonal = self.V().span(orthogonal.rows())
        # We have to make sure we get a submodule
        orthogonal = orthogonal.intersection(self.V())
        orthogonal = self.submodule(orthogonal.gens())
        return orthogonal

    @cached_method
    def normal_form(self, partial=False):
        r"""
        Return the normal form of this torsion quadratic module.

        Two torsion quadratic modules are isomorphic if and only if they have
        the same value modules and the same normal form.

        A torsion quadratic module `(T,q)` with values in `\QQ/n\ZZ` is
        in normal form if the rescaled quadratic module `(T, q/n)`
        with values in `\QQ/\ZZ` is in normal form.

        For the definition of normal form see [MirMor2009]_ IV Definition 4.6.
        Below are some of its properties.
        Let `p` be odd and `u` be the smallest non-square modulo `p`.
        The normal form is a diagonal matrix with diagonal entries either `p^n`
        or `u p^n`.

        If `p = 2` is even, then the normal form consists of
        1 x 1 blocks of the form

        .. MATH::

            (0), \quad 2^n(1),\quad 2^n(3),\quad 2^n(5) ,\quad 2^n(7)

        or of `2 \times 2` blocks of the form

        .. MATH::

            2^n
            \left(\begin{matrix}
                2 & 1\\
                1 & 2
            \end{matrix}\right), \quad
            2^n
            \left(\begin{matrix}
                0 & 1\\
                1 & 0
            \end{matrix}\right).

       The blocks are ordered by their valuation.

        INPUT:

        - partial - bool (default: ``False``) return only a partial normal form
          it is not unique but still useful to extract invariants

        OUTPUT:

        - a torsion quadratic module

        EXAMPLES::

            sage: L1=IntegralLattice(matrix([[-2,0,0],[0,1,0],[0,0,4]]))
            sage: L1.discriminant_group().normal_form()
            Finite quadratic module over Integer Ring with invariants (2, 4)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/2   0]
            [  0 1/4]
            sage: L2=IntegralLattice(matrix([[-2,0,0],[0,1,0],[0,0,-4]]))
            sage: L2.discriminant_group().normal_form()
            Finite quadratic module over Integer Ring with invariants (2, 4)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/2   0]
            [  0 1/4]

        We check that :trac:`24864` is fixed::

            sage: L1=IntegralLattice(matrix([[-4,0,0],[0,4,0],[0,0,-2]]))
            sage: AL1=L1.discriminant_group()
            sage: L2=IntegralLattice(matrix([[-4,0,0],[0,-4,0],[0,0,2]]))
            sage: AL2=L2.discriminant_group()
            sage: AL1.normal_form()
            Finite quadratic module over Integer Ring with invariants (2, 4, 4)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1/2   0   0]
            [  0 1/4   0]
            [  0   0 5/4]
            sage: AL2.normal_form()
            Finite quadratic module over Integer Ring with invariants (2, 4, 4)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1/2   0   0]
            [  0 1/4   0]
            [  0   0 5/4]

        Some exotic cases::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: D4_gram = Matrix(ZZ,4,4,[2,0,0,-1,0,2,0,-1,0,0,2,-1,-1,-1,-1,2])
            sage: D4 = FreeQuadraticModule(ZZ,4,D4_gram)
            sage: D4dual = D4.span(D4_gram.inverse())
            sage: T = TorsionQuadraticModule((1/6)*D4dual,D4)
            sage: T
            Finite quadratic module over Integer Ring with invariants (6, 6, 12, 12)
            Gram matrix of the quadratic form with values in Q/(1/3)Z:
            [1/18 5/36    0    0]
            [5/36 1/18 5/36 5/36]
            [   0 5/36 1/36 1/72]
            [   0 5/36 1/72 1/36]
            sage: T.normal_form()
            Finite quadratic module over Integer Ring with invariants (6, 6, 12, 12)
            Gram matrix of the quadratic form with values in Q/(1/3)Z:
            [ 1/6 1/12    0    0    0    0    0    0]
            [1/12  1/6    0    0    0    0    0    0]
            [   0    0 1/12 1/24    0    0    0    0]
            [   0    0 1/24 1/12    0    0    0    0]
            [   0    0    0    0  1/9    0    0    0]
            [   0    0    0    0    0  1/9    0    0]
            [   0    0    0    0    0    0  1/9    0]
            [   0    0    0    0    0    0    0  1/9]

        TESTS:

        A degenerate case::

            sage: T = TorsionQuadraticModule((1/6)*D4dual, D4, modulus=1/36)
            sage: T.normal_form()
            Finite quadratic module over Integer Ring with invariants (6, 6, 12, 12)
            Gram matrix of the quadratic form with values in Q/(1/18)Z:
            [1/36 1/72    0    0    0    0    0    0]
            [1/72 1/36    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
        """
        gens = []
        from sage.quadratic_forms.genera.normal_form import p_adic_normal_form, _normalize
        for p in self.annihilator().gen().prime_divisors():
            D_p = self.primary_part(p)
            q_p = D_p.gram_matrix_quadratic()
            q_p = q_p / D_p._modulus_qf

            # continue with the non-degenerate part
            r = q_p.rank()
            if r != q_p.ncols():
                U = q_p._clear_denom()[0].hermite_form(transformation=True)[1]
            else:
                U = q_p.parent().identity_matrix()
            kernel = U[r:, :]
            nondeg = U[:r, :]
            q_p = nondeg * q_p * nondeg.T

            # the normal form is implemented for p-adic lattices
            # so we should work with the lattice q_p --> q_p^-1
            q_p1 = q_p.inverse()
            prec = self.annihilator().gen().valuation(p) + 5
            D, U = p_adic_normal_form(q_p1, p, precision=2*prec + 5, partial=partial)
            # if we compute the inverse in the p-adics everything explodes --> go to ZZ
            U = U.change_ring(ZZ).inverse().transpose()

            # the inverse is in normal form - so to get a normal form for the original one
            # it is enough to massage each 1x1 resp. 2x2 block.
            U = U.change_ring(Zp(p, type='fixed-mod', prec=prec)).change_ring(ZZ)
            D = U * q_p * U.T * p**q_p.denominator().valuation(p)
            D = D.change_ring(Zp(p, type='fixed-mod', prec=prec))
            _, U1 = _normalize(D, normal_odd=False)
            U = U1.change_ring(ZZ) * U

            # reattach the degenerate part
            nondeg = U * nondeg
            U = nondeg.stack(kernel)

            # apply U to the generators
            n = U.ncols()
            gens_p = []
            for i in range(n):
                g = self.V().zero()
                for j in range(n):
                    g += D_p.gens()[j].lift() * U[i, j]
                gens_p.append(g)
            gens += gens_p
        return self.submodule_with_gens(gens)

    def primary_part(self, m):
        r"""
        Return the ``m``-primary part of this torsion quadratic module
        as a submodule.

        INPUT:

        - ``m`` -- an integer

        OUTPUT:

        - a submodule

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: T = TorsionQuadraticModule((1/6)*ZZ^3,ZZ^3)
            sage: T
            Finite quadratic module over Integer Ring with invariants (6, 6, 6)
            Gram matrix of the quadratic form with values in Q/(1/3)Z:
            [1/36    0    0]
            [   0 1/36    0]
            [   0    0 1/36]
            sage: T.primary_part(2)
            Finite quadratic module over Integer Ring with invariants (2, 2, 2)
            Gram matrix of the quadratic form with values in Q/(1/3)Z:
            [1/4   0   0]
            [  0 1/4   0]
            [  0   0 1/4]

        TESTS::

            sage: T == T.primary_part(T.annihilator().gen())
            True
        """
        annihilator = self.annihilator().gen()
        a = annihilator.prime_to_m_part(m)
        return self.submodule((a * self.V()).gens())

    def submodule_with_gens(self, gens):
        r"""
        Return a submodule with generators given by ``gens``.

        INPUT:

        - ``gens`` -- a list of generators that convert into ``self``

        OUTPUT:

        - a submodule with the specified generators

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ,3,matrix.identity(3)*10)
            sage: T = TorsionQuadraticModule((1/10)*V, V)
            sage: g = T.gens()
            sage: new_gens = [2*g[0], 5*g[0]]
            sage: T.submodule_with_gens(new_gens)
            Finite quadratic module over Integer Ring with invariants (10,)
            Gram matrix of the quadratic form with values in Q/2Z:
            [2/5   0]
            [  0 1/2]

        The generators do not need to be independent::

            sage: new_gens = [g[0], 2*g[1], g[0], g[1]]
            sage: T.submodule_with_gens(new_gens)
            Finite quadratic module over Integer Ring with invariants (10, 10)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1/10    0 1/10    0]
            [   0  2/5    0  1/5]
            [1/10    0 1/10    0]
            [   0  1/5    0 1/10]

        TESTS:

        Test that things work without specified gens too::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ,3,matrix.identity(3)*5)
            sage: T = TorsionQuadraticModule((1/5)*V, V)
            sage: T
            Finite quadratic module over Integer Ring with invariants (5, 5, 5)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/5   0   0]
            [  0 1/5   0]
            [  0   0 1/5]
            sage: T.submodule(T.gens()[:2])
            Finite quadratic module over Integer Ring with invariants (5, 5)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/5   0]
            [  0 1/5]
        """
        gens = tuple(self(v) for v in gens)
        V = self.V().submodule([v.lift() for v in gens]) + self._W
        W = self.W()
        return TorsionQuadraticModule(V, W, gens=gens, modulus=self._modulus,
                                      modulus_qf=self._modulus_qf, check=False)

    def twist(self, s):
        r"""
        Return the torsion quadratic module with quadratic form scaled by ``s``.

        If the old form was defined modulo `n`, then the new form is defined
        modulo `n s`.

        INPUT:

        - ``s`` - a rational number

        EXAMPLES::

            sage: q = TorsionQuadraticForm(matrix.diagonal([3/9, 1/9]))
            sage: q.twist(-1)
            Finite quadratic module over Integer Ring with invariants (3, 9)
            Gram matrix of the quadratic form with values in Q/Z:
            [2/3   0]
            [  0 8/9]

        This form is defined modulo `3`::

            sage: q.twist(3)
            Finite quadratic module over Integer Ring with invariants (3, 9)
            Gram matrix of the quadratic form with values in Q/3Z:
            [  1   0]
            [  0 1/3]

        The next form is defined modulo `4`::

            sage: q.twist(4)
            Finite quadratic module over Integer Ring with invariants (3, 9)
            Gram matrix of the quadratic form with values in Q/4Z:
            [4/3   0]
            [  0 4/9]
        """
        s = self.base_ring().fraction_field()(s)
        n = self.V().degree()
        inner_product_matrix = s * self.V().inner_product_matrix()
        ambient = FreeQuadraticModule(self.base_ring(), n, inner_product_matrix)
        V = ambient.span(self.V().basis())
        W = ambient.span(self.W().basis())
        modulus = s.abs() * self._modulus
        modulus_qf = s.abs() * self._modulus_qf
        return TorsionQuadraticModule(V, W, modulus=modulus, modulus_qf=modulus_qf)

    def value_module(self):
        r"""
        Return `\QQ / m\ZZ` with `m = (V, W)`.

        This is where the inner product takes values.

        EXAMPLES::

            sage: A2 = Matrix(ZZ, 2, 2, [2,-1,-1,2])
            sage: L = IntegralLattice(2*A2)
            sage: D = L.discriminant_group()
            sage: D
            Finite quadratic module over Integer Ring with invariants (2, 6)
            Gram matrix of the quadratic form with values in Q/2Z:
            [  1 1/2]
            [1/2 1/3]
            sage: D.value_module()
            Q/Z
        """
        return QmodnZ(self._modulus)

    def value_module_qf(self):
        r"""
        Return `\QQ / n\ZZ` with `n\ZZ = (V,W) + \ZZ \{ (w,w) | w \in W \}`.

        This is where the torsion quadratic form takes values.

        EXAMPLES::

            sage: A2 = Matrix(ZZ, 2, 2, [2,-1,-1,2])
            sage: L = IntegralLattice(2*A2)
            sage: D = L.discriminant_group()
            sage: D
            Finite quadratic module over Integer Ring with invariants (2, 6)
            Gram matrix of the quadratic form with values in Q/2Z:
            [  1 1/2]
            [1/2 1/3]
            sage: D.value_module_qf()
            Q/2Z
        """
        return QmodnZ(self._modulus_qf)

    def to_smith(self):
        r"""
        Return the transformation matrix from the user to smith form generators.

        To go in the other direction use :meth:`to_gens`.

        OUTPUT:

        - an integer matrix

        EXAMLES::

            sage: L2 = IntegralLattice(3 * matrix([[-2,0,0],[0,1,0],[0,0,-4]]))
            sage: D = L2.discriminant_group().normal_form()
            sage: D
            Finite quadratic module over Integer Ring with invariants (3, 6, 12)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/2   0   0   0   0]
            [  0 1/4   0   0   0]
            [  0   0 1/3   0   0]
            [  0   0   0 1/3   0]
            [  0   0   0   0 2/3]
            sage: D.to_smith()
            [0 3 0]
            [0 0 3]
            [0 2 0]
            [1 0 0]
            [0 0 4]
            sage: T = D.to_smith()*D.to_gens()
            sage: T
            [ 3  0 33  0  0]
            [ 0 33  0  0  3]
            [ 2  0 22  0  0]
            [ 0  0  0  1  0]
            [ 0 44  0  0  4]

        The matrix `T` now satisfies a certain congruence::

            sage: for i in range(T.nrows()):
            ....:     T[:,i] = T[:,i] % D.gens()[i].order()
            sage: T
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
        """
        to_smith = matrix(self.base_ring(), [t.vector() for t in self.gens()])
        return to_smith

    @cached_method
    def to_gens(self):
        r"""
        Return the transformation matrix from smith form to user generators.

        To go in the other direction use :meth:`to_smith`.

        OUTPUT:

        - an integer matrix

        EXAMPLES::

            sage: L2 = IntegralLattice(3 * matrix([[-2,0,0],[0,1,0],[0,0,-4]]))
            sage: D = L2.discriminant_group().normal_form()
            sage: D
            Finite quadratic module over Integer Ring with invariants (3, 6, 12)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/2   0   0   0   0]
            [  0 1/4   0   0   0]
            [  0   0 1/3   0   0]
            [  0   0   0 1/3   0]
            [  0   0   0   0 2/3]
            sage: D.to_gens()
            [ 0  0  0  1  0]
            [ 1  0 11  0  0]
            [ 0 11  0  0  1]
            sage: T = D.to_gens()*D.to_smith()
            sage: T
            [ 1  0  0]
            [ 0 25  0]
            [ 0  0 37]

        This matrix satisfies the congruence::

            sage: for i in range(T.ncols()):
            ....:     T[:, i] = T[:, i] % D.smith_form_gens()[i].order()
            sage: T
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        invs = self.invariants()
        R = IntegerModRing(invs[-1])
        E = matrix.identity(R, len(invs))
        # view self as a submodule of (ZZ/nZZ)^(len(invs))
        B = self.to_smith().change_ring(R)
        for k in range(B.ncols()):
            B[:, k] *= invs[-1] // invs[k]
            E[:, k] *= invs[-1] // invs[k]
        to_gens = B.solve_left(E)
        to_gens = to_gens.change_ring(ZZ)
        # double check that we did not mess up
        # tmp = to_gens * self.to_smith()
        # E = matrix.identity(R, len(invs))
        # for k in range(E.ncols()):
        #     assert 0 == (tmp[:, k]-E[:, k]) % invs[k]
        return to_gens

    _to_smith = to_smith
    _to_gens = to_gens

    def direct_sum(self, other, return_maps=False):
        r"""
        """
        V, fVs, fVo = self.V().direct_sum(other.V(),return_embeddings=True)
        W = self.W().direct_sum(other.W())
        n = len(self.V().gens())
        T = TorsionQuadraticModule(V, W, modulus=self._modulus)
        fs = self.hom([fVs(g.lift()) for g in self.gens()], T)
        fo = other.hom([fVo(g.lift()) for g in other.gens()], T)
        if return_maps:
            return T, fs, fo, fVs, fVo
        else:
            return T, fs, fo

    def _subgroup_to_gap(self, S):
        r"""
        """
        A = self.orthogonal_group().domain()
        gensS = [A(self(s)).gap() for s in S.gens()]
        Sgap = A.gap().Subgroup(gensS)
        return Sgap

    def _subgroup_from_gap(self, S):
        r"""
        """
        A = self.orthogonal_group().domain()
        S = self.submodule([self.linear_combination_of_smith_form_gens(
                            A(g).exponents())
                            for g in S.GeneratorsOfGroup()])
        return S

    def stab(self, G, S, backend="sage"):
        r"""

        EXAMPLES::

            sage: L = IntegralLattice("A4")
            sage: D = L.discriminant_group()
            sage: G = D.orthogonal_group()
            sage: S = D.submodule([D.0])
            sage: D.stab(G,S).order()
            2
            sage: D.stab(G,S,backend="magma").order()
            2

        NEED MORE TESTS!
        """
        if self.cardinality()==1:
            return G
        if backend=="sage":
            mu = libgap.function_factory("mu:=function(x,g) return(Image(g,x)); end;")
            Sgap = self._subgroup_to_gap(S)
            stab = libgap.Stabilizer(G.gap(), Sgap, mu)
            return G._subgroup_constructor(stab)
        elif backend=="magma":
            from sage.interfaces.magma import magma as m
            from sage.rings.all import GF
            degree = len(self.gens())
            p = self.invariants()[0]
            assert p.is_prime() and all(p==g for g in self.invariants())
            field =  GF(p)
            V = field**degree
            gensG = [g.matrix().change_ring(field) for g in G.gens()]
            GL = m.GeneralLinearGroup(degree, field)
            mG = GL.sub(gensG)
            gensS = [V(self(s).vector()) for s in S.gens()]
            mS = m(V).sub(gensS)
            stab = mG.Stabiliser(mS)
            stab = stab.Generators().SetToSequence()
            stab = [g.Matrix().sage().change_ring(ZZ) for g in stab]
            stab = G.subgroup(stab)
            assert all(s*g in S for s in S.gens() for g in stab.gens())
            return stab

    def subgroup_representatives(self, H, G, order, g=1, algorithm=None):
        r"""
        Return representatives of the subgroups of `H` modulo the action of `G`.

        INPUT:

        - ``H`` -- a submodule of `self`

        - ``G`` -- a group of automorphisms

        - ``algorithm`` -- one of the following

          * ``"hulpke"`` -- following an algorithm of A. Hulpke

          * ``"brute force""`` -- enumerates all subgroups first and takes orbits.

          * ``"elementary"`` -- H must be elementary abelian, used linear algebra in gap

          * ``"magma"`` -- use magma; H must be elementary abelian

        OUTPUT:

        - a list of submodules

        EXAMLES::

            sage: D = TorsionQuadraticForm(matrix.diagonal([1/2,1/2,1/2]))
            sage: G = D.orthogonal_group()
            sage: g = G.one()
            sage: L = D.subgroup_representatives(D, G, 4, g=g, algorithm='elementary')
            sage: L.sort(key=lambda x: x.representative().gram_matrix_quadratic())
            sage: L
            [Orbit represented by
            Finite quadratic module over Integer Ring with invariants (2, 2)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1/2   0]
            [  0 1/2], Orbit represented by
            Finite quadratic module over Integer Ring with invariants (2, 2)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1/2   0]
            [  0   1], Orbit represented by
            Finite quadratic module over Integer Ring with invariants (2, 2)
            Gram matrix of the quadratic form with values in Q/2Z:
            [  1 1/2]
            [1/2   1]]


        TESTS::

            sage: T = TorsionQuadraticForm(matrix.diagonal([2/3,2/9,4/9]))
            sage: H = T.V() & (1/3)*T.W()
            sage: H = T.submodule(H.gens())
            sage: G = T.orthogonal_group()
            sage: g = G.one()
            sage: subs1 = T.subgroup_representatives(H, G, 9, g=g, algorithm="brute force")
            sage: subs2 = T.subgroup_representatives(H, G, 9, g=g, algorithm="magma")
            sage: subs3 = T.subgroup_representatives(H, G, 9, g=g, algorithm="elementary")
            sage: len(subs1)
            4
            sage: len(subs2)
            4
            sage: len(subs3)
            4

        Larger instances only work with magma as backend::

            sage: U = matrix(ZZ,2,[0,1,1,0])/2
            sage: T = matrix.block_diagonal(5*[U])
            sage: T = TorsionQuadraticForm(T)
            sage: G = T.orthogonal_group()
            sage: subs = T.subgroup_representatives(T,G,4, g=G.one(),algorithm='magma')
            sage: subs = [x.representative().normal_form() for x in subs]
            sage: subs = [x.gram_matrix_quadratic() for x in subs]
            sage: subs.sort()

        By witts theorem the subspace orbits are given by their isomorphism class as quadratic forms::

            sage: subs
            [
            [0 0]  [  0 1/2]  [1 0]  [  1 1/2]
            [0 0], [1/2   0], [0 0], [1/2   1]
            ]

        """
        g = G(g)
        if H.cardinality() == 1:
            return [orbit_sage(rep=H,stab=G)]
        if algorithm is None:
            if not H.invariants()[0].is_prime():
                algorithm = 'hulpke'
            elif H.cardinality() < 2**8:
                algorithm = 'elementary'
            else:
                algorithm = 'magma'
        # treat the trivial case
        if algorithm == "magma":
            from sage.rings.all import GF,ZZ
            # we continue with magma
            G0 = H.orthogonal_group()
            action_hom = G.hom(G.gens(), codomain=G0,check=False)
            p = H.invariants()[0]
            gens = [i.matrix().change_ring(GF(p))
                    for i in action_hom.image(G).gens()]
            orb = OrbitsOfSpaces(gens, k=order.valuation(p),
                                 D=self, H=H, action_hom=action_hom,
                                 g=G0(g).matrix())
            return orb
        # some algorithms implemented in gap follow
        A = G.domain()
        Hgap = A.subgroup([A(self(h)) for h in H.gens()]).gap()
        from sage.libs.gap.libgap import libgap
        if algorithm=="brute force":
            mu = libgap.function_factory("mu:=function(x,g) return(Image(g,x)); end;")
            all_subgroups = Hgap.AllSubgroups()
            subgroup_reps = G.gap().ExternalOrbitsStabilizers(all_subgroups, mu)
            # take a representative in each orbit
            reps = [orbit_sage(rep=self._subgroup_from_gap(S.Representative()),
                               stab=G.subgroup(S.Stabilizer().GeneratorsOfGroup()))
                    for S in subgroup_reps]
            reps = [r for r in reps if r.representative().cardinality() == order]
            return reps
        elif algorithm == "hulpke":
            from sage.env import SAGE_EXTCODE
            gapcode = SAGE_EXTCODE + '/gap/subgroup_orbits.g'
            libgap.Read(gapcode)
            subgroup_reps = libgap.function_factory("SubgroupRepresentatives")
            subgroup_reps = subgroup_reps(Hgap, G, order)
            subgroup_reps = [S for S in subgroup_reps if S['repr'].Size()==order]
        elif algorithm =="elementary":
            if not H.invariants()[-1].is_prime():
                raise ValueError("")
            from sage.env import SAGE_EXTCODE
            gapcode = SAGE_EXTCODE + '/gap/subgroup_orbits.g'
            libgap.Read(gapcode)
            subgroup_reps = libgap.function_factory("SubgroupRepresentatives_elementary_equiv")
            subgroup_reps = subgroup_reps(Hgap, G, order, g.gap())
            #  [dict(S) for S in
        else:
            raise ValueError("not a valid algorithm")
        return [orbit_sage(rep=self._subgroup_from_gap(S['repr']), stab=G.subgroup(S['stab'].GeneratorsOfGroup())) for S in subgroup_reps]

    def all_primitive_prime_equiv(self, other, H1, H2, G1, G2, h1, h2, glue_valuation, H10=None, H20=None, target_genus=None, qlist=None,check=False):
        r"""
        Return all totally isotropic subgroups `S` of `H1 + H2` such that
        ``H1 & S == 1`` and ``H2 & S = 1`` modulo the subgroup
        G of the orthogonal group of self

        Input:

        - ``H1``, ``H2`` -- subgroups of ``self``, ``other``
        - ``H10``, ``H20`` -- (default: ``None``) subgroups of ``H10<H1``, ``H20<H2``
        - ``G1``, ``G2`` -- subgroups of the orthogonal group of ``self``
          and ``other`` preserving ``Hi`` and ``Hi0``
        - ``h1``, ``h2`` -- elements in the center of ``G1``, ``G2``
        - ``target_genus`` -- (default: ``None``) a genus

        OUTPUT:

        A list consisting of pairs (H,G) where
        H is a totally isotropic subgroup of H1 + H2 such that H \cap Hi = \empty
        and Hi0 < projection_i(H)
        and G is a group of isometries on H^\perp/H


        EXAMPLES::

            sage: q1 = TorsionQuadraticForm(matrix.diagonal([2/3,2/27]))
            sage: q2 = TorsionQuadraticForm(matrix.diagonal([4/3,2/9]))
            sage: H1 = q1.submodule([(g.order()//3)*g for g in q1.gens()])
            sage: H2 = q2.submodule([(g.order()//3)*g for g in q2.gens()])
            sage: G1 = q1.orthogonal_group()
            sage: G2 = q2.orthogonal_group()
            sage: q1.all_primitive_prime_equiv(q2, H1, H2, G1, G2, G1.one(),G2.one(),2)
            [[Finite quadratic module over Integer Ring with invariants (3, 3)
            Gram matrix of the quadratic form with values in Q/2Z:
            [0 0]
            [0 0], Group of isometries of
            Finite quadratic module over Integer Ring with invariants (27,)
            Gram matrix of the quadratic form with values in Q/2Z:
            [8/27]
            generated by 1 elements]]
        """
        if qlist is not None:
            assert qlist[0].ncols() == glue_valuation
            assert qlist[0].ncols() >= max(len(H10.invariants()), len(H20.invariants()))

        h1 = G1.subgroup([h1]).gen(0)
        h2 = G2.subgroup([h2]).gen(0)
        if not h1 in G1:#.center():
            raise ValueError()
        if not h2 in G2:#.center():
            raise ValueError()
        glue_valuation = ZZ(glue_valuation)
        if not glue_valuation >= 0:
            raise ValueError()
        D, i1, i2, iV1, iV2 = self.direct_sum(other,return_maps=True)

        from sage.libs.gap.libgap import libgap
        from sage.env import SAGE_EXTCODE
        gapcode = SAGE_EXTCODE + '/gap/subgroup_orbits.g'
        libgap.Read(gapcode)
        OnSubgroups = libgap.function_factory("OnSubgroups")
        OD =  D.orthogonal_group(gens=())
        primitive_extensions = []
        if glue_valuation == 0:
            OD = D.orthogonal_group()
            embedG1, embedG2 = direct_sum_embed(D, i1, i2, OD, G1, G2)
            gens = embedG1.Image(G1.gap()).GeneratorsOfGroup()
            gens = gens.Concatenation(embedG2.Image(G2.gap()).GeneratorsOfGroup())
            return [[D.submodule([]),OD.subgroup(gens)]]
        if H1.cardinality()==1 or H2.cardinality()==1:
            return []
        p1 = H1.invariants()[-1]
        p2 = H2.invariants()[-1]
        if p1 != p2:
            raise ValueError("invariants do not match")
        if not p1.is_prime():
            raise ValueError("not a prime number")
        p = p1

        glue_order = p**glue_valuation
        # these may not be invariant subspaces!!! ---> crap
        # uhm really?
        subs1 = self.subgroup_representatives(H1, G1, order=glue_order, g=h1)
        subs2 = other.subgroup_representatives(H2, G2, order=glue_order, g=h2)

        # make sure they contain Hi0
        # TODO: use this fact for the computation of the subspace representatives
        # compute in Hi/Hi0
        if H10 is not None:
            subs1 = [s for s in subs1 if H10<=s.representative()]
        if H20 is not None:
            subs2 = [s for s in subs2 if H20<=s.representative()]

        subs1 = [[_normalize(s.representative().twist(-1)),s] for s in subs1]
        subs2 = [[_normalize(s.representative()),s] for s in subs2]


        G12 = G1.gap().DirectProduct(G2.gap())
        embG1 = G12.Embedding(1)
        embG2 = G12.Embedding(2)
        D1 = G1.one().gap().Source()
        D2 = G2.one().gap().Source()
        D12 = D1.DirectProduct(D2)
        embD1 = D12.Embedding(1)
        embD2 = D12.Embedding(2)
        projD1 = D12.Projection(1)
        projD2 = D12.Projection(2)

        Dg = OD.domain()
        imgs = [Dg(
                  i1(self.linear_combination_of_smith_form_gens((G1.domain()(projD1.Image(d))).exponents()))
                 +i2(other.linear_combination_of_smith_form_gens((G2.domain()(projD2.Image(d))).exponents()))
                 ).gap() for d in D12.GeneratorsOfGroup()]

        D12toDg = D12.GroupHomomorphismByImages(Dg.gap(),imgs)
        InducedAut = ("""InducedAut:=function(iso,aut)
        local f,g,emb1,emb2,proj1,proj2,aut1,aut2,MyImage;
        f:=Range(iso);
        g:=Source(iso);
        emb1 :=Embedding(g,1);
        emb2 :=Embedding(g,2);
        proj1:=Projection(g,1);
        proj2:=Projection(g,2);
        aut1:=proj1*aut[1]*emb1;
        aut2:=proj2*aut[2]*emb2;
        MyImage:=function(aut,x)
              return Image(aut1,x)*Image(aut2,x);
              end;
        aut:= GroupHomomorphismByImagesNC(f,f,
           GeneratorsOfGroup(f),
           List(GeneratorsOfGroup(f), function(i)
                  return
                   Image(iso,
                     MyImage(aut,PreImagesRepresentative(iso,i)
                       ) );
                end ));
        SetIsInjective( aut, true );
        SetIsSurjective( aut, true );
        return aut;
        end;""")
        InducedAut = libgap.function_factory(InducedAut)

        for S1 in subs1:
            S1n, S1 = S1
            q1 = S1n.gram_matrix_quadratic()
            n = len(S1n.invariants())
            if qlist is not None and q1 not in qlist:
                continue
            stab1 = S1.stabiliser()
            O1 = S1.representative().orthogonal_group()
            act1 = stab1.hom([O1(x) for x in stab1.gens()],codomain=O1,check=check)
            im1 = act1.image(stab1)
            ker1 = [embG1.Image(k.gap()) for k in act1.kernel().gens()]
            for S2 in subs2:
                S2n, S2 = S2
                q2 = S2n.gram_matrix_quadratic()
                if q1 != q2:
                    continue
                stab2 = S2.stabiliser()
                # there is a glue map
                O2 = S2.representative().orthogonal_group()

                act2 = stab2.hom([O2(x) for x in stab2.gens()],codomain=O2,check=check)
                im2 = act2.image(stab2)
                ker2 = [embG2.Image(k.gap()) for k in act2.kernel().gens()]


                A1 = O1.domain()
                A2 = O2.domain()
                gens1 = [A1(g).gap() for g in S1n.gens()]
                gens2 = [A2(g).gap() for g in S2n.gens()]

                phi = A1.gap().GroupHomomorphismByImages(A2.gap(), gens1, gens2)

                h1_on_S2 = phi.InducedAutomorphism(O1(h1).gap())
                h2_on_S2 = O2(h2).gap()
                if not O2.gap().IsConjugate(h1_on_S2, h2_on_S2):
                    # this glue map cannot be modified to be equivariant
                    continue
                else:
                    # make it equivariant
                    g0 = O2.gap().RepresentativeAction(h1_on_S2, h2_on_S2)
                    phi = phi*g0
                h1_on_S2 = phi.InducedAutomorphism(O1(h1).gap())
                assert h1_on_S2 == h2_on_S2
                center = O2.gap().Centraliser(h2_on_S2)

                stab1phi= [phi.InducedAutomorphism(O1(g).gap()) for g in stab1.gens()]
                stab1phi = center.Subgroup(stab1phi)

                stab2c = O2.subgroup(stab2.gens())
                reps = center.DoubleCosetRepsAndSizes(stab1phi,stab2c)
                for g in reps:
                    g = g[0]
                    phig = phi*g
                    g0g = O2(g0*g)
                    # graph of phig
                    gens = [i1(S1n.gen(k)) + i2(S2n.gen(k)*g0g) for k in range(n)]

                    ext = D.submodule(gens)
                    perp = D.orthogonal_submodule_to(ext)
                    disc = perp/ext
                    if target_genus is not None:
                        # it is not guaranteed that W is an instance of FreeQuadraticModule
                        if hasattr(ext.W(),"overlattice"):
                            W = ext.W()
                        else:
                            W = FreeQuadraticModule_integer_symmetric(
                                    ambient=ext.W().ambient_module(),
                                    basis=ext.W().basis(),
                                    inner_product_matrix=ext.W().inner_product_matrix(),
                                    already_echelonized=False)
                        ext_genus = W.overlattice([x.lift() for x in ext.gens()]).genus()
                        if ext_genus != target_genus:
                            continue
                    ###############################
                    # we also need the stabiliser in S1 x S2
                    # of the graph of phi
                    # this is done by hand to avoid an
                    # expensive stabilizer computation
                    im2_phi = O1.subgroup([phig*g1.gap()*phig.Inverse()
                                            for g1 in im2.gens()])
                    im = im1.subgroup(im1.intersection(im2_phi).gap().SmallGeneratingSet())
                    stab = [(act1.lift(x),
                             act2.lift(im2(phig.Inverse()*x.gap()*phig)))
                            for x in im.gens()]
                    stab = [embG1.Image(x[0].gap())*embG2.Image(x[1].gap())
                            for x in stab]
                    stab += ker1
                    stab += ker2   # tends to generate a huge group

                    # convert to sage
                    stab = [InducedAut(D12toDg,s) for s in stab]
                    stab = [matrix(
                            [disc(D.linear_combination_of_smith_form_gens(Dg(s.Image(Dg(D(x)).gap())).exponents()))
                                    for x in disc.gens()])
                            for s in stab]   # slow
                    stab = disc.orthogonal_group(tuple(stab))  # less huge
                    primitive_extensions.append([ext, stab])
        return primitive_extensions

    def all_primitive(self, other, H1, H2, G1, G2, h1=None, h2=None, glue_order=None):
        r"""
        Return all totally isotropic subgroups `S` of `H1 + H2` such that
        ``H1 & S == 1`` and ``H2 & S = 1`` modulo the subgroup
        G of the orthogonal group of self.
        Does not assume ``H1`` and ``H2`` to be vector spaces.

        Input:

        - ``H1``, ``H2`` -- subgroups of ``self``, ``other``
        - ``G1``, ``G2`` -- subgroups of ``self``, ``other``
        - ``h1``, ``h2`` -- elements in the center of ``G1``, ``G2``

        EXAMPLES::

            sage: q1 = TorsionQuadraticForm(matrix.diagonal([2/3,2/27]))
            sage: q2 = TorsionQuadraticForm(matrix.diagonal([4/3,2/9]))
            sage: G1 = q1.orthogonal_group()
            sage: G2 = q2.orthogonal_group()
            sage: extensions = q1.all_primitive(q2, q1, q2, G1, G2, G1.one(),G2.one(),9)
            sage: len(extensions)
            3
        """
        if h1 is None:
            h1 = G1.one()
        if h2 is None:
            h2 = G2.one()
        if glue_order is None:
            glue_order = ZZ.gcd(H1.cardinality(), H2.cardinality())
        h1 = G1.subgroup([h1]).gen(0)
        h2 = G2.subgroup([h2]).gen(0)
        if not h1 in G1:#.center():
            raise ValueError()
        if not h2 in G2:#.center():
            raise ValueError()
        if not glue_order > 0:
            raise ValueError()
        if not (H1.is_submodule(self) and H2.is_submodule(other)):
            raise ValueError()

        D, i1, i2, iV1, iV2 = self.direct_sum(other,return_maps=True)
        OD = D.orthogonal_group()
        embedG1, embedG2 = direct_sum_embed(D, i1, i2, OD, G1, G2)

        from sage.libs.gap.libgap import libgap
        from sage.env import SAGE_EXTCODE
        gapcode = SAGE_EXTCODE + '/gap/subgroup_orbits.g'
        libgap.Read(gapcode)
        OnSubgroups = libgap.function_factory("OnSubgroups")

        primitive_extensions = []
        if glue_order == 1 or H1.cardinality()==1 or H2.cardinality()==1:
            gens = embedG1.Image(G1.gap()).GeneratorsOfGroup()
            gens = gens.Concatenation(embedG2.Image(G2.gap()).GeneratorsOfGroup())
            return [[D.submodule([]),OD.subgroup(gens)],iV1,iV2]

        # these may not be invariant subspaces!!! ---> crap
        subs1 = self.subgroup_representatives(H1, G1, algorithm="hulpke",
                                              order=glue_order, g=h1.gap())
        subs2 = other.subgroup_representatives(H2, G2, algorithm="hulpke",
                                               order=glue_order, g=h2.gap())


        for S1 in subs1:
            for S2 in subs2:
                is_isom, isom = S1.representative().is_isomorphic_to(S2.representative().twist(-1))
                if not is_isom:
                    continue
                isom = [isom[0],[S2.representative()(g) for g in isom[1]]]

                # there are glue maps
                O1 = S1.representative().orthogonal_group()
                O2 = S2.representative().orthogonal_group()

                A1 = O1.domain()
                A2 = O2.domain()
                gens1 = [A1(g).gap() for g in isom[0]]
                gens2 = [A2((g)).gap() for g in isom[1]]

                # create a glue map
                phi = A1.gap().GroupHomomorphismByImages(A2.gap(), gens1, gens2)

                h1_on_S2 = phi.InducedAutomorphism(O1(h1,False).gap())
                h2_on_S2 = O2(h2, False).gap()
                if not O2.gap().IsConjugate(h2_on_S2, h1_on_S2):
                    # this glue map cannot be modified to be equivariant
                    continue
                else:
                    # make it equivariant
                    g = O2.gap().RepresentativeAction(h2_on_S2, h1_on_S2)
                    phi = phi*g
                assert h2_on_S2==phi.InducedAutomorphism(O1(h1).gap())

                center = O2.gap().Centraliser(h2_on_S2)
                stab1phi= [phi.InducedAutomorphism(O1(g,False).gap()) for g in S1.stabilizer().gens()]
                stab1phi = center.Subgroup(stab1phi)
                stab2c = O2.subgroup(O2(g,False) for g in S2.stabilizer().gens())
                reps = center.DoubleCosetRepsAndSizes(stab1phi,stab2c)

                for g in reps:
                    g = g[0]
                    phig = phi*g
                    g = O2(g,False)
                    gens = [i1(isom[0][k]) + i2(isom[1][k]*g) for k in range(len(isom[0]))]
                    ext = D.submodule(gens)
                    # we also need the centralizer of h1 x h2 in S1 x S2
                    ###############################
                    phig_graph = OD.domain().subgroup([OD.domain()(a) for a in gens]).gap()
                    S1_times_S2 =  [embedG1.Image(s.gap()) for s in S1.stabilizer().gens()]
                    S1_times_S2 += [embedG2.Image(s.gap()) for s in S2.stabilizer().gens()]
                    S1_times_S2 = OD.gap().Subgroup(S1_times_S2)
                    stab = S1_times_S2.Stabilizer(phig_graph, OnSubgroups).GeneratorsOfGroup()
                    ##############################
                    primitive_extensions.append([ext, OD.subgroup(stab)])
        return primitive_extensions, iV1, iV2

def _brown_indecomposable(q, p):
    r"""
    Return the Brown invariant of the indecomposable form ``q``.

    The values are taken from Table 2.1 in [Shim2016]_.

    INPUT:

    - ``q`` - an indecomposable quadratic form represented by a
      rational `1 \times 1` or `2 \times 2` matrix
    - ``p`` - a prime number

    EXAMPLES::

        sage: from sage.modules.torsion_quadratic_module import _brown_indecomposable
        sage: q = Matrix(QQ, [1/3])
        sage: _brown_indecomposable(q,3)
        6
        sage: q = Matrix(QQ, [2/3])
        sage: _brown_indecomposable(q,3)
        2
        sage: q = Matrix(QQ, [5/4])
        sage: _brown_indecomposable(q,2)
        5
        sage: q = Matrix(QQ, [7/4])
        sage: _brown_indecomposable(q,2)
        7
        sage: q = Matrix(QQ, 2, [0,1,1,0])/2
        sage: _brown_indecomposable(q,2)
        0
        sage: q = Matrix(QQ, 2, [2,1,1,2])/2
        sage: _brown_indecomposable(q,2)
        4
    """
    v = q.denominator().valuation(p)
    if p == 2:
        # brown(U) = 0
        if q.ncols() == 2:
            if q[0, 0].valuation(2) > v + 1 and q[1, 1].valuation(2) > v + 1:
                # type U
                return mod(0, 8)
            else:
                # type V
                return mod(4 * v, 8)
        u = q[0, 0].numerator()
        return mod(u + v * (u**2 - 1) / 2, 8)
    if p % 4 == 1:
        e = -1
    if p % 4 == 3:
        e = 1
    if v % 2 == 1:
        u = q[0, 0].numerator() // 2
        if legendre_symbol(u, p) == 1:
            return mod(1 + e, 8)
        else:
            return mod(-3 + e, 8)
    return mod(0, 8)

def _normalize(D):
    r"""
    INPUT:

    - a possibly degenerate torsion quadratic form over a field

    OUTPUT:

    - a TorsionQuadraticModule in normal form
    """
    if D.cardinality() == 1:
        return D
    D = D.normal_form()
    if not D.is_degenerate():
        return D
    p = D.invariants()[-1]
    if not p.is_prime():
        raise ValueError("")
    gens = list(D.gens())
    kerb = D.orthogonal_submodule_to(D)
    n = len(gens)

    kergens = [gens[i] for i in range(n) if gens[i] in kerb]
    nondeg = [gens[i] for i in range(n) if not gens[i] in kerb]
    k = len(kergens)

    for i in range(k):
        if kergens[i].q()==1:
            for j in range(i+1, k):
                if kergens[j].q() == 1:
                    kergens[j] += kergens[i]
            kergens = kergens[:i] + kergens[i+1:] + [kergens[i]]
            break
    gens = kergens + nondeg
    # translate all others to the kernel of q
    return D.submodule_with_gens(gens)

def direct_sum_embed(D, i1, i2, OD, G1, G2,as_hom=True):
    r"""

    EXAMPLES::

        sage: from sage.modules.torsion_quadratic_module import direct_sum_embed
        sage: T = TorsionQuadraticForm(matrix.diagonal([2/3]*3))
        sage: T1 = TorsionQuadraticForm(matrix.diagonal([2/3]*3))
        sage: T2 = TorsionQuadraticForm(matrix.diagonal([2/3]*2))
        sage: D, i1, i2 = T1.direct_sum(T2)
        sage: OD = D.orthogonal_group()
        sage: G1 = T1.orthogonal_group()
        sage: G2 = T2.orthogonal_group()
        sage: phi1, phi2 = direct_sum_embed(D,i1,i2,OD,G1,G2)
    """
    Dgap = OD.domain()
    direct_sum_gens = [(Dgap(D(g))).gap() for g in i1.image().gens() + i2.image().gens()]

    n1 = len(i1.image().gens())
    n2 = len(i2.image().gens())
    if as_hom:
        G1g = G1.gens()
        G2g = G2.gens()
    else:
        G1g = G1
        G2g = G2

    gensG1_imgs = []
    for f in G1g:
        imgs = [Dgap(i1(a*f)).gap() for a in i1.domain().gens()]
        imgs += direct_sum_gens[n1:]
        f = Dgap.gap().GroupHomomorphismByImagesNC(Dgap.gap(), direct_sum_gens, imgs)
        gensG1_imgs.append(f)

    gensG2_imgs = []
    for f in G2g:
        imgs = [Dgap(i2(a*f)).gap() for a in i2.domain().gens()]
        imgs = direct_sum_gens[:n1] + imgs
        f = Dgap.gap().GroupHomomorphismByImagesNC(Dgap.gap(), direct_sum_gens, imgs)
        gensG2_imgs.append(f)
    if as_hom:
        embed1 = G1.gap().GroupHomomorphismByImagesNC(OD.gap(),[g.gap() for g in G1.gens()], gensG1_imgs)
        embed2 = G2.gap().GroupHomomorphismByImagesNC(OD.gap(),[g.gap() for g in G2.gens()], gensG2_imgs)
        return embed1, embed2
    return gensG1_imgs, gensG2_imgs

def OrbitsOfSpaces(gens, k, D, H, action_hom, g=None ):
    r"""
    Return orbit representatives of the ``k`` dimensional subspaces.

    INPUT:

    - ``gens`` -- a list of matrices over a finite field
    - ``k`` -- the dimension of the subspaces
    - ``g`` -- (default:``None``) if given return only subspaces preserved by ``g``

    OUTPUT:

    - a tuple ``representative, stabiliser``
    where `representative[n]` belongs to `stabiliser[n]`
    """
    from sage.interfaces.magma import magma as m
    degree = gens[0].ncols()
    assert k<= degree
    field =  gens[0].base_ring()
    GL = m.GeneralLinearGroup(degree, field)
    G = GL.sub(gens)
    orb = m.OrbitsOfSpaces(G,k)
    if g is not None:
        g = GL(g)
        # filter the invariant subspaces
        orb = [v for v in orb if v[2] == v[2]*g]
    orb = [orbit_magma(rep=v[2], G=G,D=D,H=H,
                       action_hom=action_hom)
           for v in orb]
    return orb


class orbit(object):
    r"""
    """
    def __repr__(self):
        return "Orbit represented by \n%s"%str(self.representative())

    def representative(self):
        raise NotImplementedError()

    def list(self):
        raise NotImplementedError()

    def stabilizer(self):
        raise NotImplementedError()

    stabiliser = stabilizer


class orbit_magma(orbit):
    r"""
    """
    def __init__(self, rep, G, D, H, action_hom):
        r"""
        """
        self._G = G
        self._H = H
        self._action_hom = action_hom
        self._rep_magma = rep
        rep = rep.BasisMatrix().sage().change_ring(ZZ)
        gens = [H.linear_combination_of_smith_form_gens(r)
                for r in rep]
        rep = D.submodule(gens)
        self._rep = rep


    @cached_method
    def representative(self):
        return self._rep

    @cached_method
    def stabilizer(self,algorithm='magma'):
        stab = self._G.Stabiliser(self._rep_magma)
        stab = stab.Generators().SetToSequence()
        stab = [g.Matrix().sage().change_ring(ZZ) for g in stab]
        # we continue with magma
        G0 = self._H.orthogonal_group()
        # we could speed this up a lot by computing
        # the kernel only once and lift
        S = G0.subgroup(stab)
        stab = self._action_hom.preimage(S)
        return stab

    stabiliser = stabilizer

class orbit_sage(orbit):
    def __init__(self,rep , stab):
        self._rep = rep
        self._stab = stab

    def representative(self):
        return self._rep

    def stabilizer(self):
        return self._stab

    stabiliser = stabilizer

