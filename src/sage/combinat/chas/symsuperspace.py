r"""
Symmetric functions in super space

Symmetric functions in super space is an algebra whose bases are indexed by
``SuperPartitions``.

Let `{\mathbb Q}[\Theta_N, X_N]` be the polynomial ring in two sets of variables
`\Theta_N = \theta_1, \theta_2, \ldots, \theta_N` and
`X_N = x_1, x_2, \ldots, x_N`
where the first set of variables anti-commute, `\theta_i^2 = 0`,
and the second set
of variables commute and commutes with the first set.  The symmetric group
acts diagonally on this polynomial ring and the symmetric functions in
superspace are isomorphic to the invariants in this polynomial ring.
See [DLM2006]_ for a description of this space presented here.

When `N=2`, the following are examples of symmetric polynomials in super space

.. MATH::

    x_1^2x_2^2, \qquad \theta_1 x_1^4 + \theta_2 x_2^4, \qquad
    \theta_1 x_2^2 + \theta_2 x_1^2,  \qquad
    \theta_1 \theta_2(x_1^3 x_2 - x_1 x_2^3)

In this implementation (as with the symmetric functions), the symmetric
functions in super space are developed without referencing the variables
in the ambient polynomial ring.  Instead, the ring is represented as a vector
space with distinguished bases that are indexed by combinatorial objects
called super partitions.

.. SEEALSO::

    :class:`~sage.combinat.superpartition.SuperPartitions`,
    :class:`~sage.combinat.superpartition.SuperPartition`,
    :class:`~sage.combinat.sf.sf.SymmetricFunctions`

``SummetricFunctionsinSuperSpace`` is isomorphic to the graded ring with
generators `p_{(0;)}, p_{(;1)}, p_{(1;)}, p_{(;2)}, p_{(2;)}, \ldots`
satisfying the relations

.. MATH::

    p_{(;i)} p_{(;j)} = p_{(;j)} p_{(;i)}, \qquad
    p_{(i;)} p_{(;j)} = p_{(;j)} p_{(i;)}, \qquad
    p_{(i;)} p_{(j;)} = - p_{(j;)} p_{(i;)}, \qquad p_{(i;)}^2 = 0

There are three natural sets of generators that are analogues of the power,
complete and elementary bases of the symmetric functions.  They are indexed
by the super partition `(; n)` and `(n; )`.  Define the generators

.. MATH::

    p_{(; n)} = p_n, \quad e_{(; n)} = e_n, \quad h_{(; n)} = h_n

    p_{(; n)} = \sum_{i\geq1} x_i^n

    e_{(; n)} = \sum_{1 \leq i_1<i_2<\cdots<i_n} x_{i_1} x_{i_2} \cdots x_{i_n}

    h_{(; n)} = \sum_{1 \leq i_1\leq i_2\leq\cdots\leq i_n} x_{i_1} x_{i_2} \cdots x_{i_n}

where `p_n, e_n` and `h_n` are the symmetric function elements in the `X_N`
variables.  In addition there are fermonionic generators which include variables
from `\Theta_N`.  These are given

.. MATH::

    p_{(n; )} = \sum_{i\geq1} \theta_i x_i^n

    e_{(n; )} = \sum_{i\geq1} \sum_{J:|J|=n,i \notin J} \theta_i x_{j_1} x_{j_2} \cdots x_{j_n}

The complete fermionic generator is most clearly expressed in terms of the
monomial basis (discussed below).

The bases of the symmetric functions in super space are indexed by super
partitions `\Lambda = (\Lambda_1, \ldots, \Lambda_m; \Lambda_{m+1}, \ldots, \Lambda_N)`
where `\Lambda_i > \Lambda_{i+1} \geq 0` for `1 \leq i \leq m-1` and
`\Lambda_j \geq \Lambda_{j+1} \geq 0` for `m+1 \leq j \leq N-1`.  The
degree (or bosonic degree) of `\Lambda` is `|\Lambda| = \sum_{i=1}^N \Lambda_i`
and the fermionic degree (or sector) is `m`.  The super partitions of degree
`n` and bosonic degree `m` will be denoted `SPar(n|m)`.
The power, elementary and complete basis elements are defined by the
following products of the generators

.. MATH::

    p_\Lambda = p_{(\Lambda_1; 0)} \cdots p_{(\Lambda_m; 0)}
    p_{(; \Lambda_{m+1})} \cdots p_{(; \Lambda_N)}~,

    e_\Lambda = e_{(\Lambda_1; 0)} \cdots e_{(\Lambda_m; 0)}
    e_{(; \Lambda_{m+1})} \cdots e_{(; \Lambda_N)}~,

    h_\Lambda = h_{(\Lambda_1; 0)} \cdots h_{(\Lambda_m; 0)}
    h_{(; \Lambda_{m+1})} \cdots h_{(; \Lambda_N)}~.

The monomial symmetric polynomial `m_\Lambda[X_N;\Theta_N]` indexed by the
super partition `\Lambda \in SPar(n|m)` is equal to the sum over all distinct
rearrangments of the monomial

.. MATH::

    \theta_1 \theta_2 \cdots \theta_m x_1^{\Lambda_1} \cdots x_N^{\Lambda_N}

where the indicies of both alphabets are permuted simultaneously.  The basis
element of the symmetric functions in super space will be denoted `m_\Lambda`
and will have the same coefficients in the products as `m_\Lambda[X_N;\Theta_N]`
as long as `N` is sufficiently large.

The complete homogeneous basis is then defined as

.. MATH::

    h_{(n; )} = \sum_{\Lambda \in SPar(n|1)} m_\Lambda, \qquad
    h_{(; n)} = \sum_{\lambda \vdash n} m_\lambda~.

The space of symmetric functions in super space is bigraded by the bosonic and
fermionic degree.

There are several acceptable input formats for the index set of basis
elements.  A super partition can be represented as a pair consisting
of a strict partition and a partition.  The entries can also be expressed
as a list of negative entries or zero for the fermionic parts (and hence
must all be distinct) followed by positive entries for the bosonic parts.
Alternatively, the input may also be a ``Partition``.  The default display
of basis elements is with the fermionic parts first, followed by a semicolon,
followed by the bosonic parts.  This can be changed through the options
in super partitions.::

    sage: SymmetricFunctionsinSuperSpace(QQ).inject_shorthands(verbose=False)
    sage: h(((2,1),(1,1)))
    h[2, 1; 1, 1]
    sage: p([[2,1],[1,1]])
    p[2, 1; 1, 1]
    sage: e[-1,0,2,1]
    e[1, 0; 2, 1]
    sage: e[[1,0],[2,1]]
    e[1, 0; 2, 1]
    sage: e(Partition([2,1]))
    e[; 2, 1]
    sage: SuperPartitions.options.display = "list"
    sage: h[[3,1,0],[2,2,1]]
    h[-3, -1, 0, 2, 2, 1]
    sage: SuperPartitions.options.display = "pair"
    sage: h[[3,1,0],[2,2,1]]
    h[[3, 1, 0], [2, 2, 1]]
    sage: SuperPartitions.options.display = "default"
    sage: h[[3,1,0],[2,2,1]]
    h[3, 1, 0; 2, 2, 1]

There are four bases of the symmetric functions in super space which are
analogues of the Schur functions.  These bases were introduced in [JL2016]_
and are related to specializations of an analogue of the Macdonald symmetric
functions.

The shorthand for these bases are ``s``, ``sb``, ``ss`` and ``ssb`` for the
Schur, Schur-bar, Schur-star and Schur-star-bar bases.  All of the bases
are positive in the monomial basis except for the Schur-star-bar basis.

For `\lambda \vdash n`, it is the case that `s_\lambda = s^\ast_\lambda =
{\bar s}_\lambda = {\bar s}^\ast_\lambda`, but for a general super partition
`\Lambda` the bases are different.

::

    sage: SymmetricFunctionsinSuperSpace(QQ).inject_shorthands(verbose=False)
    sage: all(s(la)==ss(la) for la in SuperPartitions(5,0))
    True
    sage: all(s(la)==sb(la) for la in SuperPartitions(5,0))
    True
    sage: all(s(la)==ssb(la) for la in SuperPartitions(5,0))
    True
    sage: m(s[-2,0,1])
    m[2, 0; 1]
    sage: m(sb[-2,0,1])
    2*m[1, 0; 1, 1] + m[1, 0; 2] + m[2, 0; 1]
    sage: m(ss[-2,0,1])
    2*m[1, 0; 1, 1] + m[1, 0; 2] + 3*m[2, 0; 1] + 2*m[2, 1; ] + 2*m[3, 0; ]
    sage: m(ssb[-2,0,1])
    -m[2, 0; 1] - 2*m[2, 1; ] - 2*m[3, 0; ]
    sage: m(ssb[-2,1])
    2*m[0; 1, 1, 1] + m[0; 2, 1] + 3*m[1; 1, 1] + m[1; 2] + 2*m[2; 1] - m[3; ]

There is a scalar product defined so that for super partitions
`\Lambda = (\Lambda^a; \Lambda^s)` and `\Gamma`,
`\left<\!\left< p_\Lambda, p_\Gamma \right>\!\right> = \delta_{\Lambda\Gamma} z_{\Lambda^s}`
where `z_{\Lambda^s}` is the usual constant `z_\lambda` which is the size of
the centralizer of an element of cycle type `\lambda`.

.. SEEALSO::

    :meth:`~sage.combinat.partition.Partition.centralizer_size`,
    :meth:`~sage.combinat.superpartition.SuperPartition.zee`

Given this definition, the complete and monomial bases of the symmetric
functions in super space are dual.

::

    sage: SymmetricFunctionsinSuperSpace(QQ).inject_shorthands(verbose=False)
    sage: matrix([[h(la).scalar(m(mu)) for mu in SuperPartitions(2)]
    ....: for la in SuperPartitions(2)]).is_one()
    True
    sage: matrix([[p(la).scalar(p(mu)) for mu in SuperPartitions(2)]
    ....: for la in SuperPartitions(2)])
    [2 0 0 0 0 0 0 0]
    [0 2 0 0 0 0 0 0]
    [0 0 1 0 0 0 0 0]
    [0 0 0 1 0 0 0 0]
    [0 0 0 0 2 0 0 0]
    [0 0 0 0 0 2 0 0]
    [0 0 0 0 0 0 1 0]
    [0 0 0 0 0 0 0 1]

The four Schur bases are related by duality with respect to this scalar product
in pairs.  For all `\Lambda, \Gamma \in SPar(n|m)` ,

.. MATH::

    \left<\!\left< s_\Lambda, s^\ast_\Gamma \right>\!\right>
    = \left<\!\left< {\bar s}_\Lambda, {\bar s}^\ast_\Gamma \right>\!\right>
    = \delta_{\Lambda\Gamma}

::

    sage: SymmetricFunctionsinSuperSpace(QQ).inject_shorthands(verbose=False)
    sage: all(s(la).scalar(ss(la)) for la in SuperPartitions(5,2))
    True
    sage: all(sb(la).scalar(ssb(la)) for la in SuperPartitions(5,2))
    True

As with the space of symmetric functions, there is an involution `\omega`
that for all super partitions `\Lambda`, satisfies `\omega(h_\Lambda) = e_\Lambda`.
This involution is an algebra morphism with
`\omega(p_{(; r)}) = (-1)^{r-1} p_{(; r)}` and
`\omega(p_{(r; )}) = (-1)^r p_{(r; )}`.  The Schur bases are related by

.. MATH::

    s^\ast_\Lambda = (-1)^{\binom{m}{2}} \omega {\bar s}_{\Lambda'} \qquad
    and \qquad {\bar s}^\ast_\Lambda = (-1)^{\binom{m}{2}} \omega s_{\Lambda'}

where `\Lambda` is a super partition with fermionic sector `m`.

::

    sage: SymmetricFunctionsinSuperSpace(QQ).inject_shorthands(verbose=False)
    sage: all(e(la).omega()==h(la) for la in SuperPartitions(5,2))
    True
    sage: [f.omega() for f in [p[r] for r in range(-4,5)]]
    [p[4; ], -p[3; ], p[2; ], -p[1; ], p[0; ], p[; 1], -p[; 2], p[; 3], -p[; 4]]
    sage: [ssb(s(la).omega()) for la in SuperPartitions(3,2)]
    [-ssb[1, 0; 2], -ssb[1, 0; 1, 1], -ssb[2, 0; 1], -ssb[2, 1; ], -ssb[3, 0; ]]
    sage: [ss(sb(la).omega()) for la in SuperPartitions(3,2)]
    [-ss[1, 0; 2], -ss[1, 0; 1, 1], -ss[2, 0; 1], -ss[2, 1; ], -ss[3, 0; ]]

The algebra of symmetric functions in super space is also a (bi-)graded
Hopf algebra and the coproduct is defined by declaring that the power sum
generators are primitive. ::

    sage: SymmetricFunctionsinSuperSpace(QQ).inject_shorthands(verbose=False)
    sage: p[3].coproduct()
    p[; ] # p[; 3] + p[; 3] # p[; ]
    sage: p[-3].coproduct()
    p[; ] # p[3; ] + p[3; ] # p[; ]
    sage: sb[-2,1].coproduct()
    sb[; ] # sb[2; 1] + sb[; 1] # sb[1; 1] + sb[; 1] # sb[2; ] + sb[; 1, 1]
    # sb[1; ] + sb[; 2] # sb[1; ] + sb[; 2, 1] # sb[0; ] + sb[0; ] # sb[; 2, 1]
    + sb[1; ] # sb[; 1, 1] + sb[1; ] # sb[; 2] + sb[1; 1] # sb[; 1] + sb[2; ]
    # sb[; 1] + sb[2; 1] # sb[; ]
    sage: [f.antipode() for f in [p[r] for r in range(-3,4)]]
    [-p[3; ], -p[2; ], -p[1; ], -p[0; ], -p[; 1], -p[; 2], -p[; 3]]
    sage: ss(sb[-2,2].antipode())
    -ss[0; 2, 2]

TESTS::

    sage: SymmetricFunctionsinSuperSpace(QQ).inject_shorthands(verbose=False)
    sage: bases = [m,p,e,h,s,ss,sb,ssb]
    sage: all(b1(b2(b1([[3,2],[1,1]]))) == b1([[3,2],[1,1]])
    ....:     for b1 in bases for b2 in bases if b1!=b2) # long time (17s)
    True

AUTHORS:

- Mike Zabrocki (May 2018) - initial implementation
"""
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.rings import Rings
from sage.categories.fields import Fields
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.permutation import Permutation
from sage.sets.family import Family
from sage.misc.bindable_class import BindableClass
from sage.rings.integer import Integer
from sage.combinat.ncsf_qsym.generic_basis_code import AlgebraMorphism
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.combinat.superpartition import SuperPartition, SuperPartitions
from sage.combinat.sf.sf import SymmetricFunctions
from sage.misc.misc_c import prod
from sage.combinat.partition import Partition, Partitions
from sage.functions.other import factorial
from sage.sets.set import Set
from sage.combinat.subset import Subsets
from sage.combinat.permutation import Permutations
from sage.misc.misc import uniq
from sage.functions.other import binomial
from sage.rings.finite_rings.integer_mod_ring import Integers

class SFSuperSpaceAlgebraMorphism(AlgebraMorphism):
    r"""
    An algebra morphism on symmetric functions in super space.

    The generators are indexed by integers, non-positive integers are
    indices of the fermionic generators and the positive integers
    are mapped to the bosonic generators.
    """
    def _on_basis(self, sp):
        r"""
        Computes the image of this morphism on the basis element indexed by
        ``sp``.

        INPUT:

        - ``sp`` -- a SuperPartition

        OUTPUT:

        - element of the codomain

        EXAMPLES::

            sage: SFSS = SymmetricFunctionsinSuperSpace(QQ)
            sage: p = SFSS.p()
            sage: h = SFSS.h()
            sage: from sage.combinat.chas.symsuperspace import SFSuperSpaceAlgebraMorphism
            sage: f = SFSuperSpaceAlgebraMorphism(p, lambda i : h[[],[abs(i)]], codomain=h)
            sage: f._on_basis([[1,0],[3,2]])
            h[; 3, 2, 1]
            sage: f = SFSuperSpaceAlgebraMorphism(p, lambda i : h[[abs(i)],[]], codomain=h)
            sage: f._on_basis([[0],[3, 1]])
            h[3, 1, 0; ]
            sage: f._on_basis([[1],[1]])
            0
            sage: f = SFSuperSpaceAlgebraMorphism(p, lambda i : h[[abs(i)],[]], codomain=h, anti=True)
            sage: f._on_basis([[0],[3, 1]])
            -h[3, 1, 0; ]
        """
        c = ([-a for a in sp[0]]+[a for a in sp[1]])[::1-2*int(self._anti)]
        return self.codomain().prod(self._on_generators(i) for i in c)

def valid_fermionic_matchings(sp1, sp2, a_sp3):
    r"""
    Generator for ways of matching the fermionic parts.

    For a fixed fermionic part ``a_sp3`` such that the length is equal to the
    sums of the lengths of the fermionic parts of ``sp1`` and ``sp2``, find
    all ways of covering ``a_sp3`` with the fermionic parts of ``sp1`` and ``sp2``.
    The output of this function will be a list of matchings, where a matching
    is a list of 4 objects:

    - a permutation of a subset of the parts of ``a_sp3``
    - the difference of the bosonic parts of ``sp1`` and this permutation
    - a permutation of the complement of the subset of parts of ``a_sp3``
    - the difference of the complement and the bosonic parts of ``sp2``

    OUTPUT:

    - a generator of matchings

    EXAMPLES::

        sage: from sage.combinat.chas.symsuperspace import valid_fermionic_matchings
        sage: [p for p in valid_fermionic_matchings([[1,0],[1]], [[0],[2,1,1]], [2,1,0])]
        [[[1, 2], [0, 2], [0], [0]],
         [[2, 1], [1, 1], [0], [0]],
         [[2, 0], [1, 0], [1], [1]]]
        sage: [p for p in valid_fermionic_matchings([[1,0],[1]], [[0],[2,1,1]], [3,1,0])]
        [[[3, 1], [2, 1], [0], [0]], [[3, 0], [2, 0], [1], [1]]]
    """
    for SS in Subsets(a_sp3,len(sp1[0])): #pick a subset of fermionic part of a_sp3
        for pi in Permutations(SS): #rearrange them in all possible ways
            vs = [pi[i]-sp1[0][i] for i in range(len(sp1[0]))]
            if min(vs+[0])>=0 and all(vs.count(v)<=sp2[1].count(v) for v in vs if v>0):
                for tau in Permutations(Set(a_sp3).difference(SS)):
                    ws = [tau[i]-sp2[0][i] for i in range(len(sp2[0]))]
                    if min(ws+[0])>=0 and all(ws.count(w)<=sp1[1].count(w) for w in ws if w>0):
                        yield [pi,vs,tau,ws]

def remove_parts(la, SS):
    r"""
    Remove multiple parts from a partition.

    OUTPUT:

    - a list of integers

    EXAMPLES::

        sage: from sage.combinat.chas.symsuperspace import remove_parts
        sage: remove_parts([4,4,3,3,3,2,2,2,1,1],[3,3,2,1,1])
        [4, 4, 3, 2, 2]
    """
    return sum(([k]*(la.count(k)-SS.count(k)) for k in uniq(la)[::-1]),[])

def m_prod_fixed_a_match(SFSS_m, sp1, sp2, a_sp3, mt):
    r"""
    The terms with fixed fermionic part and matching in a monomial product.

    The matching ``mt`` is a list consisting of permutations of the
    parts of ``a_sp3`` and lists indicating which bosonic parts will
    be used from ``sp1`` and ``sp2`` in the parts of ``a_sp3``.  There is
    one term in the sum for each term in the product of the remaining bosonic
    parts of ``sp1`` and ``sp2`` considered as monomial symmetric functions.

    INPUT:

    - ``SFSS_m`` -- a monomial basis
    - ``sp1, sp2`` -- super partitions
    - ``a_sp3`` -- a strict partition
    - ``mt`` -- a list

    OUTPUT:

    - an element in the ``SFSS_m`` basis

    Below we have calculated the terms in example (2.41) and (2.42) from
    [DLM2006]_ in detail.

    EXAMPLES::

        sage: m = SymmetricFunctionsinSuperSpace(QQ).m()
        sage: from sage.combinat.chas.symsuperspace import m_prod_fixed_a_match
        sage: m_prod_fixed_a_match(m,[[1,0],[1]],[[0],[2,1,1]],[2,1,0],[[1,2],[0,2],[0],[0]])
        3*m[2, 1, 0; 1, 1, 1] + m[2, 1, 0; 2, 1]
        sage: m_prod_fixed_a_match(m,[[1,0],[1]],[[0],[2,1,1]],[2,1,0],[[2,1],[1,1],[0],[0]])
        m[2, 1, 0; 2, 1] + m[2, 1, 0; 3]
        sage: m_prod_fixed_a_match(m,[[1,0],[1]],[[0],[2,1,1]],[2,1,0],[[2,0],[1,0],[1],[1]])
        m[2, 1, 0; 2, 1]
        sage: m_prod_fixed_a_match(m,[[1,0],[1]],[[0],[2,1,1]],[3,1,0],[[3,1],[2,1],[0],[0]])
        2*m[3, 1, 0; 1, 1] + m[3, 1, 0; 2]
        sage: m_prod_fixed_a_match(m,[[1,0],[1]],[[0],[2,1,1]],[3,1,0],[[3,0],[2,0],[1],[1]])
        m[3, 1, 0; 1, 1]
    """
    m = SFSS_m._Sym_m
    ga = remove_parts(sp2[1], mt[1])
    nu = remove_parts(sp1[1], mt[3])
    return SFSS_m.sum(SFSS_m.term(SuperPartition([a_sp3,list(la)]),c) \
        for (la,c) in m(ga)*m(nu))

def sign_matching(sp1, sp2, mt):
    r"""
    The sign for a matching between ``sp1`` and ``sp2``.

    INPUT:

    - ``sp1, sp2`` -- super partitions
    - ``mt`` -- a list

    OUTPUT:

    - the sign of the permutation

    EXAMPLES::

        sage: from sage.combinat.chas.symsuperspace import sign_matching
        sage: sign_matching([[1,0],[1]],[[0],[2,1,1]],[[1,2],[0,2],[0],[0]])
        -1
        sage: sign_matching([[1,0],[1]],[[0],[2,1,1]],[[2,1],[1,1],[0],[0]])
        1
        sage: sign_matching([[1,0],[1]],[[0],[2,1,1]],[[2,0],[1,0],[1],[1]])
        -1
        sage: sign_matching([[1,0],[1]],[[0],[2,1,1]],[[3,1],[2,1],[0],[0]])
        1
        sage: sign_matching([[1,0],[1]],[[0],[2,1,1]],[[3,0],[2,0],[1],[1]])
        -1
    """
    tp = list(mt[0])+list(mt[2])
    prm2 = Permutation([j for (i,j) in sorted([(tp[-k-1],k+1) for k in range(len(tp))])])
    return prm2.sign()

def m_m_mult_fixed_a(SFSS_m, sp1, sp2, a_sp3):
    r"""
    For a fixed a-matching, return the sum of the terms in the monomial product.

    INPUT:

    - ``SFSS_m`` -- a monomial basis
    - ``sp1, sp2`` -- super partitions
    - ``a_sp3`` -- a strict partition

    OUTPUT:

    - an element of the monomial basis

    EXAMPLES::

        sage: from sage.combinat.chas.symsuperspace import m_m_mult_fixed_a
        sage: m = SymmetricFunctionsinSuperSpace(QQ).m()
        sage: m_m_mult_fixed_a(m,[[1,0],[1]],[[0],[2,1,1]],[2,1,0])
        -3*m[2, 1, 0; 1, 1, 1] - m[2, 1, 0; 2, 1] + m[2, 1, 0; 3]
        sage: m_m_mult_fixed_a(m,[[1,0],[1]], [[0],[2,1,1]], [3,1,0])
        m[3, 1, 0; 1, 1] + m[3, 1, 0; 2]
    """
    if len(sp1[0])+len(sp2[0])!=len(a_sp3):
        return SFSS_m.zero()
    if sum(sp1[0])+sum(sp2[0])>sum(a_sp3):
        return SFSS_m.zero()
    return SFSS_m.sum(sign_matching(sp1, sp2, mt)*m_prod_fixed_a_match(SFSS_m, sp1, sp2, a_sp3, mt)
            for mt in valid_fermionic_matchings(sp1, sp2, a_sp3))

class SymmetricFunctionsinSuperSpace(UniqueRepresentation, Parent):
    def __init__(self, R):
        r"""
        TESTS::

            sage: SFSS1 = SymmetricFunctionsinSuperSpace(FiniteField(23))
            sage: SFSS2 = SymmetricFunctionsinSuperSpace(Integers(23))
            sage: TestSuite(SymmetricFunctionsinSuperSpace(QQ)).run()
        """
        assert(R in Fields() or R in Rings())
        self._base = R # Won't be needed once CategoryObject won't override base_ring
        Parent.__init__(self, category = HopfAlgebras(R).Graded().Connected().WithRealizations())
        h = self.Complete()
        e = self.Elementary()
        p = self.Power()
        m = self.Monomial()
        s = self.Schur()
        ss = self.Schur_s()
        sb = self.Schur_b()
        ssb = self.Schur_sb()

        # coercions between multiplicative bases
        for b1 in ['e','h','p']:
            for b2 in ['e','h','p']:
                if b1!=b2:
                    eval(b1+'.algebra_morphism('+b2+'.'+b1+\
                        '_to_self_on_generator, codomain = '\
                        +b2+').register_as_coercion()')
        p.module_morphism(m.power_to_self_on_basis, codomain = m).register_as_coercion()
        m.module_morphism(function=p.monomial_to_self_by_triangularity, codomain = p).register_as_coercion()
        ss.module_morphism(ss.self_to_complete_on_basis, codomain = h).register_as_coercion()
        h.module_morphism(function=ss.complete_to_self_by_triangularity, codomain = ss).register_as_coercion()
        ssb.module_morphism(ssb.self_to_complete_on_basis, codomain = h).register_as_coercion()
        h.module_morphism(function=ssb.complete_to_self_by_triangularity, codomain = ssb).register_as_coercion()
        sb.module_morphism(sb.self_to_elementary_on_basis, codomain = e).register_as_coercion()
        e.module_morphism(function=sb.elementary_to_self, codomain = sb).register_as_coercion()
        s.module_morphism(s.self_to_elementary_on_basis, codomain = e).register_as_coercion()
        e.module_morphism(function=s.elementary_to_self, codomain = s).register_as_coercion()

    def a_realization(self):
        r"""
        A basis of the symmetric functions in super space.

        This particular realization is the Power basis.

        EXAMPLES::

            sage: SymmetricFunctionsinSuperSpace(QQ).a_realization()
            Symmetric Functions in super space over the Rational Field in the Power basis
        """
        return self.Power()

    def _repr_(self):
        r"""
        The name of the symmetric functions in super space over a base ring.

        EXAMPLES::

            sage: SFSS = SymmetricFunctionsinSuperSpace(ZZ)
            sage: SFSS._repr_()
            'Symmetric Functions in super space over the Integer Ring'
        """
        return "Symmetric Functions in super space over the %s"%self.base_ring()

    _shorthands = tuple(['e', 'h', 'p', 'm', 'ss', 'sb', 'ssb', 's'])

    class Bases(Category_realization_of_parent):
        """
        Category of bases of symmetric functions in super space.

        EXAMPLES::

            sage: SFSS = SymmetricFunctionsinSuperSpace(QQ)
            sage: SFSS.Bases()
            Category of bases of Symmetric Functions in super space over the Rational Field
            sage: p = SFSS.Power()
            sage: p in SFSS.Bases()
            True
        """
        def super_categories(self):
            r"""
            Return the super categories of the category of bases of
            symmetric functions in super space.

            OUTPUT:

            - a list

            TESTS::

                sage: SFSS = SymmetricFunctionsinSuperSpace(QQ)
                sage: SFSS.Bases().super_categories()
                [Category of realizations of Symmetric Functions in super space over the Rational Field,
                 Join of Category of realizations of hopf algebras over Rational Field and Category of graded algebras over Rational Field,
                 Category of graded connected hopf algebras with basis over Rational Field]
            """
            R = self.base().base_ring()
            cat = HopfAlgebras(R).Graded().WithBasis().Graded()
            return [self.base().Realizations(),
                    HopfAlgebras(R).Graded().Realizations(), cat.Connected()]

        class ParentMethods:
            def one_basis(self):
                r"""
                The empty super-partition.

                OUTPUT:

                - The empty super-partition.

                EXAMPLES::

                    sage: p=SymmetricFunctionsinSuperSpace(QQ).Power()
                    sage: p.one_basis()
                    [; ]
                    sage: p.one()
                    p[; ]
                """
                return SuperPartitions()([[],[]])

            def __getitem__(self, c, *rest):
                r"""
                This method implements the abuses of notations::

                    sage: p = SymmetricFunctionsinSuperSpace(QQ).Power()
                    sage: p[[2],[1]]
                    p[2; 1]
                    sage: p[[[2],[1]]]
                    p[2; 1]
                    sage: p[SuperPartition([[2],[1]])]
                    p[2; 1]
                    sage: p[-2, 1, 0]
                    p[2, 0; 1]
                """
                if isinstance(c, SuperPartition):
                    assert len(rest) == 0
                else:
                    if isinstance(c, (int, Integer)):
                        c = self._indices([c])
                    elif list(c) in self._indices:
                        c = self._indices(list(c))
                    else:
                        raise ValueError("%s must be a %s"%(c,self._indices))
                return self.monomial(c)

            def _repr_(self):
                r"""
                TESTS::

                    sage: h = SymmetricFunctionsinSuperSpace(QQ).Complete()
                    sage: h._repr_()
                    'Symmetric Functions in super space over the Rational Field in the Complete basis'
                """
                return "%s in the %s basis" % (self.realization_of(), self._realization_name())

            def degree_on_basis(self, sp):
                r"""
                The degree of a basis element is the degree of the super partition.

                EXAMPLES::

                    sage: h = SymmetricFunctionsinSuperSpace(QQ).Complete()
                    sage: h.degree_on_basis(SuperPartition([[2,1],[3,1,1]]))
                    8
                    sage: h.degree_on_basis(SuperPartition([[],[]]))
                    0
                """
                return sp.degree()

            def bi_degree_on_basis(self, sp):
                r"""
                The bi-degree of the super partition is the bosonic and fermionic degree.

                EXAMPLES::

                    sage: h = SymmetricFunctionsinSuperSpace(QQ).Complete()
                    sage: h.bi_degree_on_basis(SuperPartition([[2,1],[3,1,1]]))
                    (8, 2)
                    sage: h.bi_degree_on_basis(SuperPartition([[],[]]))
                    (0, 0)
                """
                return sp.bi_degree()

        class ElementMethods:
            def is_homogeneous_bi_degree(self):
                r"""
                Test if the support has elements of the same bi-degree.

                EXAMPLES::

                    sage: h = SymmetricFunctionsinSuperSpace(QQ).Complete()
                    sage: (h[[1],[2,1]]+h[[1,0],[3]]).is_homogeneous_bi_degree()
                    False
                    sage: (h[[1],[2,1]]+h[[0],[4]]).is_homogeneous_bi_degree()
                    True
                """
                bi_degree_on_basis = self.parent().bi_degree_on_basis
                bi_degree = None
                for m in self.support():
                    if bi_degree is None:
                        bi_degree = bi_degree_on_basis(m)
                    else:
                        if bi_degree != bi_degree_on_basis(m):
                            return False
                return True

            def bi_degree(self):
                r"""
                The bi-degree element if the element is homogeneous bi_degree.

                EXAMPLES::

                    sage: h = SymmetricFunctionsinSuperSpace(QQ).Complete()
                    sage: (h[[1],[2,1]]+h[[1,0],[3]]).bi_degree()
                    Traceback (most recent call last):
                    ...
                    ValueError: element is not of homogeneous bi-degree
                    sage: (h[[1],[2,1]]+h[[0],[4]]).bi_degree()
                    (4, 1)
                """
                if not self.support():
                    raise ValueError("the zero element does not have a well-defined bi-degree")
                if not self.is_homogeneous_bi_degree():
                    raise ValueError("element is not of homogeneous bi-degree")
                return self.parent().bi_degree_on_basis(self.leading_support())

            def maximal_bi_degree(self):
                r"""
                The maximal bi-degree over the support of an element.

                EXAMPLES::

                    sage: h = SymmetricFunctionsinSuperSpace(QQ).Complete()
                    sage: (h[[1],[2,1]]+h[[1,0],[3]]).maximal_bi_degree()
                    (4, 2)
                    sage: (h[[1],[2,1]]+h[[0],[4]]).maximal_bi_degree()
                    (4, 1)
                """
                if self.is_zero():
                    raise ValueError("the zero element does not have a well-defined degree")
                bi_degree_on_basis = self.parent().bi_degree_on_basis
                return max(bi_degree_on_basis(m) for m in self.support())

            def all_bi_degrees(self):
                r"""
                A list of all the bi-degrees in the support of the element.

                OUTPUT:

                - a list of tuples of integers

                EXAMPLES::

                    sage: h = SymmetricFunctionsinSuperSpace(QQ).h()
                    sage: (h[0,3,2]+h[-1,1]).all_bi_degrees()
                    [(2, 1), (5, 1)]
                    sage: (h[[]]+h[0]+h[-1]).all_bi_degrees()
                    [(0, 0), (0, 1), (1, 1)]
                """
                bi_degree_on_basis = self.parent().bi_degree_on_basis
                return uniq(bi_degree_on_basis(m) for m in self.support())

            def homogeneous_bi_degree_component(self, bd):
                r"""
                Take the component of fixed bi_degree.

                INPUT:

                - ``bd`` -- a tuple of integers corresponding to a bidegree

                OUTPUT:

                - a supersymmetric function

                EXAMPLES::

                    sage: h = SymmetricFunctionsinSuperSpace(QQ).h()
                    sage: (h[1,1]+h[-1,2]+h[-2,1]).homogeneous_bi_degree_component((3,1))
                    h[1; 2] + h[2; 1]
                    sage: (h[1,1]+h[-1,2]+h[-2,1]).homogeneous_bi_degree_component((2,0))
                    h[; 1, 1]
                    sage: (h[1,1]+h[-1,2]+h[-2,1]).homogeneous_bi_degree_component((1,0))
                    0
                """
                b = self.parent()
                return b.sum(c*b(sp) for (sp,c) in self if sp.bi_degree()==bd)

            def scalar(self, other):
                r"""
                The scalar product between ``self`` and ``other``.

                The power sum elements are orthonormal.

                INPUT:

                - ``other`` -- an element of symmetric functions in super space

                OUTPUT:

                - an element of the base ring

                EXAMAPLES::

                    sage: h = SymmetricFunctionsinSuperSpace(QQ).h()
                    sage: h[[],[3,1]].scalar(h[[1],[2,1]])
                    0
                    sage: h[[1],[2,1]].scalar(h[[1],[2,1]])
                    9
                """
                p = self.parent().realization_of().Power()
                px = p(self)
                py = p(other)
                return sum(sp.zee()*c*py.coefficient(sp) for (sp,c) in px)

            def skew_by(self, other, side='right'):
                r"""
                Compute the ``skew_by`` operation by coersion to the power basis.

                INPUT:

                - ``other`` -- a supersymmetric function

                - ``side`` -- either 'right' or 'left' (default: 'right')

                OUTPUT:

                - a supersymmetric function

                EXAMPLES::

                    sage: h = SymmetricFunctionsinSuperSpace(QQ).h()
                    sage: (h[2]+h[1]+h[[]]).skew_by(h[1]+h[[]])
                    2*h[; ] + 2*h[; 1] + h[; 2]
                    sage: (h[2]+h[1]+h[[]]).skew_by(h[2])
                    h[; ]
                    sage: (h[2]+h[1]+h[[]]).skew_by(h[0])
                    0
                    sage: (h[-1,1]+h[-1]+h[-2]).skew_by(h[0])
                    h[; 1] + h[; 1, 1] + h[; 2]
                """
                p = self.parent().realization_of().Power()
                return self.parent()(p(self).skew_by(other, side))

            def omega(self):
                r"""
                The involution omega.

                The involution omega is computed by coercion to the power
                sum basis.

                EXAMPLES::

                    sage: SFSS = SymmetricFunctionsinSuperSpace(QQ)
                    sage: (h,e) = (SFSS.h(), SFSS.e())
                    sage: all(h(sp).omega()==e(sp) for sp in SuperPartitions(3))
                    True
                    sage: e[[1,0],[2,2,1]].omega()
                    -e[1, 0; 1, 1, 1, 1, 1] + 2*e[1, 0; 2, 1, 1, 1] - e[1, 0; 2, 2, 1]
                """
                p = self.parent().realization_of().Power()
                return self.parent()(p(self).omega())

    class MultiplicativeBases(Category_realization_of_parent):
        def super_categories(self):
            r"""
            Return the super categories of the category of multiplicative
            bases of the summetric functions in super space.

            OUTPUT:

            - a list

            EXAMPLES::

                sage: SFSS = SymmetricFunctionsinSuperSpace(QQ)
                sage: SFSS.MultiplicativeBases().super_categories()
                [Category of bases of Symmetric Functions in super space over the Rational Field]
            """
            return [self.base().Bases()]

        class ParentMethods:
            def algebra_generators(self):
                r"""
                A family of generators associated to a multiplicative basis.

                The generators are ordered so that the non-positive integers
                are mapped to a fermionic generator and the positive
                integers are mapped to the bosonic generator.

                EXAMPLES::

                    sage: e = SymmetricFunctionsinSuperSpace(QQ).Elementary()
                    sage: e.algebra_generators()
                    Lazy family (pos_int_to_gen(i))_{i in Integer Ring}
                    sage: it = e.algebra_generators().__iter__()
                    sage: [next(it) for i in range(6)]
                    [e[0; ], e[; 1], e[1; ], e[; 2], e[2; ], e[; 3]]
                """
                def pos_int_to_gen(n):
                    if n<=0:
                        return self.monomial(self._indices([[abs(n)],[]]))
                    else:
                        return self.monomial(self._indices([[],[n]]))
                return Family(Integers(), pos_int_to_gen)

            def algebra_morphism(self, on_generators, **keywords):
                """
                Extend a map on the generators of a multiplicative basis to
                the whole basis.

                Given a map defined on the generators of the multiplicative
                basis ``self``, return the algebra morphism that extends
                this map to the whole algebra of symmetric functions in super
                space.

                INPUT:

                - ``on_generators`` -- a function defined on the index set of
                  the generators (that is, on the positive integers)
                - ``codomain`` -- the codomain of ``on_generators``
                - ``anti`` -- a boolean; defaults to ``False``
                - ``category`` -- a category; defaults to ``None``

                OUTPUT:

                - The algebra morphism of ``self`` which is defined by
                  ``on_generators`` in the basis ``self``. When ``anti``
                  is set to ``True``, an algebra anti-morphism is
                  computed instead of an algebra morphism.

                EXAMPLES::

                    sage: SFSS = SymmetricFunctionsinSuperSpace(QQ)
                    sage: p = SFSS.p()
                    sage: neg = lambda i: p.algebra_generators()[-i]
                    sage: f = p.algebra_morphism(neg, codomain = p)
                    sage: f
                    Generic endomorphism of Symmetric Functions in super space over the Rational Field in the Power basis
                    sage: f(2*p[[1],[2]] + 3 * p[[2,1],[]] + p[[2,0],[]])
                    3*p[; 2, 1] + p[0; 2] + 2*p[2; 1]
                    sage: f.category()
                    Category of endsets of unital magmas and right modules over Rational Field and left modules over Rational Field

                When extra properties about the morphism are known, one
                can specify the category of which it is a morphism::

                    sage: alpha = lambda i: -p.algebra_generators()[i]
                    sage: f = p.algebra_morphism(alpha, codomain = p, category = GradedHopfAlgebrasWithBasis(QQ))
                    sage: f
                    Generic endomorphism of Symmetric Functions in super space over the Rational Field in the Power basis
                    sage: f(2*p[[],[1]] + 3 * p[[1],[]] + p[[0],[1]] )
                    -2*p[; 1] + p[0; 1] - 3*p[1; ]
                    sage: f.category()
                    Category of endsets of hopf algebras over Rational Field and graded modules over Rational Field

                If ``anti`` is true, this returns an anti-algebra morphism::

                    sage: f = p.algebra_morphism(alpha, codomain = p, anti=True)
                    sage: f
                    Generic endomorphism of Symmetric Functions in super space over the Rational Field in the Power basis
                    sage: f(2*p[[3,1],[]] + 3 * p[[2,1],[1]])
                    3*p[2, 1; 1] - 2*p[3, 1; ]
                    sage: f.category()
                    Category of endsets of modules with basis over Rational Field
                """
                return SFSuperSpaceAlgebraMorphism(self, on_generators, **keywords)

            def product_on_basis(self, spartition1, spartition2):
                r"""
                Product of two multiplicative elements indexed by super partitions.

                Generators from the symmetric component commute while
                generators from the antisymmetric component satisfy
                `a_i a_j = - a_j a_i`

                INPUT:

                - ``spartition1, spartition2`` -- super partitions or lists

                OUTPUT:

                - a monomial in the basis

                EXAMPLES::

                    sage: e = SymmetricFunctionsinSuperSpace(QQ).e()
                    sage: e.product_on_basis([[1,0],[1,1]],[[2],[2]])
                    e[2, 1, 0; 2, 1, 1]
                    sage: e.product_on_basis([[1],[1,1]],[[2],[2]])
                    -e[2, 1; 2, 1, 1]
                    sage: e.product_on_basis([[2],[1,1]],[[2],[2]])
                    0
                """
                pa = spartition1[0]+spartition2[0]
                pb = sorted([(pa[i],i+1) for i in range(len(pa))],reverse=True)
                sign = Permutation([v[1] for v in pb]).sign()
                ps = sorted(spartition1[1]+spartition2[1],reverse=True)
                if len(uniq(pa))==len(pa):
                    return sign*self.monomial(SuperPartition([[v[0] for v in pb], ps]))
                else:
                    return self.zero()

    class Monomial(CombinatorialFreeModule, BindableClass):
        def __init__(self, SFSS):
            r"""
            TESTS::

                sage: m = SymmetricFunctionsinSuperSpace(QQ).Monomial()
                sage: TestSuite(m).run()
            """
            CombinatorialFreeModule.__init__(self, SFSS.base_ring(), SuperPartitions(),
                                             prefix='m', bracket=False, category=SFSS.Bases())
            self._Sym_m = SymmetricFunctions(self.base_ring()).m()
            self._Sym_p = SymmetricFunctions(self.base_ring()).p()

        @cached_method
        def power_to_self_on_basis(self, sp):
            r"""
            Expand a power basis element in the monomial basis.

            INPUT:

            - ``sp`` -- a super partition

            EXAMPLES::

                sage: m = SymmetricFunctionsinSuperSpace(QQ).Monomial()
                sage: m.power_to_self_on_basis(SuperPartition([[2,0],[1,1]]))
                2*m[2, 0; 1, 1] + m[2, 0; 2] + 2*m[2, 1; 1] + 2*m[3, 0; 1] + 2*m[3, 1; ] + m[4, 0; ]
            """
            return self.sum(self.term(SuperPartition([sp[0],[]]))*\
                self.term(SuperPartition([[],list(la)]),c) \
                for (la,c) in self._Sym_m(self._Sym_p(sp[1])))

        @cached_method
        def product_on_basis(self, sp1, sp2):
            r"""
            Product of two monomial basis elements.

            The terms in the product of monomial basis elements are calculated
            by fixing a possible fermionic partition and determining the
            possible ways that two possible monomials in super space have a
            product with the fermionic partition.

            INPUT:

            - ``sp1``, ``sp2`` -- super partitions

            EXAMPLES::

                sage: m = SymmetricFunctionsinSuperSpace(QQ).Monomial()
                sage: m.product_on_basis(SuperPartition([[1],[1]]),SuperPartition([[0],[2]]))
                m[1, 0; 2, 1] + m[1, 0; 3] + m[3, 0; 1] + m[3, 1; ]
                sage: m.product_on_basis(SuperPartition([[1],[]]),SuperPartition([[1],[]]))
                0
                sage: m.product_on_basis(SuperPartition([[0],[]]),SuperPartition([[1],[]]))
                -m[1, 0; ]
            """
            min_a_size = sum(sp1[0]+sp2[0])
            out = self.zero()
            for r in range(min_a_size, min_a_size+sum(sp1[1]+sp2[1])+1):
                pars = Partitions(r, max_slope=-1, length=len(sp1[0]+sp2[0]))
                pars2 = Partitions(r, max_slope=-1, length=len(sp1[0]+sp2[0])-1)
                out += self.sum(m_m_mult_fixed_a(self, sp1, sp2, list(a_sp3)) for a_sp3 in pars)
                out += self.sum(m_m_mult_fixed_a(self, sp1, sp2, list(a_sp3)+[0]) for a_sp3 in pars2)
            return out

    m = Monomial

    class Schur_sb(CombinatorialFreeModule, BindableClass):
        r"""
        The Schur-star-bar basis.

        The Pieri rule for this basis is described in [JL2016]_ .
        The Schur-star-bar basis elements are (up to sign) an application
        of omega applied to the Schur basis.
        The basis is also dual to the Schur-bar basis.

        EXAMPLES::

            sage: ssb = SymmetricFunctionsinSuperSpace(QQ).Schur_sb()
            sage: ssb[[1,0],[]]*ssb[[],[3]]
            ssb[1, 0; 3] + ssb[3, 0; 1] + ssb[4, 0; ]
            sage: ssb[[1],[2]]*ssb[[],[3]]
            ssb[1; 3, 2] + ssb[1; 4, 1] + ssb[1; 5]
            sage: ssb[[1],[2]]*ssb[[3],[]]
            -ssb[3, 1; 2] - ssb[4, 1; 1] - ssb[5, 1; ]
            sage: ssb[[0],[1]].coproduct()
            ssb[; ] # ssb[0; 1] + ssb[0; 1] # ssb[; ]
            sage: ssb[[1],[2]].omega()
            -ssb[0; 1, 1, 1] - 2*ssb[0; 2, 1] - ssb[0; 3]
            sage: sb = SymmetricFunctionsinSuperSpace(QQ).Schur_b()
            sage: sb(ssb[[1],[2]])
            -sb[0; 2, 1] + sb[0; 3] + 2*sb[1; 2] - sb[2; 1] - sb[3; ]
            sage: ssb[[2],[1]].antipode()
            -ssb[0; 1, 1, 1] + 2*ssb[0; 2, 1] + ssb[0; 3] + ssb[1; 2] + ssb[2; 1]
            sage: [ssb[-2,0,2,1].scalar(sb(la)) for la in SuperPartitions(5,2)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
            sage: s = SymmetricFunctionsinSuperSpace(QQ).Schur()
            sage: ssb(s[-2,0,2,1].omega())
            -ssb[3, 0; 2]
        """
        def __init__(self, SFSS):
            r"""
            TESTS::

                sage: ssb = SymmetricFunctionsinSuperSpace(QQ).Schur_sb()
                sage: TestSuite(ssb).run()
            """
            CombinatorialFreeModule.__init__(self, SFSS.base_ring(), \
                SuperPartitions(), prefix='ssb', bracket=False, \
                category=SFSS.Bases())
            self._h = SFSS.Complete()

        @cached_method
        def self_to_complete_on_basis(self, sp):
            r"""
            Convert the Schur basis to complete basis using the Pieri rule.

            INPUT:

            - ``sp`` -- a super-partition

            OUTPUT:

            - an expression in the complete basis

            EXAMPLES::

                sage: ssb = SymmetricFunctionsinSuperSpace(QQ).Schur_sb()
                sage: ssb.self_to_complete_on_basis(SuperPartition([[1],[1]]))
                h[1; 1] - h[2; ]
                sage: ssb.self_to_complete_on_basis(SuperPartition([[1],[1,1]]))
                h[1; 1, 1] - h[1; 2] - h[2; 1] + h[3; ]
            """
            if SuperPartition(sp).bosonic_length()>0:
                sp2 = SuperPartition([sp[0],sp[1][1:]])
                return self.self_to_complete_on_basis(sp2)*self._h([[],[sp[1][0]]])\
                    - sum(self.self_to_complete_on_basis(spa) for spa in \
                          sp2.add_horizontal_border_strip_star_bar(sp[1][0]) if spa!=sp)
            else:
                return self._h.prod(self._h([[a],[]]) for a in sp[0])

        def complete_to_self_by_triangularity(self, h_expr):
            r"""
            Convert an expression in the complete basis to the Schur basis.

            INPUT:

            - ``h_expr`` -- an expression in the complete basis

            OUTPUT:

            - an expression in the Schur-star basis

            EXAMPLES::

                sage: ssb = SymmetricFunctionsinSuperSpace(QQ).Schur_sb()
                sage: h = SymmetricFunctionsinSuperSpace(QQ).Complete()
                sage: ssb.complete_to_self_by_triangularity(-h[[0],[2]]+h[[1],[1]])
                -ssb[0; 2] + ssb[1; 1]
                sage: ssb.self_to_complete_on_basis(SuperPartition([[1],[1,1]]))
                h[1; 1, 1] - h[1; 2] - h[2; 1] + h[3; ]
                sage: ssb(ssb.self_to_complete_on_basis(SuperPartition([[1],[1,1]])))
                ssb[1; 1, 1]
            """
            out = self.zero()
            while h_expr!=0:
                (a, sp, c) = min([[-sum(sp[0]),sp,c] for (sp,c) in h_expr])
                htrm = self.self_to_complete_on_basis(sp)
                h_expr -= c*htrm.coefficient(sp)*htrm
                out += self.term(sp, c*htrm.coefficient(sp))
            return out

    ssb = Schur_sb

    class Schur_s(CombinatorialFreeModule, BindableClass):
        r"""
        The Schur-star basis.

        The basis is calculated by the Pieri rule in order to obtain the
        the expansion in the complete basis.  It is calculated from the
        complete basis by identifying a triangularity relation.

        The Pieri rule for this basis is described in [JL2016]_ .
        The Schur-star basis elements are (up to sign) an element
        of the Schur-bar basis.  The Schur-star basis is dual to the
        Schur basis.

        EXAMPLES::

            sage: ss = SymmetricFunctionsinSuperSpace(QQ).Schur_s()
            sage: ss[[1,0],[]]*ss[[],[3]]
            ss[1, 0; 3]
            sage: ss[[1],[2]]*ss[[],[3]]
            ss[0; 4, 2] + ss[1; 3, 2] + ss[1; 4, 1] + ss[1; 5]
            sage: ss[[1],[2]]*ss[[3],[]]
            ss[1, 0; 5] - ss[3, 1; 2]
            sage: ss[[0],[1]].coproduct()
            ss[; ] # ss[0; 1] + ss[; 1] # ss[0; ] + ss[0; ] # ss[; 1] + ss[0; 1] # ss[; ]
            sage: ss[[1],[2]].omega()
            ss[0; 1, 1, 1] + ss[0; 2, 1] - ss[1; 1, 1]
            sage: sb = SymmetricFunctionsinSuperSpace(QQ).Schur_b()
            sage: sb(ss[[1],[2]].omega())
            sb[1; 2]
            sage: ss[[2],[1]].antipode()
            2*ss[0; 1, 1, 1] + 2*ss[0; 2, 1] - 2*ss[1; 1, 1] - 2*ss[1; 2] + ss[2; 1]
            sage: s = SymmetricFunctionsinSuperSpace(QQ).Schur()
            sage: [ss[-2,0,2,1].scalar(s(la)) for la in SuperPartitions(5,2)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
        """
        def __init__(self, SFSS):
            r"""
            TESTS::

                sage: ss = SymmetricFunctionsinSuperSpace(QQ).Schur_s()
                sage: TestSuite(ss).run()
            """
            CombinatorialFreeModule.__init__(self, SFSS.base_ring(), \
                SuperPartitions(), prefix='ss', bracket=False, \
                category=SFSS.Bases())
            self._h = SFSS.Complete()

        @cached_method
        def self_to_complete_on_basis(self, sp):
            r"""
            Convert the Schur-star basis to complete basis using the Pieri rule.

            INPUT:

            - ``sp`` -- a super-partition

            OUTPUT:

            - an expression in the complete basis

            EXAMPLES::

                sage: ss = SymmetricFunctionsinSuperSpace(QQ).Schur_s()
                sage: ss.self_to_complete_on_basis(SuperPartition([[1],[1]]))
                -h[0; 2] + h[1; 1]
                sage: ss.self_to_complete_on_basis(SuperPartition([[1],[1,1]]))
                -h[0; 2, 1] + h[0; 3] + h[1; 1, 1] - h[1; 2]
            """
            if SuperPartition(sp).bosonic_length()>0:
                sp2 = SuperPartition([sp[0],sp[1][1:]])
                return self.self_to_complete_on_basis(sp2)*self._h([[],[sp[1][0]]])\
                    - sum(self.self_to_complete_on_basis(spa) for spa in \
                          sp2.add_horizontal_border_strip_star(sp[1][0]) if spa!=sp)
            else:
                return self._h.prod(self._h([[a],[]]) for a in sp[0])

        def complete_to_self_by_triangularity(self, h_expr):
            r"""
            Convert an expression in the complete basis to the Schur-star basis.

            INPUT:

            - ``h_expr`` -- an expression in the complete basis

            OUTPUT:

            - an expression in the Schur-star basis

            EXAMPLES::

                sage: ss = SymmetricFunctionsinSuperSpace(QQ).Schur_s()
                sage: h = SymmetricFunctionsinSuperSpace(QQ).Complete()
                sage: ss.complete_to_self_by_triangularity(-h[[0],[2]]+h[[1],[1]])
                ss[1; 1]
                sage: ss.self_to_complete_on_basis(SuperPartition([[1],[1,1]]))
                -h[0; 2, 1] + h[0; 3] + h[1; 1, 1] - h[1; 2]
                sage: ss(ss.self_to_complete_on_basis(SuperPartition([[1],[1,1]])))
                ss[1; 1, 1]
            """
            out = self.zero()
            while h_expr!=0:
                (a, sp, c) = min([[-sum(sp[0]),sp,c] for (sp,c) in h_expr])
                htrm = self.self_to_complete_on_basis(sp)
                h_expr -= c*htrm.coefficient(sp)*htrm
                out += self.term(sp, c*htrm.coefficient(sp))
            return out

    ss = Schur_s

    class Schur_b(CombinatorialFreeModule, BindableClass):
        r"""
        The Schur-bar basis as the `q=t=\infty` specialization of Macdonald.

        The basis is calculated from the complete expansion of the Schur-star
        basis and then the same coefficients are used for the elementary
        basis for the Schur-bar.  It is calculated from the
        elementary basis by identifying a triangularity relation.

        This basis is dual to the Schur-star-bar basis.
        It is also (up to sign) the involution omega applied to the
        Schur-star basis.

        The Pieri rule for this basis is proven in [JL2016]_ and it was
        conjectured in [BFM2015]_ .

        EXAMPLES::

            sage: sb = SymmetricFunctionsinSuperSpace(QQ).Schur_b()
            sage: sb[[1,0],[]]*sb[[],[3]]
            sb[4, 0; ]
            sage: ss = SymmetricFunctionsinSuperSpace(QQ).Schur_s()
            sage: sb(ss[-2,0,2,1].omega())
            -sb[3, 0; 2]
            sage: ssb = SymmetricFunctionsinSuperSpace(QQ).Schur_sb()
            sage: [sb[-2,0,2,1].scalar(ssb(la)) for la in SuperPartitions(5,2)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
        """
        def __init__(self, SFSS):
            r"""
            TESTS::

                sage: sb = SymmetricFunctionsinSuperSpace(QQ).Schur_b()
                sage: TestSuite(sb).run()
            """
            CombinatorialFreeModule.__init__(self, SFSS.base_ring(), \
                SuperPartitions(), prefix='sb', bracket=False, \
                category=SFSS.Bases())
            self._e = SFSS.Elementary()
            self._ss = SFSS.Schur_s()

        @cached_method
        def self_to_elementary_on_basis(self, sp):
            r"""
            Convert the Schur-bar basis to elementary basis using omega.

            INPUT:

            - ``sp`` -- a super-partition

            OUTPUT:

            - an expression in the elementary basis

            EXAMPLES::

                sage: sb = SymmetricFunctionsinSuperSpace(QQ).Schur_b()
                sage: sb.self_to_elementary_on_basis(SuperPartition([[1],[1]]))
                e[0; 2]
                sage: sb.self_to_elementary_on_basis(SuperPartition([[1],[1,1]]))
                e[0; 3]
            """
            mc2 = binomial(sp.fermionic_degree(),2)
            return (-1)**mc2*self._e.sum(self._e.term(sp,c) for (sp,c) in \
                self._ss.self_to_complete_on_basis(sp.conjugate()))

        def elementary_to_self(self, e_expr):
            r"""
            Convert an expression in the elementary basis to the Schur-bar
            basis.

            The conversion applies ``omega`` to the expression, converts to the
            Schur-star basis and then applies the relation that ``omega``
            applied to a Schur-star is a Schur-bar up to sign.

            INPUT:

            - ``e_expr`` -- an expression in the elementary basis

            OUTPUT:

            - an expression in the Schur-bar basis

            EXAMPLES::

                sage: sb = SymmetricFunctionsinSuperSpace(QQ).Schur_b()
                sage: e = SymmetricFunctionsinSuperSpace(QQ).Elementary()
                sage: sb.elementary_to_self(-e[[0],[2]]+e[[1],[1]])
                sb[0; 2]
                sage: sb.self_to_elementary_on_basis(SuperPartition([[1],[1,1]]))
                e[0; 3]
            """
            ss_e_expr = self._ss(e_expr.omega())
            return self.sum_of_terms((sp.conjugate(), \
                (-1)**binomial(sp.fermionic_sector(),2)*c) \
                for (sp,c) in ss_e_expr)

    sb = Schur_b

    class Schur(CombinatorialFreeModule, BindableClass):
        r"""
        The Schur basis as the `q=t=0` specialization of Macdonald.

        The basis is calculated from the complete expansion of the Schur-star-bar
        basis and then the same coefficients are used for the elementary
        basis for the Schur.  It is calculated from the
        elementary basis by identifying a triangularity relation.

        This basis is the one-true Schur basis (except that there are 4).
        It is dual to the Schur-star basis.
        It is also the involution omega (up to sign) applied to the
        Schur-star-bar basis.

        The Pieri rule for this basis is proven in [JL2016]_ and conjectured
        in [BFM2015]_ .

        EXAMPLES::

            sage: s = SymmetricFunctionsinSuperSpace(QQ).Schur()
            sage: s[[1,0],[]]*s[[],[3]]
            s[1, 0; 3] + s[3, 0; 1] + s[4, 0; ]
            sage: ssb = SymmetricFunctionsinSuperSpace(QQ).Schur_sb()
            sage: s(ssb[-2,0,2,1].omega())
            -s[3, 0; 2]
            sage: ss = SymmetricFunctionsinSuperSpace(QQ).Schur_s()
            sage: [s[-2,0,2,1].scalar(ss(la)) for la in SuperPartitions(5,2)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]

        This basis appears in [BFM2015]_ . The following computations
        are examples 3.3, 3.4, 4.1, 4.2 and 4.3 from that reference. ::

            sage: m = SymmetricFunctionsinSuperSpace(QQ).m()
            sage: m(s[-3,0,4,1])
            4*m[3, 0; 1, 1, 1, 1, 1] + 3*m[3, 0; 2, 1, 1, 1] + 2*m[3, 0; 2, 2, 1]
             + 2*m[3, 0; 3, 1, 1] + m[3, 0; 3, 2] + m[3, 0; 4, 1]
             + 3*m[3, 1; 1, 1, 1, 1] + 2*m[3, 1; 2, 1, 1] + m[3, 1; 2, 2]
             + m[3, 1; 3, 1] + 2*m[3, 2; 1, 1, 1] + m[3, 2; 2, 1]
            sage: m(s[-1,4])
            m[1; 1, 1, 1, 1] + m[1; 2, 1, 1] + m[1; 2, 2] + m[1; 3, 1] + m[1; 4]
             + m[2; 1, 1, 1] + m[2; 2, 1] + m[2; 3] + m[3; 1, 1] + m[3; 2]
            sage: s[-4,0,3]*s[-3] # long time (14s)
            -s[4, 3, 0; 3] + s[5, 4, 0; 1] + s[5, 4, 1; ] + s[6, 4, 0; ]
            sage: s[-1,2,1]*s[0,1,1,1]
            s[1, 0; 2, 1, 1, 1, 1] + s[1, 0; 2, 2, 1, 1] + s[1, 0; 3, 1, 1, 1]
             + s[1, 0; 3, 2, 1] + s[2, 0; 3, 1, 1] + s[2, 0; 3, 2]
            sage: s[-2,0,1]*s[3]
            s[2, 0; 3, 1] + s[2, 0; 4] + s[3, 0; 2, 1] + s[4, 0; 1, 1]
             + s[4, 0; 2] + s[5, 0; 1]
            sage: s[-2,0,1]*s[-3]
            s[3, 2, 0; 1] + s[4, 2, 0; ]
            sage: s[-2,0,1]*s[1,1]
            s[2, 0; 1, 1, 1] + s[2, 0; 2, 1] + s[2, 1; 2] + s[3, 0; 1, 1]
             + s[3, 0; 2]
            sage: s[-2,0,1]*s[0,1,1]
            s[2, 1, 0; 2]
        """
        def __init__(self, SFSS):
            r"""
            TESTS::

                sage: s = SymmetricFunctionsinSuperSpace(QQ).Schur()
                sage: TestSuite(s).run()
            """
            CombinatorialFreeModule.__init__(self, SFSS.base_ring(), \
                SuperPartitions(), prefix='s', bracket=False, \
                category=SFSS.Bases())
            self._e = SFSS.Elementary()
            self._ssb = SFSS.Schur_sb()

        @cached_method
        def self_to_elementary_on_basis(self, sp):
            r"""
            Convert the Schur basis to elementary basis using omega.

            INPUT:

            - ``sp`` -- a super-partition

            OUTPUT:

            - an expression in the elementary basis

            EXAMPLES::

                sage: s = SymmetricFunctionsinSuperSpace(QQ).Schur()
                sage: s.self_to_elementary_on_basis(SuperPartition([[1],[1]]))
                e[0; 2] - e[2; ]
                sage: s.self_to_elementary_on_basis(SuperPartition([[1],[1,1]]))
                e[0; 3] - e[3; ]
            """
            mc2 = binomial(sp.fermionic_degree(),2)
            return (-1)**mc2*self._e.sum(self._e.term(sp,c) for (sp,c) in \
                self._ssb.self_to_complete_on_basis(sp.conjugate()))

        def elementary_to_self(self, e_expr):
            r"""
            Convert an expression in the elementary basis to the Schur basis

            The conversion applies ``omega`` to the expression, converts to the
            Schur-star-bar basis and then applies the relation that ``omega``
            applied to a Schur-star-bar is a Schur up to sign.

            INPUT:

            - ``e_expr`` -- an expression in the elementary basis

            OUTPUT:

            - an expression in the Schur basis

            EXAMPLES::

                sage: s = SymmetricFunctionsinSuperSpace(QQ).Schur()
                sage: e = SymmetricFunctionsinSuperSpace(QQ).Elementary()
                sage: s.elementary_to_self(-e[[0],[2]]+e[[1],[1]])
                s[0; 2] - s[1; 1]
                sage: s.elementary_to_self(e[0,2]-e[-2])
                s[1; 1]
                sage: s.self_to_elementary_on_basis(SuperPartition([[1],[1,1]]))
                e[0; 3] - e[3; ]
            """
            ssb_e_expr = self._ssb(e_expr.omega())
            return self.sum_of_terms((sp.conjugate(), \
                (-1)**binomial(sp.fermionic_sector(),2)*c) \
                for (sp,c) in ssb_e_expr)

    s = Schur

    class Elementary(CombinatorialFreeModule, BindableClass):
        def __init__(self, SFSS):
            r"""
            TESTS::

                sage: e = SymmetricFunctionsinSuperSpace(QQ).Elementary()
                sage: TestSuite(e).run()
            """
            CombinatorialFreeModule.__init__(self, SFSS.base_ring(), \
                SuperPartitions(), prefix='e', bracket=False, \
                category=SFSS.MultiplicativeBases())

        def p_to_self_on_generator(self, n):
            r"""
            A Power generator in the Complete basis.

            INPUT:

            - ``n`` -- an integer

            The generators are ordered so that the non-positive integers
            are mapped to a fermionic generator and the positive
            integers are mapped to the bosonic generator.

            EXAMPLES::

                sage: e = SymmetricFunctionsinSuperSpace(QQ).e()
                sage: e.p_to_self_on_generator(2)
                e[; 1, 1] - 2*e[; 2]
                sage: e.p_to_self_on_generator(-2)
                e[0; 1, 1] - e[0; 2] - e[1; 1] + e[2; ]
            """
            if n<=0:
                return self.agen_p_to_self(abs(n))
            else:
                return self.sgen_p_to_self(n)

        @cached_method
        def sgen_p_to_self(self, n):
            r"""
            A symmetric Power generator in the Elementary basis.

            INPUT:

            - ``n`` -- a positive integer

            EXAMPLES::

                sage: e = SymmetricFunctionsinSuperSpace(QQ).e()
                sage: e.sgen_p_to_self(4)
                e[; 1, 1, 1, 1] - 4*e[; 2, 1, 1] + 2*e[; 2, 2] + 4*e[; 3, 1] - 4*e[; 4]
            """
            n = Integer(n)
            if n==0:
                return self.term(SuperPartition([[],[]]))
            return self.term(SuperPartition([[],[n]]),(-1)**(n+1)*n) - \
                self.sum(self.sgen_p_to_self(r)*self.term(SuperPartition([[],[n-r]]),\
                    (-1)**(n+r)) for r in range(1,n))

        @cached_method
        def agen_p_to_self(self, n):
            r"""
            A fermionic Power generator in the Elementary basis.

            INPUT:

            - ``n`` -- a non-negative integer

            This is equation (3.44) of [DLM2006]_.

            EXAMPLES::

                sage: e = SymmetricFunctionsinSuperSpace(QQ).e()
                sage: e.agen_p_to_self(3)
                e[0; 1, 1, 1] - 2*e[0; 2, 1] + e[0; 3] - e[1; 1, 1] + e[1; 2] + e[2; 1] - e[3; ]
            """
            n = Integer(n)
            if n==0:
                return self.term(SuperPartition([[0],[]]))
            return (-1)**n*self.term(SuperPartition([[n],[]])) + \
               self.sum(self.sgen_p_to_self(r)*self.term(\
                   SuperPartition([[n-r],[]]),(-1)**(r+n)/(n+1)) \
                   for r in range(1,n+1)) - \
               self.sum(self.agen_p_to_self(r)*self.term(\
                   SuperPartition([[],[n-r]]),(-1)**(n+r)*(r+1)/(n+1)) \
                   for r in range(n))

        def h_to_self_on_generator(self, n):
            r"""
            A Complete generator in the Elementary basis.

            INPUT:

            - ``n`` -- an integer

            The generators are ordered so that the non-positive integers
            are mapped to a fermionic generator and the positive
            integers are mapped to the bosonic generator.

            EXAMPLES::

                sage: e = SymmetricFunctionsinSuperSpace(QQ).e()
                sage: e.h_to_self_on_generator(2)
                e[; 1, 1] - e[; 2]
                sage: e.h_to_self_on_generator(-2)
                3*e[0; 1, 1] - 2*e[0; 2] - 2*e[1; 1] + e[2; ]
            """
            if n<=0:
                return self.agen_h_to_self(abs(n))
            else:
                return self.sgen_h_to_self(n)

        @cached_method
        def sgen_h_to_self(self, n):
            r"""
            A symmetric Complete generator expanded in the Elementary basis.

            INPUT:

            - ``n`` -- a positive integer

            EXAMPLES::

                sage: e = SymmetricFunctionsinSuperSpace(QQ).e()
                sage: e.sgen_h_to_self(4)
                e[; 1, 1, 1, 1] - 3*e[; 2, 1, 1] + e[; 2, 2] + 2*e[; 3, 1] - e[; 4]
            """
            if n==0:
                return self.term(SuperPartition([[],[]]))
            return self.sum(self.term(SuperPartition([[],la.to_list()]), \
                (-1)**(n-len(la))*factorial(len(la))/prod([factorial(a) \
                    for a in la.to_exp()])) for la in Partitions(n))

        @cached_method
        def agen_h_to_self(self, n):
            r"""
            A fermionic Complete generator expanded in the Elementary basis.

            INPUT:

            - ``n`` -- a non-negative integer

            This is equation (3.20) of [DLM2006]_.

            EXAMPLES::

                sage: e = SymmetricFunctionsinSuperSpace(QQ).e()
                sage: e.agen_h_to_self(2)
                3*e[0; 1, 1] - 2*e[0; 2] - 2*e[1; 1] + e[2; ]
            """
            if n==0:
                return self.term(SuperPartition([[0],[]]))
            return -self.sum((-1)**r*self.term(SuperPartition([[],[r]]))*\
                self.agen_h_to_self(n-r) for r in range(1,n+1)) +\
                self.sum((-1)**r*self.term(SuperPartition([[r],[]]))*\
                self.sgen_h_to_self(n-r) for r in range(n+1))

    e = Elementary

    class Complete(CombinatorialFreeModule, BindableClass):
        def __init__(self, SFSS):
            r"""
            TESTS::

                sage: h = SymmetricFunctionsinSuperSpace(QQ).Complete()
                sage: TestSuite(h).run()
            """
            CombinatorialFreeModule.__init__(self, \
                SFSS.base_ring(), SuperPartitions(), prefix='h', \
                bracket=False, category=SFSS.MultiplicativeBases())

        def p_to_self_on_generator(self, n):
            r"""
            A Power generator in the Complete basis.

            INPUT:

            - ``n`` -- an integer

            The generators are ordered so that the non-positive integers
            are mapped to a fermionic generator and the positive
            integers are mapped to the bosonic generator.

            EXAMPLES::

                sage: h = SymmetricFunctionsinSuperSpace(QQ).h()
                sage: h.p_to_self_on_generator(2)
                -h[; 1, 1] + 2*h[; 2]
                sage: h.p_to_self_on_generator(-2)
                h[0; 1, 1] - h[0; 2] - h[1; 1] + h[2; ]
            """
            if n<=0:
                return self.agen_p_to_self(abs(n))
            else:
                return self.sgen_p_to_self(n)

        @cached_method
        def sgen_p_to_self(self, n):
            r"""
            A symmetric Power generator in the Complete basis.

            INPUT:

            - ``n`` -- a positive integer

            EXAMPLES::

                sage: h = SymmetricFunctionsinSuperSpace(QQ).h()
                sage: h.sgen_p_to_self(4)
                -h[; 1, 1, 1, 1] + 4*h[; 2, 1, 1] - 2*h[; 2, 2] - 4*h[; 3, 1] + 4*h[; 4]
            """
            n = Integer(n)
            if n==0:
                return self.term(SuperPartition([[],[]]))
            return self.term(SuperPartition([[],[n]]),n) - \
                self.sum(self.sgen_p_to_self(r)*self.term(SuperPartition([[],[n-r]])) \
                    for r in range(1,n))

        @cached_method
        def agen_p_to_self(self, n):
            r"""
            A fermionic Power generator in the Complete basis.

            INPUT:

            - ``n`` -- a non-negative integer

            This is equation (3.43) of [DLM2006]_.

            EXAMPLES::

                sage: h = SymmetricFunctionsinSuperSpace(QQ).h()
                sage: h.agen_p_to_self(3)
                -h[0; 1, 1, 1] + 2*h[0; 2, 1] - h[0; 3] + h[1; 1, 1] - h[1; 2] - h[2; 1] + h[3; ]
            """
            n = Integer(n)
            if n==0:
                return self.term(SuperPartition([[0],[]]))
            return self.term(SuperPartition([[n],[]])) - \
               self.sum(self.sgen_p_to_self(r)*self.term(SuperPartition([[n-r],[]]),1/(n+1)) \
                        for r in range(1,n+1)) - \
               self.sum(self.agen_p_to_self(r)*self.term(SuperPartition([[],[n-r]]),(r+1)/(n+1)) \
                        for r in range(n))

        def e_to_self_on_generator(self, n):
            r"""
            An Elementary generator in the Complete basis.

            INPUT:

            - ``n`` -- an integer

            The generators are ordered so that the non-positive integers
            are mapped to a fermionic generator and the positive
            integers are mapped to the bosonic generator.

            EXAMPLES::

                sage: h = SymmetricFunctionsinSuperSpace(QQ).h()
                sage: h.e_to_self_on_generator(2)
                h[; 1, 1] - h[; 2]
                sage: h.e_to_self_on_generator(-2)
                3*h[0; 1, 1] - 2*h[0; 2] - 2*h[1; 1] + h[2; ]
            """
            if n<=0:
                return self.agen_e_to_self(abs(n))
            else:
                return self.sgen_e_to_self(n)

        @cached_method
        def sgen_e_to_self(self, n):
            r"""
            A symmetric Elementary generator expanded in the Complete basis.

            INPUT:

            - ``n`` -- a positive integer

            EXAMPLES::

                sage: h = SymmetricFunctionsinSuperSpace(QQ).h()
                sage: h.sgen_e_to_self(4)
                h[; 1, 1, 1, 1] - 3*h[; 2, 1, 1] + h[; 2, 2] + 2*h[; 3, 1] - h[; 4]
            """
            if n==0:
                return self.term(SuperPartition([[],[]]))
            return self.sum(self.term(SuperPartition([[],la.to_list()]), \
                (-1)**(n-len(la))*factorial(len(la))/prod([factorial(a) for a in \
                    la.to_exp()])) for la in Partitions(n))

        @cached_method
        def agen_e_to_self(self, n):
            r"""
            A fermionic Elementary generator expanded in the Complete basis.

            INPUT:

            - ``n`` -- a non-negative integer

            This is equation (3.20) of [DLM2006]_.

            EXAMPLES::

                sage: h = SymmetricFunctionsinSuperSpace(QQ).h()
                sage: h.agen_e_to_self(2)
                3*h[0; 1, 1] - 2*h[0; 2] - 2*h[1; 1] + h[2; ]
            """
            n = Integer(n)
            if n==0:
                return self.term(SuperPartition([[0],[]]))
            return self.sum((-1)**(n+r)*self.sgen_e_to_self(r)*\
                self.term(SuperPartition([[n-r],[]])) for r in range(n+1)) - \
                self.sum((-1)**(n+r)*self.agen_e_to_self(r)*\
                self.term(SuperPartition([[],[n-r]])) for r in range(n))

    h = Complete
    Homogeneous = Complete

    class Power(CombinatorialFreeModule, BindableClass):
        def __init__(self, SFSS):
            r"""
            TESTS::

                sage: p = SymmetricFunctionsinSuperSpace(QQ).Power()
                sage: TestSuite(p).run()
            """
            CombinatorialFreeModule.__init__(self, SFSS.base_ring(), \
                SuperPartitions(), prefix='p', bracket=False, \
                category=SFSS.MultiplicativeBases())

        def monomial_to_self_by_triangularity(self, m_expr):
            r"""
            Change a monomial expression to the power sum basis.

            The change of basis between the monomial and power sum basis
            elements have coefficients.

            INPUT:

            - ``m_expr`` -- an expression in the monomial basis

            EXAMPLES::

                sage: m = SymmetricFunctionsinSuperSpace(QQ).Monomial()
                sage: p = SymmetricFunctionsinSuperSpace(QQ).Power()
                sage: p.monomial_to_self_by_triangularity(m[[2],[1,1]])
                1/2*p[2; 1, 1] - 1/2*p[2; 2] - p[3; 1] + p[4; ]
            """
            out = self.zero()
            SFSS_m = m_expr.parent() # monomial basis
            while m_expr!=self.zero():
                (sp, c) = m_expr.trailing_item()
                ft = prod(factorial(p) for p in Partition(sp[1]).to_exp())
                out += self.term(sp, c/ft)
                m_expr -= c/ft*SFSS_m.power_to_self_on_basis(sp)
            return out

        @cached_method
        def h_to_self_on_generator(self, n):
            r"""
            A Complete generator expanded in the Power basis.

            INPUT:

            - ``n`` -- an integer

            The generators are ordered so that the non-positive integers
            are mapped to a fermionic generator and the positive
            integers are mapped to the bosonic generator.

            This is equation (3.61) and equation (3.62) of [DLM2006]_.

            EXAMPLES::

                sage: p = SymmetricFunctionsinSuperSpace(QQ).p()
                sage: p.h_to_self_on_generator(-2)
                1/2*p[0; 1, 1] + 1/2*p[0; 2] + p[1; 1] + p[2; ]
                sage: p.h_to_self_on_generator(3)
                1/6*p[; 1, 1, 1] + 1/2*p[; 2, 1] + 1/3*p[; 3]
            """
            return self.sum(self.term(sp,1/sp.zee()) for sp \
                in SuperPartitions(abs(n),Integer(n<=0)))

        @cached_method
        def e_to_self_on_generator(self, n):
            r"""
            An Elementary generator expanded in the Power basis.

            INPUT:

            - ``n`` -- an integer

            The generators are ordered so that the non-positive integers
            are mapped to a fermionic generator and the positive
            integers are mapped to the bosonic generator.

            This is equation (3.61) and equation (3.62) of [DLM2006]_.

            EXAMPLES::

                sage: p = SymmetricFunctionsinSuperSpace(QQ).p()
                sage: p.e_to_self_on_generator(-2)
                1/2*p[0; 1, 1] - 1/2*p[0; 2] - p[1; 1] + p[2; ]
                sage: p.e_to_self_on_generator(3)
                1/6*p[; 1, 1, 1] - 1/2*p[; 2, 1] + 1/3*p[; 3]
            """
            return self.sum(self.term(sp,sp.sign()/sp.zee()) for sp \
                in SuperPartitions(abs(n),Integer(n<=0)))

        def coproduct_on_generators(self, n):
            r"""
            The coproduct on the power sum generators.

            INPUT:

            - ``n`` -- an integer

            OUTPUT:

            - an element of the tensor square of the power basis

            EXAMPLES::

                sage: p = SymmetricFunctionsinSuperSpace(QQ).p()
                sage: p.coproduct_on_generators(-2)
                p[; ] # p[2; ] + p[2; ] # p[; ]
                sage: p.coproduct_on_generators(3)
                p[; ] # p[; 3] + p[; 3] # p[; ]
            """
            from sage.categories.all import tensor
            x = self.algebra_generators()[n]
            return tensor([self.one(),x])+tensor([x,self.one()])

        @lazy_attribute
        def coproduct(self):
            r"""
            Return the coproduct morphism in the basis ``self``.

            OUTPUT:

            - an element of the tensor square of the power basis

            EXAMPLES::

                sage: p = SymmetricFunctionsinSuperSpace(QQ).p()
                sage: p([[2],[1]]).coproduct()
                p[; ] # p[2; 1] + p[; 1] # p[2; ] + p[2; ] # p[; 1] + p[2; 1] # p[; ]
            """
            return self.algebra_morphism(self.coproduct_on_generators, \
                codomain = self.tensor_square())

        def antipode_on_basis(self, sp):
            r"""
            The antipode on a power basis element indexed by ``sp``.

            INPUT:

            - ``sp`` -- a super-partition

            OUPUT:

            - an element of the power sum basis

            EXAMPLES::

                sage: p = SymmetricFunctionsinSuperSpace(QQ).p()
                sage: p.antipode_on_basis(SuperPartition([[2,1],[2]]))
                p[2, 1; 2]
                sage: p.antipode_on_basis(SuperPartition([[2],[1]]))
                p[2; 1]
            """
            return self.term(sp,(-1)**(sp.length()+binomial(sp.fermionic_degree(),2)))

        class Element(CombinatorialFreeModule.Element):

            def scalar(self, other):
                r"""
                The scalar product between ``self`` and ``other``.

                The power sum elements are orthonormal.

                INPUT:

                - ``other`` -- an element of symmetric functions in super space

                EXAMAPLES::

                    sage: p = SymmetricFunctionsinSuperSpace(QQ).p()
                    sage: p[[],[3,1]].scalar(p[[1],[2,1]])
                    0
                    sage: p[[1],[2,1]].scalar(p[[1],[2,1]])
                    2
                """
                p = self.parent()
                py = p(other)
                return sum(sp.zee()*c*py.coefficient(sp) for (sp,c) in self)

            def omega(self):
                r"""
                The involution omega.

                The involution omega is computed by coercion to the power
                sum basis.

                EXAMPLES::

                    sage: p = SymmetricFunctionsinSuperSpace(QQ).p()
                    sage: p[[1,0],[2,2,1]].omega()
                    -p[1, 0; 2, 2, 1]
                    sage: all(p[[n],[]].omega()==(-1)**n*p[[n],[]] for n in range(5))
                    True
                    sage: all(p[[],[n]].omega()==(-1)**(n+1)*p[[],[n]] for n in range(1,6))
                    True
                """
                p = self.parent()
                return p.sum(sp.sign()*c*p(sp) for (sp, c) in self)

            def skew_by(self, other, side='right'):
                r"""
                Operation that is dual to multiplication in the scalar product.

                INPUT:

                - ``other`` -- a supersymmetric function

                - ``side`` -- either 'right' or 'left' (default: 'right')

                OUTPUT:

                - a supersymmetric function

                EXAMPLES::

                    sage: p = SymmetricFunctionsinSuperSpace(QQ).p()
                    sage: p[-2,0,3,2].skew_by(p[-2])
                    -p[0; 3, 2]
                    sage: p[-2,0,3,2].skew_by(p[-2],'left')
                    p[0; 3, 2]
                    sage: p[-2,0,3,2].skew_by(p[0])
                    p[2; 3, 2]
                    sage: p[-2,0,3,2].skew_by(p[3])
                    3*p[2, 0; 2]
                """
                p = self.parent()
                if side=='right':
                    f = lambda sp,bd: p(sp)*other.homogeneous_bi_degree_component(bd)
                else:
                    f = lambda sp,bd: p(other).homogeneous_bi_degree_component(bd)*p(sp)
                return sum(self.scalar(f(sp,bd2))*p(sp)/sp.zee() \
                    for bd1 in self.all_bi_degrees() \
                    for bd2 in other.all_bi_degrees() \
                    for sp in SuperPartitions(bd1[0]-bd2[0], bd1[1]-bd2[1]))

    p = Power
