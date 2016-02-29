r"""
Kac Polynomials and DT-Invariants
"""
#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>,
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
#*****************************************************************************

from sage.misc.all import prod
from sage.rings.all import ZZ
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.partition_tuple import PartitionTuple
from sage.combinat.partition import Partitions
from sage.arith.misc import moebius
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
import itertools


def cohomology_polynomial(Q, mu, var=None):
    r"""
    Return the cohomology polynomial of the quiver `Q`.

    INPUT:

    - `Q` -- a quiver as a digraph
    - ``mu`` -- a partition tuple
    - `var` -- (optional) the variable `q`

    The cohomology polynomial is given by:

    .. MATH::

        \sum_i \dim\left( H_c^{2i}(\mathcal{Q}_{\tilde{v}} ; \CC)_{\epsilon
        \chi^{\mu}} \right) q^{i - d_{\tilde{v}}} = \mathbb{H}_{\mu}^s(q)

    where `\widetilde{v}` is the extended dimension vector given
    in [HLRV2012]_, `\epsilon` is the sign character, and

    .. MATH::

        \mathbb{H}_{\mu}^s(q) := \bigl\langle \mathbb{H}(x_1, \ldots, x_r; q),
        s_{\mu} \bigr\rangle,

    where `\mathbb{H}` is given by Equation (1.4) in [HLRV2012]_.

    EXAMPLES::

        sage: from sage.combinat.kac_polynomial import cohomology_polynomial
        sage: cohomology_polynomial(DiGraph([[1],[]]), PartitionTuple([[1]]))
        1?

    REFERENCES:

    .. [HLRV2012] Tamas Hausel, Emmanuel Letellier, and Fernando
       Rodriguez-Villegas. *Positivity of Kac polynomials and DT-invariants
       for quivers*. (2012). :arxiv:`1204.2375v1`.
    """
    R = PolynomialRing(ZZ, 'q')
    q = R.gen()
    HLP = SymmetricFunctions(R).hall_littlewood(q).P()
    n = ZZ(mu.size())

    ret = R.zero()
    for d in n.divisors():  # We must have d dividing n
        coeff = {}
        for pt in itertools.product(*[Partitions(p.size() // d)
                                      for p in mu.components()]):
            c = prod(q ** pt[i].pairing(pt[j])
                     for i, j in Q.edges(False))
            c /= prod(q ** la.pairing(la) *
                      prod(1 - q ** (-j)
                           for mk in la.to_exp_dict().values()
                           for j in range(1, mk + 1))
                      for la in pt)
            coeff[pt] = c
        for D in d.divisors():
            # here the product *should* live in the tensor product of Sym ?
            power = sum(coeff[p](q=q ** D) *
                        prod(HLP[p].frobenius(D) for pi in p)
                        for p in coeff)
            ret += moebius(D) * power ** (d / D).coefficient(mu) / D

    ret *= (q - 1)

    if var is None:
        return ret
    else:
        return ret(q=var)


def kac_polynomial(Q, v, q=None):
    r"""
    Return the Kac polynomial for the quiver `Q` and the dimension vector `v`.

    INPUT:

    - `Q` -- a quiver as a digraph
    - `v` -- the dimension vector
    - `q` -- (optional) the variable `q`

    The vertices of the quiver must be the integers `0,\dots,n`.

    EXAMPLES::

        sage: from sage.combinat.kac_polynomial import kac_polynomial
        sage: kac_polynomial(DiGraph([[1],[]]), (1,))
        1?
        sage: kac_polynomial(DiGraph({0:[1,2,3]}), (2,1,1,1))
        1+q?
    """
    return cohomology_polynomial(Q, PartitionTuple([[x] for x in v]), q)


def DT_invariant(Q, v, q=None):
    r"""
    Return the DT-invariant for the quiver `Q` and the dimension vector `v`.

    DT is a shorthand for Donaldson-Thomas.

    The vertices of the quiver must be the integers `0,\dots,n`.

    INPUT:

    - `Q` -- a quiver as a digraph
    - `v` -- the dimension vector
    - `q` -- (optional) the variable `q`

    EXAMPLES::

        sage: from sage.combinat.kac_polynomial import DT_invariant
        sage: DT_invariant(DiGraph([[1],[]]), (1,))
        1?
        sage: kac_polynomial(DiGraph({0:[1,2,3]}), (2,1,1,1))
        1+q?
    """
    return cohomology_polynomial(Q, PartitionTuple([[1] * x for x in v]), q)
