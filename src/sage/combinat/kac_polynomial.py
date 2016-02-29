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
from sage.rings.arith import moebius
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
import itertools


def cohomology_polynomial(Q, mu, var=None):
    r"""
    Return the cohomology polynomial of the quiver `Q` .

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
    r = mu.level()
    fac_r = (-1) ** (r - 1) / r

    ret = R.zero()
    for d in n.divisors():  # We must have d dividing n
        qd = q ** d
        terms = R.zero()
        # We only have a chance if each partition's size times d is equal to
        # the size of the corresponding partition in mu
        for pt in itertools.product(*[Partitions(p.size() // d)
                                      for p in mu.components()]):
            Z = zip(mu.components(), pt)
            t = prod(HLP[p].frobenius(d).coefficient(m) for m, p in Z)
            if t != 0:  # If there's something to do
                t *= prod(qd ** pt[i].pairing(pt[j])
                          for i, j in Q.edges(False))
                t /= prod(qd ** la.pairing(la) *
                          prod(1 - qd ** -j for mk in la.to_exp_dict().values()
                               for j in range(1, mk + 1))
                          for la in pt)
                terms += t
        ret += moebius(d) * terms / d

    ret *= fac_r * (q - 1)

    if var is None:
        return ret
    else:
        return ret(q=var)


def kac_polynomial(Q, v, q=None):
    """
    Return the Kac polynomial for the quiver `Q` and the dimension vector `v`.

    INPUT:

    - `Q` -- a quiver as a digraph
    - `v` -- the dimension vector
    - `q` -- (optional) the variable `q`

    EXAMPLES::

        sage: from sage.combinat.kac_polynomial import kac_polynomial
        sage: kac_polynomial(DiGraph([[1],[]]), (1,))
        1?
    """
    return cohomology_polynomial(Q, PartitionTuple([[x] for x in v]), q)


def DT_invariant(Q, v, q=None):
    """
    Return the DT-invariant for the quiver `Q` and the dimension vector `v`.

    INPUT:

    - `Q` -- a quiver as a digraph
    - `v` -- the dimension vector
    - `q` -- (optional) the variable `q`

    EXAMPLES::

        sage: from sage.combinat.kac_polynomial import DT_invariant
        sage: DT_invariant(DiGraph([[1],[]]), (1,))
        1?
    """
    return cohomology_polynomial(Q, PartitionTuple([[1] * x for x in v]), q)
