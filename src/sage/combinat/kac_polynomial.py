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

from sage.misc.misc import prod
from sage.rings.all import ZZ
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.partition_tuple import PartitionTuple, PartitionTuples
from sage.combinat.partition import Partitions
from sage.combinat.cartesian_product import CartesianProduct
from sage.rings.arith import moebius

def cohomology_polynomial(Q, mu, q):
    r"""
    Return the cohomology polynomial.

    The cohomology polynomial is given by:

    .. MATH::

        \sum_i \dim\left( H_c^{2i}(\mathcal{Q}_{\tilde{v}} ; \CC)_{\epsilon
        \chi^{\mu}} \right) q^{i - d_{\tilde{v}}} = \mathbb{H}_{\mu}^s(q)

    where `\widetilde{v}` is the extend dimension vector given in [HLRV2012]_,
    `\epsilon` is the sign character, and

    .. MATH::

        \mathbb{H}_{\mu}^s(q) := \bigl\langle \mathbb{H}(x_1, \ldots, x_r; q),
        s_{\mu} \bigr\rangle,

    where `\mathbb{H}` is given by Equation (1.4) in [HLRV2012]_.

    REFERENCES:

    .. [HLRV2012] Tamas Hausel, Emmanuel Letellier, and Fernando
       Rodriguez-Villegas. *Positivity of Kac polynomials and DT-invariants
       for quivers*. (2012). :arxiv:`1204.2375v1`.
    """
    q_poly = ZZ['q'].gen(0)
    FF = q_poly.parent().fraction_field()
    qp = FF.gen(0)
    Sym = SymmetricFunctions(FF)
    s = Sym.s()
    HLP = Sym.hall_littlewood(qp).P()
    n = ZZ(mu.size())
    r = mu.level()

    ret = FF.zero()
    for d in n.divisors(): # We must have d dividing n
        qd = qp**d
        terms = FF.zero()
        # We only have a chance if each partition's size times d is equal to
        #   the size of the corresponding partition in mu
        PT = CartesianProduct(*[Partitions(p.size() // d) for p in mu.components()])

        for pt in PT:
            Z = zip(mu.components(), pt)
            t = FF.prod(s(HLP[p]).adams_operation(d).coefficient(m) for m,p in Z)
            if t != 0: # If there's something to do
                t *= FF.prod(qd**(pt[i].pairing(pt[j])) for i,j in Q.edges(False))
                t /= FF.prod(qd**(la.pairing(la)) * FF.prod(1 - qd**-j for mk in la.to_exp_dict().values()
                                                          for j in range(1, mk+1))
                             for la in pt)
                terms += t
        ret += moebius(d) * terms / d

    ret = (-1)**(r-1) / r * (qp - 1) * ret
    return ret
    if q is None:
        q = q_poly
    P = q.parent()
    return P(ret.substitute(qp=q))

def kac_polynomial(Q, v, q=None):
    """
    Return the Kac polynomial.

    INPUT:

    - ``Q`` -- a quiver as a digraph
    - ``v`` -- the dimension vector
    - ``q`` -- (optional) the variable `q`
    """
    return cohomology_polynomial(Q, PartitionTuple([[x] for x in v]), q)

def DT_invariant(Q, v, q=None):
    """
    Return the DT-invariant.

    INPUT:

    - ``Q`` -- a quiver as a digraph
    - ``v`` -- the dimension vector
    - ``q`` -- (optional) the variable `q`
    """
    return cohomology_polynomial(Q, PartitionTuple([[1]*x for x in v]), q)

