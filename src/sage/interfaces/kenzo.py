r"""
Library interface to Kenzo

Kenzo is a set of lisp functions to compute homology and
homotopy groups of topological spaces.

AUTHORS:

- Miguel Marco, Ana Romero (2019-01): Initial version


For this interface, Kenzo is loaded into ECL which is itself loaded
as a C library in Sage. Kenzo objects in this interface are nothing
but wrappers around ECL objects.



#*****************************************************************************
#       Copyright (C) 2019 Miguel Marco <mmarco@unizar.es>
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
"""


from __future__ import print_function
from __future__ import absolute_import
from six import string_types

from sage.structure.sage_object import SageObject
from sage.homology.homology_group import HomologyGroup
from sage.rings.integer_ring import ZZ
from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup

from sage.libs.ecl import EclObject, ecl_eval, EclListIterator


# Redirection of ECL and Maxima stdout to /dev/null
ecl_eval(r"""(defparameter *dev-null* (make-two-way-stream
              (make-concatenated-stream) (make-broadcast-stream)))""")
ecl_eval("(setf original-standard-output *standard-output*)")
ecl_eval("(setf *standard-output* *dev-null*)")

# Loading and initialization of Kenzo
# Note that it will load kenzo.fas file from $SAGE_LOCAL/lib/ecl/
ecl_eval("(require :kenzo)")
ecl_eval("(in-package :cat)")
ecl_eval("(setf *HOMOLOGY-VERBOSE* nil)")


# defining the auxiliary functions as wrappers over the kenzo ones
chcm_mat = EclObject("chcm-mat")
homologie = EclObject("homologie")
sphere = EclObject("sphere")
crts_prdc = EclObject("crts-prdc")
moore = EclObject("moore")
k_z = EclObject("k-z")
k_z2 = EclObject("k-z2")
echcm = EclObject("echcm")
loop_space = EclObject("loop-space")
tnsr_prdc = EclObject("tnsr-prdc")
typep = EclObject("typep")
classifying_space = EclObject("classifying-space")
suspension = EclObject("suspension")


def Sphere(n):
    r"""
    Return the ``n`` dimensional sphere as a Kenzo simplicial set.

    INPUT:

    - ``n`` -- The dimension of the sphere

    OUTPUT:

    - A KenzoSimplicialSetEXAMPLES::

        sage: from sage.interfaces.kenzo import Sphere # optional - kenzo
        sage: s2 = Sphere(2)                           # optional - kenzo
        sage: s2                                       # optional - kenzo
        [K1 Simplicial-Set]
        sage: [s2.homology(i) for i in range(8)]       # optional - kenzo
        [Z, 0, Z, 0, 0, 0, 0, 0]
    """
    kenzosphere = sphere(n)
    return KenzoSimplicialSet(kenzosphere)


def MooreSpace(m, n):
    r"""
    Return the Moore space ``M(m, n)`` as a Kenzo simplicial set.

    The Moore space ``M(m, n)`` is the space whose n'th homology group
    is isomorphic to the cyclic group of order ``m``, and the rest of the
    homology groups are trivial.

    INPUT:

    - ``m```- A positive integer. The order of the nontrivial homology group.

    - ``n```- The dimension in which the homology is not trivial


    OUTPUT:

    - A KenzoSimplicialSet

    EXAMPLES::

        sage: from sage.interfaces.kenzo import MooreSpace   # optional - kenzo
        sage: m24 = MooreSpace(2,4)                          # optional - kenzo
        sage: m24                                            # optional - kenzo
        [K10 Simplicial-Set]
        sage: [m24.homology(i) for i in range(8)]            # optional - kenzo
        [Z, 0, 0, 0, C2, 0, 0, 0]
    """

    kenzomoore = moore(m, n)
    return KenzoSimplicialSet(kenzomoore)


def EilenbergMacLaneSpace(G, n):
    r"""Return the Eilenberg-MacLane space ``k(G, n)`` as a Kenzo simplicial group

    The Eilenberg-MacLane space ``k(G, n)`` is the space whose has n'th homotopy
    group isomorphic to ``G``, and the rest of the homotopy groups are trivial.

    INPUT:

    - ``G`` -- A group. Currently only ``ZZ`` and the group of two elements are
    supported.

    - ``n`` -- The dimension in which the homotopy is not trivial


    OUTPUT:

    - A KenzoSimplicialGroup

    EXAMPLES::

        sage: from sage.interfaces.kenzo import EilenbergMacLaneSpace
        sage: e3 = EilenbergMacLaneSpace(ZZ, 3)                          # optional - kenzo
        sage: [e3.homology(i) for i in range(8)]                         # optional - kenzo
        [Z, 0, 0, Z, 0, C2, 0, C3]
        sage: f3 = EilenbergMacLaneSpace(AdditiveAbelianGroup([2]), 3)   # optional - kenzo
        sage: [f3.homology(i) for i in range(8)]                         # optional - kenzo
        [Z, 0, 0, C2, 0, C2, C2, C2]
    """
    if G == ZZ:
        kenzospace = k_z(n)
        return KenzoSimplicialGroup(kenzospace)
    elif G == AdditiveAbelianGroup([2]):
        kenzospace = k_z2(n)
        return KenzoSimplicialGroup(kenzospace)
    else:
        raise NotImplementedError("Eilenberg-MacLane spaces are only supported over ZZ and ZZ_2")


class KenzoChainComplex(SageObject):
    r"""
    Wraper to kenzo chain complexes

    INPUT:

        - ``kenzoobject``- a wrapper around a kenzo chain complex
        (which is an ecl object).
    """

    def __init__(self, kenzoobject):
        r"""
        Construct the chain complex.

        TESTS::

            sage: from sage.interfaces.kenzo import KenzoChainComplex  # optional -kenzo
            sage: from sage.interfaces.kenzo import sphere             # optional -kenzo
            sage: ks = sphere(2)                                       # optional -kenzo
            sage: ks                                                   # optional -kenzo
            <ECL: [K1 Simplicial-Set]>
            sage: s2 = KenzoChainComplex(ks)                           # optional -kenzo
            sage: s2                                                   # optional -kenzo
            [K1 Simplicial-Set]
            sage: TestSuite(s2).run(skip='_test_pickling')             # optional -kenzo

        """
        self._kenzo = kenzoobject

    def _repr_(self):
        r"""
        Represent the object. It just uses the ecl representation, removing the
        ecl decoration

        EXAMPLES::

            sage: from sage.interfaces.kenzo import MooreSpace    # optional - kenzo
            sage: m2 = MooreSpace(2,4)                            # optional - kenzo
            sage: m2._repr_()                                     # optional - kenzo
            '[K10 Simplicial-Set]'
        """
        kenzo_string = self._kenzo.__repr__()
        return kenzo_string[6:-1]

    def homology(self, n):
        r"""
        Return the ``n``'th homology group of the kenzo chain complex

        INPUT:

        - ``n`` -- The dimension in which compute the homology

        EXAMPLES::

                sage: from sage.interfaces.kenzo import Sphere   # optional - kenzo
                sage: s2 = Sphere(2)                             # optional - kenzo
                sage: s2                                         # optional - kenzo
                [K1 Simplicial-Set]
                sage: s2.homology(2)                             # optional - kenzo
                Z
        """
        echcm1 = echcm(self._kenzo)
        m1 = chcm_mat(echcm1, n)
        m2 = chcm_mat(echcm1, n+1)
        homology = homologie(m1, m2)
        lhomomology = [i for i in EclListIterator(homology)]
        res = []
        for component in lhomomology:
            pair = [i for i in EclListIterator(component)]
            res.append(pair[0].python())
        return HomologyGroup(len(res), ZZ, res)

    def tensor_product(self, other):
        r"""
        Return the tensor product of ``self`` and ``other``.

        INPUT:

        - ``other`` --  The Kenzo object with which to compute the tensor product

        OUTPUT:

        - A KenzoChainComplex

        sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
        sage: s2 = Sphere(2)                              # optional - kenzo
        sage: s3 = Sphere(3)                              # optional - kenzo
        sage: p = s2.tensor_product(s3)                   # optional - kenzo
        sage: type(p)                                     # optional - kenzo
        <class 'sage.interfaces.kenzo.KenzoChainComplex'>
        sage: [p.homology(i) for i in range(8)]           # optional - kenzo
        [Z, 0, Z, Z, 0, Z, 0, 0]
        """
        return KenzoChainComplex(tnsr_prdc(self._kenzo, other._kenzo))


class KenzoSimplicialSet(KenzoChainComplex):
    r"""
    Wrapper to Kenzo simplicial sets.
    """

    def loop_space(self, n=1):
        r"""
        Return the ``n`` th iterated loop space.

        INPUT:

        - ``n`` -- (default: 1) the number of times to iterate the loop space
        construction

        OUTPUT:

        - A KenzoSimplicialGroup

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: l2 = s2.loop_space()                        # optional - kenzo
            sage: type(l2)                                    # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialGroup'>
            sage: l2 = s2.loop_space()                        # optional - kenzo
            sage: [l2.homology(i) for i in range(8)]          # optional - kenzo
            [Z, Z, Z, Z, Z, Z, Z, Z]
        """
        return KenzoSimplicialGroup(loop_space(self._kenzo, n))

    def cartesian_product(self, other):
        r"""
        Return the cartesian product of ``self`` and ``other``.

        INPUT:

        - ``other`` -- The Kenzo simplicial set with which the product is made

        OUTPUT:

        - A KenzoSimplicialSet

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: s3 = Sphere(3)                              # optional - kenzo
            sage: p = s2.cartesian_product(s3)                # optional - kenzo
            sage: type(p)                                     # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [p.homology(i) for i in range(6)]           # optional - kenzo
            [Z, 0, Z, Z, 0, Z]
        """
        prod_kenzo = crts_prdc(self._kenzo, other._kenzo)
        return KenzoSimplicialSet(prod_kenzo)

    def suspension(self):
        r"""
        Return the suspension of the simplicial set.

        OUTPUT:

        - A KenzoSimplicialSet

        EXAMPLES::

            sage: from sage.interfaces.kenzo import EilenbergMacLaneSpace    # optional - kenzo
            sage: e3 = EilenbergMacLaneSpace(ZZ, 3)                          # optional - kenzo
            sage: s = e3.suspension()                                        # optional - kenzo
            sage: type(s)                                                    # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [s.homology(i) for i in range(6)]                          # optional - kenzo
            [Z, 0, 0, 0, Z, 0]
        """
        return KenzoSimplicialSet(suspension(self._kenzo))


class KenzoSimplicialGroup(KenzoSimplicialSet):
    r"""
    Wrapper around Kenzo simplicial groups.
    """

    def classifying_space(self):
        r"""
        Return the classifying space.

        OUTPUT:

        - A KenzoSimplicialGroup

        EXAMPLES::
                                                                    # optional - kenzo
            sage: from sage.interfaces.kenzo import MooreSpace      # optional - kenzo
            sage: m2 = MooreSpace(2,4)                              # optional - kenzo
            sage: l2 = m2.loop_space()                              # optional - kenzo
            sage: c = l2.classifying_space()                        # optional - kenzo
            sage: type(c)                                           # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialGroup'>
            sage: [c.homology(i) for i in range(8)]                  # optional - kenzo
            [Z, 0, 0, 0, C2, 0, 0, 0]
        """
        return KenzoSimplicialGroup(classifying_space(self._kenzo))
