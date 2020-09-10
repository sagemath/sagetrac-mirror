r"""
Library interface to Kenzo

Kenzo is a set of lisp functions to compute homology and
homotopy groups of topological spaces.

AUTHORS:

- Miguel Marco, Ana Romero (2019-01): Initial version


For this interface, Kenzo is loaded into ECL which is itself loaded
as a C library in Sage. Kenzo objects in this interface are nothing
but wrappers around ECL objects.
"""
# ****************************************************************************
#       Copyright (C) 2019 Miguel Marco <mmarco@unizar.es>
#                      and Ana Romero <ana.romero@unirioja.es>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import print_function, absolute_import


from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject
from sage.homology.homology_group import HomologyGroup
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup
from sage.groups.abelian_gps.abelian_group import AbelianGroup
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups

from sage.matrix.all import matrix

from sage.libs.ecl import EclObject, ecl_eval, EclListIterator
from sage.features.kenzo import Kenzo


# defining the auxiliary functions as wrappers over the kenzo ones
kenzo_names = ['2absm-acrpr',
               '2cmbn-add',
               '2cmbn-sbtr',
               'absm1',
               'absm2',
               'absm-aux',
               'add',
               'array-dimensions',
               'assoc-to-function',
               'basis-aux',
               'bicomplex-spectral-sequence',
               'build-finite-ss2',
               'build-mrph-aux',
               'change-sorc-trgt-aux',
               'chcm-mat',
               'chcm-mat2',
               'classifying-space',
               'cmbn-aux',
               'cmbn-degr',
               'cmbn-list',
               'cmbn-opps',
               'cmpr-aux',
               'cmps',
               'cone2',
               'constant-smmr',
               'convertmatrice',
               'crpr-absms-aux',
               'crts-prdc',
               'degenerate-p',
               'degr-aux',
               'delta',
               'dffr-aux',
               'dffr-aux1',
               'dgop',
               'dgop-int-ext',
               'dstr-change-sorc-trgt-aux',
               'echcm',
               'eilenberg-moore-spectral-sequence',
               'evaluate-cmbn',
               'evaluate-simplex',
               'evaluation-aux',
               'gmsm',
               'homologie',
               'homotopy-list',
               'idnm',
               'idnt-mrph',
               'idnt-smmr',
               'join',
               'k',
               'k-z',
               'k-z2',
               'k-zp',
               'kabstractsimplex-aux',
               'kchaincomplex-aux',
               'kchaincomplexmorphism-aux',
               'ksimplicialsetmorphism-aux',
               'listp',
               'loop-space',
               'make-array-from-lists',
               'make-array-to-lists',
               'moore',
               'n-cmbn',
               'n-mrph',
               'ncol',
               'nlig',
               'non-degenerate-p',
               'nreverse',
               'nth',
               'opps',
               'orgn-aux',
               'pushout',
               'r-proj-space',
               'sbtr',
               'serre-spectral-sequence-product',
               'serre-whitehead-spectral-sequence',
               'sfinitesimplicialset-aux',
               'smash-product',
               'sorc-aux',
               'spectral-sequence-differential-matrix',
               'spectral-sequence-group',
               'sphere',
               'suspension',
               'tnsr-prdc',
               'trgt-aux',
               'type-of',
               'wedge',
               'zero-cmbn',
               'zero-mrph']


# Now initialize Kenzo. For each string s in kenzo_names, the
# following defines __s__, a wrapper for a Kenzo function. For
# example __sphere__ is defined as EclObject("sphere"). Hyphens
# are replaced with underscores to get valid Python identifiers.
if Kenzo().is_present():
    ecl_eval("(require :kenzo)")
    ecl_eval("(in-package :cat)")
    ecl_eval("(setf *HOMOLOGY-VERBOSE* nil)")
    for s in kenzo_names:
        name = '__{}__'.format(s.replace('-', '_'))
        exec('{} = EclObject("{}")'.format(name, s))


def Sphere(n):
    r"""
    Return the ``n`` dimensional sphere as a Kenzo simplicial set.

    INPUT:

    - ``n`` -- the dimension of the sphere

    OUTPUT:

    - A :class:`KenzoSimplicialSet`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import Sphere # optional - kenzo
        sage: s2 = Sphere(2)                           # optional - kenzo
        sage: s2                                       # optional - kenzo
        [K... Simplicial-Set]
        sage: [s2.homology(i) for i in range(8)]       # optional - kenzo
        [Z, 0, Z, 0, 0, 0, 0, 0]
    """
    kenzosphere = __sphere__(n)
    return KenzoSimplicialSet(kenzosphere)


def MooreSpace(m, n):
    r"""
    Return the Moore space ``M(m, n)`` as a Kenzo simplicial set.

    The Moore space ``M(m, n)`` is the space whose n'th homology group
    is isomorphic to the cyclic group of order ``m``, and the rest of the
    homology groups are trivial.

    INPUT:

    - ``m`` -- A positive integer. The order of the nontrivial homology group.

    - ``n`` -- The dimension in which the homology is not trivial

    OUTPUT:

    - A KenzoSimplicialSet

    EXAMPLES::

        sage: from sage.interfaces.kenzo import MooreSpace   # optional - kenzo
        sage: m24 = MooreSpace(2,4)                          # optional - kenzo
        sage: m24                                            # optional - kenzo
        [K... Simplicial-Set]
        sage: [m24.homology(i) for i in range(8)]            # optional - kenzo
        [Z, 0, 0, 0, C2, 0, 0, 0]
    """
    kenzomoore = __moore__(m, n)
    return KenzoSimplicialSet(kenzomoore)


def EilenbergMacLaneSpace(G, n):
    r"""
    Return the Eilenberg-MacLane space ``K(G, n)`` as a Kenzo simplicial group.

    The Eilenberg-MacLane space ``K(G, n)`` is the space whose has n'th homotopy
    group isomorphic to ``G``, and the rest of the homotopy groups are trivial.

    INPUT:

    - ``G`` -- group. Currently only ``ZZ`` and ``ZZ_n`` are supported.

    - ``n`` -- the dimension in which the homotopy is not trivial

    OUTPUT:

    - A :class:`KenzoSimplicialGroup`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import EilenbergMacLaneSpace    # optional - kenzo
        sage: e3 = EilenbergMacLaneSpace(ZZ, 3)                          # optional - kenzo
        sage: [e3.homology(i) for i in range(8)]                         # optional - kenzo
        [Z, 0, 0, Z, 0, C2, 0, C3]
        sage: f3 = EilenbergMacLaneSpace(AdditiveAbelianGroup([2]), 3)   # optional - kenzo
        sage: [f3.homology(i) for i in range(8)]                         # optional - kenzo
        [Z, 0, 0, C2, 0, C2, C2, C2]
    """
    if G == ZZ:
        kenzospace = __k_z__(n)
        return KenzoSimplicialGroup(kenzospace)
    elif G == AdditiveAbelianGroup([2]):
        kenzospace = __k_z2__(n)
        return KenzoSimplicialGroup(kenzospace)
    elif G in CommutativeAdditiveGroups() and G.is_cyclic():
        kenzospace = __k_zp__(G.cardinality(), n)
        return KenzoSimplicialGroup(kenzospace)
    else:
        raise NotImplementedError("Eilenberg-MacLane spaces are only supported over ZZ and ZZ_n")

def KenzoId(n):
    r"""
    Return the KenzoObject whose identification number is ``n``.

    INPUT:

    - ``n`` -- A positive integer.

    EXAMPLES::

        sage: from sage.interfaces.kenzo import Sphere, KenzoId, IdNumber    # optional -kenzo
        sage: s2 = Sphere(2)                                                 # optional -kenzo
        sage: s2 == KenzoId(IdNumber(s2))                                    # optional -kenzo
        True
    """
    kenzoobj = __k__(n)
    if str(kenzoobj) == 'NIL':
        raise AssertionError("There is no KenzoObject with identification number {}".format(n))
    return translate_type(kenzoobj)

def IdNumber(kobject):
    r"""
    Return the slot :idnm of the KenzoObject ``kobject`` if it has such a slot.

    INPUT:

    - ``kobject`` -- A KenzoObject

    OUTPUT:

    -  A positive integer, number plate for this object. This is generated by the
       Kenzo system in a sequential way, each time a new Kenzo object is created.
       Not all objects in Kenzo has :idnm slot (for example abstract simplexes absm).

    EXAMPLES::

        sage: from sage.interfaces.kenzo import MooreSpace, KenzoId, IdNumber  # optional -kenzo
        sage: m34 = MooreSpace(3,4)                                            # optional -kenzo
        sage: MooreSpace(3,4) == KenzoId(IdNumber(m34))                        # optional -kenzo
        True
    """
    return Integer(str(__idnm__(kobject._kenzo)))


class KenzoObject(SageObject):
    r"""
    Wrapper to Kenzo objects

    INPUT:

    - ``kenzo_object`` -- a wrapper around a Kenzo object
      (which is an ecl object).
    """

    def __init__(self, kenzo_object):
        r"""
        Construct the Kenzo object.

        TESTS::

            sage: from sage.interfaces.kenzo import KenzoObject  # optional -kenzo
            sage: from sage.interfaces.kenzo import __sphere__   # optional -kenzo
            sage: ks = __sphere__(2)                             # optional -kenzo
            sage: ks                                             # optional -kenzo
            <ECL: [K... Simplicial-Set]>
            sage: s2 = KenzoObject(ks)                           # optional -kenzo
            sage: s2                                             # optional -kenzo
            [K... Simplicial-Set]
            sage: TestSuite(s2).run(skip='_test_pickling')       # optional -kenzo

        """
        self._kenzo = kenzo_object

    def _repr_(self):
        r"""
        Represent the object.

        It just uses the ecl representation, removing the
        ecl decoration.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import MooreSpace     # optional - kenzo
            sage: m24 = MooreSpace(2,4)                            # optional - kenzo
            sage: m24._repr_()                                     # optional - kenzo
            '[K... Simplicial-Set]'
        """
        kenzo_string = repr(self._kenzo)
        return kenzo_string[6:-1]


class KenzoSpectralSequence(KenzoObject):
    r"""
    Wrapper around Kenzo spectral sequences
    """

    def group(self, p, i, j):
        r"""
        Return the ``i,j``'th group of the ``p`` page.

        INPUT:

        - ``p`` -- the page to take the group from.

        - ``i`` -- the column where the group is taken from.

        - ``j`` -- the row where the group is taken from.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere # optional - kenzo
            sage: S2 = Sphere(2)                           # optional - kenzo
            sage: EMS = S2.em_spectral_sequence()          # optional - kenzo
            sage: EMS.group(0, -1, 2)                      # optional - kenzo
            Additive abelian group isomorphic to Z
            sage: EMS.group(0, -1, 3)                      # optional - kenzo
            Trivial group
        """
        invs = __spectral_sequence_group__(self._kenzo, p, i, j).python()
        if not invs:
            invs = []
        return AdditiveAbelianGroup(invs)

    def matrix(self, p, i, j):
        r"""
        Return the matrix that determines the differential from the
        ``i,j``'th group of the ``p``'th page.

        INPUT:

        - ``p`` -- the page.

        - ``i`` -- the column of the differential domain.

        - ``j`` -- the row of the differential domain.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere   # optional - kenzo
            sage: S3 = Sphere(3)                             # optional - kenzo
            sage: L = S3.loop_space()                        # optional - kenzo
            sage: EMS = L.em_spectral_sequence()             # optional - kenzo
            sage: EMS.table(1, -5, -2, 5, 8)                 # optional - kenzo
              0   Z   Z + Z + Z   Z + Z + Z
              0   0   0           0
              0   0   Z           Z + Z
              0   0   0           0
            sage: EMS.matrix(1, -2 ,8)                       # optional - kenzo
            [ 3 -2  0]
            [ 3  0 -3]
            [ 0  2 -3]
        """
        klist = __spectral_sequence_differential_matrix__(self._kenzo, p, i, j)
        plist = klist.python()
        if plist is None or plist == [None]:
            i = len(self.group(p, i, j).invariants())
            j = len(self.group(p, i - p, j + p - 1).invariants())
            return matrix(i, j)
        return matrix(plist)

    def differential(self, p, i, j):
        r"""
        Return the ``(p, i, j)`` differential morphism of the spectral sequence.

        INPUT:

        - ``p`` -- the page.

        - ``i`` -- the column of the differential domain.

        - ``j`` -- the row of the differential domain.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere   # optional - kenzo
            sage: S3 = Sphere(3)                             # optional - kenzo
            sage: L = S3.loop_space()                        # optional - kenzo
            sage: EMS = L.em_spectral_sequence()             # optional - kenzo
            sage: EMS.table(1,-5,-2,5,8)                     # optional - kenzo
              0   Z   Z + Z + Z   Z + Z + Z
              0   0   0           0
              0   0   Z           Z + Z
              0   0   0           0
            sage: EMS.matrix(1, -3, 8)                       # optional - kenzo
            [ 2 -2  2]
            sage: EMS.differential(1, -3, 8)                 # optional - kenzo
            Morphism from module over Integer Ring with invariants (0, 0, 0) to module with invariants (0,) that sends the generators to [(2), (-2), (2)]
        """
        domain = self.group(p, i, j)
        codomain = self.group(p, i - p, j + p - 1)
        M = self.matrix(p, i, j)
        images = [codomain(r) for r in M.columns()]
        return domain.hom(images, codomain=codomain)

    def table(self, p, i1, i2, j1, j2):
        r"""
        Return a table printing the groups in the ``p`` page.

        INPUT:

        - ``p`` -- the page to print.

        -- ``i1`` -- the first column to print.

        -- ``i2`` -- the last column to print.

        -- ``j1`` -- the first row to print.

        -- ``j2`` -- the last row to print.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere # optional - kenzo
            sage: S2 = Sphere(2)                           # optional - kenzo
            sage: EMS = S2.em_spectral_sequence()          # optional - kenzo
            sage: EMS.table(0, -2, 2, -2, 2)               # optional - kenzo
              0   Z   0   0   0
              0   0   0   0   0
              0   0   Z   0   0
              0   0   0   0   0
              0   0   0   0   0
        """
        from sage.misc.table import table
        groups = []
        for j in range(j2 - j1 + 1):
            row = []
            for i in range(i1, i2 + 1):
                group = self.group(p, i, j2 - j)
                if group.invariants():
                    row.append(group.short_name())
                else:
                    row.append('0')
            groups.append(row)
        return table(groups)


class KenzoAbstractSimplex(KenzoObject):
    r"""
    Wrapper to the kenzo objects of type :absm. In Kenzo, an abstract simplex
    is a pair consisting of:

    - A (possibly iterated) degeneracy operator.

    - A 'geometric' simplex, i.e. a non-degenerate simplex.

    An abstract simplex is represented internally in the Kenzo system by the lisp
    object <AbSm ext-dgop gmsm> where:

    - ext-dgop is a strictly decreasing sequence of integer numbers representing
       a list of degeneracy operators.

    - gmsm is a geometric simplex, i.e. any kind lisp object modelizing a
      non-degenerate simplex, to which is applied the degeneracy operators ext-dgop.
    """
    def degeneracies(self):
        r"""
        Returns the list representing the degeneracy operators of ``self``.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(3, (0,1,3), name = 'v')
            sage: kv = v._kenzo_repr                                           # optional - kenzo
            sage: kv                                                           # optional - kenzo
            <AbSm 5-2-0 S...>
            sage: kv.degeneracies()                                            # optional - kenzo
            [5, 2, 0]
            sage: kv.degeneracies() == v.degeneracies()                        # optional - kenzo
            True
        """
        rslt = __dgop_int_ext__(__dgop__(self._kenzo)).python()
        if rslt == None:
            return []
        return rslt

    def nondegenerate(self):
        r"""
        Returns the basic non-degenerate simplex part of ``self``.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(3, (0,1,3), name = 'v')
            sage: w = AbstractSimplex(2, (0,1), name = 'w')
            sage: kv = v._kenzo_repr                                           # optional - kenzo
            sage: kv                                                           # optional - kenzo
            <AbSm 5-2-0 S...>
            sage: kw = w._kenzo_repr                                           # optional - kenzo
            sage: kw                                                           # optional - kenzo
            <AbSm 2-0 S...>
            sage: pr = kv.product(kw)                                          # optional - kenzo
            sage: pr                                                           # optional - kenzo
            <AbSm 2-0 <CrPr 3 S... - S...>>
            sage: pr.nondegenerate()                                           # optional - kenzo
            <CrPr 3 S... - S...>
            sage: pr.nondegenerate().factors()                                 # optional - kenzo
            (<AbSm 3 S...>, <AbSm - S...>)
        """
        return translate_type(__gmsm__(self._kenzo))

    def is_degenerate(self):
        r"""
        Return True if ``self`` is degenerate, False otherwise.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(3, (0,1,3), name = 'v')
            sage: kv = v._kenzo_repr                                           # optional - kenzo
            sage: kv.degeneracies()                                            # optional - kenzo
            [5, 2, 0]
            sage: kv.is_degenerate()                                           # optional - kenzo
            True
        """
        rslt = __degenerate_p__(self._kenzo).python()
        if rslt == None:
            return False
        return rslt

    def is_nondegenerate(self):
        r"""
        Return True if ``self`` is non-degenerate, False otherwise.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(3, name = 'v')
            sage: kv = v._kenzo_repr                                           # optional - kenzo
            sage: kv.is_nondegenerate()                                        # optional - kenzo
            True
        """
        rslt = __non_degenerate_p__(self._kenzo).python()
        if rslt == None:
            return False
        return rslt

    def product(self, other):
        r"""
        Returns the cartesian product of both KenzoAbstractSimplexes ``self`` and ``other``.

        INPUT:

        - ``other`` - A :class:`KenzoAbstractSimplex`.

        OUTPUT:

        - If ``self`` and ``other`` are degenerate and the degeneracy operators
            have no common `\eta_j`, a KenzoCRPRSimplex is returned, otherwise a
            KenzoAbstractSimplex is returned.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(3, (0,1,3), name = 'v')
            sage: w = AbstractSimplex(2, (0,1), name = 'w')
            sage: x = AbstractSimplex(5, name = 'x')
            sage: kv = v._kenzo_repr
            sage: kv
            <AbSm 5-2-0 S...>
            sage: kw = w._kenzo_repr
            sage: kw
            <AbSm 2-0 S...>
            sage: kx = x._kenzo_repr
            sage: kx
            <AbSm - S...>
            sage: prvw = kv.product(kw)
            sage: prvw
            <AbSm 2-0 <CrPr 3 S... - S...>>
            sage: prvx = kv.product(kx)
            sage: prvx
            <CrPr 5-2-0 S... - S...>
        """
        if not isinstance(other, KenzoAbstractSimplex):
            raise TypeError("{} must be a KenzoAbstractSimplex".format(other))
        deg1 = self.degeneracies()
        deg2 = other.degeneracies()
        if list(set(deg1) & set(deg2)) == []:
            return KenzoCRPRSimplex(__crpr_absms_aux__(self._kenzo, other._kenzo))
        else:
            return KenzoAbstractSimplex(__2absm_acrpr__(self._kenzo, other._kenzo))


class KenzoCRPRSimplex(KenzoObject):
    r"""
    Wrapper to the Kenzo objects of type :crpr. The
    non-degenerate simplices of the cartesian product `X\times Y` of simplcial
    sets are represented internally in the Kenzo system by a lisp object of the
    form <CrPr ext-dgop1 gmsm1 ext-dgop2 gmsm2> where:

    - ext-dgop1 is a sequence of degeneracy operators.

    - gmsm1 is a non-degenerate simplex of `X`, to which is applied the degeneracy
      operators ext-dgop1.

    - ext-dgop2 is a sequence of degeneracy operators.

    - gmsm2 is a non-degenerate simplex of `Y`, to which is applied the degeneracy
      operators ext-dgop2.

    The degeneracy operators in ext-dgop1 and ext-dgop2 must have a void intersection.
    """
    def factors(self, n=None):
        r"""
        Returns the factors of ``self``, that is, the KenzoAbstractSimplexes
        `A_1` and `A_2` such that ``self`` = `A_1\times A_2`.

        OUTPUT:

        - A tuple of two elements: the factors of ``self``.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: K = T._kenzo_repr                                  # optional - kenzo
            sage: crpr = K.basis(1)[1]                               # optional - kenzo
            sage: crpr                                               # optional - kenzo
            <CrPr - S1 0 *>
            sage: A1 = crpr.factors(0)                               # optional - kenzo
            sage: A1                                                 # optional - kenzo
            <AbSm - S1>
            sage: A2 = crpr.factors(1)                               # optional - kenzo
            sage: A2                                                 # optional - kenzo
            <AbSm 0 *>
            sage: A1.product(A2)                                     # optional - kenzo
            <CrPr - S1 0 *>
        """
        factor1 = KenzoAbstractSimplex(__absm1__(self._kenzo))
        factor2 = KenzoAbstractSimplex(__absm2__(self._kenzo))
        rslt = (factor1, factor2)
        if n==0:
            return factor1
        elif n==1:
            return factor2
        else:
            return rslt


class KenzoCombination():
    r"""
    Wrapper to Kenzo combinations (type :cmbn) and functions handling combinations.
    An element `\sigma_i` of any basis corresponding to one of the groups `C_p`
    (`p \in \ZZ`) of a Kenzo chain complex is called a *generator*. A product
    such as `\lambda_i\sigma_i` where `\lambda_i\in\ZZ` is called a *term* and a
    sum of terms, a *combination*.

    A Kenzo combination contains:

    - The degree of the combination corresponding to the index `p\in\ZZ` of the
      group `C_p` to which the combination belongs.

    - The list of the internal representation of the terms, namely the list of
      pairs `(\lambda_i\ .\ \sigma_i)`.

    This choice of representation in Kenzo implies that only homogenous combinations
    will be considered.
    """
    def __init__(self, degree, terms=None, Kchcm=None, kenzo=None):
        r"""
        Construct the combination of degree ``degree`` sum of the terms
        determined by ``terms`` (in the KenzoChainComplex ``Kchcm`` when given).

        INPUT:

        - ``degree`` -- An integer number.

        - ``terms`` -- A list of even length (default None). The even indexes
          determine the coefficients of the generators, which are determined by
          the odd indexes. For example, the ``terms`` parameter of a combination
          like 3A - 5B + C is the list [3, A, -5, B, 1, C].

        - ``Kchcm`` -- A KenzoChainComplex or None (default). In case it is not
          None, the system checks if each one of the generators in ``terms`` is
          indeed in the basis of `C_p`, where `p=` ``degree``.

        - ``kenzo`` -- The kenzo representation of the combination (default None).

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KenzoCombination             # optional - kenzo
            sage: D = ChainComplex({6: matrix(ZZ, 2, 2, [1,0,0,2])}, degree = -1)
            sage: K = D._kenzo_repr                                              # optional - kenzo
            sage: KenzoCombination(5, [-4, 'G5G0'], K)                           # optional - kenzo
            <BLANKLINE>
            {K... CMBN 5}
            --------------------------------------
            <-4 * G5G0>
            --------------------------------------
            <BLANKLINE>
            sage: KenzoCombination(19, Kchcm=K)                                        # optional - kenzo
            <BLANKLINE>
            {K... CMBN 19}
            --------------------------------------
            <BLANKLINE>
            --------------------------------------
            <BLANKLINE>
        """
        if not degree.is_integer():
            raise AssertionError("degree parameter must be an integer number")
        if Kchcm and not isinstance(Kchcm, KenzoChainComplex):
            raise AssertionError("Kchcm parameter must be a KenzoChainComplex or None")
        self._degree = degree
        self._terms = terms
        self._chcm = Kchcm
        self._kenzo = kenzo
        if terms:
            if not kenzo:
                self._kenzo = __cmbn_aux__(degree, pairing(terms))
            if Kchcm:
                strbasis = [str(x) for x in Kchcm.basis(degree)]
                if strbasis==[]:
                    raise AssertionError("There are no generators in dimension {}"
                                          .format(degree))
                for k in range(1, len(terms), 2):
                    if not str(terms[k]) in strbasis:
                        raise AssertionError("Generator {} is not of degree {} in {}"
                                             .format(terms[k], degree, Kchcm))
        else:
            if not kenzo:
                self._kenzo = __zero_cmbn__(degree)

    def __repr__(self):
        r"""
        Print representation.

        EXAMPLES::

            sage: C = ChainComplex({2: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])}, degree=-1)
            sage: K = C._kenzo_repr                                              # optional - kenzo
            sage: K.basis(2)                                                     # optional - kenzo
            ['G2G0', 'G2G1', 'G2G2']
            sage: c1 = K(2, [1, 'G2G0', -9, 'G2G1', 6, 'G2G2']); c1              # optional - kenzo
            <BLANKLINE>
            {K... CMBN 2}
            --------------------------------------
            <1 * G2G0>
            <-9 * G2G1>
            <6 * G2G2>
            --------------------------------------
            <BLANKLINE>
        """
        if self._chcm:
            s = "{{K{} CMBN {}}}".format(__idnm__(self._chcm._kenzo).python(),
                                         self._degree)
        else:
            s = "{{CMBN {}}}".format(self._degree)
        sep = "\n" + "-"*38
        s += sep
        if self._terms:
            for k in range(0, len(self._terms), 2):
                s += "\n<" + str(self._terms[k]) + " * " + str(self._terms[k+1]) + ">"
            return s + sep
        return s + "\n" + sep

    def __eq__(self, other):
        r"""
        Return ``True`` if and only if ``self == other``.

        EXAMPLES::

            sage: C = ChainComplex({3: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])}, degree=-1)
            sage: K1 = C._kenzo_repr                                             # optional - kenzo
            sage: null6 = K1(6)                                                  # optional - kenzo
            sage: null6 == K1(6)                                                 # optional - kenzo
            True
            sage: null6 == K1(5)                                                 # optional - kenzo
            False
            sage: c1 = K1(3, [-4, 'G3G1'])                                       # optional - kenzo
            sage: c1 == c1 + K1(3)                                               # optional - kenzo
            True
            sage: D = ChainComplex({3: matrix(ZZ, 2, 2, [3, 0, 0, 0])}, degree=-1)
            sage: K2 = D._kenzo_repr                                             # optional - kenzo
            sage: c2 = K2(3, [-4, 'G3G1'])                                       # optional - kenzo
            sage: c1 == c2                                                       # optional - kenzo
            False
        """
        return bool(self._chcm==other._chcm \
                    and self._degree==other._degree \
                    and self._terms==other._terms)

    def __neg__(self):
        r"""
        Return ``-self``.

        EXAMPLES::

            sage: E = ChainComplex({0: matrix(ZZ, 2, 2, [1,0,0,2])}, degree=-1)
            sage: K = E._kenzo_repr                                              # optional - kenzo
            sage: c = K(-1, [7, 'G-1G0', -3, 'G-1G1']); c                        # optional - kenzo
            <BLANKLINE>
            {K... CMBN -1}
            --------------------------------------
            <7 * G-1G0>
            <-3 * G-1G1>
            --------------------------------------
            <BLANKLINE>
            sage: -c                                                             # optional - kenzo
            <BLANKLINE>
            {K... CMBN -1}
            --------------------------------------
            <-7 * G-1G0>
            <3 * G-1G1>
            --------------------------------------
            <BLANKLINE>
            sage: -c == K(-1) - c                                                # optional - kenzo
            True
        """
        return Kenzocmbn(__cmbn_opps__(self._kenzo), self._chcm)

    def __add__(self, other):
        r"""
        Return ``self + other``.

        EXAMPLES::

            sage: C = ChainComplex({13: matrix(ZZ, 2, 3, [5, 0, 0, 0, 0, 0]), 14: matrix(ZZ, [[0, 1]])})
            sage: K = C.dual()._kenzo_repr                                       # optional - kenzo
            sage: c1 = K(13, [1, 'G13G0', -2, 'G13G1', 3, 'G13G2'])              # optional - kenzo
            sage: c2 = K(13, [2, 'G13G1', -2, 'G13G2'])                          # optional - kenzo
            sage: c1 + c2                                                        # optional - kenzo
            <BLANKLINE>
            {K... CMBN 13}
            --------------------------------------
            <1 * G13G0>
            <1 * G13G2>
            --------------------------------------
            <BLANKLINE>
            sage: -(-c1 + -c2) == c1 + c2                                        # optional - kenzo
            True
        """
        if not self._chcm == other._chcm:
            raise AssertionError("Combinations must live in the same complex")
        if self._chcm:
            cmpr = __cmpr_aux__(self._chcm._kenzo)
        else:
            cmpr = ecl_eval("#'s-cmpr")
        return Kenzocmbn(__2cmbn_add__(cmpr, self._kenzo, other._kenzo), self._chcm)

    def __sub__(self, other):
        r"""
        Return ``self - x``.

        EXAMPLES::

            sage: A = simplicial_sets.KleinBottle()
            sage: K = A._kenzo_repr                                              # optional - kenzo
            sage: c1 = K(1, [7, 'CELL_1_0', -3, 'CELL_1_2'])                     # optional - kenzo
            sage: c2 = K(1, [-3, 'CELL_1_0', -5, 'CELL_1_1'])                    # optional - kenzo
            sage: c1 - c2                                                        # optional - kenzo
            <BLANKLINE>
            {K... CMBN 1}
            --------------------------------------
            <10 * CELL_1_0>
            <5 * CELL_1_1>
            <-3 * CELL_1_2>
            --------------------------------------
            <BLANKLINE>
            sage: c2 - c1 == -c1 + c2                                            # optional - kenzo
            True
        """
        if not self._chcm == other._chcm:
            raise AssertionError("Combinations must live in the same complex")
        if self._chcm:
            cmpr = __cmpr_aux__(self._chcm._kenzo)
        else:
            cmpr = ecl_eval("#'s-cmpr")
        return Kenzocmbn(__2cmbn_sbtr__(cmpr, self._kenzo, other._kenzo), self._chcm)

    def __rmul__(self, x):
        r"""
        Return ``x * self`` if ``x`` is an integer number.

        EXAMPLES::

            sage: B = simplicial_sets.RealProjectiveSpace(4)
            sage: K = B._kenzo_repr                                              # optional - kenzo
            sage: c = K(2, [-4, 2]); c                                           # optional - kenzo
            <BLANKLINE>
            {K... CMBN 2}
            --------------------------------------
            <-4 * 2>
            --------------------------------------
            <BLANKLINE>
            sage: 3*c                                                            # optional - kenzo
            <BLANKLINE>
            {K... CMBN 2}
            --------------------------------------
            <-12 * 2>
            --------------------------------------
            <BLANKLINE>
            sage: 3*c == 7*c - c*4                                               # optional - kenzo
            True
        """
        assert(x.is_integer(), "Coefficient must be an integer number")
        if x == 0:
            return KenzoCombination(self._degree, Kchcm=self._chcm)
        return Kenzocmbn(__n_cmbn__(x, self._kenzo), self._chcm)

    def __mul__(self, x):
        r"""
        Return ``self * x`` if ``x`` is an integer number.

        EXAMPLES::

            sage: C = ChainComplex({-3: matrix(ZZ, 2, 2, [7, 0, 0, 0]), -2: matrix(ZZ, [[0, 1]])})
            sage: K = C.dual()._kenzo_repr                                       # optional - kenzo
            sage: c = K(-3, [1, 'G-3G0', -2, 'G-3G1']); c                        # optional - kenzo
            <BLANKLINE>
            {K... CMBN -3}
            --------------------------------------
            <1 * G-3G0>
            <-2 * G-3G1>
            --------------------------------------
            <BLANKLINE>
            sage: c * -7                                                         # optional - kenzo
            <BLANKLINE>
            {K... CMBN -3}
            --------------------------------------
            <-7 * G-3G0>
            <14 * G-3G1>
            --------------------------------------
            <BLANKLINE>
            sage: c*1 == c == 1*c                                                # optional - kenzo
            True
        """
        return self.__rmul__(x)

def Kenzocmbn(cmbn, Kchcm=None):
    r"""
    Auxiliary function for wrapping Kenzo combinations.

    INPUT:

    - ``cmbn`` -- An ECL object representing a Kenzo combination (type :cmbn).

    - ``Kchcm`` -- A KenzoChainComplex or None (default).

    OUTPUT:

    - The corresponding KenzoCombination instance of ``cmbn`` in ``Kchcm``.

    EXAMPLES::

        sage: from sage.interfaces.kenzo import __cmbn_aux__, Kenzocmbn          # optional - kenzo
        sage: from sage.libs.ecl import ecl_eval
        sage: C = ChainComplex({5: matrix(ZZ, 2, 2, [7, 0, 0, 0]), 6: matrix(ZZ, [[0, 1]])})
        sage: K = C.dual()._kenzo_repr                                           # optional - kenzo
        sage: c = __cmbn_aux__(5, [(4, 'G5G0'), (-7, 'G5G1')])                   # optional - kenzo
        sage: d = ecl_eval("(cmbn 5 4 'G5G0 -7 'G5G1)")                          # optional - kenzo
        sage: Kenzocmbn(c, K) == Kenzocmbn(d, K)                                 # optional - kenzo
        True
        sage: type(Kenzocmbn(c, K))                                              # optional - kenzo
        <class 'sage.interfaces.kenzo.KenzoCombination'>
    """
    cmbndegr = __cmbn_degr__(cmbn).python()
    #cmbnlist = __cmbn_list__(cmbn).python()
    cmbnlist = translate_type(__cmbn_list__(cmbn))
    return KenzoCombination(degree=cmbndegr, terms=unpairing(cmbnlist), Kchcm=Kchcm, kenzo=cmbn)


class KenzoChainComplex(KenzoObject):
    r"""
    Wrapper to Kenzo chain complexes. Kenzo simplicial sets are a particular case
    of Kenzo chain complexes.
    """

    def __call__(self, degree, terms=None):
        r"""
        INPUT:

        - ``degree`` -- An integer number.

        - ``terms`` -- A list of even length (default None). The even indexes
          determine the coefficients of the generators, which are determined by
          the odd indexes. For example, the ``terms`` parameter of a combination
          like 3A - 5B + C is the list [3, A, -5, B, 1, C].

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KenzoCombination             # optional - kenzo
            sage: D = ChainComplex({6: matrix(ZZ, 2, 2, [1,0,0,2])}, degree = -1)
            sage: K = D._kenzo_repr                                              # optional - kenzo
            sage: KenzoCombination(5, [-4, 'G5G0'], K) == K(5, [-4, 'G5G0'])     # optional - kenzo
            True
        """
        return KenzoCombination(degree=degree, terms=terms, Kchcm=self)

    def homology(self, n):
        r"""
        Return the ``n``'th homology group of the chain complex associated to this
        kenzo object.

        INPUT:

        - ``n`` -- the dimension in which compute the homology

        OUTPUT:

        - An homology group.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere   # optional - kenzo
            sage: s2 = Sphere(2)                             # optional - kenzo
            sage: s2                                         # optional - kenzo
            [K... Simplicial-Set]
            sage: s2.homology(2)                             # optional - kenzo
            Z
        """
        echcm1 = __echcm__(self._kenzo)
        m1 = __chcm_mat__(echcm1, n)
        m2 = __chcm_mat__(echcm1, n + 1)
        homology = __homologie__(m1, m2)
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

        - A :class:`KenzoChainComplex`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: s3 = Sphere(3)                              # optional - kenzo
            sage: p = s2.tensor_product(s3)                   # optional - kenzo
            sage: type(p)                                     # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoChainComplex'>
            sage: [p.homology(i) for i in range(8)]           # optional - kenzo
            [Z, 0, Z, Z, 0, Z, 0, 0]
        """
        return KenzoChainComplex(__tnsr_prdc__(self._kenzo, other._kenzo))

    def basis(self, dim):
        r"""
        Return the list of generators of the chain complex associated to the kenzo
        object ``self`` in dimension ``dim``.

        INPUT:

        - ``dim`` -- An integer number

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex   # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)         # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                                # optional - kenzo
            sage: kenzo_chcm                                                           # optional - kenzo
            [K... Chain-Complex]
            sage: for i in range(6):                                                   # optional - kenzo
            ....:     print("Basis in dimension %i: %s" % (i, kenzo_chcm.basis(i)))    # optional - kenzo
            Basis in dimension 0: ['G0G0', 'G0G1', 'G0G2']
            Basis in dimension 1: ['G1G0', 'G1G1']
            Basis in dimension 2: []
            Basis in dimension 3: ['G3G0', 'G3G1']
            Basis in dimension 4: ['G4G0', 'G4G1']
            Basis in dimension 5: ['G5G0', 'G5G1', 'G5G2']

        """
        return translate_type(__basis_aux__(self._kenzo, dim))

    def identity_morphism(self):
        r"""
        Return the identity morphism (degree 0) between ``self`` and itself.

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere                   # optional - kenzo
            sage: s2 = Sphere(2)                                             # optional - kenzo
            sage: tp = s2.tensor_product(s2)                                 # optional - kenzo
            sage: idnt = tp.identity_morphism()                              # optional - kenzo
            sage: type(idnt)                                                 # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoChainComplexMorphism'>
        """
        return KenzoChainComplexMorphism(__idnt_mrph__(self._kenzo))

    def null_morphism(self, target=None, degree=0):
        r"""
        Return the null morphism between the chain complexes ``self`` and ``target``
        of degree ``degree``.

        INPUT:

        - ``target`` -- A KenzoChainComplex or None (default).

        - ``degree`` -- An integer number (default 0).

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism` representing the null morphism
          between ``self`` and ``target`` of degree ``degree``. If ``target``
          takes None value, ``self`` is assumed as the target chain complex.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere                   # optional - kenzo
            sage: s2 = Sphere(2)                                             # optional - kenzo
            sage: s3 = Sphere(3)                                             # optional - kenzo
            sage: tp22 = s2.tensor_product(s2)                               # optional - kenzo
            sage: tp22                                                       # optional - kenzo
            [K... Chain-Complex]
            sage: tp23 = s2.tensor_product(s3)                               # optional - kenzo
            sage: tp23                                                       # optional - kenzo
            [K... Chain-Complex]
            sage: null1 = tp22.null_morphism()                               # optional - kenzo
            sage: null1                                                      # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null2 = tp22.null_morphism(target = tp23, degree = -3)     # optional - kenzo
            sage: null2                                                      # optional - kenzo
            [K... Morphism (degree -3): K... -> K...]
        """
        if target is None:
            target = self
        if not isinstance(target, KenzoChainComplex):
            raise ValueError("'target' parameter must be a KenzoChainComplex instance")
        elif (not degree == 0) and (not degree.is_integer()):
            raise ValueError("'degree' parameter must be an Integer number")
        else:
            return KenzoChainComplexMorphism(__zero_mrph__(self._kenzo, target._kenzo, degree))

    def differential(self, comb=None):
        r"""
        Return the differential of a combination or, if ``comb`` is not provided,
        return the differential :class:`KenzoMorphismChainComplex` of ``self``.

        INPUT:

        - ``comb`` -- A :class:`KenzoCombination` in ``self``.

        OUTPUT:

        - If ``comb`` is not None, it returns a KenzoCombination
          representing the differential of ``comb`` in ``self``. On the other
          hand, if ``comb`` take None value, the differential
          :class:`KenzoMorphismChainComplex` of ``self`` is returned.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                 # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)
            sage: kenzo_chcm = KChainComplex(sage_chcm)                           # optional - kenzo
            sage: kenzo_chcm                                                      # optional - kenzo
            [K... Chain-Complex]
            sage: kenzo_chcm.basis(4)                                             # optional - kenzo
            ['G4G0', 'G4G1']
            sage: c1 = kenzo_chcm(4, [1, 'G4G0'])                                 # optional - kenzo
            sage: kenzo_chcm.differential(c1)                                     # optional - kenzo
            <BLANKLINE>
            {K... CMBN 3}
            --------------------------------------
            <1 * G3G0>
            <3 * G3G1>
            --------------------------------------
            <BLANKLINE>
            sage: kenzo_chcm.basis(5)                                             # optional - kenzo
            ['G5G0', 'G5G1', 'G5G2']
            sage: c2 = kenzo_chcm(5, [1, 'G5G0', 2, 'G5G2'])                      # optional - kenzo
            sage: kenzo_chcm.differential(c2)                                     # optional - kenzo
            <BLANKLINE>
            {K... CMBN 4}
            --------------------------------------
            <6 * G4G0>
            <-3 * G4G1>
            --------------------------------------
            <BLANKLINE>
        """
        if not comb:
            return KenzoChainComplexMorphism(__dffr_aux__(self._kenzo))
        else:
            if not isinstance(comb, KenzoCombination) or not comb._chcm == self:
                raise AssertionError("The parameter must be a KenzoCombination in {}".format(self))
            return Kenzocmbn(__dffr_aux1__(self._kenzo, comb._kenzo), Kchcm=self)

    def orgn(self):
        r"""
        Return the :orgn slot of Kenzo, which stores as a list the origin of the object.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere                # optional - kenzo
            sage: s2 = Sphere(2)                                          # optional - kenzo
            sage: l2 = s2.loop_space()                                    # optional - kenzo
            sage: l2.orgn()                                               # optional - kenzo
            '(LOOP-SPACE [K... Simplicial-Set])'
            sage: A = l2.cartesian_product(s2)                            # optional - kenzo
            sage: A.orgn()                                                # optional - kenzo
            '(CRTS-PRDC [K... Simplicial-Group] [K... Simplicial-Set])'
        """
        return str(__orgn_aux__(self._kenzo))

    def __eq__(self, other):
        r"""
        Return ``True`` if and only if ``self == other`` as KenzoObjects.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import MooreSpace   # optional - kenzo
            sage: m24 = MooreSpace(2,4)                          # optional - kenzo
            sage: m24 == MooreSpace(2,4)                         # optional - kenzo
            True
            sage: m24 == MooreSpace(3,4)                         # optional - kenzo
            False
        """
        return IdNumber(self) == IdNumber(other)


class KenzoSimplicialSet(KenzoChainComplex):
    r"""
    Wrapper to Kenzo simplicial sets.

    In Kenzo, the homology of a simplicial set is computed from its associated
    chain complex. Hence, this class inherits from `KenzoChainComplex`.
    """
    def identity_simplicial_morphism(self):
        r"""
        Return the identity morphism of ``self``.

        EXAMPLES::

            sage:
            sage:
        """
        return KenzoSimplicialSetMorphism(__idnt_smmr__(self._kenzo))

    def loop_space(self, n=1):
        r"""
        Return the ``n`` th iterated loop space.

        INPUT:

        - ``n`` -- (default: 1) the number of times to iterate the loop space
          construction

        OUTPUT:

        - A :class:`KenzoSimplicialGroup`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: l2 = s2.loop_space()                        # optional - kenzo
            sage: type(l2)                                    # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialGroup'>
            sage: [l2.homology(i) for i in range(8)]          # optional - kenzo
            [Z, Z, Z, Z, Z, Z, Z, Z]
        """
        return KenzoSimplicialGroup(__loop_space__(self._kenzo, n))

    def cartesian_product(self, other):
        r"""
        Return the cartesian product of ``self`` and ``other``.

        INPUT:

        - ``other`` -- the Kenzo simplicial set with which the product is made

        OUTPUT:

        - A :class:`KenzoSimplicialSet`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: s3 = Sphere(3)                              # optional - kenzo
            sage: p = s2.cartesian_product(s3)                # optional - kenzo
            sage: type(p)                                     # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [p.homology(i) for i in range(6)]           # optional - kenzo
            [Z, 0, Z, Z, 0, Z]
            sage: p._factors                                  # optional - kenzo
            ([K... Simplicial-Set], [K... Simplicial-Set])
            sage: p._factors == (s2, s3)                         # optional - kenzo
            True
        """
        prod_kenzo = __crts_prdc__(self._kenzo, other._kenzo)
        result = KenzoSimplicialSet(prod_kenzo)
        setattr(result, '_factors', (self, other))
        return result

    def suspension(self):
        r"""
        Return the suspension of the simplicial set.

        OUTPUT:

        - A :class:`KenzoSimplicialSet`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import EilenbergMacLaneSpace    # optional - kenzo
            sage: e3 = EilenbergMacLaneSpace(ZZ, 3)                          # optional - kenzo
            sage: s = e3.suspension()                                        # optional - kenzo
            sage: type(s)                                                    # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [s.homology(i) for i in range(6)]                          # optional - kenzo
            [Z, 0, 0, 0, Z, 0]
        """
        result = KenzoSimplicialSet(__suspension__(self._kenzo))
        setattr(result, '_base', self)
        return result

    def homotopy_group(self, n):
        """
        Return the n'th homotopy group of ``self``

        INPUT:

        - ``n`` -- the dimension of the homotopy group to be computed

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere      # optional - kenzo
            sage: s2 = Sphere(2)                                # optional - kenzo
            sage: p = s2.cartesian_product(s2)                  # optional - kenzo
            sage: p.homotopy_group(3)                           # optional - kenzo
            Multiplicative Abelian group isomorphic to Z x Z


        .. WARNING::

            This method assumes that the underlying space is simply connected.
            You might get wrong answers if it is not.
        """
        if n not in ZZ or n < 2:
            raise ValueError("""homotopy groups can only be computed
                for dimensions greater than 1""")
        lgens = __homotopy_list__(self._kenzo, n).python()
        if lgens is not None:
            trgens = [0 if i == 1 else i for i in sorted(lgens)]
            return AbelianGroup(trgens)
        else:
            return AbelianGroup([])

    def em_spectral_sequence(self):
        r"""
        Return the Eilenberg-Moore spectral sequence of ``self``.

        OUTPUT:

        - A :class:`KenzoSpectralSequence`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere # optional - kenzo
            sage: S2 = Sphere(2)                           # optional - kenzo
            sage: EMS = S2.em_spectral_sequence()          # optional - kenzo
            sage: EMS.table(0, -2, 2, -2, 2)               # optional - kenzo
              0   Z   0   0   0
              0   0   0   0   0
              0   0   Z   0   0
              0   0   0   0   0
              0   0   0   0   0


        .. WARNING::

            This method assumes that the underlying space is simply connected.
            You might get wrong answers if it is not.
        """
        if self.homology(1).invariants():
            raise ValueError("""Eilenberg-Moore spectral sequence implemented
                only for 1-reduced simplicial sets""")
        return KenzoSpectralSequence(__eilenberg_moore_spectral_sequence__(self._kenzo))

    def sw_spectral_sequence(self):
        r"""
        Return the Serre sequence of the first step of the Whitehead tower.

        OUTPUT:

        - A :class:`KenzoSpectralSequence`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere  # optional - kenzo
            sage: S3 = Sphere(3)                            # optional - kenzo
            sage: E = S3.sw_spectral_sequence()             # optional - kenzo
            sage: T = E.table(0, 0, 4, 0, 4)                # optional - kenzo
            sage: T                                         # optional - kenzo
              Z   0   0   Z   0
              0   0   0   0   0
              Z   0   0   Z   0
              0   0   0   0   0
              Z   0   0   Z   0
        """
        if self.homology(1).invariants():
            raise ValueError("""Eilenberg-Moore spectral sequence implemented
                only for 1-reduced simplicial sets""")
        return KenzoSpectralSequence(__serre_whitehead_spectral_sequence__(self._kenzo))

    def serre_spectral_sequence(self):
        r"""
        Return the spectral sequence of ``self``.

        The object self must be created as a cartesian product (twisted or not).

        OUTPUT:

        - A :class:`KenzoSpectralSequence`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere
            sage: S2 = Sphere(2)                            # optional - kenzo
            sage: S3 = Sphere(3)                            # optional - kenzo
            sage: P = S2.cartesian_product(S3)              # optional - kenzo
            sage: E = P.serre_spectral_sequence()           # optional - kenzo
            sage: E.table(0, 0, 2, 0, 3)                    # optional - kenzo
              Z   0   Z
              0   0   0
              0   0   0
              Z   0   Z

        .. WARNING::

            This method assumes that the underlying space is simply connected.
            You might get wrong answers if it is not.
        """
        if self.homology(1).invariants():
            raise ValueError("""Eilenberg-Moore spectral sequence implemented
                only for 1-reduced simplicial sets""")
        return KenzoSpectralSequence(__serre_spectral_sequence_product__(self._kenzo))

    def wedge(self, other):
        r"""
        Return the wedge of ``self`` and ``other``.

        INPUT:

        - ``other`` -- the Kenzo simplicial set with which the wedge is made

        OUTPUT:

        - A :class:`KenzoSimplicialSet`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: s3 = Sphere(3)                              # optional - kenzo
            sage: w = s2.wedge(s3)                            # optional - kenzo
            sage: type(w)                                     # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [w.homology(i) for i in range(6)]           # optional - kenzo
            [Z, 0, Z, Z, 0, 0]
            sage: w._factors == (s2, s3)                      # optional - kenzo
            True
        """
        wedge_kenzo = __wedge__(self._kenzo, other._kenzo)
        result = KenzoSimplicialSet(wedge_kenzo)
        setattr(result, '_factors', (self, other))
        return result

    def join(self, other):
        r"""
        Return the join of ``self`` and ``other``.

        INPUT:

        - ``other`` -- the Kenzo simplicial set with which the join is made

        OUTPUT:

        - A :class:`KenzoSimplicialSet`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: s3 = Sphere(3)                              # optional - kenzo
            sage: j = s2.join(s3)                             # optional - kenzo
            sage: type(j)                                     # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [j.homology(i) for i in range(6)]           # optional - kenzo
            [Z, 0, 0, 0, 0, 0]
        """
        join_kenzo = __join__(self._kenzo, other._kenzo)
        return KenzoSimplicialSet(join_kenzo)

    def smash_product(self, other):
        r"""
        Return the smash product of ``self`` and ``other``.

        INPUT:

        - ``other`` -- the Kenzo simplicial set with which the smash product is made

        OUTPUT:

        - A :class:`KenzoSimplicialSet`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: s3 = Sphere(3)                              # optional - kenzo
            sage: s = s2.smash_product(s3)                    # optional - kenzo
            sage: type(s)                                     # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [s.homology(i) for i in range(6)]           # optional - kenzo
            [Z, 0, 0, 0, 0, Z]
            sage: s._factors == (s2, s3)                      # optional - kenzo
            True
        """
        smash_kenzo = __smash_product__(self._kenzo, other._kenzo)
        result = KenzoSimplicialSet(smash_kenzo)
        setattr(result, '_factors', (self, other))
        return result


class KenzoSimplicialGroup(KenzoSimplicialSet):
    r"""
    Wrapper around Kenzo simplicial groups.
    """

    def classifying_space(self):
        r"""
        Return the classifying space.

        OUTPUT:

        - A :class:`KenzoSimplicialGroup`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import MooreSpace      # optional - kenzo
            sage: m2 = MooreSpace(2,4)                              # optional - kenzo
            sage: l2 = m2.loop_space()                              # optional - kenzo
            sage: c = l2.classifying_space()                        # optional - kenzo
            sage: type(c)                                           # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialGroup'>
            sage: [c.homology(i) for i in range(8)]                 # optional - kenzo
            [Z, 0, 0, 0, C2, 0, 0, 0]
        """
        return KenzoSimplicialGroup(__classifying_space__(self._kenzo))


def k2s_matrix(kmatrix):
    r"""
    Convert an array of ECL to a matrix of Sage.

    INPUT:

    - ``kmatrix`` -- An array in ECL

    EXAMPLES::

        sage: from sage.interfaces.kenzo import k2s_matrix         # optional - kenzo
        sage: from sage.libs.ecl import EclObject
        sage: M = EclObject("#2A((1 2 3) (3 2 1) (1 1 1))")
        [1 2 3]
        [3 2 1]
        [1 1 1]
    """
    dimensions = __array_dimensions__(kmatrix).python()
    kmatrix_list = __make_array_to_lists__(kmatrix).python()
    return matrix(dimensions[0], dimensions[1], kmatrix_list)


def s2k_matrix(smatrix):
    r"""
    Convert a matrix of Sage to an array of ECL.

    INPUT:

    - ``smatrix`` -- A matrix in Sage

    OUTPUT:

    - A :class:`EclObject`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import s2k_matrix      # optional - kenzo
        sage: A = Matrix([[1,2,3],[3,2,1],[1,1,1]])
        sage: s2k_matrix(A)                                     # optional - kenzo
        <ECL: #2A((1 2 3) (3 2 1) (1 1 1))>
    """
    initcontents = []
    dimensions = smatrix.dimensions()
    for i in smatrix.rows():
        initcontents.append(i.list())
    return __make_array_from_lists__(dimensions[0], dimensions[1], initcontents)


def s2k_dictmat(sdictmat):
    r"""
    Convert a dictionary in Sage, whose values are matrices, to an assoc list
    in ECL.

    INPUT:

    - ``sdictmat`` -- A dictionary in Sage

    OUTPUT:

    - A :class:`EclObject`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import s2k_dictmat   # optional - kenzo
        sage: A = Matrix([[1,2,3],[3,2,1],[1,1,1]])
        sage: B = Matrix([[1,2],[2,1],[1,1]])
        sage: d = {1 : A, 2 : B}
        sage: s2k_dictmat(d)                                  # optional - kenzo
        <ECL: ((2 . #2A((1 2) (2 1) (1 1))) (1 . #2A((1 2 3) (3 2 1) (1 1 1))))>
    """
    rslt = EclObject([])
    for k in sdictmat.keys():
        rslt = EclObject(k).cons(s2k_matrix(sdictmat[k])).cons(rslt)
    return rslt

def pairing(slist):
    r"""
    Convert a list of Sage (which has an even length) to an assoc list in ECL.

    INPUT:

    - ``slist`` -- A list in Sage.

    OUTPUT:

    - An :class:`EclObject`.

    EXAMPLES::

        sage: from sage.interfaces.kenzo import pairing   # optional - kenzo
        sage: l = [1, 2, 3, 4]
        sage: pairing(l)                                  # optional - kenzo
        <ECL: ((1 . 2) (3 . 4))>
    """
    l = len(slist)
    if l % 2 == 1:
        raise ValueError("The parameter must be a list of even length")
    rslt = EclObject([])
    for k in range(l - 1, 0, -2):
        if hasattr(slist[k], '_kenzo'):
            even = slist[k]._kenzo
        else:
            even = EclObject(slist[k])
        rslt = EclObject(slist[k - 1]).cons(even).cons(rslt)
    return rslt

def unpairing(klist):
    r"""
    Convert a list of tuples in a plain list.

    EXAMPLES::

        sage: from sage.interfaces.kenzo import unpairing, __cmbn_list__         # optional - kenzo
        sage: from sage.libs.ecl import ecl_eval                                 # optional - kenzo
        sage: c = ecl_eval("(cmbn 4 -2 'A 3 'B)")                                # optional - kenzo
        sage: m = __cmbn_list__(c).python(); m                                   # optional - kenzo
        [(-2, 'A'), (3, 'B')]
        sage: unpairing(m)                                                       # optional - kenzo
        [-2, 'A', 3, 'B']
    """
    if klist==None:
        return None
    result = []
    for x in klist:
        result += list(x)
    return result

def KChainComplex(chain_complex):
    r"""
    Construct a KenzoChainComplex from a ChainComplex of degree = -1 in Sage.

    INPUT:

    - ``chain_complex`` -- A ChainComplex of degree = -1

    OUTPUT:

    - A :class:`KenzoChainComplex`.

    EXAMPLES::

        sage: from sage.interfaces.kenzo import KChainComplex                 # optional - kenzo
        sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
        sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
        sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
        sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)    # optional - kenzo
        sage: kenzo_chcm = KChainComplex(sage_chcm)                           # optional - kenzo
        sage: kenzo_chcm                                                      # optional - kenzo
        [K... Chain-Complex]
        sage: kenzo_chcm.homology(5)                                          # optional - kenzo
        Z x Z
    """
    if chain_complex.degree_of_differential() != -1:
        raise AssertionError('Chain complexes in Kenzo has degree -1.')
    diff = chain_complex.differential()
    chcm = s2k_dictmat(diff)
    str_orgn = 'sage.' + str(diff)[1:-1].translate({ord(i): None for i in ["\n", ",", ":"]}) \
               .replace("[ ", "[").replace("  ", " ").replace(" ", ".")
    return KenzoChainComplex(__kchaincomplex_aux__(chcm, str_orgn))


def SChainComplex(kchaincomplex, start=0, end=15):
    r"""
    Convert the KenzoChainComplex ``kchcm`` (between dimensions ``start`` and
    ``end``) to a ChainComplex.

    INPUT:

    - ``kchaincomplex`` -- A KenzoChainComplex

    - ``start`` -- An integer number (optional, default 0)

    - ``end`` -- An integer number greater than or equal to ``start`` (optional, default 15)

    OUTPUT:

    - A ChainComplex

    EXAMPLES::

        sage: from sage.interfaces.kenzo import KChainComplex, SChainComplex   # optional - kenzo
        sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
        sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
        sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
        sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)     # optional - kenzo
        sage: SChainComplex(KChainComplex(sage_chcm)) == sage_chcm             # optional - kenzo
        True

    ::

        sage: from sage.interfaces.kenzo import SChainComplex, Sphere     # optional - kenzo
        sage: S4 = Sphere(4)                       # optional - kenzo
        sage: C = SChainComplex(S4)                # optional - kenzo
        sage: C                                    # optional - kenzo
        Chain complex with at most 3 nonzero terms over Integer Ring
        sage: C._ascii_art_()                      # optional - kenzo
        0 <-- C_4 <-- 0  ...  0 <-- C_0 <-- 0
        sage: [C.homology(i) for i in range(6)]    # optional - kenzo
        [Z, 0, 0, 0, Z, 0]
    """
    from sage.homology.chain_complex import ChainComplex
    matrices = {}
    for i in range(start, end):
        dffr_i = __chcm_mat2__(kchaincomplex._kenzo, i)
        nlig = __nlig__(dffr_i).python()
        ncol = __ncol__(dffr_i).python()
        if ((nlig != 0) and (ncol != 0)):
            matrices[i] = k2s_matrix(__convertmatrice__(dffr_i))
        else:
            matrices[i] = matrix(nlig, ncol)
    return ChainComplex(matrices, degree=-1)


def SAbstractSimplex(simplex, dim):
    r"""
    Convert an abstract simplex of Kenzo to an AbstractSimplex.

    INPUT:

    - ``simplex`` -- An abstract simplex of Kenzo.

    - ``dim``-- The dimension of ``simplex``.

    OUTPUT:

    - An AbstractSimplex.

    EXAMPLES::

        sage: from sage.libs.ecl import EclObject, ecl_eval
        sage: from sage.interfaces.kenzo import KenzoObject,\
        ....: SAbstractSimplex                                  # optional - kenzo
        sage: KAbSm = KenzoObject(ecl_eval("(ABSM 15 'K)"))     # optional - kenzo
        sage: SAbSm1 = SAbstractSimplex(KAbSm, 2)               # optional - kenzo
        sage: SAbSm2 = SAbstractSimplex(KAbSm, 7)               # optional - kenzo
        sage: SAbSm1.degeneracies()                             # optional - kenzo
        [3, 2, 1, 0]
        sage: SAbSm1.dimension()                                # optional - kenzo
        6
        sage: SAbSm2.dimension()                                # optional - kenzo
        11
    """
    from sage.homology.simplicial_set import AbstractSimplex
    degeneracies = __dgop_int_ext__(__dgop__(simplex._kenzo)).python()
    if degeneracies is None:
        degeneracies = []
    else:
        degeneracies = tuple(degeneracies)
    name = __gmsm__(simplex._kenzo).python()
    return AbstractSimplex(dim, degeneracies, name=name)


def KAbstractSimplex(simplex):
    r"""
    Convert an AbstractSimplex in Sage to an abstract simplex of Kenzo.

    INPUT:

    - ``simplex`` -- An AbstractSimplex.

    OUTPUT:

    - An abstract simplex of Kenzo.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import AbstractSimplex
        sage: from sage.interfaces.kenzo import KAbstractSimplex,\
        ....: SAbstractSimplex                                          # optional - kenzo
        sage: SAbSm = AbstractSimplex(1, (2,0,3,2,1), name = 'SAbSm')   # optional - kenzo
        sage: KAbSm = KAbstractSimplex(SAbSm)                           # optional - kenzo
        sage: SAbSm2 = SAbstractSimplex(KAbSm, 1)                       # optional - kenzo
        sage: SAbSm.degeneracies() == SAbSm2.degeneracies()             # optional - kenzo
        True
        sage: SAbSm.dimension() == SAbSm2.dimension()                   # optional - kenzo
        True
    """
    nondegenerate = simplex.nondegenerate()
    result = KenzoAbstractSimplex(__kabstractsimplex_aux__(simplex.degeneracies(),
                                                's' + str(hash(nondegenerate))))
    setattr(result, '_dim', simplex.dimension())
    return result


def KFiniteSimplicialSet(sset):
    r"""
    Convert a finite SimplicialSet in Sage to a finite simplicial set of Kenzo.

    INPUT:

    - ``sset`` -- A finite SimplicialSet.

    OUTPUT:

    - A finite simplicial set of Kenzo.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
        sage: from sage.interfaces.kenzo import KFiniteSimplicialSet    # optional - kenzo
        sage: s0 = AbstractSimplex(0, name='s0')
        sage: s1 = AbstractSimplex(0, name='s1')
        sage: s2 = AbstractSimplex(0, name='s2')
        sage: s01 = AbstractSimplex(1, name='s01')
        sage: s02 = AbstractSimplex(1, name='s02')
        sage: s12 = AbstractSimplex(1, name='s12')
        sage: s012 = AbstractSimplex(2, name='s012')
        sage: Triangle = SimplicialSet({s01: (s1, s0),\
        ....: s02: (s2, s0), s12: (s2, s1)}, base_point = s0)
        sage: KTriangle = KFiniteSimplicialSet(Triangle)                # optional - kenzo
        sage: [KTriangle.homology(i) for i in range(3)]                 # optional - kenzo
        [Z, Z, 0]
        sage: KTriangle.basis(1)                                        # optional - kenzo
        ['CELL_1_0', 'CELL_1_1', 'CELL_1_2']
        sage: S1 = simplicial_sets.Sphere(1)
        sage: S3 = simplicial_sets.Sphere(3)
        sage: KS1vS3 = KFiniteSimplicialSet(S1.wedge(S3))               # optional - kenzo
        sage: [KS1vS3.homology(i) for i in range(4)]                    # optional - kenzo
        [Z, Z, 0, Z]
    """
    if hasattr(sset, '_kenzo_repr'):
        return sset._kenzo_repr
    else:
        allcells = sset.cells()
        namecells = {c: 'cell_{}_{}'.format(d, allcells[d].index(c))
                     for d in allcells for c in allcells[d]}
        dim = sset.dimension()
        list_rslt = [namecells[i] for i in sset.n_cells(0)]
        if (dim > 0):
            for k in range(1, dim + 1):
                k_cells = sset.n_cells(k)
                if k_cells:
                    list_rslt.append(k)
                    for x in k_cells:
                        list_rslt.append(namecells[x])
                        auxiliar_list = []
                        for z in sset.faces(x):
                            degen_z = z.degeneracies()
                            name = namecells[z.nondegenerate()]
                            degen_z.append(name)
                            auxiliar_list.append(degen_z)
                        list_rslt.append(auxiliar_list)
        return KenzoSimplicialSet(__build-finite-ss2__(list_rslt))


def SFiniteSimplicialSet(ksimpset, limit):
    r"""
    Convert the ``limit``-skeleton of a finite simplicial set in Kenzo to a
    finite SimplicialSet in Sage.

    INPUT:

    - ``ksimpset`` -- A finite simplicial set in Kenzo.

    - ``limit`` -- A natural number.

    OUTPUT:

    - A finite SimplicialSet.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import SimplicialSet
        sage: from sage.interfaces.kenzo import AbstractSimplex,\
        ....:  KFiniteSimplicialSet, SFiniteSimplicialSet, Sphere   # optional - kenzo
        sage: s0 = AbstractSimplex(0, name='s0')                    # optional - kenzo
        sage: s1 = AbstractSimplex(0, name='s1')                    # optional - kenzo
        sage: s2 = AbstractSimplex(0, name='s2')                    # optional - kenzo
        sage: s01 = AbstractSimplex(1, name='s01')                  # optional - kenzo
        sage: s02 = AbstractSimplex(1, name='s02')                  # optional - kenzo
        sage: s12 = AbstractSimplex(1, name='s12')                  # optional - kenzo
        sage: s012 = AbstractSimplex(2, name='s012')                # optional - kenzo
        sage: Triangle = SimplicialSet({s01: (s1, s0),\
        ....: s02: (s2, s0), s12: (s2, s1)}, base_point = s0)       # optional - kenzo
        sage: KTriangle = KFiniteSimplicialSet(Triangle)            # optional - kenzo
        sage: STriangle = SFiniteSimplicialSet(KTriangle, 1)        # optional - kenzo
        sage: STriangle.homology()                                  # optional - kenzo
        {0: 0, 1: Z}
        sage: S1 = simplicial_sets.Sphere(1)                        # optional - kenzo
        sage: S3 = simplicial_sets.Sphere(3)                        # optional - kenzo
        sage: KS1vS3 = KFiniteSimplicialSet(S1.wedge(S3))           # optional - kenzo
        sage: SS1vS3 = SFiniteSimplicialSet(KS1vS3, 3)              # optional - kenzo
        sage: SS1vS3.homology()                                     # optional - kenzo
        {0: 0, 1: Z, 2: 0, 3: Z}
    """
    from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
    list_orgn = __orgn_aux__(ksimpset._kenzo).python()
    if __nth__(0, list_orgn).python()[0] == 'CRTS-PRDC':
        return SFiniteSimplicialSet(
            KenzoSimplicialSet(__nth__(1, list_orgn)), limit).cartesian_product(
                SFiniteSimplicialSet(KenzoSimplicialSet(__nth__(2, list_orgn)), limit))
    rslt = {}
    simplices = []
    faces = []
    bases = []
    names = []
    for k in range(limit + 1):
        basis_k = __basis_aux__(ksimpset._kenzo, k)
        names_k = ksimpset.basis(k)
        lbasis_k = [AbstractSimplex(k, name=i) for i in EclListIterator(basis_k)]
        bases.append(lbasis_k)
        names.append(names_k)
    all_simplices = __sfinitesimplicialset_aux__(ksimpset._kenzo, limit)
    lall_simplices = [i for i in EclListIterator(all_simplices)]
    dim = 1
    for Kdim in lall_simplices:
        for simp in Kdim:
            index1 = names[dim].index(str(simp.car()))
            lKdim_cdr = []
            for i in EclListIterator(simp.cdr()):
                degenop = __dgop_int_ext__(__dgop__(i)).python()
                if degenop is None:
                    degenop = []
                index2 = names[dim - len(degenop) - 1].index(str(__gmsm__(i)))
                lKdim_cdr.append(bases[dim - len(degenop) - 1][index2].apply_degeneracies(*degenop))
            simplices.append(bases[dim][index1])
            faces.append(tuple(lKdim_cdr))
        dim += 1
    for i in range(len(simplices)):
        rslt[simplices[i]] = faces[i]
    return SimplicialSet(rslt)


class KenzoChainComplexMorphism(KenzoObject):
    r"""
    Wrapper to Kenzo morphisms between chain complexes.
    """

    def __call__(self, kcmbn):
        r"""
        Returns the image of ``kcmbn`` under the morphism.

        INPUT:

        - ``kcmbn`` -- A :class:`KenzoCombination` in the source complex of the
          morphism.

        OUTPUT:

        - A :class:`KenzoCombination` in the target complex of the morphism.

        EXAMPLES::

            sage: C = ChainComplex({5: matrix(ZZ, 2, 2, [7, 0, 0, 0]), 6: matrix(ZZ, [[0, 1]])})
            sage: K = C.dual()._kenzo_repr                                       # optional - kenzo
            sage: c = K(5, [4, 'G5G0', 6, 'G5G1'])                               # optional - kenzo
            sage: one = K.identity_morphism()                                    # optional - kenzo
            sage: one(c)                                                         # optional - kenzo
            <BLANKLINE>
            {K... CMBN 5}
            --------------------------------------
            <4 * G5G0>
            <6 * G5G1>
            --------------------------------------
            <BLANKLINE>
            sage: D = ChainComplex({8: matrix(ZZ, 2, 2, [3, 0, 0, 0]), 9: matrix(ZZ, [[0, 1]])})
            sage: L = D.dual()._kenzo_repr                                       # optional - kenzo
            sage: zero = K.null_morphism(L)                                      # optional - kenzo
            sage: zero(c)                                                        # optional - kenzo
            <BLANKLINE>
            {K... CMBN 5}
            --------------------------------------
            <BLANKLINE>
            --------------------------------------
            <BLANKLINE>
            sage: zero(c)._chcm == L                                             # optional - kenzo
            True
            sage: zero(c)._chcm == K                                             # optional - kenzo
            False
        """
        if not isinstance(kcmbn, KenzoCombination):
            raise TypeError("The parameter must be a KenzoCombination instance")
        if kcmbn._chcm != self.source_complex():
            raise AssertionError("{} must be a combination in {}") \
                                  .format(kcmbn, self.source_complex())
        image = __evaluation_aux__(self._kenzo, kcmbn._kenzo)
        return Kenzocmbn(image, self.target_complex())

    def source_complex(self):
        r"""
        Return the source chain complex of the morphism.

        OUTPUT:

        - A :class:`KenzoChainComplex`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)   # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                          # optional - kenzo
            sage: kenzo_chcm                                                     # optional - kenzo
            [K... Chain-Complex]
            sage: differential_morphism = kenzo_chcm.differential()              # optional - kenzo
            sage: differential_morphism                                          # optional - kenzo
            [K... Morphism (degree -1): K... -> K...]
            sage: differential_morphism.source_complex() == kenzo_chcm           # optional - kenzo
            True
        """
        return KenzoChainComplex(__sorc_aux__(self._kenzo))

    def target_complex(self):
        r"""
        Return the target chain complex of the morphism.

        OUTPUT:

        - A :class:`KenzoChainComplex`.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)   # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                          # optional - kenzo
            sage: kenzo_chcm                                                     # optional - kenzo
            [K... Chain-Complex]
            sage: differential_morphism = kenzo_chcm.differential()              # optional - kenzo
            sage: differential_morphism                                          # optional - kenzo
            [K... Morphism (degree -1): K... -> K...]
            sage: differential_morphism.target_complex() == kenzo_chcm           # optional - kenzo
            True
        """
        return KenzoChainComplex(__trgt_aux__(self._kenzo))

    def degree(self):
        r"""
        Return the degree of the morphism.

        OUTPUT:

        - An integer number, the degree of the morphism.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)   # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                          # optional - kenzo
            sage: kenzo_chcm                                                     # optional - kenzo
            [K... Chain-Complex]
            sage: differential_morphism = kenzo_chcm.differential()              # optional - kenzo
            sage: differential_morphism                                          # optional - kenzo
            [K... Morphism (degree -1): K... -> K...]
            sage: differential_morphism.degree()                                 # optional - kenzo
            -1
            sage: differential_morphism.composite(differential_morphism).degree() # optional - kenzo
            -2
            sage: kenzo_chcm.null_morphism().degree()                            # optional - kenzo
            0
        """
        return __degr_aux__(self._kenzo).python()

    def evaluation(self, comb):
        r"""
        Apply the morphism on a Kenzocombination ``comb``.

        INPUT:

        - ``comb`` -- A :class:`KenzoCombination`.

        OUTPUT:

        - A KenzoCombination representing the result of applying the morphism on
            ``comb``.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                 # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)
            sage: kenzo_chcm = KChainComplex(sage_chcm)                           # optional - kenzo
            sage: kenzo_chcm                                                      # optional - kenzo
            [K... Chain-Complex]
            sage: differential_morphism = kenzo_chcm.differential()               # optional - kenzo
            sage: differential_morphism                                           # optional - kenzo
            [K... Morphism (degree -1): K... -> K...]
            sage: dif_squared = differential_morphism.composite(differential_morphism)  # optional - kenzo
            sage: dif_squared                                                     # optional - kenzo
            [K... Morphism (degree -2): K... -> K...]
            sage: kenzo_chcm.basis(5)                                             # optional - kenzo
            ['G5G0', 'G5G1', 'G5G2']
            sage: c1 = kenzo_chcm(5, [1, 'G5G0', 2, 'G5G2'])                      # optional - kenzo
            sage: kenzo_chcm.differential(c1)                                     # optional - kenzo
            <BLANKLINE>
            {K... CMBN 4}
            --------------------------------------
            <6 * G4G0>
            <-3 * G4G1>
            --------------------------------------
            <BLANKLINE>
            sage: differential_morphism.evaluation(c1)                            # optional - kenzo
            <BLANKLINE>
            {K... CMBN 4}
            --------------------------------------
            <6 * G4G0>
            <-3 * G4G1>
            --------------------------------------
            <BLANKLINE>
            sage: dif_squared.evaluation(c1)                                      # optional - kenzo
            <BLANKLINE>
            {K... CMBN 3}
            --------------------------------------
            <BLANKLINE>
            --------------------------------------
            <BLANKLINE>
            sage: idnt = kenzo_chcm.identity_morphism()                           # optional - kenzo
            sage: idx2 = idnt.sum(idnt)                                           # optional - kenzo
            sage: idnt.evaluation(c1)                                             # optional - kenzo
            <BLANKLINE>
            {K... CMBN 5}
            --------------------------------------
            <1 * G5G0>
            <2 * G5G2>
            --------------------------------------
            <BLANKLINE>
            sage: idx2.evaluation(c1)                                             # optional - kenzo
            <BLANKLINE>
            {K... CMBN 5}
            --------------------------------------
            <2 * G5G0>
            <4 * G5G2>
            --------------------------------------
            <BLANKLINE>
        """
        if not isinstance(comb, KenzoCombination):
            raise ValueError("The parameter must be a KenzoCombination")
        return Kenzocmbn(__evaluation_aux__(self._kenzo, comb._kenzo), self.source_complex())

    def opposite(self):
        r"""
        Return the opposite morphism of ``self``, i.e., -1 x ``self``.

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)   # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                          # optional - kenzo
            sage: kenzo_chcm                                                     # optional - kenzo
            [K... Chain-Complex]
            sage: idnt = kenzo_chcm.identity_morphism()                          # optional - kenzo
            sage: idnt                                                           # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: opps_id = idnt.opposite()                                      # optional - kenzo
            sage: opps_id                                                        # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: kenzo_chcm.basis(4)                                            # optional - kenzo
            ['G4G0', 'G4G1']
            sage: c1 = kenzo_chcm(4, [2, 'G4G0', -5, 'G4G1'])                    # optional - kenzo
            sage: idnt.evaluation(c1)                                            # optional - kenzo
            <BLANKLINE>
            {K... CMBN 4}
            --------------------------------------
            <2 * G4G0>
            <-5 * G4G1>
            --------------------------------------
            <BLANKLINE>
            sage: opps_id.evaluation(c1)                                         # optional - kenzo
            <BLANKLINE>
            {K... CMBN 4}
            --------------------------------------
            <-2 * G4G0>
            <5 * G4G1>
            --------------------------------------
            <BLANKLINE>
        """
        return KenzoChainComplexMorphism(__opps__(self._kenzo))

    def composite(self, *morphisms):
        r"""
        Return the composite of ``self`` and the morphism(s) given by the parameter
        ``morphisms``.

        INPUT:

        - ``morphisms`` -- A KenzoChainComplexMorphism instance, a KenzoChainComplex
          instance or a tuple of KenzoChainComplexMorphism and KenzoChainComplex
          instances.

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism`: if ``morphisms`` is a KenzoChainComplexMorphism,
          the composite of ``self`` and ``morphisms`` is returned. If ``morphisms``
          is a KenzoChainComplex, the composite of ``self`` and the differential
          morphism of ``morphisms`` is returned. If ``morphisms`` is a tuple,
          the composite of ``self`` and the morphisms or the differential morphisms
          of the given chain complexes in ``morphisms`` is returned (if ``morphisms``
          is None, ``self`` morphism is returned).

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere                   # optional - kenzo
            sage: s2 = Sphere(2)                                             # optional - kenzo
            sage: s3 = Sphere(3)                                             # optional - kenzo
            sage: tp22 = s2.tensor_product(s2)                               # optional - kenzo
            sage: tp23 = s2.tensor_product(s3)                               # optional - kenzo
            sage: idnt = tp22.identity_morphism()                            # optional - kenzo
            sage: idnt                                                       # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null = tp23.null_morphism(target = tp22, degree = 4)       # optional - kenzo
            sage: null                                                       # optional - kenzo
            [K... Morphism (degree 4): K... -> K...]
            sage: idnt.composite(tp22, null)                                 # optional - kenzo
            [K... Morphism (degree 3): K... -> K...]
        """
        if morphisms is None:
            return self
        if isinstance(morphisms, KenzoChainComplexMorphism):
            return KenzoChainComplexMorphism(__cmps__(self._kenzo, morphisms._kenzo))
        elif isinstance(morphisms, KenzoChainComplex):
            return KenzoChainComplexMorphism(__cmps__(self._kenzo, __dffr_aux__(morphisms._kenzo)))
        elif isinstance(morphisms, tuple):
            rslt = self._kenzo
            for mrph in morphisms:
                rslt = __cmps__(rslt, mrph._kenzo)
            return KenzoChainComplexMorphism(rslt)

    def sum(self, *morphisms):
        r"""
        Return a morphism, sum of the morphism ``self`` and the morphism(s) given
        by the parameter ``morphisms``.

        INPUT:

        - ``morphisms`` -- A KenzoChainComplexMorphism instance or a tuple of
          KenzoChainComplexMorphism instances.

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism`, sum of the morphism ``self`` and
          the morphism(s) given by ``morphisms`` (if ``morphisms`` is None, ``self``
          morphism is returned).

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)   # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                          # optional - kenzo
            sage: kenzo_chcm                                                     # optional - kenzo
            [K... Chain-Complex]
            sage: idnt = kenzo_chcm.identity_morphism()                          # optional - kenzo
            sage: idnt                                                           # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: opps_id = idnt.opposite()                                      # optional - kenzo
            sage: opps_id                                                        # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null = kenzo_chcm.null_morphism()                              # optional - kenzo
            sage: null                                                           # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: idx2 = idnt.sum(idnt)                                          # optional - kenzo
            sage: idx5 = idx2.sum(\
            ....: opps_id, idnt, idnt, null, idx2.sum(idnt), opps_id)            # optional - kenzo
            sage: kenzo_chcm.basis(4)                                            # optional - kenzo
            ['G4G0', 'G4G1']
            sage: c1 = kenzo_chcm(4, [2, 'G4G0', -5, 'G4G1'])                    # optional - kenzo
            sage: idx2.evaluation(c1)                                            # optional - kenzo
            <BLANKLINE>
            {K... CMBN 4}
            --------------------------------------
            <4 * G4G0>
            <-10 * G4G1>
            --------------------------------------
            <BLANKLINE>
            sage: idx5.evaluation(c1)                                            # optional - kenzo
            <BLANKLINE>
            {K... CMBN 4}
            --------------------------------------
            <10 * G4G0>
            <-25 * G4G1>
            --------------------------------------
            <BLANKLINE>
        """
        if morphisms is None:
            return self
        if isinstance(morphisms, KenzoChainComplexMorphism):
            return KenzoChainComplexMorphism(__add__(self._kenzo, morphisms._kenzo))
        elif isinstance(morphisms, tuple):
            rslt = self._kenzo
            for mrph in morphisms:
                rslt = __add__(rslt, mrph._kenzo)
            return KenzoChainComplexMorphism(rslt)

    def substract(self, *morphisms):
        r"""
        Return a morphism, difference of the morphism ``self`` and the morphism(s)
        given by the parameter ``morphisms``.

        INPUT:

        - ``morphisms`` -- A KenzoChainComplexMorphism instance or a tuple of
          KenzoChainComplexMorphism instances.

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism`, difference of the morphism ``self``
          and the morphism(s) given by ``morphisms`` (if ``morphisms`` is None,
          ``self`` morphism is returned).
          For example, if ``morphisms`` = (mrph1, mrph2, mrph3) the result is
          ``self`` - mrph1 - mrph2 - mrph3.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)   # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                          # optional - kenzo
            sage: kenzo_chcm                                                     # optional - kenzo
            [K... Chain-Complex]
            sage: idnt = kenzo_chcm.identity_morphism()                          # optional - kenzo
            sage: idnt                                                           # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: opps_id = idnt.opposite()                                      # optional - kenzo
            sage: opps_id                                                        # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null = kenzo_chcm.null_morphism()                              # optional - kenzo
            sage: null                                                           # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: idx2 = idnt.substract(opps_id)                                 # optional - kenzo
            sage: opps_idx2 = idx2.substract\
            ....: (opps_id, idnt, idnt, null, idx2.substract(opps_id))           # optional - kenzo
            sage: kenzo_chcm.basis(4)                                            # optional - kenzo
            ['G4G0', 'G4G1']
            sage: c1 = kenzo_chcm(4, [2, 'G4G0', -5, 'G4G1'])                    # optional - kenzo
            sage: idx2.evaluation(c1)                                            # optional - kenzo
            <BLANKLINE>
            {K... CMBN 4}
            --------------------------------------
            <4 * G4G0>
            <-10 * G4G1>
            --------------------------------------
            <BLANKLINE>
            sage: opps_idx2.evaluation(c1)                                       # optional - kenzo
            <BLANKLINE>
            {K... CMBN 4}
            --------------------------------------
            <-4 * G4G0>
            <10 * G4G1>
            --------------------------------------
            <BLANKLINE>
        """
        if morphisms is None:
            return self
        if isinstance(morphisms, KenzoChainComplexMorphism):
            return KenzoChainComplexMorphism(__sbtr__(self._kenzo, morphisms._kenzo))
        elif isinstance(morphisms, tuple):
            rslt = self._kenzo
            for mrph in morphisms:
                rslt = __sbtr__(rslt, mrph._kenzo)
            return KenzoChainComplexMorphism(rslt)

    def mul(self, n):
        r"""
        Return the morphism ``n * self``.

        INPUT:

        - ``n`` -- An integer number.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)   # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                          # optional - kenzo
            sage: d = kenzo_chcm.differential()                                  # optional - kenzo
            sage: dx6 = d.mul(6)                                                 # optional - kenzo
            sage: kenzo_chcm.basis(4)                                            # optional - kenzo
            ['G4G0', 'G4G1']
            sage: c1 = kenzo_chcm(4, [-1, 'G4G0', 2, 'G4G1'])                    # optional - kenzo
            sage: d.evaluation(c1)                                               # optional - kenzo
            <BLANKLINE>
            {K... CMBN 3}
            --------------------------------------
            <3 * G3G0>
            <9 * G3G1>
            --------------------------------------
            <BLANKLINE>
            sage: dx6.evaluation(c1)                                             # optional - kenzo
            <BLANKLINE>
            {K... CMBN 3}
            --------------------------------------
            <18 * G3G0>
            <54 * G3G1>
            --------------------------------------
            <BLANKLINE>
        """
        if not n.is_integer():
            raise AssertionError("Parameter 'n' must be an integer number")
        if n == 1:
            return self
        if n == 0:
            return self.source_complex().null_morphism(target=self.target_complex(), \
                                                       degree=self.degree())
        if n == -1:
            return self.opposite()
        return KenzoChainComplexMorphism(__n_mrph__(n, self._kenzo))

    def change_source_target_complex(self, source=None, target=None):
        r"""
        Build, from the morphism ``self``, a new morphism with ``source``
        and ``target`` as source and target KenzoChainComplexes.

        INPUT:

        - ``source`` -- A KenzoChainComplex instance or None (default).

        - ``target`` -- A KenzoChainComplex instance or None (default).

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism` inheriting from ``self`` the
          degree (:degr slot in Kenzo), the algorithm (:intr slot in Kenzo)
          and the strategy (:strt slot in Kenzo). The source and target slots
          of this new morphism are given by the parameters ``source`` and
          ``target`` respectively; if any parameter is ommited, the corresponding
          slot is inherited from ``self``.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere, KenzoChainComplex # optional - kenzo
            sage: from sage.libs.ecl import ecl_eval
            sage: ZCC = KenzoChainComplex(ecl_eval("(z-chcm)"))               # optional - kenzo
            sage: ZCC                                                         # optional - kenzo
            [K... Chain-Complex]
            sage: s2 = Sphere(2)                                              # optional - kenzo
            sage: s3 = Sphere(3)                                              # optional - kenzo
            sage: tp = s2.tensor_product(s3)                                  # optional - kenzo
            sage: tp                                                          # optional - kenzo
            [K... Chain-Complex]
            sage: null = ZCC.null_morphism(tp)                                # optional - kenzo
            sage: null                                                        # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null2 = null.change_source_target_complex(source = tp)      # optional - kenzo
            sage: null2                                                       # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null2.source_complex() == null.source_complex()             # optional - kenzo
            False
        """
        source = source or self.source_complex()
        target = target or self.target_complex()
        return KenzoChainComplexMorphism(
            __change_sorc_trgt_aux__(self._kenzo, source._kenzo, target._kenzo))

    def destructive_change_source_target_complex(self, source=None, target=None):
        r"""
        Modify destructively the morphism ``self`` taking ``source`` and ``target``
        as source and target KenzoChainComplexes of ``self``.

        INPUT:

        - ``source`` -- A KenzoChainComplex instance or None (default).

        - ``target`` -- A KenzoChainComplex instance or None (default).

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism`. The source and target slots of
          ``self`` are replaced respectively by the parameters ``source`` and
          ``target``; if any parameter is ommited, the corresponding slot is
          inherited from ``self``.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere, KenzoChainComplex # optional - kenzo
            sage: from sage.libs.ecl import ecl_eval
            sage: ZCC = KenzoChainComplex(ecl_eval("(z-chcm)"))               # optional - kenzo
            sage: ZCC                                                         # optional - kenzo
            [K... Chain-Complex]
            sage: s2 = Sphere(2)                                              # optional - kenzo
            sage: s3 = Sphere(3)                                              # optional - kenzo
            sage: tp = s2.tensor_product(s3)                                  # optional - kenzo
            sage: tp                                                          # optional - kenzo
            [K... Chain-Complex]
            sage: null = ZCC.null_morphism(tp)                                # optional - kenzo
            sage: null                                                        # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null.target_complex()                                       # optional - kenzo
            [K... Chain-Complex]
            sage: null.destructive_change_source_target_complex(target = ZCC) # optional - kenzo
            [K... Cohomology-Class on K... of degree 0]
            sage: null.target_complex()                                       # optional - kenzo
            [K... Chain-Complex]
        """
        source = source or self.source_complex()
        target = target or self.target_complex()
        return KenzoChainComplexMorphism(
            __dstr_change_sorc_trgt_aux__(self._kenzo, source._kenzo, target._kenzo))

    def orgn(self):
        r"""
        Return the :orgn slot of Kenzo, which stores as a list the origin of the object.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere                # optional - kenzo
            sage: s2 = Sphere(2)                                          # optional - kenzo
            sage: i = s2.identity_morphism()                              # optional - kenzo
            sage: i.orgn()                                                # optional - kenzo
            '(IDNT-MRPH [K... Simplicial-Set])'
            sage: d = s2.differential()                                   # optional - kenzo
            sage: d.orgn()                                                # optional - kenzo
            '(ZERO-MRPH [K... Simplicial-Set] [K... Simplicial-Set] -1)'
        """
        return str(__orgn_aux__(self._kenzo))

    def __eq__(self, other):
        r"""
        Return ``True`` if and only if ``self == other`` as KenzoObjects.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import MooreSpace   # optional - kenzo
            sage: m24 = MooreSpace(2,4)                          # optional - kenzo
            sage: z = m24.null_morphism()                        # optional - kenzo
            sage: z == MooreSpace(2,4).null_morphism()           # optional - kenzo
            True
        """
        return Integer(str(__idnm__(self._kenzo))) == Integer(str(__idnm__(other._kenzo)))

def translate_type (kenzoobj):
    r"""
    Assign a class in Sage to the Kenzo objects.

    INPUT::

    - ``kenzoobj`` -- An ECL object in Kenzo

    EXAMPLES::

        sage: from sage.libs.ecl import ecl_eval
        sage: from sage.interfaces.kenzo import translate_type                 # optional - kenzo
        sage: M = ecl_eval("(moore 4 5)")                                      # optional - kenzo
        sage: M                                                                # optional - kenzo
        <ECL: [K... Simplicial-Set]>
        sage: TM = translate_type(M) ; TM                                      # optional - kenzo
        [K... Simplicial-Set]
        sage: type(TM)                                                         # optional - kenzo
        <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
        sage: L = ecl_eval("(list (z-chcm) (absm 7 'g) 1906)")                 # optional - kenzo
        sage: translate_type(L)                                                # optional - kenzo
        [[K... Chain-Complex], <AbSm 2-1-0 G>, 1906]
    """
    if __listp__(kenzoobj).python():
        return [translate_type(i) for i in EclListIterator(kenzoobj)]
    kenzotype = __type_of__(kenzoobj).python()
    if kenzotype == 'CHAIN-COMPLEX':
        return KenzoChainComplex(kenzoobj)
    elif kenzotype == 'MORPHISM':
        return KenzoChainComplexMorphism(kenzoobj)
    elif kenzotype == 'SIMPLICIAL-SET':
        return KenzoSimplicialSet(kenzoobj)
    elif kenzotype == 'SIMPLICIAL-MRPH':
        return KenzoSimplicialSetMorphism(kenzoobj)
    elif kenzotype == 'SPECTRAL-SEQUENCE':
        return KenzoSpectralSequence(kenzoobj)
    elif kenzotype == 'CRPR':
        return KenzoCRPRSimplex(kenzoobj)
    elif kenzotype == 'ABSM':
        return KenzoAbstractSimplex(kenzoobj)
    elif kenzotype == 'CMBN':
        return Kenzocmbn(kenzoobj)
    elif kenzotype == 'SYMBOL':
        return kenzoobj.python()
    else:
        return KenzoObject(kenzoobj)


class KenzoSimplicialSetMorphism(KenzoChainComplexMorphism):
    r"""
    Wrapper to simplicial set morphisms in Kenzo.

    EXAMPLES::

        sage: from sage.homology.simplicial_set_morphism import SimplicialSetMorphism
        sage: K = simplicial_sets.Simplex(1)
        sage: S1 = simplicial_sets.Sphere(1)
        sage: v0 = K.n_cells(0)[0]
        sage: v1 = K.n_cells(0)[1]
        sage: e01 = K.n_cells(1)[0]
        sage: w = S1.n_cells(0)[0]
        sage: sigma = S1.n_cells(1)[0]
        sage: f = {v0: w, v1: w, e01: sigma}
        sage: Sf = SimplicialSetMorphism(f, K, S1)
        sage: KSf = Sf._kenzo_repr ; KSf                                               # optional - kenzo
        [K... Simplicial-Morphism K... -> K...]
        sage: type(KSf)                                                                # optional - kenzo
        <class 'sage.interfaces.kenzo.KenzoSimplicialSetMorphism'>
    """

    def __call__(self, x):
        r"""
        Returns the image of ``x`` under the morphism.

        INPUT::

        - ``x`` -- A :class:`KenzoCombination` in the domain of ``self`` or a
        simplex of the domain.

        EXAMPLES::

            sage: from sage.homology.simplicial_set_morphism import SimplicialSetMorphism
            sage: K = simplicial_sets.Simplex(1)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: v0 = K.n_cells(0)[0]
            sage: v1 = K.n_cells(0)[1]
            sage: e01 = K.n_cells(1)[0]
            sage: w = S1.n_cells(0)[0]
            sage: sigma = S1.n_cells(1)[0]
            sage: f = {v0: w, v1: w, e01: sigma}
            sage: Sf = SimplicialSetMorphism(f, K, S1)
            sage: KSf = Sf._kenzo_repr
            sage: KSf(e01.apply_degeneracies(2,1))
            <AbSm 2-1 S...>
        """
        if isinstance(x, KenzoCombination):
            return Kenzocmbn(__evaluate_cmbn__(self._kenzo, x._kenzo)) #, self.target_simplicial_set())
        return translate_type(__evaluate_simplex__(self._kenzo, x.dimension(),
                                                   x._kenzo_repr._kenzo))

    def source_simplicial_set(self):
        r"""
        Returns the domain of the morphism ``self``.

        EXAMPLES::

            sage: S3 = simplicial_sets.Sphere(3)
            sage: f = S3.base_point_map()
            sage: Kf = f._kenzo_repr
            sage: Kf.source_simplicial_set() == f.domain()._kenzo_repr
            True
        """
        return KenzoSimplicialSet(__sorc_aux__(self._kenzo))

    def target_simplicial_set(self):
        r"""
        Returns the codomain of the morphism ``self``.

        EXAMPLES::

            sage: S = simplicial_sets.Simplex(7)
            sage: f = S.constant_map()
            sage: Kf = f._kenzo_repr
            sage: Kf.target_simplicial_set() == f.codomain()._kenzo_repr
            True
        """
        return KenzoSimplicialSet(__trgt_aux__(self._kenzo))

    def suspension(self, n=1):
        r"""
        Return the `n`-th suspension of this morphism of KenzoSimplicialSets.

        INPUT:

        - ``n`` (optional) -- non-negative integer, default 1

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KSimplicialSetMorphism      # optional - kenzo
            sage: T = simplicial_sets.Torus()
            sage: KT = T._kenzo_repr                                            # optional - kenzo
            sage: f = T.identity().suspension()
            sage: Kf = KSimplicialSetMorphism(f)                                # optional - kenzo
            sage: Kf.orgn()                                                     # optional - kenzo
            '(IDNT-SMMR [K... Simplicial-Set])'
            sage: Kf.source_simplicial_set() == KT.suspension()                 # optional - kenzo
            True
        """
        result = KenzoSimplicialSetMorphism(__suspension__(self._kenzo, n))
        setattr(result, '_base', self)
        return result

    def cone(self):
        r"""
        Return the cone of ``self``.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KSimplicialSetMorphism      # optional - kenzo
            sage: T = simplicial_sets.Torus()
            sage: f = T.identity()
            sage: Kf = KSimplicialSetMorphism(f)                                # optional - kenz
            sage: Kf.cone().orgn()     # description as a cone in Kenzo
            '(CONE2 [K... Simplicial-Morphism K... -> K...])'
        """
        return KenzoChainComplex(__cone2__(self._kenzo))

    def pushout(self, other):
        r"""
        Return the pushout of ``self`` and ``other``.

        INPUT:

        - ``other`` -- A :class:`KenzoSimplicialSetMorphism` with the same domain
          as ``self``.

        OUTPUT::

        - A :class:`KenzoSimplicialSet`.

        EXAMPLES::

            sage: S4 = simplicial_sets.Sphere(4)
            sage: T = simplicial_sets.Torus()
            sage: f = S4.base_point_map()
            sage: g = T.base_point_map()
            sage: p = f.pushout(g)
            sage: Kf = f._kenzo_repr
            sage: Kg = g._kenzo_repr
            sage: Kp = Kf.pushout(Kg)
            sage: Kp._maps == (Kf, Kg)
            True
        """
        if not self.source_simplicial_set() == other.source_simplicial_set():
            raise AssertionError('The domain of the morphisms must be the same')
        result = KenzoSimplicialSet(__pushout__(self._kenzo, other._kenzo))
        setattr(result, '_maps', (self, other))
        return result


def KSimplicialSetMorphism(simplicialmrph):
    r"""
    Construct a KenzoSimplicialSetMorphism from a SimplicialSetMorphism in Sage.

    INPUT:

    - ``simplicialmrph`` -- A morphism between simplicial sets.

    OUTPUT:

    - A :class:`KenzoSimplicialSetMorphism`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import KSimplicialSetMorphism        # optional - kenzo
        sage: S2 = simplicial_sets.Sphere(2)
        sage: f = S2.base_point_map()
        sage: Kf = KSimplicialSetMorphism(f)                                  # optional - kenzo
        sage: Kf                                                              # optional - kenzo
        [K... Simplicial-Morphism K... -> K...]
        sage: A = f.domain().n_cells(0)[0]
        sage: Kf(A)                                                           # optional - kenzo
        <AbSm - S...>
        sage: Kf(A).nondegenerate() == f(A)._kenzo_repr.nondegenerate()       # optional - kenzo
        True
        sage: Kf(A).degeneracies() == f(A)._kenzo_repr.degeneracies()         # optional - kenzo
        True
    """
    if simplicialmrph.is_identity():
        kdomain = simplicialmrph.domain()._kenzo_repr
        return KenzoSimplicialSetMorphism(__idnt_smmr__(kdomain._kenzo))
    sdict = simplicialmrph._dictionary
    data = {sigma: sdict[sigma] for sigma in sdict if sigma.is_nondegenerate()}
    smmrdict = []
    for x in sdict:
        if hasattr(x._kenzo_repr, 'nondegenerate'):
            eclx = EclObject(x._kenzo_repr.nondegenerate())
        else:
            eclx = x._kenzo_repr._kenzo
        smmrdict += [EclObject((eclx, sdict[x]._kenzo_repr._kenzo))]
    smmrdict = EclObject(smmrdict)
    source = KFiniteSimplicialSet(simplicialmrph.domain())
    target = KFiniteSimplicialSet(simplicialmrph.codomain())
    return KenzoSimplicialSetMorphism(
       __ksimplicialsetmorphism_aux__(smmrdict, source._kenzo, target._kenzo))


def build_morphism(source_complex, target_complex, degree, algorithm, strategy, orgn):
    r"""
    Build a morphism of chain complexes by means of the corresponding build-mrph
    Kenzo function.

    INPUT:

    - ``source_complex`` -- The source object as a KenzoChainComplex instance

    - ``target_complex`` -- The target object as a KenzoChainComplex instance

    - ``degree`` -- An integer number representing the degree of the morphism

    - ``algorithm`` -- A Lisp function defining the mapping (:intr slot in Kenzo)

    - ``strategy`` -- The strategy (:strt slot in Kenzo), which must be one of
      the two strings ``gnrt`` or ``cmbn``, depending if the ``algorithm`` (a Lisp
      function) uses as arguments a degree and a generator or a combination,
      respectively.

    - ``orgn`` -- A list containing a description about the origin of the morphism

    OUTPUT:

    - A :class:`KenzoChainComplexMorphism`.

    EXAMPLES::

        sage: from sage.interfaces.kenzo import KenzoChainComplex,\
        ....: build_morphism                                            # optional - kenzo
        sage: from sage.libs.ecl import ecl_eval
        sage: ZCC = KenzoChainComplex(ecl_eval("(z-chcm)"))             # optional - kenzo
        sage: A = build_morphism(ZCC, ZCC, -1,\
        ....: ecl_eval("#'(lambda (comb) (cmbn (1- (degr comb))))"),\
        ....: "cmbn", ["zero morphism on ZCC"])                         # optional - kenzo
        sage: A.target_complex()                                        # optional - kenzo
        [K... Chain-Complex]
        sage: A.degree()                                                # optional - kenzo
        -1
        sage: type(A)                                                   # optional - kenzo
        <class 'sage.interfaces.kenzo.KenzoChainComplexMorphism'>
    """
    return KenzoChainComplexMorphism(
        __build_mrph_aux__(source_complex._kenzo, target_complex._kenzo,
                       degree, algorithm, ":"+strategy, orgn))


def morphism_dictmat(morphism):
    r"""
    Computes a list of matrices in ECL associated to a morphism in Sage.

    INPUT:

    - ``morphism`` -- A morphism of chain complexes

    OUTPUT:

    - A :class:`EclObject`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import morphism_dictmat    # optional - kenzo
        sage: X = simplicial_complexes.Simplex(1)
        sage: Y = simplicial_complexes.Simplex(0)
        sage: g = Hom(X,Y)({0:0, 1:0})
        sage: f = g.associated_chain_complex_morphism()
        sage: morphism_dictmat(f)                                   # optional - kenzo
        <ECL: ((2 . #2A()) (1 . #2A()) (0 . #2A((1 1))))>
    """
    rslt = EclObject([])
    source = morphism.domain()
    d = source.differential()
    for k in d.keys():
        rslt = EclObject(k).cons(s2k_matrix(morphism.in_degree(k))).cons(rslt)
    return rslt


def KChainComplexMorphism(morphism):
    r"""
    Construct a KenzoChainComplexMorphism from a ChainComplexMorphism in Sage.

    INPUT:

    - ``morphism`` -- A morphism of chain complexes

    OUTPUT:

    - A :class:`KenzoChainComplexMorphism`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import KChainComplexMorphism           # optional - kenzo
        sage: C = ChainComplex({0: identity_matrix(ZZ, 1)}).dual()
        sage: D = ChainComplex({0: zero_matrix(ZZ, 1), 1: zero_matrix(ZZ, 1)}).dual()
        sage: f = Hom(D,C)({0: identity_matrix(ZZ, 1), 1: zero_matrix(ZZ, 1)})
        sage: g = KChainComplexMorphism(f)                                      # optional - kenzo
        sage: g                                                                 # optional - kenzo
        [K... Morphism (degree 0): K... -> K...]
        sage: g.source_complex()                                                # optional - kenzo
        [K... Chain-Complex]
        sage: g.target_complex()                                                # optional - kenzo
        [K... Chain-Complex]
    """
    if morphism.is_identity() == True:
        return KIdentity(morphism)
    source = KChainComplex(morphism.domain())
    target = KChainComplex(morphism.codomain())
    matrix_list = morphism_dictmat(morphism)
    return KenzoChainComplexMorphism(
        __kchaincomplexmorphism_aux__(matrix_list, source._kenzo, target._kenzo))

def KIdentity(identity):
    r"""
    Return the equivalent Kenzo identity morphism.

    INPUT:

    - ``identity`` -- An IdentityMorphism instance between ChainComplexes.

    OUTPUT:

    - The :meth:`identity_morphism` of ``identity`` domain.

    EXAMPLES::

        sage: from sage.interfaces.kenzo import KIdentity             # optional - kenzo
        sage: C = ChainComplex({0: identity_matrix(ZZ, 1)}, degree = -1)
        sage: iden = Hom(C,C).identity()
        sage: Kiden = KIdentity(iden)                                 # optional - kenzo
        sage: Kiden.composite(Kiden) == Kiden                         # optional - kenzo
        True
    """
    from sage.homology.chain_complex import ChainComplex_class
    domain = identity.domain()
    if isinstance(domain, ChainComplex_class) and hasattr(domain, '_kenzo_repr'):
        return domain._kenzo_repr.identity_morphism()
    else:
        raise AssertionError("The parameter is not an IdentityMorphism of valid \
                              Kenzo chain complexes")

def s2k_listofmorphisms(l):
    r"""
    Computes a list of morphisms of chain complexes in Kenzo from a list of morphisms in Sage.

    INPUT:

    - ``l`` -- A list of morphisms of chain complexes

    OUTPUT:

    - A :class:`EclObject`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import s2k_listofmorphisms  # optional - kenzo
        sage: C1 = ChainComplex({1: matrix(ZZ, 0, 2, [])}, degree_of_differential=-1)
        sage: C2 = ChainComplex({1: matrix(ZZ, 1, 2, [1, 0])},degree_of_differential=-1)
        sage: C3 = ChainComplex({0: matrix(ZZ, 0,2 , [])},degree_of_differential=-1)
        sage: M1 = Hom(C2,C1)({1: matrix(ZZ, 2, 2, [2, 0, 0, 2])})
        sage: M2 = Hom(C3,C2)({0: matrix(ZZ, 1, 2, [2, 0])})
        sage: l = [M1, M2]
        sage: s2k_listofmorphisms(l)                                 # optional - kenzo
        <ECL: ([K... Morphism (degree 0): K... -> K...] [K... Morphism (degree 0): K... -> K...])>
    """
    rslt = EclObject([])
    for m in l:
        rslt = EclObject(KChainComplexMorphism(m)._kenzo).cons(rslt)
    return __nreverse__(rslt)


def BicomplexSpectralSequence(l):
    r"""
    Construct the spectral sequence associated to the bicomplex given by a list of morphisms.

    INPUT:

    - ``l`` -- A list of morphisms of chain complexes

    OUTPUT:

    - A :class:`KenzoSpectralSequence`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import BicomplexSpectralSequence        # optional - kenzo
        sage: C1 = ChainComplex({1: matrix(ZZ, 0, 2, [])}, degree_of_differential=-1)
        sage: C2 = ChainComplex({1: matrix(ZZ, 1, 2, [1, 0])},degree_of_differential=-1)
        sage: C3 = ChainComplex({0: matrix(ZZ, 0,2 , [])},degree_of_differential=-1)
        sage: M1 = Hom(C2,C1)({1: matrix(ZZ, 2, 2, [2, 0, 0, 2])})
        sage: M2 = Hom(C3,C2)({0: matrix(ZZ, 1, 2, [2, 0])})
        sage: l = [M1, M2]
        sage: E = BicomplexSpectralSequence(l)                                   # optional - kenzo
        sage: E.group(2,0,1)                                                     # optional - kenzo
        Additive abelian group isomorphic to Z/2 + Z
        sage: E.table(3,0,2,0,2)                                                 # optional - kenzo
        0           0   0
        Z/2 + Z/4   0   0
        0           0   Z
        sage: E.matrix(2,2,0)                                                    # optional - kenzo
        [ 0  0]
        [-4  0]
    """
    return KenzoSpectralSequence(__bicomplex_spectral_sequence__(s2k_listofmorphisms(l)))
