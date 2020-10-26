# -*- coding: utf-8 -*-
r"""
Access to the KnotInfo database

This module contains the class :class:`KnotInfoBase` which is derived from :class:`Enum`
and provides knots and links listed in the databases at `KnotInfo <https://knotinfo.math.indiana.edu/>`__
and `LinkInfo <https://linkinfo.sitehost.iu.edu/>`__ as its items.

This interface contains a set of about twenty knots and links statically as demonstration cases. The complete
database can be installed as an optional Sage package using ``sage -i knotinfo``. This will be necessary to
have access to all the properties recorded in the databases, as well.

Be aware that there are a couple of conventions used differently on KnotInfo as in Sage, especially
concerning the selection of the symmetry version of the link. In our transitions to Sage objects
these are translated (by default) in order to avoid confusion about exchanged mirror versions.

Briefly, these differences are:

   - ``pd_notation`` --        KnotInfo: counter clockwise Sage: clockwise, see note in
     :meth:`KnotInfoBase.link`

   - ``homfly_polynomial`` --  KnotInfo: ``v``  Sage: `1/a`, see note in :meth:`KnotInfoBase.homfly_polynomial`.

   - ``braid_notation``    --  This is used accordingly: The crossing of the braid generators are positive
     in both systems. Here it is listed because there could arise confusion from the source where they are
     taken from. There, the braid generators are assumed to have a negative crossing
     (see definition 3  of `Gittings, T., "Minimum Braids: A Complete Invariant of Knots and Links <https://arxiv.org/abs/math/0401051>`__).


EXAMPLES::

    sage: from sage.knots.knotinfo import KnotInfo
    sage: L = KnotInfo.L4a1_0
    sage: L.pd_notation()
    [[6, 1, 7, 2], [8, 3, 5, 4], [2, 5, 3, 6], [4, 7, 1, 8]]
    sage: L.pd_notation(original=True)
    '{{6, 1, 7, 2}, {8, 3, 5, 4}, {2, 5, 3, 6}, {4, 7, 1, 8}}'
    sage: L.is_knot()
    False
    sage: L.num_components()
    2

Obtaining an instance of :class:`~sage.groups.braid.Braid`::

    sage: L.braid()
    s0*s1^-1*s2*s1^-1*s0^-1*s1^-1*s2^-1*s1^-1
    sage: type(_)
    <class 'sage.groups.braid.BraidGroup_class_with_category.element_class'>

Obtaining an instance of :class:`Link`::

    sage: l = L.link(); l
    Link with 2 components represented by 4 crossings
    sage: type(l)
    <class 'sage.knots.link.Link'>

If you have `SnapPy <https://snappy.math.uic.edu/index.html>`__ installed inside Sage
you can obtain an instance of :class:`~spherogram.links.links_base.Link`, too::

    sage: L6 = KnotInfo.L6a2_0
    sage: l6s = L6.link(snappy=True); l6s                              # optional - snappy
    Plink failed to import tkinter.
    <Link: 2 comp; 6 cross>
    sage: type(l6s)                                                    # optional - snappy
    <class 'spherogram.links.invariants.Link'>
    sage: l6  = L6.link()
    sage: l6 == l6s.sage_link()                                        # optional - snappy
    True
    sage: l6sn = L6.link(use_item=L6.items.name, snappy=True); l6sn    # optional - snappy
    <Link L6a2: 2 comp; 6 cross>
    sage: l6s == l6sn                                                  # optional - snappy
    False
    sage: l6sn.sage_link().is_isotopic(l6)                             # optional - snappy
    True

But observe that the name conversion to SnapPy does not distingiush orientation types::

    sage: L6b = KnotInfo.L6a2_1
    sage: l6bsn = L6b.link(use_item=L6b.items.name, snappy=True)       # optional - snappy
    sage: l6bsn.PD_code() == l6sn.PD_code()                            # optional - snappy
    True

Obtaining the HOMFLY-PT polynomial::

    sage: L.homfly_polynomial()
    -v^-1*z - v^-3*z - v^-3*z^-1 + v^-5*z^-1
    sage: L.homfly_polynomial(sage_convention=True)
    L^5*M^-1 - L^3*M - L^3*M^-1 - L*M
    sage: _ == l.homfly_polynomial(normalization='az')
    True

Items for knots need a leading ``K`` for technical reason::

    sage: K = KnotInfo.K4_1
    sage: K.is_knot()
    True
    sage: K.crossing_number()
    4
    sage: K.gauss_notation()
    [-1, 2, -3, 1, -4, 3, -2, 4]
    sage: K.dt_notation()
    [4, 6, 8, 2]
    sage: K.determinant()
    5
    sage: K.symmetry_type()
    'fully amphicheiral'
    sage: _ == K[K.items.symmetry_type]
    True
    sage: K.is_reversible()
    True
    sage: K.is_amphicheiral()
    True

Obtaining the original string from the database for an arbitrary property::

    sage: K[K.items.classical_conway_name]         # optional - database_knotinfo
    '4_1'

Using the ``column_type`` of a property::

    sage: [i.column_name() for i in K.items if i.column_type() != i.types.OnlyLinks and K[i] == 'Y']     # optional - database_knotinfo
    ['Alternating', 'Fibered', 'Quasialternating', 'Adequate']

You can launch web-pages attached to the links::

    sage: K.diagram()                 # not tested
    True
    sage: L.diagram(single=True)      # not tested
    True
    sage: L.knot_atlas_webpage()      # not tested
    True
    sage: K.knotilus_webpage()        # not tested
    True

and the description web-pages of the properties::

    sage: K.items.positive.description_webpage()  # not tested
    True

To see all the properties available in this interface you can use "tab-completion".
For example type ``K.items.`` and than hit the "tab-key". You can select the item
you want from the list. If you know some first letters type them first to obtain a
reduced selection list.

In a similar way you may select the knots and links. Here you have to type ``KnotInfo.``
or ``KnotInfo.L7`` before stroking the "tab-key". In the latter case  the selection list
will be reduced to proper links with 7 crossings.

Finally there is a method :meth:`Link.identify_knotinfo` of class :class:`Link` to find an instance
in the KnotInfo database::

    sage: L = Link([[3,1,2,4], [8,9,1,7], [5,6,7,3], [4,18,6,5],
    ....:           [17,19,8,18], [9,10,11,14], [10,12,13,11],
    ....:           [12,19,15,13], [20,16,14,15], [16,20,17,2]])
    sage: L.identify_knotinfo()
    (<KnotInfo.K0_1: [0, 1]>, False)


REFERENCES:

- `KnotInfo <https://knotinfo.math.indiana.edu/>`__
- `LinkInfo <https://linkinfo.sitehost.iu.edu/>`__



AUTHORS:

- Sebastian Oehms August 2020: initial version
"""


##############################################################################
#       Copyright (C) 2020 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################



from enum import Enum
from sage.misc.cachefunc import cached_method, cached_function
from sage.misc.sage_eval import sage_eval
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.groups.braid import BraidGroup
from sage.knots.knot import Knots
from sage.databases.knotinfo_db import KnotInfoColumnTypes, KnotInfoColumns, db




def is_knotinfo_available(raise_error=False):
    r"""
    Return whether the KnotInfo databases are installed or not.

    INPUT:

    - ``raise_error`` -- boolean (default ``False``) if set to ``True``
      an import error is raised in case KnotInfo is not installed

    EXAMPLES::

        sage: from sage.knots.knotinfo import is_knotinfo_available
        sage: is_knotinfo_available()     # optional - database_knotinfo
        True
    """
    res = db.is_available()
    if not res and raise_error:
        raise ImportError('This functionality needs KnotInfo to be installed! Type `sage -i knotinfo` to have this done')
    return res

@cached_function
def knotinfo_matching_list(number_of_crossings, num_components, homfly_polynomial=None):
    r"""
    Return a list of links from the KontInfo and LinkInfo tables with given properties.

    INPUT:

    - ``number_of_crossings``  -- Python ``int`` giving the (not necessarily minimal)
      number of crossings to be matched
    - ``num_components``   -- Python ``int`` giving the number of components
      to be matched
    - ``homfly_polynomial``  -- instance of :class:`~sage.rings.polynomial.laurent_polynomial_ring.LaurentPolynomial_mpair`
      giving the HOMFLY-PT polynomial to be matched

    EXAMPLES::

        sage: from sage.knots.knotinfo import KnotInfo, knotinfo_matching_list
        sage: knotinfo_matching_list(3,1)
        [<KnotInfo.K0_1: [0, 1]>, <KnotInfo.K3_1: [1, 1]>]
        sage: [l.name for l in knotinfo_matching_list(4,2)]
        ['L2a1_0', 'L2a1_1', 'L4a1_0', 'L4a1_1']
        sage: L = KnotInfo.L6a3_0                            # optional - database_knotinfo
        sage: h = L.homfly_polynomial(sage_convention=True)  # optional - database_knotinfo
        sage: l = knotinfo_matching_list(L.crossing_number(), L.num_components(), h)  # optional - database_knotinfo
        sage: len(l) == 1 and l[0] == L                      # optional - database_knotinfo
        True
    """
    res = []
    if homfly_polynomial:
        l = knotinfo_matching_list(number_of_crossings, num_components)
        for L in l:
            if homfly_polynomial:
                if L.homfly_polynomial(sage_convention=True) != homfly_polynomial:
                    continue
            res.append(L)
        return res

    for L in KnotInfo:
        if L.crossing_number() > number_of_crossings:
            continue
        if L.num_components() != num_components:
            continue
        res.append(L)

    return res


def eval_knotinfo(string, locals={}, to_tuple=True):
    r"""
    Preparse a string from the KnotInfo database and evaluate it by ``sage_eval``.

    INPUT:

    - ``string``  -- string that gives a value of some database entry
    - ``locals``      -- dictionary of locals passed to ``sage_eval``

    EXAMPLES::

        sage: from sage.knots.knotinfo import KnotInfo, eval_knotinfo
        sage: L = KnotInfo.L4a1_0
        sage: L.braid_notation(original=True)
        '{4, {1, -2, 3, -2, -1, -2, -3, -2}}'
        sage: eval_knotinfo(_)
        (4, (1, -2, 3, -2, -1, -2, -3, -2))
    """
    if to_tuple:
        new_string = string.replace('{', '(')
        new_string = new_string.replace('}', ')')
    else:
        new_string = string.replace('{', '[')
        new_string = new_string.replace('}', ']')
    new_string = new_string.replace(';', ',')
    return sage_eval(new_string, locals=locals)

def knotinfo_bool(string):
    r"""
    Preparse a string from the KnotInfo database representing a boolean.

    INPUT:

    - ``string``  -- string that gives a value of some database entry

    EXAMPLES::

        sage: from sage.knots.knotinfo import knotinfo_bool
        sage: knotinfo_bool('Y')
        True
    """
    if string == 'Y':
        return True
    elif string == 'N':
        return False
    raise ValueError('%s is not a KnotInfo boolean')



# ---------------------------------------------------------------------------------
# KnotInfoBase
# ---------------------------------------------------------------------------------
class KnotInfoBase(Enum):
    r"""
    Enum class to select the knots and links listed in the databases at
    `KnotInfo <https://knotinfo.math.indiana.edu/>`__ and `LinkInfo <https://linkinfo.sitehost.iu.edu/>`__.

    EXAMPLES::

        sage: from sage.knots.knotinfo import KnotInfo
        sage: [knot.name for knot in KnotInfo if knot.crossing_number() < 5]
        ['K0_1', 'K3_1', 'K4_1', 'L2a1_0', 'L2a1_1', 'L4a1_0', 'L4a1_1']

    TESTS:

        sage: KnotInfo.K7_1.inject()
        Defining K7_1
        sage: TestSuite(K7_1).run()
    """
    @property
    def items(self):
        r"""
        Return an Enum class to select a column item of the KnotInfo database.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: it = L.items
            sage: [i.name for i in it if i.name.endswith('notation')]   # optional - database_knotinfo
            ['dt_notation',
             'conway_notation',
             'two_bridge_notation',
             'gauss_notation',
             'enhanced_gauss_notation',
             'pd_notation',
             'braid_notation',
             'positive_braid_notation',
             'positive_pd_notation',
             'strongly_quasipositive_braid_notation',
             'quasipositive_braid_notation',
             'arc_notation']
            sage: L.items.dt_notation.column_name()
            'DT Notation'

        To check if the item is available for proper links or only knots type::

            sage: it.gauss_notation.column_type()
            <KnotInfoColumnTypes.KnotsAndLinks: 'B'>
            sage: it.dt_notation.column_type()
            <KnotInfoColumnTypes.OnlyKnots: 'K'>

        To see the description of the item in your web browser type::

            sage: it.gauss_notation.description_webpage()    # not tested
            True
        """
        return db.columns()

    @cached_method
    def __getitem__(self, item):
        r"""
        sage: from sage.knots.knotinfo import KnotInfo
        sage: L = KnotInfo.L4a1_0
        sage: L[L.items.alternating]
        'Y'
        sage: L[L.items.arc_notation]
        '{{6, 4}, {3, 5}, {4, 2}, {1, 3}, {2, 6}, {5, 1}}'
        sage: L[L.items.braid_notation]
        '{4, {1, -2, 3, -2, -1, -2, -3, -2}}'
        sage: L[0]
        Traceback (most recent call last):
        ...
        KeyError: "Item must be an instance of <enum 'KnotInfoColumns'>"
        """
        if not isinstance(item, KnotInfoColumns):
            raise KeyError('Item must be an instance of %s' %(KnotInfoColumns))
        if item.column_type() == item.types.OnlyLinks and self.is_knot():
            raise KeyError('Item not available for knots' %(KnotInfoColumns))
        if item.column_type() == item.types.OnlyKnots and not self.is_knot():
            raise KeyError('Item not available for links' %(KnotInfoColumns))

        l = db.read(item)
        offset = 0
        if item.column_type() == item.types.OnlyLinks:
            offset = self._offset_knots()

        return l[self.value[0]-offset]

    def _offset_knots(self):
        r"""
        Return the list index of the first proper link in a combined
        list containing knots and proper links together which is the
        case for columns used for KnotInfo and LinkInfo in common.
        This index is exactly the total number of knots recorded
        in KnotInfo.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L._offset_knots()          # optional - database_knotinfo
            2978
        """
        return db.read_num_knots()

    @cached_method
    def _braid_group(self):
        r"""
        Return the braid group corresponding to the braid index
        of ``self``.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L._braid_group()
            Braid group on 4 strands
        """
        n = self.braid_index()
        if n == 1:
            return BraidGroup(2)
        else:
            return BraidGroup(n)


    @cached_method
    def _homfly_pol_ring(self, var1, var2):
        r"""
        Return the parent Laurent polynomial ring for the HOMFLY-PT
        polynomial according to Sage's internal one.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_1
            sage: L._homfly_pol_ring('u', 'v')
            Multivariate Laurent Polynomial Ring in u, v over Integer Ring
        """
        K3_1 = Knots().from_table(3,1)
        return K3_1.homfly_polynomial(var1=var1, var2=var2).parent()

    @cached_method
    def pd_notation(self, original=False):
        r"""
        Return the value of column ``pd_notation`` for this
        link as a Python list of Python lists.

        INPUT:

        - ``original`` -- boolean (optional, default ``False``) if set to
          ``True`` the original table entry is returned as a string

        OUTPUT:

        Python list of python lists each entry of the outer list
        representing a crossing.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L.pd_notation()
            [[6, 1, 7, 2], [8, 3, 5, 4], [2, 5, 3, 6], [4, 7, 1, 8]]
            sage: L.pd_notation(original=True)
            '{{6, 1, 7, 2}, {8, 3, 5, 4}, {2, 5, 3, 6}, {4, 7, 1, 8}}'
            sage: K = KnotInfo.K4_1
            sage: K.pd_notation()
            [[4, 2, 5, 1], [8, 6, 1, 5], [6, 3, 7, 4], [2, 7, 3, 8]]
        """
        if self.is_knot():
            pd_notation = self[self.items.pd_notation]
        else:
            pd_notation = self[self.items.pd_notation_vector]

        if original:
            return pd_notation

        if not pd_notation:
            # don't forget the unknot
            return []

        return eval_knotinfo(pd_notation, to_tuple=False)

    @cached_method
    def dt_notation(self, original=False):
        r"""
        Return the value of column ``dt_notation`` for this
        link as a Python list of Python lists.

        INPUT:

        - ``original`` -- boolean (optional, default ``False``) if set to
          ``True`` the original table entry is returned as a string

        OUTPUT:

        Python list of python lists each entry of the outer list
        representing a crossing.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L.dt_notation()
            [[6, 8], [2, 4]]
            sage: L.dt_notation(original=True)
            '[{6, 8}, {2, 4}]'
            sage: L = KnotInfo.L4a1_0
            sage: K = KnotInfo.K4_1
            sage: K.dt_notation()
            [4, 6, 8, 2]
        """
        if self.is_knot():
            dt_notation = self[self.items.dt_notation]
        else:
            dt_notation = self[self.items.dt_code]

        if original:
            return dt_notation

        if not dt_notation:
            # don't forget the unknot
            return []

        return eval_knotinfo(dt_notation, to_tuple=False)

    @cached_method
    def gauss_notation(self, original=False):
        r"""
        Return the value of column ``gauss_notation`` for this
        link as a Python list of Python lists.

        INPUT:

        - ``original`` -- boolean (optional, default ``False``) if set to
          ``True`` the original table entry is returned as a string

        OUTPUT:

        Python list of

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L.gauss_notation()
            [[1, -3, 2, -4], [3, -1, 4, -2]]
            sage: L.gauss_notation(original=True)
            '{{1, -3, 2, -4}, {3, -1, 4, -2}}'
        """
        gauss_notation = self[self.items.gauss_notation]
        if original:
            return gauss_notation

        if not gauss_notation:
            # don't forget the unknot
            return []

        return eval_knotinfo(gauss_notation, to_tuple=False)

    @cached_method
    def braid_notation(self, original=False):
        r"""
        Return the value of column ``braid_notation`` for this
        link as a Python tuple (Tietze form).

        INPUT:

        - ``original`` -- boolean (optional, default ``False``) if set to
          ``True`` the original table entry is returned as a string

        OUTPUT:

        Python tuple representing the braid whose closure is ``self``
        in Tietze form.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L.braid_notation()
            (1, -2, 3, -2, -1, -2, -3, -2)
            sage: L.braid_notation(original=True)
            '{4, {1, -2, 3, -2, -1, -2, -3, -2}}'
        """
        braid_notation = self[self.items.braid_notation]
        if original:
            return braid_notation

        if not braid_notation:
            # don't forget the unknot
            return (1, -1)

        braid_notation = eval_knotinfo(braid_notation)
        if type(braid_notation) is list:
            # in some cases there are a pair of braid representations
            # in the database. If this is the case we select the
            # corresponding to the braid index.
            if type(braid_notation[0]) is tuple:
                i = self.braid_index()
                for b in braid_notation:
                    if -i < min(b) and max(b) < i:
                        braid_notation = b
                        break

        if not self.is_knot():
            # in linkinfo the braid_notation includes the braid_index as first item of a pair
            braid_notation = braid_notation[1]
        return braid_notation

    @cached_method
    def braid_index(self):
        r"""
        Return the value of column ``braid_index`` for this
        link as a Python int.

        OUTPUT:

        Python int giving the minimum of strands needed to
        represent ``self`` as closure of a braid.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L4a1_0
            sage: L.braid_index()
            4
        """
        if self.is_knot():
            return int(self[self.items.braid_index])
        else:
            braid_notation = self[self.items.braid_notation]
            braid_notation = eval_knotinfo(braid_notation)
            return int(braid_notation[0])

    @cached_method
    def braid_length(self):
        r"""
        Return the value of column ``braid_length`` for this
        link as a Python int.

        OUTPUT:

        Python int giving the minimum length of a braid word
        needed to represent ``self`` as closure of a braid.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K = KnotInfo.K3_1
            sage: K.braid_length()
            3
        """
        return int(self[self.items.braid_length])

    @cached_method
    def braid(self):
        r"""
        Return the braid notation of self as an instance of :class:`~sage.groups.braid.Braid`.


        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K = KnotInfo.K3_1
            sage: K.braid()
            s^3
            sage: K.braid_notation()
            (1, 1, 1)
        """
        return self._braid_group()(self.braid_notation())

    @cached_method
    def num_components(self):
        r"""
        Return the number of components of ``self``.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.L6a1_0.num_components()
            2
        """
        return self.value[1]

    @cached_method
    def crossing_number(self):
        r"""
        Return the minimal number of crossings of ``self``.

        .. NOTE::

           In contrast to the number of crossings displayed for instances
           of :class:`Link` this number is the minimum over all possible
           diagrams of the link. The number of crossings displayed in
           the representation string of :class:`Link` refers to the
           special diagram which could be larger.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.L4a1_0.crossing_number()
            4
            sage: KnotInfo.K3_1.crossing_number()
            3
            sage: Link(KnotInfo.L4a1_0.braid())
            Link with 2 components represented by 8 crossings
        """
        return int(self[self.items.crossing_number])

    @cached_method
    def determinant(self):
        r"""
        Return the determinant of ``self``.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.L4a1_0.determinant()
            4
            sage: KnotInfo.K3_1.determinant()
            3
        """
        return int(self[self.items.determinant])

    @cached_method
    def is_knot(self):
        r"""
        Return whether ``self`` is a knot or a proper link.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.L7a1_0.is_knot()      # optional - database_knotinfo
            False
            sage: KnotInfo.K6_3.is_knot()
            True
        """
        return self.num_components() == 1

    @cached_method
    def name_unoriented(self):
        r"""
        Return the the part of the name of ``self`` which is independent on the orientation.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.L10a122_1_0.name_unoriented()  # optional - database_knotinfo
            'L10a122'
        """
        return self[self.items.name_unoriented]

    @cached_method
    def symmetry_type(self):
        r"""
        Return the symmetry type of ``self``.

        From the KnotInfo description page:

            If a knot is viewed as the oriented diffeomorphism
            class of an oriented pair, `K = (S_3, S_1)`, with `S_i`
            diffeomorphic to `S^i`, there are four oriented knots
            associated to any particular knot `K`. In addition to
            `K` itself, there is the reverse, `K^r = (S_3, -S_1)`,
            the concordance inverse, `-K = (-S_3, -S_1)`, and the
            mirror image, `K^m = (-S_3, S_1)`. A knot is called
            reversible if `K = K^r`, negative amphicheiral if
            `K = -K`, and positive amphicheiral if `K = K^m`.

            A knot possessing any two of these types of symmetry
            has all three. Thus, in the table, a knot is called
            reversible if that is the only type of symmetry it has,
            and likewise for negative amphicheiral. If it has none
            of these types of symmetry it is called chiral, and if
            it has all three it is called fully amphicheiral.

            For prime knots with fewer than 12 crossings, all
            amphicheiral knots are negative amphicheiral.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: [(L.name, L.symmetry_type()) for L in KnotInfo if L.is_knot() and L.crossing_number() < 6]
            [('K0_1', 'fully amphicheiral'),
            ('K3_1', 'reversible'),
            ('K4_1', 'fully amphicheiral'),
            ('K5_1', 'reversible'),
            ('K5_2', 'reversible')]
        """
        if not self.is_knot():
            raise NotImplementedError('This is only available for knots')
        if not self[self.items.symmetry_type] and self.crossing_number() == 0:
            return 'fully amphicheiral'
        return self[self.items.symmetry_type]

    @cached_method
    def is_reversible(self):
        r"""
        Return whether ``self`` is reversible.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.K6_3.is_reversible()
            True
        """
        if self.symmetry_type() == 'reversible':
            return True
        if self.symmetry_type() == 'fully amphicheiral':
            return True
        return False

    @cached_method
    def is_amphicheiral(self, positive=False):
        r"""
        Return whether ``self`` is amphicheiral.

        INPUT:

        - ``positive`` -- boolean (default False) whether to check
          if ``self`` is positive or negative amphicheiral (see
          doctest of :meth:`symmetry_type`)

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K = KnotInfo.K12a_427                 # optional - database_knotinfo
            sage: K.is_amphicheiral()                   # optional - database_knotinfo
            False
            sage: K.is_amphicheiral(positive=True)      # optional - database_knotinfo
            True
        """
        if positive:
            if self.symmetry_type() == 'positive amphicheiral':
                return True
        else:
            if self.symmetry_type() == 'negative amphicheiral':
                return True
        if self.symmetry_type() == 'fully amphicheiral':
            return True
        return False

    @cached_method
    def is_alternating(self):
        r"""
        Return whether ``self`` is alternating.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.K5_2.is_alternating()
            True
        """
        return knotinfo_bool(self[self.items.alternating])

    @cached_method
    def is_almost_alternating(self):
        r"""
        Return whether ``self`` is almost alternating.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.K5_2.is_almost_alternating()        # optional - database_knotinfo
            False
        """
        is_knotinfo_available(raise_error=True) # column not available in demo-version
        return knotinfo_bool(self[self.items.almost_alternating])

    @cached_method
    def is_quasi_alternating(self):
        r"""
        Return whether ``self`` is quasi alternating.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.K5_2.is_quasi_alternating()         # optional - database_knotinfo
            True
        """
        is_knotinfo_available(raise_error=True) # column not available in demo-version
        return knotinfo_bool(self[self.items.quasi_alternating])

    @cached_method
    def is_adequate(self):
        r"""
        Return whether ``self`` is adequate.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.K5_2.is_adequate()                  # optional - database_knotinfo
            True
        """
        is_knotinfo_available(raise_error=True) # column not available in demo-version
        return knotinfo_bool(self[self.items.adequate])

    @cached_method
    def is_positive(self):
        r"""
        Return whether ``self`` is positive.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.K5_2.is_positive()
            True
        """
        return knotinfo_bool(self[self.items.positive])

    @cached_method
    def is_quasipositive(self):
        r"""
        Return whether ``self`` is quasi-positive.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.K5_2.is_quasipositive()              # optional - database_knotinfo
            True
        """
        is_knotinfo_available(raise_error=True) # column not available in demo-version
        return knotinfo_bool(self[self.items.quasipositive])

    @cached_method
    def is_strongly_quasipositive(self):
        r"""
        Return whether ``self`` is strongly quasi-positive.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.K5_2.is_strongly_quasipositive()     # optional - database_knotinfo
            True
        """
        is_knotinfo_available(raise_error=True) # column not available in demo-version
        return knotinfo_bool(self[self.items.strongly_quasipositive])

    @cached_method
    def is_positive_braid(self):
        r"""
        Return whether ``self`` is a positive braid.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.K5_2.is_positive_braid()             # optional - database_knotinfo
            False
        """
        is_knotinfo_available(raise_error=True) # column not available in demo-version
        return knotinfo_bool(self[self.items.positive_braid])

    @cached_method
    def is_fibered(self):
        r"""
        Return whether ``self`` is fibered.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.K6_3.is_fibered()
            True
        """
        return knotinfo_bool(self[self.items.fibered])

    @cached_method
    def is_oriented(self):
        r"""
        Return whether ``self`` is oriented.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.L6a2_1.is_oriented()
            True
        """
        return not knotinfo_bool(self[self.items.unoriented])


    @cached_method
    def homfly_polynomial(self, var1=None, var2=None, original=False, sage_convention=False):
        r"""
        Return the HOMFLY-PT polynomial according to the value of column
        ``homfly_polynomial`` for this knot or link (in this case the
        column ``homflypt_polynomial`` is used) as an instance of the
        element class according to the output of :meth:`Link.homfly_polynomial`
        of :class:`Link`.

        The HOMFLY-PT polynomial `P(L)` of a link `L` satisfies the following skein relation
        (see the corresponding `KnotInfo description page <https://knotinfo.math.indiana.edu/descriptions/jones_homfly_kauffman_description/polynomial_defn.html)>`__):

        .. MATH::

            P(O) = 1,\,\,\,   v^{-1} P(L_+) -  v P(L_-) = z P(L_0)

        INPUT:

        - ``var1`` -- string for the name of the first variable (default depending
          on keyword ``sage_convention``: ``'v'`` or ``'L'`` if ``sage_convention == True``)
        - ``var2`` -- string for the name of the second variable (default depending
          on keyword ``sage_convention``: ``'z'`` or ``'M'`` if ``sage_convention == True``)
        - ``original`` -- boolean (default ``False``) if set to
          ``True`` the original table entry is returned as a string
        - ``sage_convention`` -- boolean (default ``False``) if set to ``True`` the conversion
          to Sage's conventions (see the note below) is performed

        OUTPUT:

        A Laurent polynomial over the integers, more precisely an instance of
        :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial_mpair`.
        If ``original`` is set to ``False`` then a string is returned.

        .. NOTE::

            The skein-relation for the HOMFLY-PT polynomial given on KnotInfo
            differs from the ones used in Sage.

            Using Sage's HOMFLY-PT polynomial with ``normalization='az'``
            the corresponding skein-relation is (see :meth:`Link.homfly_polynomial`
            of :class:`Link`):

            .. MATH::

                P(O) = 1,\,\,\,    a P(L_+) - a^{-1} P(L_-) = z P(L_0)

            Thus, the HOMFLY-PT polynomial of KnotInfo compares to the one of Sage
            by replacing ``v`` by ``~a``. To keep them comparable this translation
            can be performed by setting the keyword ``sage_convention`` to ``True``.


        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K3_1 = KnotInfo.K3_1
            sage: PK3_1 = K3_1.homfly_polynomial(); PK3_1
            -v^4 + v^2*z^2 + 2*v^2
            sage: K3_1.homfly_polynomial(original=True)
            '(2*v^2-v^4)+ (v^2)*z^2'
            sage: PK3_1s = K3_1.homfly_polynomial(sage_convention=True); PK3_1s
            L^-2*M^2 + 2*L^-2 - L^-4
            sage: PK3_1s == K3_1.link().homfly_polynomial(normalization='az')
            True

        for proper links::

            sage: L4a1_1 = KnotInfo.L4a1_1
            sage: PL4a1_1 = L4a1_1.homfly_polynomial(var1='x', var2='y'); PL4a1_1
            -x^5*y + x^3*y^3 - x^5*y^-1 + 3*x^3*y + x^3*y^-1
            sage: PL4a1_1s = L4a1_1.homfly_polynomial(var1='x', var2='y', sage_convention=True); PL4a1_1s
            x^-3*y^3 + 3*x^-3*y + x^-3*y^-1 - x^-5*y - x^-5*y^-1
            sage: PL4a1_1s == L4a1_1.link().homfly_polynomial(var1='x', var2='y', normalization='az')
            True

        check the skein-relation from the KnotInfo description page (applied to one of
        the positive crossings of the right-handed trefoil)::

            sage: R = PK3_1.parent()
            sage: PO = R.one()
            sage: L2a1_1 = KnotInfo.L2a1_1
            sage: PL2a1_1 = L2a1_1.homfly_polynomial()
            sage: v, z = R.gens()
            sage: ~v*PK3_1 -v*PO == z*PL2a1_1
            True

        check the skein-relation given in the doc string of :meth:`Link.homfly_polynomial` of
        :class:`Link` (applied to one of the positive crossings of the right-handed trefoil)::

            sage: Rs = PK3_1s.parent()
            sage: POs = Rs.one()
            sage: PL2a1_1s = L2a1_1.homfly_polynomial(sage_convention=True)
            sage: a, z = Rs.gens()
            sage: a*PK3_1s - ~a*POs == z*PL2a1_1s
            True


        TESTS::

            all(L.homfly_polynomial() == L.link().homfly_polynomial(normalization='az') for L in KnotInfo if L.crossing_number() > 0 and L.crossing_number() < 7)
            True

        REFERENCES:

        - :wikipedia:`HOMFLY_polynomial`
        """
        if self.is_knot():
            homfly_polynomial = self[self.items.homfly_polynomial]
        else:
            homfly_polynomial = self[self.items.homflypt_polynomial]

        if original:
            return homfly_polynomial

        if sage_convention:
            if not var1:
                var1='L'
            if not var2:
                var2='M'
        else:
            if not var1:
                var1='v'
            if not var2:
                var2='z'

        R = self._homfly_pol_ring(var1, var2)
        if not homfly_polynomial and self.crossing_number() == 0:
            return R.one()

        L, M = R.gens()
        if sage_convention:
            lc = {'v': ~L, 'z':M}  # see note above
        else:
            lc = {'v': L, 'z':M}
        return eval_knotinfo(homfly_polynomial, locals=lc)

    @cached_method
    def kauffman_polynomial(self, var1='a', var2='z', original=False):
        r"""
        Return the Kauffman polynomial according to the value of column
        ``kauffman_polynomial`` for this knot or link as an instance of
        :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial_mpair`.

        The Kauffman polynomial `F(L)` respectivlely its corresponding invariant under
        regular isotopy `\Delta (L) = a^{w(L)} F(L)` where `w(L)` is the writhe of
        the link `L` satisfies the following skein relation
        (see the corresponding `KnotInfo description page <https://knotinfo.math.indiana.edu/descriptions/jones_homfly_kauffman_description/polynomial_defn.html)>`__):

        .. MATH::

            \Delta(O) = 1,\,\,\,   \Delta(L_+) -  \Delta(L_-) = z (\Delta(L_0 + \Delta(L_{\infty}))

        Furthermore, removing a curl of sign `\epsilon` leads to a multiplication of `\Delta(L)`
        with `a^{\epsilon}`.

        INPUT:

        - ``var1`` -- (default: ``'a'``) the first variable
        - ``var2`` -- (default: ``'z'``) the second variable
        - ``original`` -- boolean (optional, default ``False``) if set to
          ``True`` the original table entry is returned as a string

        OUTPUT:

        A Laurent polynomial over the integers, more precisely an instance of
        :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial_mpair`.
        If ``original`` is set to ``False`` then a string is returned.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L2a1_1
            sage: L.kauffman_polynomial()
            a^-1*z - a^-1*z^-1 + a^-2 + a^-3*z - a^-3*z^-1
            sage: K = KnotInfo.K4_1
            sage: K.kauffman_polynomial()
            a^2*z^2 + a*z^3 - a^2 - a*z + 2*z^2 + a^-1*z^3 - 1 - a^-1*z + a^-2*z^2 - a^-2

        Comparison with Jones polynomial::

            sage: k    = _
            sage: a, z = k.variables()
            sage: j    = K.jones_polynomial(skein_normalization=True)
            sage: t,   = j.variables()
            sage: k.subs(a=-t^3, z=~t+t) == j.subs(t=t^4)
            True

        Check the skein relation::

            sage: K3_1    = KnotInfo.K3_1
            sage: FK3_1   = K3_1.kauffman_polynomial()
            sage: FL2a1_1 = L.kauffman_polynomial()
            sage: z, a    = FK3_1.variables()
            sage: ΔK3_1   = FK3_1   * a**K3_1.link().writhe()
            sage: ΔL2a1_1 = FL2a1_1 * a**L.link().writhe()
            sage: ΔO1p    = a          # unknot with one positive curl
            sage: ΔO2n    = a**-2      # unknot with two negative curls
            sage: ΔK3_1 + ΔO1p == z*(ΔL2a1_1 + ΔO2n)
            True

        REFERENCES:

        - :wikipedia:`Kauffman_polynomial`
        """
        kauffman_polynomial = self[self.items.kauffman_polynomial]

        if original:
            return kauffman_polynomial

        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
        R = LaurentPolynomialRing(ZZ, (var1, var2))
        if not kauffman_polynomial and self.crossing_number() == 0:
            return R.one()

        a, z = R.gens()
        lc = {'a':  a, 'z': z}
        return R(eval_knotinfo(kauffman_polynomial, locals=lc))


    @cached_method
    def jones_polynomial(self, variab=None, skein_normalization=False, puiseux=False, original=False, sage_convention=False):
        r"""
        Return the Jones polynomial according to the value of column ``jones_polynomial``
        for this knot or link as an element of the symbolic ring :class:`~sage.symbolic.ring.SR`
        or an instance of :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial`
        depending on the keyword ``skein_normalization``. Using the keyword ``puiseux`` instead
        of an element of the symbolic ring an instance of :class:`~sage.rings.puiseux_series_ring_element.PuiseuxSeries`
        can be returned.

        The Jones polynomial `V(L)` of a link `L` satisfies the following skein relation
        (see the corresponding `KnotInfo description page <https://knotinfo.math.indiana.edu/descriptions/jones_homfly_kauffman_description/polynomial_defn.html)>`__):

        .. MATH::

            V(O) = 1,\,\,\,   t^{-1} V(L_+) -  t V(L_-) = (t^{\frac{1}{2}} - t^{-\frac{1}{2}}) V(L_0)

        INPUT:

        - ``variab`` -- variable (default: ``None``) used according to :meth:`Link.jones_polynomial`
        - ``skein_normalization`` -- boolean (default: ``False``) used according to
          :meth:`Link.jones_polynomial`
        - ``puiseux`` -- boolean (default ``True``) only uses in case ``skein_normalization=False``.
          if set to ``True`` instead of an element of the symbolic ring an instance of
          :class:`~sage.rings.puiseux_series_ring_element.PuiseuxSeries` is returned
        - ``original`` -- boolean (default ``False``) if set to
          ``True`` the original table entry is returned as a string
        - ``sage_convention`` -- boolean (default ``False``) if set to ``True`` the conversion
          to Sage's conventions (see the note below) is performed


        OUTPUT:

        Depends on the keywords (in excluding order):

        - ``original=True`` a string according to the original value from the database
        - ``skein_normalization=True`` a Laurent polynomial over the integers, more precisely
          an instance of :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial`
        - ``puiseux=True`` a puiseux series over the integers, more precisely an instance of
          :class:`~sage.rings.puiseux_series_ring_element.PuiseuxSeries`

        In all other cases an element of the symbolic ring :class:`~sage.symbolic.ring.SR`.

        .. NOTE::

            The only difference of conventions concerning the Jones polynomial is its representation
            in the case of proper links. KnotInfo does not display these polynomials in the indeterminate
            `t` used in the skein relation. Instead a variable `x` is used defined by `x^2 = t`.
            Sage uses `t` in both cases, knots and proper links. Thus, to obtain the Jones polynomial
            for a proper link in `t` you have to set the keyword ``sage_convention`` to ``True``.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K = KnotInfo.K4_1
            sage: Kj = K.jones_polynomial(); Kj
            t^2 - t - 1/t + 1/t^2 + 1
            sage: Kjs = K.jones_polynomial(skein_normalization=True); Kjs
            A^-8 - A^-4 + 1 - A^4 + A^8
            sage: Kjp = K.jones_polynomial(puiseux=True); Kjp
            t^-2 - t^-1 + 1 - t + t^2

        for proper links::

            sage: L = KnotInfo.L2a1_1
            sage: Lj = L.jones_polynomial(); Lj
            -x^5 - x
            sage: Ljt = L.jones_polynomial(sage_convention=True); Ljt
            -t^(5/2) - sqrt(t)
            sage: Ljp = L.jones_polynomial(puiseux=True, sage_convention=True); Ljp
            -t^(1/2) - t^(5/2)
            sage: Ljs = L.jones_polynomial(skein_normalization=True); Ljs
            -A^2 - A^10
            sage: Lj.parent()
            Symbolic Ring
            sage: Ljt.parent()
            Symbolic Ring
            sage: Ljp.parent()
            Puiseux Series Ring in t over Integer Ring
            sage: Ljs.parent()
            Univariate Laurent Polynomial Ring in A over Integer Ring

        Comparison with Sage's results::

            sage: k = K.link()
            sage: kj = k.jones_polynomial()
            sage: bool(Kj == kj)
            True
            sage: kjs = k.jones_polynomial(skein_normalization=True)
            sage: Kjs == kjs
            True
            sage: l = L.link()
            sage: lj = l.jones_polynomial()
            sage: bool(Lj == lj)
            False
            sage: bool(Ljt == lj)   # see note above
            True
            sage: ljs = l.jones_polynomial(skein_normalization=True)
            sage: Ljs == ljs
            True

        Check the skein-relation from the KnotInfo description page (applied to one of
        the positive crossings of the right-handed trefoil)::

            sage: K3_1  = KnotInfo.K3_1
            sage: K3_1j = K3_1.jones_polynomial()
            sage: L2a1_1j = Ljt     # see note above
            sage: R = L2a1_1j.parent()
            sage: Oj = R(1)
            sage: t = R('t')
            sage: lhs = expand(~t*K3_1j - t*Oj)
            sage: rhs = expand((sqrt(t) - ~sqrt(t))*L2a1_1j)
            sage: bool(lhs == rhs)
            True

        The same with the Puiseux series version::

            sage: K3_1jp = K3_1.jones_polynomial(puiseux=True)
            sage: L2a1_1jp = Ljp
            sage: R = L2a1_1jp.parent()
            sage: Ojp = R(1)
            sage: t = R('t')
            sage: ~t*K3_1jp - t*Ojp == (t^(1/2)-~t^(1/2))*L2a1_1jp
            True

        The same in the case of skein normalization (using `t = A^4`)::

            sage: K3_1js = K3_1.jones_polynomial(skein_normalization=True)
            sage: L2a1_1js = L.jones_polynomial(skein_normalization=True)
            sage: Rs = K3_1js.parent()
            sage: Ojs = Rs.one()
            sage: A, = Rs.gens()
            sage: ~A^4*K3_1js - A^4*Ojs == (A^2-~A^2)*L2a1_1js
            True

        REFERENCES:

        - :wikipedia:`Jones_polynomial`
        """
        jones_polynomial = self[self.items.jones_polynomial]

        if original:
            return jones_polynomial

        if skein_normalization:
            if not variab:
                variab='A'
            from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
            R = LaurentPolynomialRing(ZZ, variab)
        else:
            if not variab:
                if sage_convention or self.is_knot():
                    variab='t'
                else:
                    variab='x'
            if puiseux:
                from sage.rings.puiseux_series_ring import PuiseuxSeriesRing  # since PuiseuxPolynomial is not available, so far
                R = PuiseuxSeriesRing(ZZ, variab)
            else:
                from sage.symbolic.ring import SR
                R = SR

        if not jones_polynomial and self.crossing_number() == 0:
            return R(1)

        t = R(variab)
        if skein_normalization:
            if self.is_knot():
                lc = {'t':  t**4}
            else:
                lc = {'x':  t**2}
        else:
            if self.is_knot():
                lc = {'t':  t}
            elif puiseux:
                lc = {'x':  t**(1/2)}
            elif sage_convention:
                from sage.functions.other import sqrt
                lc = {'x':  sqrt(t)}
            else:
                lc = {'x':  t}


        return R(eval_knotinfo(jones_polynomial, locals=lc))


    @cached_method
    def alexander_polynomial(self, var='t', original=False, sage_convention=False):
        r"""
        Return the Alexander polynomial according to the value of column
        ``alexander_polynomial`` for this knot as an instance of
        :class:`~sage.rings.polynomial.polynomial_element.Polynomial`.

        It is obtained from the Seifert matrix `V` of ``self`` by the following
        formula (see the KnotInfo description web-page; to launch it see the
        example below):

        .. MATH::

            A(L) = \det(V -t V^t)

        Here `V^t` stands for the transpose of `V`.


        INPUT:

        - ``var`` -- (default: ``'t'``) the variable
        - ``original`` -- boolean (optional, default ``False``) if set to
          ``True`` the original table entry is returned as a string
        - ``sage_convention`` -- boolean (default ``False``) if set to ``True`` the
          conversion to Sage's conventions (see the note below) is performed

        OUTPUT:

        A polynomial over the integers, more precisely an instance of
        :class:`~sage.rings.polynomial.polynomial_element.Polynomial`.
        If ``sage_convention`` is set to ``True`` a Laurent polynomial
        over the integers, more precisely an instance of
        :class:`~sage.rings.polynomial.laurent_polynomial.LaurentPolynomial`
        is returned. If ``original`` is set to ``True`` then a string
        is returned.

        .. NOTE::

            As an invariant the Alexander polynomial is only unique up to
            a unit factor in the Laurent polynomial ring over the integers
            in the indeterminate `t`. While the normalization of the exponents
            in KnotInfo guarantees it to be a proper polynomial, this is
            not the case for the implementation in Sage. The transition
            can be made using the keyword ``sage_convention``. But still
            there may be a difference in sign (see the example below).

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K = KnotInfo.K4_1
            sage: Ka = K.alexander_polynomial(); Ka
            t^2 - 3*t + 1

        Comparison with Sage's results::

            sage: k = K.link()
            sage: ka = k.alexander_polynomial(); ka
            -t^-1 + 3 - t
            sage: K.alexander_polynomial(sage_convention=True)
            t^-1 - 3 + t
            sage: _ == -ka
            True

        Launch the KnotInfo description web-page::

            sage: K.items.alexander_polynomial.description_webpage()    # not tested
            True
        """
        alexander_polynomial = self[self.items.alexander_polynomial]

        if original:
            return alexander_polynomial

        if sage_convention:
            from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
            R = LaurentPolynomialRing(ZZ, var)
        else:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(ZZ, var)

        if not alexander_polynomial and self.crossing_number() == 0:
            return R.one()

        t, = R.gens()
        lc = {'t':  t}
        ap = R(eval_knotinfo(alexander_polynomial, locals=lc))
        if not sage_convention or ap.is_constant():
            return ap

        exp = ap.exponents()
        return t ** ((-max(exp) - min(exp)) // 2) * ap


    @cached_method
    def link(self, use_item=db.columns().pd_notation, snappy=False):
        r"""
        Return ``self`` as an instance of :class:`Link` or optional
        :class:`~spherogram.links.invariants.Link`  (SnapPy).

        INPUT:

        - ``use_item`` -- (optional default ``self.items.pd_notation``)
          instance of :class:`KnotInfoColumns` to choose the column
          that should be used to construct the link. Allowed values
          are:
          - ``self.items.pd_notation``
          - ``self.items.braid_notation``
          - ``self.items.name``           (only for ``snappy=True``)
          - ``self.items.dt_notation``    (only for knots and ``snappy=False``)
          - ``self.items.gauss_notation`` (only for knots and ``snappy=False``)
        - ``snappy`` boolean (default ``False``) if set to ``True``
          the target of the conversion is the ``pip`` installable
          package ``SnapPy`` (explicitely, :class:`~spherogram.links.invariants.Link`).
          If SnapPy is not installed an ``ImportError`` is raised. To
          install SnapPy use ``sage -pip install snappy``.

        .. NOTE::

            We use the PD-notation to construct ``self`` as
            default. This ensures that the number of crossings
            displayed in representation string of the link
            coincides with the crossing number as a topological
            invariant.

            But attention: The convention on how the edges are
            listed are opposite to each other

            KnotInfo: counter clockwise
            Sage:     clockwise

            Therefore, we take the mirror version of the ``pd_notation``!

            Furthermore, note that the mirror version may depend
            on the used KnotInfo-notation. For example for the
            knot ``5_1`` the Gauss- and the DT-notation refer to
            the mirror image (see example below).

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K = KnotInfo.K3_1
            sage: K.link()
            Knot represented by 3 crossings
            sage: _.braid()
            s^3
            sage: _ == K.braid()
            True

        using ``dt_notation``::

            sage: K.link(use_item=K.items.dt_notation)
            Knot represented by 3 crossings
            sage: _.braid()
            s^3

            sage: L = KnotInfo.L4a1_0
            sage: L.link()
            Link with 2 components represented by 4 crossings

            sage: L.link(use_item=L.items.dt_notation)
            Traceback (most recent call last):
            ...
            ValueError: Link construction using Columns.dt_notation not possible

        using ``snappy``::

            sage: K7   = KnotInfo.K7_2
            sage: k7s  = K7.link(snappy=True); k7s                           # optional - snappy
            <Link: 1 comp; 7 cross>
            sage: k7sn = K7.link(use_item=K7.items.name, snappy=True); k7sn  # optional - snappy
            <Link 7_2: 1 comp; 7 cross>
            sage: k7s == k7sn                                                # optional - snappy
            False
            sage: k7s.sage_link().is_isotopic(k7sn.sage_link())              # optional - snappy
            True

        but observe::

            sage: L2  = KnotInfo.L2a1_1
            sage: l2  = L2.link()
            sage: l2s = L2.link(snappy=True).sage_link()  # optional -  snappy
            sage: l2 == l2s                               # optional -  snappy
            False
            sage: l2 == l2s.mirror_image()                # optional -  snappy
            True

        using ``braid_notation``::

            sage: L2.link(use_item=L.items.braid_notation) == l2
            True

        observe::

            sage: L.link(use_item=L.items.braid_notation)
            Link with 2 components represented by 8 crossings

            sage: K6_1 = KnotInfo.K6_1
            sage: K6_1.link().braid() == K6_1.braid()
            False

        also observe::

            sage: K4_1 = KnotInfo.K4_1
            sage: K4_1.link().pd_code()
            [[4, 1, 5, 2], [8, 5, 1, 6], [6, 4, 7, 3], [2, 8, 3, 7]]
            sage: K4_1.pd_notation()
            [[4, 2, 5, 1], [8, 6, 1, 5], [6, 3, 7, 4], [2, 7, 3, 8]]

            sage: K5_1 = KnotInfo.K5_1
            sage: K5_1.link().braid()
            s^5
            sage: K5_1.link(K5_1.items.dt_notation).braid()
            s^-5
            sage: K5_1.link(K5_1.items.gauss_notation).braid()
            s^-5
        """
        if not isinstance(use_item, KnotInfoColumns):
            raise TypeError('%s must be an instance of %s' %(use_item, KnotInfoColumns))

        if snappy:
            try:
                from snappy import Link
            except ImportError:
                raise ImportError('This option demands snappy to be installed')
        elif self.is_knot():
            from sage.knots.knot import Knot as Link
        else:
            from sage.knots.link import Link

        if   use_item == self.items.pd_notation:
            pd_code = [[a[0], a[3], a[2], a[1]] for a in self.pd_notation()] # take mirror version, see note above
            return Link(pd_code)
        elif use_item == self.items.braid_notation:
            return Link(self.braid())
        elif use_item == self.items.name and snappy:
            if not self.is_knot():
                use_item = self.items.name_unoriented
            return Link(self[use_item])
        elif self.is_knot() and not snappy:
            # Construction via Gauss and DT-Code only possible for knots
            from sage.knots.knot import Knots
            if use_item == self.items.dt_notation:
                return Knots().from_dowker_code(self.dt_notation())
            elif use_item == self.items.gauss_notation:
                return Knots().from_gauss_code(self.gauss_notation())

        raise ValueError('Link construction using %s not possible' %use_item)

    def diagram(self, single=False, new=0, autoraise=True):
        r"""
        Launch the diagram of ``self`` given on the KnotInfo web-page.

        INPUT:

        - ``single`` -- boolean (default ``False``) if set to ``True`` only one
          diagram is shown.
        - ``new`` -- ``int`` according to :func:`open` of :mod:`webbrowser`
          (``0`` default, ``1`` new window, ``2`` new tab)
        - ``autoraise`` -- boolean (default ``True``)

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K = KnotInfo.K3_1
            sage: K.diagram()            # not tested
            True
            sage: K.diagram(single=True) # not tested
            True
        """
        import webbrowser
        if self.is_knot():
            filename = db.filename.knots
        else:
            filename = db.filename.links

        if single:
            return webbrowser.open(filename.diagram_url(self[self.items.diagram], single=single), new=new, autoraise=autoraise)
        else:
            return webbrowser.open(filename.diagram_url(self[self.items.name]), new=new, autoraise=autoraise)


    def knot_atlas_webpage(self, new=0, autoraise=True):
        r"""
        Launch the Knot Atlas web-page for ``self``.

        INPUT:

        - ``new`` -- ``int`` according to :func:`open` of :mod:`webbrowser`
          (``0`` default, ``1`` new window, ``2`` new tab)
        - ``autoraise`` -- boolean (default ``True``)

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K = KnotInfo.K3_1
            sage: K.knot_atlas_webpage()        # not tested
            True
        """
        import webbrowser
        return webbrowser.open(self[self.items.knot_atlas_anon], new=new, autoraise=autoraise)

    def knotilus_webpage(self, new=0, autoraise=True):
        r"""
        Launch the Knotilus web-page for ``self``.

        INPUT:

        - ``new`` -- ``int`` according to :func:`open` of :mod:`webbrowser`
          (``0`` default, ``1`` new window, ``2`` new tab)
        - ``autoraise`` -- boolean (default ``True``)

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: K = KnotInfo.K3_1
            sage: K.knotilus_webpage(new=1)   # not tested
            True
        """
        import webbrowser
        return webbrowser.open(self[self.items.knotilus_page_anon], new=new, autoraise=autoraise)

    def inject(self, verbose=True):
        """
        Inject ``self`` with its name into the namespace of the
        Python code from which this function is called.

        INPUT:

        - ``verbose`` -- boolean (optional default ``True``) to supress
          the message printed on the invocation

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.K5_2.inject()
            Defining K5_2
            sage: K5_2.is_alternating()
            True
        """
        name = self.name
        if verbose:
            print("Defining %s" % (name))
        from sage.repl.user_globals import set_global
        set_global(name, self)



# --------------------------------------------------------------------------------------------
# KnotInfoSeries
# --------------------------------------------------------------------------------------------
class KnotInfoSeries(UniqueRepresentation):
    r"""
    This class can be used to access knots and links via their index
    accordung to the series they belong to.

    INPUT:

    - ``crossing_number`` -- integer giving the crossing numer of this series of links
    - ``is_knot``         -- boolean wether this series is a series of knots or proper links
    - ``is_alternating``  -- boolean wether this series is restriced to alternatimg links or not.
      This is not relevant for knots with less than 11 crossings
    - ``name_unoriented`` -- string restricting the series to all links with that ``name_unoriented``

    EXAMPLES::

        sage: from sage.knots.knotinfo import KnotInfoSeries
        sage: K6 = KnotInfoSeries(6, True, True); K6
        Series of knots K6
        sage: K6(3)
        <KnotInfo.K6_3: [7, 1]>
        sage: list(K6)
        [<KnotInfo.K6_1: [5, 1]>, <KnotInfo.K6_2: [6, 1]>, <KnotInfo.K6_3: [7, 1]>]
        sage: L6a = KnotInfoSeries(6, False, True); L6a
        Series of links L6a
        sage: L6a(2)
        Series of links L6a2
        sage: _.inject()
        Defining L6a2
        sage: list(L6a2)
        [<KnotInfo.L6a2_0: [2986, 2]>, <KnotInfo.L6a2_1: [2987, 2]>]
    """
    def __init__(self, crossing_number, is_knot, is_alternating, name_unoriented=None):
        r"""
        Python constructor.
        
        EXAMPLES::
            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: L6a = KnotInfoSeries(6, False, True); L6a
            Series of links L6a
            sage: TestSuite(L6a).run()
        """
        self._crossing_number   = crossing_number
        self._is_knot           = is_knot
        self._is_alternating    = is_alternating
        self._name_unoriented   = name_unoriented
        self._list              = None

    def list(self):
        r"""
        Return this series as a Python list.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: K6 = KnotInfoSeries(6, True, True); K6
            Series of knots K6
            sage: K6(3)
            <KnotInfo.K6_3: [7, 1]>
        """
        if self._list:
            return self._list

        is_knot  = self._is_knot
        cross_nr = self._crossing_number
        is_alt   = self._is_alternating
        n_unori  = self._name_unoriented

        self._list = []
        curr_n_unori = None
        for K in KnotInfo:
            if K.is_knot() != is_knot:
                continue
            if K.crossing_number() != cross_nr:
                continue
            if not is_knot or cross_nr > 10:
                if K.is_alternating() !=  is_alt:
                    continue
            if is_knot:
                self._list.append(K)
            else:
                this_n_unori = K.name_unoriented()
                if n_unori:
                    if this_n_unori  != n_unori:
                        continue
                    self._list.append(K)
                elif this_n_unori  != curr_n_unori:
                    if curr_n_unori:
                        self._list.append(KnotInfoSeries(cross_nr, is_knot, is_alt, curr_n_unori))
                    curr_n_unori = this_n_unori
                else:
                    continue

        if curr_n_unori:
            self._list.append(KnotInfoSeries(cross_nr, is_knot, is_alt, curr_n_unori))
        return self._list


    def __repr__(self):
        r"""
        Return the representation string of ``self``.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: KnotInfoSeries(6, True, True)
            Series of knots K6
            sage: _.__repr__()
            'Series of knots K6'
        """
        if self._is_knot:
            return 'Series of knots %s' %(self._name())
        else:
            return 'Series of links %s' %(self._name())


    def __getitem__(self, item):
        r"""
        Return the given ``item`` from the list of ``self``
        (making the Python build-in ``list`` work).

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: KnotInfoSeries(6, True, True).inject()
            Defining K6
            sage: list(K6)                      # indirect doctest
            [<KnotInfo.K6_1: [5, 1]>, <KnotInfo.K6_2: [6, 1]>, <KnotInfo.K6_3: [7, 1]>]
        """
        from sage.rings.integer import Integer
        if  not type(item) in (int, Integer):
            raise ValueError('Item must be an integer')
        l =self.list()
        max_item = len(l)
        if item < 0 or item  > max_item:
            raise ValueError('Item must be non negative and smaller than %s' %(max_item))

        return l[item]

    def __call__(self, item):
        r"""
        Return the given ``item`` from the list of ``self``
        (making the function call for ``self`` work).
        In contrast to ``__getitem__`` the first ``item``
        has to be ``1`` (not ``0``).

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: KnotInfoSeries(6, True, True).inject()
            Defining K6
            sage: K6(2)                         # indirect doctest
            <KnotInfo.K6_2: [6, 1]>
        """
        from sage.rings.integer import Integer
        if  not type(item) in (int, Integer):
            raise ValueError('Item must be an integer')
        l =self.list()
        max_item = len(l)+1
        if item < 1 or item  > max_item:
            raise ValueError('Item must be positive and smaller than %s' %(max_item))

        return l[item-1]

    def _name(self):
        r"""
        Return the name of the series.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: KnotInfoSeries(6, True, True)._name()
            'K6'
        """
        is_knot  = self._is_knot
        cross_nr = self._crossing_number
        is_alt   =  self._is_alternating
        n_unori  = self._name_unoriented

        alt = 'a'
        if not is_alt:
            alt = 'n'

        if is_knot:
            if cross_nr > 10:
                res = 'K%s%s' %(cross_nr, alt)
            else:
                res = 'K%s' %(cross_nr)
        elif n_unori:
            res = '%s' %(n_unori)
        else:
            res = 'L%s%s' %(cross_nr, alt)
        return res


    def inject(self, verbose=True):
        """
        Inject ``self`` with its name into the namespace of the
        Python code from which this function is called.

        INPUT:

        - ``verbose`` -- boolean (optional default ``True``) to supress
          the message printed on the invocation

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: KnotInfoSeries(6, True, True).inject()
            Defining K6
            sage: K6(2)
            <KnotInfo.K6_2: [6, 1]>
        """
        name = self._name()
        if verbose:
            print("Defining %s" % (name))
        from sage.repl.user_globals import set_global
        set_global(name, self)




KnotInfo = KnotInfoBase('KnotInfo', db.read_row_dict())
