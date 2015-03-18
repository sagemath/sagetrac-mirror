r"""
Rigged Configurations of `\mathcal{B}(\infty)`

AUTHORS:

- Travis Scrimshaw (2013-04-16): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
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

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.rigged_configurations.rigged_configuration_element import (
     RiggedConfigurationElement, RCNonSimplyLacedElement)
from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurationOptions

# Note on implementation, this class is used for simply-laced types only
class InfinityCrystalOfRiggedConfigurations(Parent, UniqueRepresentation):
    r"""
    Rigged configuration model for `\mathcal{B}(\infty)`.

    The crystal is generated by the empty rigged configuration with the same
    crystal structure given by the :class:`highest weight model
    <sage.combinat.rigged_configurations.rc_crystal.CrystalOfRiggedConfigurations>`
    except we remove the condition that the resulting rigged configuration
    needs to be valid when applying `f_a`.

    INPUT:

    - ``cartan_type`` -- a Cartan type

    EXAMPLES:

    For simplicity, we display all of the rigged configurations
    horizontally::

        sage: RiggedConfigurations.global_options(display='horizontal')

    We begin with a simply-laced finite type::

        sage: RC = crystals.infinity.RiggedConfigurations(['A', 3]); RC
        The infinity crystal of rigged configurations of type ['A', 3]

        sage: RC.global_options(display='horizontal')

        sage: mg = RC.highest_weight_vector(); mg
        (/)  (/)  (/)
        sage: elt = mg.f_string([2,1,3,2]); elt
        0[ ]0   -2[ ]-1   0[ ]0
                -2[ ]-1
        sage: elt.e(1)
        sage: elt.e(3)
        sage: mg.f_string([2,1,3,2]).e(2)
        -1[ ]-1  0[ ]1  -1[ ]-1
        sage: mg.f_string([2,3,2,1,3,2])
        0[ ]0  -3[ ][ ]-1  -1[ ][ ]-1
               -2[ ]-1

    Next we consider a non-simply-laced finite type::

        sage: RC = crystals.infinity.RiggedConfigurations(['C', 3])
        sage: mg = RC.highest_weight_vector()
        sage: mg.f_string([2,1,3,2])
        0[ ]0   -2[ ]-1   0[ ]0
                -2[ ]-1
        sage: mg.f_string([2,3,2,1,3,2])
        -1[ ]-1  -2[ ][ ][ ]-1  -1[ ][ ]0

    We can also construct `B(\infty)` using rigged configurations in
    affine types::

        sage: RC = crystals.infinity.RiggedConfigurations(['A', 3, 1])
        sage: mg = RC.highest_weight_vector()
        sage: mg.f_string([0,1,2,3,0,1,3])
        -1[ ]0  -1[ ]-1  1[ ]1  -1[ ][ ]-1
        -1[ ]0  -1[ ]-1

        sage: RC = crystals.infinity.RiggedConfigurations(['C', 3, 1])
        sage: mg = RC.highest_weight_vector()
        sage: mg.f_string([1,2,3,0,1,2,3,3,0])
        -2[ ][ ]-1  -1[ ]0  -1[ ]-1  -4[ ][ ][ ]-2
                    -1[ ]0  -1[ ]-1

        sage: RC = crystals.infinity.RiggedConfigurations(['A', 6, 2])
        sage: mg = RC.highest_weight_vector()
        sage: mg.f_string([1,2,3,0,1,2,3,3,0])
        0[ ]-1  0[ ]1  -1[ ]-1  -4[ ][ ][ ]-2
        0[ ]-1  0[ ]1  -1[ ]-1

    We reset the global options::

        sage: RiggedConfigurations.global_options.reset()
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type):
        r"""
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::

            sage: RC1 = crystals.infinity.RiggedConfigurations(CartanType(['A',3]))
            sage: RC2 = crystals.infinity.RiggedConfigurations(['A',3])
            sage: RC2 is RC1
            True
        """
        cartan_type = CartanType(cartan_type)
        if not cartan_type.is_simply_laced():
            vct = cartan_type.as_folding()
            return InfinityCrystalOfNonSimplyLacedRC(vct)

        return super(InfinityCrystalOfRiggedConfigurations, cls).__classcall__(cls, cartan_type)

    def __init__(self, cartan_type):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['A', 2])
            sage: TestSuite(RC).run()
            sage: RC = crystals.infinity.RiggedConfigurations(['A', 2, 1])
            sage: TestSuite(RC).run() # long time
        """
        self._cartan_type = cartan_type
        Parent.__init__(self, category=HighestWeightCrystals().Infinite())
        self._rc_index = self._cartan_type.index_set()
        # We store the cartan matrix for the vacancy number calculations for speed
        self._cartan_matrix = self._cartan_type.cartan_matrix()
        self.module_generators = (self.element_class(self, rigging_list=[[]]*cartan_type.rank()),)

    global_options = RiggedConfigurationOptions

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: crystals.infinity.RiggedConfigurations(['A', 3])
            The infinity crystal of rigged configurations of type ['A', 3]
        """
        return "The infinity crystal of rigged configurations of type {}".format(self._cartan_type)

    def _element_constructor_(self, lst=None, **options):
        """
        Construct an element of ``self`` from ``lst``.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['A', 3, 1])
            sage: ascii_art(RC(partition_list=[[1,1]]*4, rigging_list=[[1,1], [0,0], [0,0], [-1,-1]]))
            0[ ]1  0[ ]0  0[ ]0  0[ ]-1
            0[ ]1  0[ ]0  0[ ]0  0[ ]-1

            sage: RC = crystals.infinity.RiggedConfigurations(['C', 3])
            sage: ascii_art(RC(partition_list=[[1],[1,1],[1]], rigging_list=[[0],[-1,-1],[0]]))
            0[ ]0  -2[ ]-1  0[ ]0
                   -2[ ]-1
        """
        return self.element_class(self, lst, **options)

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number of the `i`-th row of the `a`-th rigged
        partition.

        This assumes that `\gamma_a = 1` for all `a` and `(\alpha_a \mid
        \alpha_b ) = A_{ab}`.

        INPUT:

        - ``partitions`` -- the list of rigged partitions we are using

        - ``a`` -- the rigged partition index

        - ``i`` -- the row index of the `a`-th rigged partition

        TESTS::

            sage: RC = crystals.infinity.RiggedConfigurations(['A', 4])
            sage: elt = RC(partition_list=[[1], [1], [], []])
            sage: RC._calc_vacancy_number(elt.nu(), 0, 0)
            -1
        """
        vac_num = 0

        if i is None:
            row_len = float("inf")
        else:
            row_len = partitions[a][i]

        for b, value in enumerate(self._cartan_matrix.row(a)):
            vac_num -= value * partitions[b].get_num_cells_to_column(row_len)

        return vac_num

    class Element(RiggedConfigurationElement):
        """
        A rigged configuration in `\mathcal{B}(\infty)` in simply-laced types.

        TESTS::

            sage: RC = crystals.infinity.RiggedConfigurations(['A', 3, 1])
            sage: elt = RC(partition_list=[[1,1]]*4, rigging_list=[[1,1], [0,0], [0,0], [-1,-1]])
            sage: TestSuite(elt).run()
        """
        def weight(self):
            """
            Return the weight of ``self``.

            EXAMPLES::

                sage: RC = crystals.infinity.RiggedConfigurations(['A', 3, 1])
                sage: elt = RC(partition_list=[[1,1]]*4, rigging_list=[[1,1], [0,0], [0,0], [-1,-1]])
                sage: elt.weight()
                0
            """
            P = self.parent().weight_lattice_realization()
            alpha = list(P.simple_roots())
            return sum(sum(x) * alpha[i] for i,x in enumerate(self))

class InfinityCrystalOfNonSimplyLacedRC(InfinityCrystalOfRiggedConfigurations):
    r"""
    Rigged configurations for `\mathcal{B}(\infty)` in non-simply-laced types.
    """
    def __init__(self, vct):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['C', 2])
            sage: TestSuite(RC).run() # long time
            sage: RC = crystals.infinity.RiggedConfigurations(['C', 2, 1])
            sage: TestSuite(RC).run() # long time
         """
        self._folded_ct = vct
        InfinityCrystalOfRiggedConfigurations.__init__(self, vct._cartan_type)

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number of the `i`-th row of the `a`-th rigged
        partition.

        INPUT:

        - ``partitions`` -- the list of rigged partitions we are using

        - ``a`` -- the rigged partition index

        - ``i`` -- the row index of the `a`-th rigged partition

        TESTS::

            sage: La = RootSystem(['C',2]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(La[1])
            sage: elt = RC(partition_list=[[1], [1]])
            sage: RC._calc_vacancy_number(elt.nu(), 0, 0)
            0
            sage: RC._calc_vacancy_number(elt.nu(), 1, 0)
            -1
        """
        I = self.index_set()
        ia = I[a]
        vac_num = 0

        if i is None:
            row_len = float("inf")
        else:
            row_len = partitions[a][i]

        gamma = self._folded_ct.scaling_factors()
        for b, value in enumerate(self._cartan_matrix.row(a)):
            ib = I[b]
            q = partitions[b].get_num_cells_to_column(gamma[ia]*row_len, gamma[ib])
            vac_num -= value * q / gamma[ib]

        return vac_num

    @lazy_attribute
    def virtual(self):
        """
        Return the corresponding virtual crystal.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['C', 3])
            sage: RC
            The infinity crystal of rigged configurations of type ['C', 3]
            sage: RC.virtual
            The infinity crystal of rigged configurations of type ['A', 5]
        """
        return InfinityCrystalOfRiggedConfigurations(self._folded_ct._folding)

    def to_virtual(self, rc):
        """
        Convert ``rc`` into a rigged configuration in the virtual crystal.

        INPUT:

        - ``rc`` -- a rigged configuration element

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['C', 2])
            sage: mg = RC.highest_weight_vector()
            sage: elt = mg.f_string([1,2,2,1,1]); elt
            <BLANKLINE>
            -3[ ][ ][ ]-2
            <BLANKLINE>
            -1[ ][ ]0
            <BLANKLINE>
            sage: velt = RC.to_virtual(elt); velt
            <BLANKLINE>
            -3[ ][ ][ ]-2
            <BLANKLINE>
            -2[ ][ ][ ][ ]0
            <BLANKLINE>
            -3[ ][ ][ ]-2
            <BLANKLINE>
            sage: velt.parent()
            The infinity crystal of rigged configurations of type ['A', 3]
        """
        gamma = map(int, self._folded_ct.scaling_factors())
        sigma = self._folded_ct._orbit
        n = self._folded_ct._folding.rank()
        vindex = self._folded_ct._folding.index_set()
        partitions = [None] * n
        riggings = [None] * n
        for a, rp in enumerate(rc):
            for i in sigma[a]:
                k = vindex.index(i)
                partitions[k] = [row_len*gamma[a] for row_len in rp._list]
                riggings[k] = [rig_val*gamma[a] for rig_val in rp.rigging]
        return self.virtual.element_class(self.virtual, partition_list=partitions,
                            rigging_list=riggings)

    def from_virtual(self, vrc):
        """
        Convert ``vrc`` in the virtual crystal into a rigged configution of
        the original Cartan type.

        INPUT:

        - ``vrc`` -- a virtual rigged configuration

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['C', 2])
            sage: elt = RC(partition_list=[[3],[2]], rigging_list=[[-2],[0]])
            sage: vrc_elt = RC.to_virtual(elt)
            sage: ret = RC.from_virtual(vrc_elt); ret
            <BLANKLINE>
            -3[ ][ ][ ]-2
            <BLANKLINE>
            -1[ ][ ]0
            <BLANKLINE>
            sage: ret == elt
            True
        """
        gamma = list(self._folded_ct.scaling_factors()) #map(int, self._folded_ct.scaling_factors())
        sigma = self._folded_ct._orbit
        n = self._cartan_type.rank()
        partitions = [None] * n
        riggings = [None] * n
        vindex = self._folded_ct._folding.index_set()
        # TODO: Handle special cases for A^{(2)} even and its dual?
        for a in range(n):
            index = vindex.index(sigma[a][0])
            partitions[a] = [row_len // gamma[a] for row_len in vrc[index]._list]
            riggings[a] = [rig_val / gamma[a] for rig_val in vrc[index].rigging]
        return self.element_class(self, partition_list=partitions,
                                  rigging_list=riggings)

    class Element(RCNonSimplyLacedElement):
        """
        A rigged configuration in `\mathcal{B}(\infty)` in
        non-simply-laced types.

        TESTS::

            sage: RC = crystals.infinity.RiggedConfigurations(['C', 3])
            sage: elt = RC(partition_list=[[1],[1,1],[1]])
            sage: TestSuite(elt).run()
        """
        def weight(self):
            """
            Return the weight of ``self``.

            EXAMPLES::

                sage: RC = crystals.infinity.RiggedConfigurations(['C', 3])
                sage: elt = RC(partition_list=[[1],[1,1],[1]], rigging_list=[[0],[-1,-1],[0]])
                sage: elt.weight()
                (1, 1, 0)

                sage: RC = crystals.infinity.RiggedConfigurations(['F', 4, 1])
                sage: mg = RC.highest_weight_vector()
                sage: elt = mg.f_string([1,0,3,4,2,2]); ascii_art(elt)
                -1[ ]-1  0[ ]1  -2[ ][ ]-2  0[ ]1  -1[ ]-1
                sage: wt = elt.weight(); wt
                Lambda[0] - Lambda[1] + 2*Lambda[2] - 3*Lambda[3] + Lambda[4]
                sage: al = RC.weight_lattice_realization().simple_roots()
                sage: wt == al[0] + al[1] + 2*al[2] + al[3] + al[4]
                True
            """
            P = self.parent().weight_lattice_realization()
            alpha = list(P.simple_roots())
            return sum(sum(x) * alpha[i] for i,x in enumerate(self))

