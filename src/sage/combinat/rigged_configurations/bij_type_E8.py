r"""
Bijection classes for type `E_8^{(1)}`

Part of the (internal) classes which runs the bijection between rigged
configurations and KR tableaux of type `E_8^{(1)}`.

AUTHORS:

- Travis Scrimshaw (2016-02): Initial version

TESTS::

    sage: from sage.combinat.rigged_configurations.bij_type_E8 import KRTToRCBijectionTypeE8
    sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['E', 8, 1], [[8,1]])
    sage: bijection = KRTToRCBijectionTypeE8(KRT.module_generators[0])
    sage: TestSuite(bijection).run()
    sage: RC = RiggedConfigurations(['E', 8, 1], [[8, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_E8 import RCToKRTBijectionTypeE8
    sage: bijection = RCToKRTBijectionTypeE8(RC(partition_list=[.....]))
    sage: TestSuite(bijection).run()
"""

#*****************************************************************************
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh umn.edu>
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

from sage.combinat.rigged_configurations.bij_abstract_class import KRTToRCBijectionAbstract
from sage.combinat.rigged_configurations.bij_abstract_class import RCToKRTBijectionAbstract
from sage.combinat.rigged_configurations.bij_type_adjoint import next_state_RC_KRT
from sage.combinat.crystals.letters import CrystalOfLetters
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method

class KRTToRCBijectionTypeE8(KRTToRCBijectionAbstract):
    r"""
    Specific implementation of the bijection from KR tableaux to rigged
    configurations for type `E_8^{(1)}`.
    """
    def next_state(self, val):
        r"""
        Build the next state for type `E_{6,7}^{(1)}`.

        TESTS::

            sage: from sage.combinat.rigged_configurations.bij_type_E8 import KRTToRCBijectionTypeE8
            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['E', 8, 1], [[8,1]])
            sage: bijection = KRTToRCBijectionTypeE8(KRT.module_generators[0])
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [1, 1])
            sage: bijection.cur_path[0].insert(0, [(1,)])
            sage: bijection.next_state((1,))
        """

    def _next_index(self, r, target):
        """
        Return the next index after ``r`` when performing a step
        in the bijection going towards ``target``.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['E', 6, 1], [[5,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_E67 import KRTToRCBijectionTypeE67
            sage: bijection = KRTToRCBijectionTypeE67(KRT.module_generators[0])
            sage: bijection._next_index(3, 5)
            2
            sage: bijection._next_index(2, 5)
            5
            sage: bijection._next_index(3, 4)
            4
            sage: bijection._next_index(1, 5)
            3
            sage: bijection._next_index(1, 4)
            3
            sage: bijection._next_index(1, 6)
            6
        """
        #     1-2-3
        #    /
        # 0-8-7-6-5-4
        if r == 0:
            return 8
        if r == 8:
            if target <= 3:
                return 1
            return 6
        if r <= 3:
            return r + 1
        # r = 7,6,5
        return r - 1

class RCToKRTBijectionTypeE8(RCToKRTBijectionAbstract):
    r"""
    Specific implementation of the bijection from rigged configurations
    to tensor products of KR tableaux for type `E_8^{(1)}`.
    """
    def next_state(self, r):
        r"""
        Build the next state for type `E_8^{(1)}`.

        TESTS::

            sage: RC = RiggedConfigurations(['E', 8, 1], [[6, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_E8 import RCToKRTBijectionTypeE8
            sage: bijection = RCToKRTBijectionTypeE8(RC(partition_list=[[1],[1,1],[1,1],[1,1,1],[1,1],[1]]))
            sage: bijection.next_state(5)
            (-2, 1)

        This returns the wrong results for the lower half::

            sage: RC = RiggedConfigurations(['E',8,1], [[6,1]])
            sage: mg = RC(partition_list=[[2,2,1,1],[2,2,2,1,1,1],[2,2,2,2,1,1,1,1],[2,2,2,2,2,2,1,1,1,1,1,1],[2,2,2,2,2,1,1,1,1,1],[2,2,2,2,1,1,1,1],[2,2,2,1,1],[2,2]])
            sage: mg.to_tensor_product_of_kirillov_reshetikhin_tableaux()
            [[(8,)], [E], [(-8,)]]

[[(8,)], [(2, 3, -4)], [(1, -2, -3, 4, -8)]]
0[ ]0  0[ ][ ]0  0[ ][ ]0  0[ ][ ]0  0[ ][ ]0  0[ ][ ]0  0[ ][ ]0  0[ ][ ]0
0[ ]0  0[ ]0     0[ ]0     0[ ][ ]0  0[ ][ ]0  0[ ][ ]0  0[ ][ ]0
       0[ ]0     0[ ]0     0[ ]0     0[ ]0     0[ ]0             
                 0[ ]0     0[ ]0     0[ ]0     0[ ]0             
                           0[ ]0     0[ ]0                       
                           0[ ]0                                 
-2*Lambda[0] + Lambda[1]

            sage: mg = RC(partition_list=[[2,2],[2,2,1,1],[2,2,2,1,1],[2,2,2,2,1,1,1,1],[2,2,2,1,1,1,1],[2,2,1,1,1,1],[2,1,1,1],[1,1]])

        """
        # Perform left_column_box
        b = self._endpoint(r)
        I = b.parent().index_set()
        P = self.rigged_con.parent()
        path = [I.index(a) for a in b.to_highest_weight()[1]]
        for ii in path:
            self.cur_partitions[ii]._list.append(1)
        self.cur_dims.insert(0, [8,1])
        for ii in path:
            vac_num = P._calc_vacancy_number(self.cur_partitions, ii, 1,
                                             dims=self.cur_dims)
            pos = self.cur_partitions[ii]._list.index(1)
            self.cur_partitions[ii].vacancy_numbers.insert(pos, vac_num)
            self.cur_partitions[ii].rigging.insert(pos, vac_num)

        # Perform \delta
        self.cur_dims.pop(0)
        return next_state_RC_KRT(self, b.parent())

    def _next_index(self, r):
        """
        Return the next index after ``r`` when performing a step
        in the bijection.

        TESTS::

            sage: RC = RiggedConfigurations(['E', 8, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_E8 import RCToKRTBijectionTypeE8
            sage: bijection = RCToKRTBijectionTypeE67(RC(partition_list=[[1],[1,1],[1,1],[1,1,1], [1,1],[1]]))
            sage: bijection._next_index(2)
            3
        """
        #     1-2-3
        #    /
        # 0-8-7-6-5-4
        if r == 1:
            return 8
        if r == 8:
            return 0
        if r <= 3:
            return r - 1
        # r = 4,5,6,7
        return r + 1

    @cached_method
    def _endpoint(self, r):
        r"""
        Return the endpoint for the bijection in type `E_8^{(1)}`.

        EXAMPLES::

            sage: from sage.combinat.rigged_configurations.bij_type_E8 import RCToKRTBijectionTypeE8, endpoint8
            sage: RC = RiggedConfigurations(['E', 8, 1], [[2, 1]])
            sage: bijection = RCToKRTBijectionTypeE8(RC(partition_list=[[1],[1,1],[1,1],[1,1,1], [1,1],[1]]))
            sage: all(bijection._endpoint(r) == endpoint8(r) for r in range(1,7))
            True
        """
        return endpoint8(r)

def endpoint8(r):
    r"""
    Return the endpoint for `B^{r,1}` in type `E_8^{(1)}`.

    EXAMPLES::

        sage: from sage.combinat.rigged_configurations.bij_type_E8 import endpoint8
        sage: endpoint8(1)
        (1, -8)
        sage: endpoint8(2)
        (-1, 2)
        sage: endpoint8(3)
        (-2, 3)
        sage: endpoint8(4)
        (4, -5)
        sage: endpoint8(5)
        (5, -6)
        sage: endpoint8(6)
        (6, -7)
        sage: endpoint8(7)
        (7, -8)
        sage: endpoint8(8)
        (8,)
    """
    C = CrystalOfLetters(['E',8])
    if r == 1:
        return C((1, -8))
    elif r == 2:
        return C((-1, 2))
    elif r == 3:
        return C((-2, 3))
    elif r == 4:
        return C((4, -5))
    elif r == 5:
        return C((5, -6))
    elif r == 6:
        return C((6, -7))
    elif r == 7:
        return C((7, -8))
    elif r == 8:
        return C.module_generators[0]  # C((8,))

