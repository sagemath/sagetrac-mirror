r"""
Rigged Configuration Bijection By Guessing

We guess what the bijection is by matching highest weight elements and the
classical crystal they generate.

AUTHORS:

- Travis Scrimshaw (2013-01-10): Initial version
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

from sage.combinat.crystals.tensor_product import FullTensorProductOfRegularCrystals

# Technically there should be an interface which the bijection abstract classes
# in :mod:`sage.combinat.rigged_configurations.bij_abstract_class` and these
# classes should inherit from. However since the only public method is
# ``run()``, we don't really care about type-checking, and python is a
# weakly-typed language, we don't need it.

def compare_classical_crystals(Lx, Ly):
    """
    Heuristic comparison of the classical crystals ``Lx`` and ``Ly``.

    ALGORITHM:

    #. Compare the size of the crystals, if not equal, then they are
       not equal.
    #. If the component has length at most 1, then they are trivially equal.
    #. Loop through ``x in Lx``

       #. Final all elements ``y in Ly`` which has the same classical weight
          as ``x`` and store in a list ``L``.
       #. See if there is an element in ``L`` which has the same `\epsilon_i`
          and `\phi_i` for all `i \in \overline{I}`. If not, then they are not
          equal. If so, remove ``y`` from ``Ly``.

    #. If completed, then we will believe that ``Lx == Ly``
       as classical crystals.

    EXAMPLES::

        sage: from sage.combinat.rigged_configurations.bij_guess import compare_classical_crystals
        sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
        sage: KR = KirillovReshetikhinTableaux(['D', 4, 1], 2, 1)
        sage: mg_rc = RC(partition_list=[[],[],[],[]])
        sage: mg_kr = KR(2,1)
        sage: S = mg_rc.subcrystal(index_set=[1,2,3,4], direction='lower')
        sage: T = mg_kr.subcrystal(index_set=[1,2,3,4], direction='lower')
        sage: compare_classical_crystals(S, T)
        True
    """
    # Make sure they are lists
    Lx = list(Lx)
    Ly = list(Ly)

    if len(Lx) != len(Ly):
        return False
    if len(Lx) <= 1:
        return True

    I = Lx[0].parent()._cartan_type.classical().index_set()
    for x in Lx:
        L = [y for y in Ly if y.classical_weight() == x.classical_weight()]
        match = None
        for y in L:
            for i in I:
                if x.epsilon(i) == y.epsilon(i) and x.phi(i) == y.phi(i):
                    match = y
                    break
            if match is not None:
                break
        if match is None:
            return False
        Ly.remove(y)
    return True

class KRTToRCBijectionGuess:
    """
    The bijection from KR tableaux to rigged configurations by guessing
    based upon highest weight and classical crystal structure.

    This class holds the state of the bijection and generates the next state.
    This class should never be created directly.
    """
    def __init__(self, tp_krt):
        """
        Initialize the bijection.

        INPUT:

        - ``parent`` -- the parent of tensor product of KR tableaux

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['E', 6, 1], [[1,1]])
            sage: from sage.combinat.rigged_configurations.bij_guess import KRTToRCBijectionGuess
            sage: krt = KRT.module_generators[0].f_string([1,3,4])
            sage: bijection = KRTToRCBijectionGuess(krt)
            sage: TestSuite(bijection).run()
        """
        self.tp_krt = tp_krt
        self.rc = tp_krt.parent().rigged_configurations()
        self.n = tp_krt.cartan_type().classical().rank()

    def __eq__(self, rhs):
        r"""
        Check equality.

        This is only here for pickling check. This is a temporary placeholder
        class, and as such, should never be compared.

        TESTS::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['E', 6, 1], [[1,1]])
            sage: from sage.combinat.rigged_configurations.bij_guess import KRTToRCBijectionGuess
            sage: krt = KRT.module_generators[0].f_string([1,3,4])
            sage: bijection = KRTToRCBijectionGuess(krt)
            sage: bijection2 = KRTToRCBijectionGuess(krt)
            sage: bijection == bijection2
            True
        """
        return isinstance(rhs, KRTToRCBijectionGuess)

    def run(self, verbose=False):
        """
        Run the bijection from a tensor product of KR tableaux to a rigged
        configuration.

        INPUT:

        - ``tp_krt`` -- a tensor product of KR tableaux

        - ``verbose`` -- (default: ``False``) display each step in the
          bijection

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['E', 6, 1], [[1, 1]])
            sage: from sage.combinat.rigged_configurations.bij_guess import KRTToRCBijectionGuess
            sage: krt = KRT.module_generators[0].f_string([1,3,4])
            sage: KRTToRCBijectionGuess(krt).run()
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
        """
        hw = self.tp_krt.to_highest_weight(index_set=self.tp_krt.cartan_type().classical().index_set())
        wt = hw[0].classical_weight()
        L = [x for x in self.rc.module_generators if x.classical_weight() == wt]
        # If there is exactly one other highest weight element with
        #   the same classical weight, we must map there
        if len(L) == 1:
            return L[0].f_string(reversed(hw[1]))

        index_set = hw[0].cartan_type().classical().index_set()
        classical_crystal = hw[0].subcrystal(index_set=index_set, direction='lower')
        for x in L:
            other = x.subcrystal(index_set=index_set, direction='lower')
            if compare_classical_crystals(other, classical_crystal):
                return x.f_string(reversed(hw[1]))
        return None

class RCToKRTBijectionGuess:
    """
    The bijection from KR tableaux to rigged configurations by guessing
    based upon highest weight and classical crystal structure.

    This class holds the state of the bijection and generates the next state.
    This class should never be created directly.
    """
    def __init__(self, RC_element):
        """
        Initialize the bijection helper.

        INPUT:

        - ``RC_element`` -- the rigged configuration

        EXAMPLES::

            sage: RC = RiggedConfigurations(['E', 6, 1], [[1, 1]])
            sage: from sage.combinat.rigged_configurations.bij_guess import RCToKRTBijectionGuess
            sage: bijection = RCToKRTBijectionGuess(RC(partition_list=[[1],[],[1],[1],[],[]]))
            sage: TestSuite(bijection).run()
        """
        self.rigged_con = RC_element.__copy__()
        self.n = RC_element.parent().cartan_type().n

    def __eq__(self, rhs):
        r"""
        Check equality.

        This is only here for pickling check. This is a temporary placeholder
        class, and as such, should never be compared.

        TESTS::

            sage: RC = RiggedConfigurations(['E', 6, 1], [[1, 1]])
            sage: from sage.combinat.rigged_configurations.bij_guess import RCToKRTBijectionGuess
            sage: bijection = RCToKRTBijectionGuess(RC(partition_list=[[1],[],[1],[1],[],[]]))
            sage: bijection2 = RCToKRTBijectionGuess(RC(partition_list=[[1],[],[1],[1],[],[]]))
            sage: bijection == bijection2
            True
        """
        return isinstance(rhs, RCToKRTBijectionGuess)

    def run(self, verbose=False):
        """
        Run the bijection from rigged configurations to tensor product of KR
        tableaux.

        INPUT:

        - ``verbose`` -- (default: ``False``) display each step in the
          bijection

        EXAMPLES::

            sage: RC = RiggedConfigurations(['E', 6, 1], [[1, 1]])
            sage: from sage.combinat.rigged_configurations.bij_guess import RCToKRTBijectionGuess
            sage: RCToKRTBijectionGuess(RC(partition_list=[[1],[],[1],[1],[],[]])).run()
            [[(-4, 2, 5)]]
        """
        index_set = self.rigged_con.cartan_type().classical().index_set()
        hw = self.rigged_con.to_highest_weight(index_set)
        KRT = self.rigged_con.parent().tensor_product_of_Kirillov_Reshetikhin_tableaux()

        classical_weight = lambda x : sum([factor.classical_weight() for factor in x])
        wt = hw[0].classical_weight()
        L = filter(lambda x: classical_weight(x) == wt, KRT._module_generators_brute_force())
        # If there is exactly one other highest weight element with
        #   the same classical weight, we must map there
        if len(L) == 1:
            return L[0].f_string(reversed(hw[1]))

        # Now check to make sure the usual weights agree
        wt = hw[0].weight()
        L = filter(lambda x: x.weight() == wt, L)
        if len(L) == 1:
            return L[0].f_string(reversed(hw[1]))

        classical_crystal = hw[0].subcrystal(index_set=index_set, direction='lower')
        for x in L:
            other = x.subcrystal(index_set=index_set, direction='lower')
            if compare_classical_crystals(other, classical_crystal):
                return x.f_string(reversed(hw[1]))
        return None

