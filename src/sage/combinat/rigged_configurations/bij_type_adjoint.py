r"""
Bijection functions for the adjoint node.

Part of the (internal) functions which runs the bijection between rigged
configurations and KR tableaux. This adds/removes the factor `B^{r,1}`,
where `r` is the node bonded to `0` in the Dynkin diagram.

AUTHORS:

- Travis Scrimshaw (2015-07-21): Initial version

TESTS::

    sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 2], [[2,1]])
    sage: from sage.combinat.rigged_configurations.bij_type_D_twisted import KRTToRCBijectionTypeDTwisted
    sage: bijection = KRTToRCBijectionTypeDTwisted(KRT(pathlist=[[-1,2]]))
    sage: TestSuite(bijection).run()
    sage: RC = RiggedConfigurations(['D', 4, 2], [[2, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_D_twisted import RCToKRTBijectionTypeDTwisted
    sage: bijection = RCToKRTBijectionTypeDTwisted(RC())
    sage: TestSuite(bijection).run()
"""

#*****************************************************************************
#       Copyright (C) 2015 Travis Scrimshaw <tscrim@ucdavis.edu>
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

def insert_cell_case_S(partition):
    """
    Insert a cell when case `(S)` holds.

    TESTS::

        sage: RC = RiggedConfigurations(['C', 2, 1], [[2, 2]])
        sage: RP = RC(partition_list=[[2],[2,2]])[1]
        sage: RP
        -4[ ][ ]-4
        -4[ ][ ]-4
        <BLANKLINE>
        sage: RP.rigging[0] = None
        sage: from sage.combinat.rigged_configurations.bij_type_adjoint import insert_cell_case_S
        sage: insert_cell_case_S(RP)
        sage: RP
        -4[ ][ ][ ]None
        -4[ ][ ]-4
        <BLANKLINE>
    """
    # Special case when adding twice to the first row
    if partition.rigging[0] is None:
        partition._list[0] += 1
        return

    num_rows = len(partition)
    for i in reversed(range(1, num_rows)):
        if partition.rigging[i] is None:
            j = i - 1
            while j >= 0 and partition._list[j] == partition._list[i]:
                partition.rigging[j+1] = partition.rigging[j] # Shuffle it along
                j -= 1
            partition._list[j+1] += 1
            partition.rigging[j+1] = None
            return
    assert False, "Bug in the bijection"

def next_state_KRT_RC(bij, letters):
    r"""
    Build the next state for the left factor `B^{r,1}`, where `r`
    is the node ajoint to node `0`.

    TESTS::

        sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 2], [[2,1]])
        sage: from sage.combinat.rigged_configurations.bij_type_D_twisted import KRTToRCBijectionTypeDTwisted
        sage: bijection = KRTToRCBijectionTypeDTwisted(KRT(pathlist=[[-1,2]]))
        sage: bijection.cur_path.insert(0, [])
        sage: bijection.cur_dims.insert(0, [0, 1])
        sage: bijection.cur_path[0].insert(0, [2])
        sage: next_state_KRT_RC(2)
    """
    n = bij.n
    tableau_height = len(bij.cur_path[0]) - 1

    if val == 'E':
        raise NotImplementedError("TODO!!!!")
        return
    if val > 0:
        # If it is a regular value, we follow the A_n type rules
        return

    pos_val = -val

    if pos_val == 0:
        if bij.ret_rig_con[pos_val - 1]:
            max_width = bij.ret_rig_con[n-1][0]
        else:
            max_width = 1
        max_width = bij.ret_rig_con[n-1].insert_cell(max_width)
        width_n = max_width + 1

        # Follow regular A_n rules
        for a in reversed(range(tableau_height, n-1)):
            max_width = bij.ret_rig_con[a].insert_cell(max_width)
            bij._update_vacancy_nums(a + 1)
            bij._update_partition_values(a + 1)
        bij._update_vacancy_nums(tableau_height)
        if tableau_height >= n - 1:
            bij._correct_vacancy_nums()
        bij._update_partition_values(tableau_height)
        if tableau_height > 0:
            bij._update_vacancy_nums(tableau_height-1)
            bij._update_partition_values(tableau_height-1)

        # Make the new string at n quasi-singular
        p = bij.ret_rig_con[n-1]
        num_rows = len(p)
        for i in range(num_rows):
            if p._list[i] == width_n:
                j = i+1
                while j < num_rows and p._list[j] == width_n \
                  and p.vacancy_numbers[j] == p.rigging[j]:
                    j += 1
                p.rigging[j-1] -= 1
                break
        return

    case_S = [None] * n
    pos_val = -val

    # Always add a cell to the first singular value in the first
    #   tableau we are updating.
    if bij.ret_rig_con[pos_val - 1]:
        max_width = bij.ret_rig_con[pos_val - 1][0]
    else:
        max_width = 1

    # Add cells similar to type A_n but we move to the right until we
    #   reach the value of n-1
    for a in range(pos_val - 1, n-1):
        max_width = bij.ret_rig_con[a].insert_cell(max_width)
        case_S[a] = max_width

    # Special case for n
    # If we find a quasi-singular string first, then we are in case (Q, S)
    #   otherwise we will find a singular string and insert 2 cells
    partition = bij.ret_rig_con[n-1]
    num_rows = len(partition)
    case_QS = False
    for i in range(num_rows + 1):
        if i == num_rows:
            max_width = 0
            if case_QS:
                partition._list.append(1)
                partition.vacancy_numbers.append(None)
                # Go through our partition until we find a length of greater than 1
                j = len(partition._list) - 1
                while j >= 0 and partition._list[j] == 1:
                    j -= 1
                partition.rigging.insert(j + 1, None)
                width_n = 1
            else:
                # Go through our partition until we find a length of greater than 2
                j = len(partition._list) - 1
                while j >= 0 and partition._list[j] <= 2:
                    j -= 1
                partition._list.insert(j+1, 2)
                partition.vacancy_numbers.insert(j+1, None)
                partition.rigging.insert(j+1, None)
            break
        elif partition._list[i] <= max_width:
            if partition.vacancy_numbers[i] == partition.rigging[i]:
                max_width = partition._list[i]
                if case_QS:
                    partition._list[i] += 1
                    width_n = partition._list[i]
                    partition.rigging[i] = None
                else:
                    j = i - 1
                    while j >= 0 and partition._list[j] <= max_width + 2:
                        partition.rigging[j+1] = partition.rigging[j] # Shuffle it along
                        j -= 1
                    partition._list.pop(i)
                    partition._list.insert(j+1, max_width + 2)
                    partition.rigging[j+1] = None
                break
            elif partition.vacancy_numbers[i] - 1 == partition.rigging[i] and not case_QS:
                case_QS = True
                partition._list[i] += 1
                partition.rigging[i] = None
                # No need to set max_width here since we will find a singular string

    # Now go back following the regular C_n (ish) rules
    for a in reversed(range(tableau_height, n-1)):
        if case_S[a] == max_width:
            insert_cell_case_S(bij.ret_rig_con[a])
        else:
            max_width = bij.ret_rig_con[a].insert_cell(max_width)
        bij._update_vacancy_nums(a + 1)
        bij._update_partition_values(a + 1)

    # Update the final rigged partitions
    bij._update_vacancy_nums(tableau_height)
    if tableau_height >= n - 1:
        bij._correct_vacancy_nums()
    bij._update_partition_values(tableau_height)

    if pos_val <= tableau_height:
        for a in range(pos_val-1, tableau_height):
            bij._update_vacancy_nums(a)
            bij._update_partition_values(a)
        if pos_val > 1:
            bij._update_vacancy_nums(pos_val - 2)
            bij._update_partition_values(pos_val - 2)
    elif tableau_height > 0:
        bij._update_vacancy_nums(tableau_height - 1)
        bij._update_partition_values(tableau_height - 1)

    if case_QS:
        # Make the new string quasi-singular
        num_rows = len(partition)
        for i in range(num_rows):
            if partition._list[i] == width_n:
                j = i+1
                while j < num_rows and partition._list[j] == width_n \
                  and partition.vacancy_numbers[j] == partition.rigging[j]:
                    j += 1
                partition.rigging[j-1] -= 1
                break

def is_simple_positive(x):
    r"""
    Helper function to see if ``x`` is a simple positive element.

    We say `b \in B` is *simple positive* if `\varphi_i(b) = 2 \delta_{ia}`
    for all `a \in I`. This means it corresponds to a simple root in
    the adjoint representation.

    EXAMPLES::

        sage: B = crystals.Letters(['E',8])
        sage: from sage.combinat.rigged_configurations.bij_type_adjoint import is_simple_positive
        sage: is_simple_positive(B((-7,8,8)))
        True
        sage: is_simple_positive(B((-8,8)))
        False 
    """
    total = 0
    for i in x.parent().index_set():
        phi = x.phi(i)
        if phi != 2 and phi != 0:
            return False
        total += phi
    return total == 2

def find_singular_string(partition, last_size, last_row=None):
    r"""
    Return the index of the singular string or ``None`` if not found.

    Helper method to find a singular string at least as long as
    ``last_size`` that has not been previously selected.

    INPUT:

    - ``partition`` -- the partition
    - ``last_size`` -- the last size found
    - ``last_row`` -- the last row found
    """
    if last_row is None:
        last_row = len(partition)
    for i in reversed(range(last_row)):
        if (partition[i] >= last_size
                and partition.vacancy_numbers[i] == partition.rigging[i]):
            return i
    return None

def next_state_RC_KRT(bij, letters):
    r"""
    Build the next state for the left factor `B^{r,1}`, where `r` is
    the node ajoint to node `0`.

    TESTS::

        sage: RC = RiggedConfigurations(['E', 8, 1], [[8, 1]])
        sage: from sage.combinat.rigged_configurations.bij_type_adjoint import RCToKRTBijectionTypeDTwisted
        sage: bijection = RCToKRTBijectionTypeDTwisted(RC(partition_list=[[2],[2,2],[2,2]]))
        sage: next_state_RC_KRT(0)
        -1
    """
    n = bij.n
    I = letters.index_set()
    # The None is a placeholder so we do not have to special case for when
    #   ell/case_S is []
    ell = {a: [None] for a in I}
    # A None signifies that we do not know if case_S holds
    case_S = {a: [None] for a in I}
    case_QS = None
    count = 0

    # Calculate the bijection

    #from sage.typeset.ascii_art import ascii_art
    #print ascii_art(bij.cur_partitions)

    # Start with the positive roots
    last_size = 0
    terminated = False
    b = letters.module_generators[0]
    while not is_simple_positive(b):
        data = [(find_singular_string(bij.cur_partitions[ii], last_size, ell[a][-1]), a, ii)
                for ii,a in enumerate(I) if b.phi(a) > 0]
        data = [(bij.cur_partitions[ii][index], index, a)
                for index,a,ii in data if index is not None]
        if not data:
            terminated = True
            break

        min_length = min(l for l,index,a in data)
        for l,index,a in data:
            if l == min_length:
                last_size = l
                ell[a].append(index)
                case_S[a].append(None)
                b = b.f(a)
                break

    # Do the path through the weight 0 element
    if not terminated:
        a = None
        partition = None
        for ii,j in enumerate(I):
            if b.phi(j) == 2:
                a = j
                partition = bij.cur_partitions[ii]
                break
        assert a is not None
        found = False

        # Modified version of find_singular_string()
        case_Q = False
        ellaset = set(ell[a][1:]) # The first element is a placeholder None
        for i in reversed(range(len(partition))):
            if partition[i] < last_size:
                continue
            # Quasi-singular condition
            if partition.vacancy_numbers[i] - 1 == partition.rigging[i] and not case_Q:
                case_Q = True
                # Check if there are no selected singular rows
                block_size = partition[i]
                for j in reversed(range(i)):
                    if partition[j] != block_size:
                        break
                    elif (partition.vacancy_numbers[j] == partition.rigging[j]
                          and j not in ellaset):
                        case_Q = False
                        break
                if case_Q:
                    b = b.f(a)
                    last_size = partition[i] + 1
                    ellaset.add(i)
                    ell[a] = [None] + sorted(ellaset, reverse=True)
                    case_S[a].append(False)
                    continue
            # Singular condition
            if partition.vacancy_numbers[i] == partition.rigging[i] and i not in ellaset:
                ell[a].append(i)
                found = True
                if partition[i] == 1:
                    b = 'E'
                    case_S[a].append(False)
                    terminated = True
                else:
                    case_S[a].append(None)
                    last_size = partition[i]
                    if case_Q:
                        b = b.f(a)
                        case_S[a][-1] = False
                        case_QS = (a, partition[i])
                    else:
                        # Find the first entry in case_S[a] for this block and make it True
                        # The first element is a placeholder None
                        for ii,index in enumerate(ell[a][1:]):
                            if partition[index] == last_size:
                                case_S[a][ii+1] = True
                                break
                        b = b.f(a).f(a)
                break

        if not found:
            terminated = True

    # Do the negative roots
    while not terminated:
        def last_row(ii):
            a = I[ii]
            partition = bij.cur_partitions[ii]
            last_index = ell[a][-1]
            row_len = partition[last_index]
            for case,index in reversed(zip(case_S[a], ell[a])):
                if case:
                    return last_index
                if case is None:
                    if index is None: # We have reached the end of the list
                        return None
                    # We moved to a different block and haven't found a None
                    if partition[index] < row_len:
                        return last_index
                    last_index = index + 1
                # At this point, ``case is False``, so we only need
                #   to check that we are still in the same sized block
                if partition[index] < row_len:
                    return last_index

        data = [(find_singular_string(bij.cur_partitions[ii], last_size, last_row(ii)), a, ii)
                for ii,a in enumerate(I) if b.phi(a) > 0]
        data = [(bij.cur_partitions[ii][index], index, a)
                for index,a,ii in data if index is not None]
        if not data:
            terminated = True
            break

        min_length = min(l for l,index,a in data)
        for l,index,a in data:
            if l == min_length:
                last_size = l
                try:
                    case_S[a][ell[a].index(index)] = True
                except ValueError: # not a case_S selected string
                    ell[a].append(index)
                    case_S[a].append(False)
                b = b.f(a)
                break

    #for a in I:
    #    print a, ":", ell[a]
    #    print a, ":", case_S[a]

    # Determine the new rigged configuration by removing boxes from the
    #   selected string and then making the new string singular
    for ii,a in enumerate(I):
        ell[a].pop(0)
        case_S[a].pop(0)
        partition = bij.cur_partitions[ii]
        for i,index in enumerate(ell[a]):
            if case_S[a][i]:
                partition.remove_cell(index, 2)
            else:
                partition.remove_cell(index, 1)

    for ii, partition in enumerate(bij.cur_partitions):
        bij._update_vacancy_numbers(ii)
        for i in range(len(partition)):
            if partition.rigging[i] is None:
                partition.rigging[i] = partition.vacancy_numbers[i]

    if case_QS is not None:
        partition = bij.cur_partitions[I.index(case_QS[0])]
        l = case_QS[1] - 1
        for i in reversed(range(len(partition))):
            if partition[i] == l and partition.rigging[i] == partition.vacancy_numbers[i]:
                partition.rigging[i] -= 1
                break

    return(b)

