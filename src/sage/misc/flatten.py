"Flatten nested lists"

import sys
def flatten(in_list, ltypes=(list, tuple), max_level=sys.maxsize):
    """
    Flattens a nested list.

    INPUT:

    - in_list -- a list or tuple
    - ltypes -- optional list of particular types to flatten
    - max_level -- the maximum level to flatten

    OUTPUT:

    a flat list of the entries of in_list

    EXAMPLES::

       sage: flatten([[1,1],[1],2])
       [1, 1, 1, 2]
       sage: flatten([[1,2,3], (4,5), [[[1],[2]]]])
       [1, 2, 3, 4, 5, 1, 2]
       sage: flatten([[1,2,3], (4,5), [[[1],[2]]]],max_level=1)
       [1, 2, 3, 4, 5, [[1], [2]]]
       sage: flatten([[[3],[]]],max_level=0)
       [[[3], []]]
       sage: flatten([[[3],[]]],max_level=1)
       [[3], []]
       sage: flatten([[[3],[]]],max_level=2)
       [3]


    In the following example, the vector isn't flattened because
    it is not given in the ltypes input. ::

       sage: flatten((['Hi',2,vector(QQ,[1,2,3])],(4,5,6)))
       ['Hi', 2, (1, 2, 3), 4, 5, 6]

    We give the vector type and then even the vector gets flattened::

       sage: flatten((['Hi',2,vector(QQ,[1,2,3])], (4,5,6)), ltypes=(list, tuple,sage.modules.vector_rational_dense.Vector_rational_dense))
       ['Hi', 2, 1, 2, 3, 4, 5, 6]

    We flatten a finite field. ::

       sage: flatten(GF(5))
       [0, 1, 2, 3, 4]
       sage: flatten([GF(5)])
       [Finite Field of size 5]
       sage: flatten([GF(5)], ltypes = (list, tuple, sage.rings.finite_rings.finite_field_prime_modn.FiniteField_prime_modn))
       [0, 1, 2, 3, 4]

    Degenerate cases::

        sage: flatten([[],[]])
        []
        sage: flatten([[[]]])
        []
    """
    index = 0
    current_level = 0
    new_list = [x for x in in_list]
    level_list = [0]*len(in_list)

    while index < len(new_list):
        len_v=None
        while isinstance(new_list[index], ltypes) and current_level<max_level:
            v = list(new_list[index])
            len_v = len(v)
            new_list[index : index + 1] = v
            old_level = level_list[index]
            level_list[index : index + 1] = [0]*len_v
            if len_v != 0:
                current_level += 1
                level_list[index + len_v - 1] = old_level + 1
            else:
                current_level -= old_level
                index -=1
                break

        # If len_v==0, then index points to a previous element, so we
        # don't need to do anything.
        if len_v!=0:
            current_level -= level_list[index]
        index += 1
    return new_list
