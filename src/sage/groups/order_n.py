r"""
<Short one-line summary that ends with no period>

<Paragraph description>

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Alexis Newton (2021-06-28): initial version


"""

# ****************************************************************************
#       Copyright (C) 2021 Alexis Newton alexis.newton@emory.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

def order_n(*k):
    import sage.libs.gap.all
    import sage.groups.all
    r"""
    Return the list of groups of order less than or equal to `n` in the GAP library.

    INPUT:

    - ``n`` -- integer; maximum order of groups

    OUTPUT: groups as a list of group descriptions

    EXAMPLES:

    This example illustrates ... ::

        sage: A = ModuliSpace()
        sage: A.point(2,3)
        xxx

    We now ... ::

        sage: B = A.point(5,6)
        sage: xxx

    It is an error to ... ::

        sage: C = A.point('x',7)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 'r' to an integer

    .. NOTE::

        This function uses the algorithm of [BCDT2001]_ to determine
        whether an elliptic curve `E` over `Q` is modular.

    ...

    .. SEEALSO::

        :func:`line`

    TESTS::

        sage: A.point(42, 0)  # Check for corner case y=0
        xxx
    """
    


    possibleparameter = ["IsAbelian","IsCyclic","IsSolvable","IsNilpotent"]
    parameter = [];
    below_order_n = gap('EmptyPlist(1)') 
    for j in [2..len(k)]:
        if k[j-1] in possibleparameter:
            parameter.append(k[j-1])
    if len(parameter) == 0:
        for i in [1..k[0]]:
            g=gap(i).AllSmallGroups()
            for x in [1..len(g)]:  
                r=gap(g[x]).StructureDescription()
                gap.Add(below_order_n,r)
    if len(parameter) == 1:
        for i in [1..k[0]]:
            g=gap(i).AllSmallGroups(parameter[0])
            for x in [1..len(g)]:  
                r=gap(g[x]).StructureDescription()
                gap.Add(below_order_n,r)
    if len(parameter) == 2:
        for i in [1..k[0]]:
            g=gap(i).AllSmallGroups(parameter[0],parameter[1])
            for x in [1..len(g)]:  
                r=gap(g[x]).StructureDescription()
                gap.Add(below_order_n,r)
    if len(parameter) == 3:
        for i in [1..k[0]]:
            g=gap(i).AllSmallGroups(parameter[0],parameter[1],parameter[2])
            for x in [1..len(g)]:  
                r=gap(g[x]).StructureDescription()
                gap.Add(below_order_n,r)
    if len(parameter) == 4:
        for i in [1..k[0]]:
            g=gap(i).AllSmallGroups(parameter[0],parameter[1],parameter[2],parameter[3])
            for x in [1..len(g)]:  
                r=gap(g[x]).StructureDescription()
                gap.Add(below_order_n,r)
    return below_order_n

        






