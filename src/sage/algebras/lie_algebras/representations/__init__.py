r"""
Lie Algebra representations

AUTHORS:

- Michael Walter (2018-07-22): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Michael Walter <m.walter@uva.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import absolute_import

from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.cartan_type import CartanType_abstract
from sage.combinat.root_system.root_system import RootSystem
from sage.combinat.partition import Partition
from sage.misc.cachefunc import cached_method
from sage.rings.all import QQ
from sage.structure.element import Element
import sage.combinat.root_system.type_A as type_A


def lie_algebra_representation(highest_weight, cartan_type=None, K=QQ):
    r"""
    Return irreducible Lie algebra representation with given highest weight.

    INPUT:

    - ``highest_weight`` -- the highest weight (either an element of a weight lattice or
      a tuple, list, :class:`~sage.combinat.partition.Partition`, etc.)
    - ``cartan_type`` -- the Cartan type of the Lie algebra (optional, only necessary
      when the highest weight is not an element of a weight lattice)
    - ``K`` -- the base field (default: $\QQ$)

    OUTPUT:

    The Lie algebra representation.

    EXAMPLES::

        sage: Lambda = RootSystem(['A', 2]).weight_lattice().fundamental_weights()
        sage: lie_algebra_representation(Lambda[1] + Lambda[2])
        Lie algebra representation of ['A', 2] with highest weight Lambda[1] + Lambda[2]

        sage: L = RootSystem(['A', 2]).ambient_lattice()
        sage: lie_algebra_representation(L([2, 1, 0]))
        Lie algebra representation of ['A', 2] with highest weight (2, 1, 0)

        sage: lie_algebra_representation(L([2, 1, 0]), cartan_type="A2")
        Lie algebra representation of ['A', 2] with highest weight (2, 1, 0)
        sage: lie_algebra_representation((2, 1, 0), cartan_type="A2")
        Lie algebra representation of ['A', 2] with highest weight (2, 1, 0)
        sage: lie_algebra_representation(Partition([2, 1]), cartan_type="A2")
        Lie algebra representation of ['A', 2] with highest weight (2, 1, 0)
        sage: lie_algebra_representation((2, 1), cartan_type="A2")
        Traceback (most recent call last):
        ...
        ValueError: highest weight (2, 1) has 2 components, expected 3 components
        sage: lie_algebra_representation(Partition([4, 3, 2, 1]), cartan_type="A2")
        Traceback (most recent call last):
        ...
        ValueError: partition [4, 3, 2, 1] has 4 parts, expected no more than 3 parts
    """
    # determine Cartan type
    if cartan_type is None:
        if highest_weight is None:
            raise ValueError("need to provide Cartan type or highest weight")
        cartan_type = highest_weight.parent().root_system.cartan_type()
    elif not isinstance(cartan_type, CartanType_abstract):
        cartan_type = CartanType(cartan_type)
    L = RootSystem(cartan_type).ambient_lattice()

    # type A
    if isinstance(cartan_type, type_A.CartanType):
        from sage.algebras.lie_algebras.representations.type_A_gelfand_tsetlin import (
            IrreducibleRepresentation
        )

        # determine highest weight
        if isinstance(highest_weight, Partition):
            if len(highest_weight) > L.dimension():
                raise ValueError(
                    "partition %s has %d parts, expected no more than %d parts"
                    % (highest_weight, len(highest_weight), L.dimension())
                )
            highest_weight = list(highest_weight)
            highest_weight += [0] * (L.dimension() - len(highest_weight))
            highest_weight = L(highest_weight)
        elif isinstance(highest_weight, (tuple, list)):
            if len(highest_weight) != L.dimension():
                raise ValueError(
                    "highest weight %s has %d components, expected %d components"
                    % (highest_weight, len(highest_weight), L.dimension())
                )
            highest_weight = L(highest_weight)

        return IrreducibleRepresentation(highest_weight, K)

    raise NotImplementedError(
        "Lie algebra representations of %s not implemented yet" % cartan_type
    )
