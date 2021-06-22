r"""
Monoterm tensor symmetry groups
"""

# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import operator

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.categories.action import Action
from sage.groups.group import FiniteGroup
from sage.rings.integer import Integer


class TensorSymmetryGroup(metaclass=ClasscallMetaclass):
    r"""
    Constructor for tensor symmetry groups

    INPUT:

    - ``nb_indices`` -- number of indices of the tensor
    - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries among
      the indices: each symmetry is described by a tuple containing the
      positions of the involved indices, with the convention ``position=0``
      for the first slot; for instance:

        * ``sym = (0, 1)`` for a symmetry between the 1st and 2nd indices
        * ``sym = [(0,2), (1,3,4)]`` for a symmetry between the 1st and 3rd
          indices and a symmetry between the 2nd, 4th and 5th indices.

    - ``antisym`` -- (default: ``None``) antisymmetry or list of antisymmetries
      among the indices, with the same convention as for ``sym``

    EXAMPLES::

        sage: from sage.groups.tensor.monoterm_symmetry import TensorSymmetryGroup
        sage: G = TensorSymmetryGroup(2, sym=(0, 1)); G
        Symmetry group of 2-index tensors, with symmetry on the index positions (0, 1)
        sage: V = FiniteRankFreeModule(QQ, 3)
        sage: G.get_action(V.tensor_module(0, 2))
        Left action
         by Symmetry group of 2-index tensors, with symmetry on the index positions (0, 1)
         on Free module of type-(0,2) tensors on the 3-dimensional vector space over the Rational Field

    Antisymmetry on 2 indices::

        sage: G = TensorSymmetryGroup(2, antisym=(0, 1)); G

    Symmetry group of tensors with 6 indices, symmetric among 3 indices (at position
    `(0, 1, 5)`) and antisymmetric among 2 indices (at position `(2, 4)`)::

        sage: G = TensorSymmetryGroup(6, sym=(0, 1, 5), antisym=(2, 4)); G

    Components with 4 indices, antisymmetric with respect to the first pair of
    indices as well as with the second pair of indices::

        sage: G = TensorSymmetryGroup(4, antisym=[(0, 1), (2, 3)]); G

    """
    @staticmethod
    def __classcall__(cls, nb_indices, sym=None, antisym=None):

        if sym:
            if isinstance(sym[0], (int, Integer)):
                # a single symmetry is provided as a tuple or a range object;
                # it is converted to a 1-item list:
                sym = [tuple(sym)]
            sym_list = []
            for isym in sym:
                if len(isym) < 2:
                    raise IndexError("at least two index positions must be " +
                                     "provided to define a symmetry")
                for i in isym:
                    if i < 0 or i > nb_indices:
                        raise IndexError("invalid index position: " + str(i) +
                                         " not in [0," + str(nb_indices) + "]")
                sym_list.append(tuple(sorted(isym)))
            sym = tuple(sorted(sym_list))
        else:
            sym = ()

        if antisym:
            if isinstance(antisym[0], (int, Integer)):
                # a single antisymmetry is provided as a tuple or a range
                # object; it is converted to a 1-item list:
                antisym = [tuple(antisym)]
            antisym_list = []
            for isym in antisym:
                if len(isym) < 2:
                    raise IndexError("at least two index positions must be " +
                                     "provided to define an antisymmetry")
                for i in isym:
                    if i < 0 or i > nb_indices:
                        raise IndexError("invalid index position: " + str(i) +
                                         " not in [0," + str(self._nid-1) + "]")
                antisym_list.append(tuple(sorted(isym)))
            antisym = tuple(sorted(antisym_list))
        else:
            antisym = ()

        # Final consistency check:
        index_list = []
        for isym in sym:
            index_list += list(isym)
        for isym in antisym:
            index_list += list(isym)
        if len(index_list) != len(set(index_list)):
            # There is a repeated index position:
            raise IndexError("incompatible lists of symmetries: the same " +
                             "index position appears more then once")

        return TensorSymmetryGroup_monoterm_permgp(nb_indices, sym, antisym)


class TensorSymmetryGroup_monoterm_permgp(FiniteGroup, UniqueRepresentation):

    def __init__(self, nb_indices, sym, antisym):
        self._nb_indices = nb_indices
        self._sym = sym
        self._antisym = antisym

    def _repr_(self):
        description = f"Symmetry group of {self._nb_indices}-index tensors"
        for isym in self._sym:
            description += ", with symmetry on the index positions " + \
                           str(tuple(isym))
        for isym in self._antisym:
            description += ", with antisymmetry on the index positions " + \
                           str(tuple(isym))
        return description

    def _get_action_(self, X, op, self_on_left):
        if op in (operator.mul, operator.matmul):
            return TensorSymmetryAction_monoterm_permgp(self, X)

    class Element(ElementWrapper):
        pass


class TensorSymmetryAction_monoterm_permgp(Action):
    r"""
    Action of a :class:`TensorSymmetryGroup_monoterm_permgp` on a tensor module.
    """

    def _act_(self, permutation, tensor):
        raise NotImplementedError
