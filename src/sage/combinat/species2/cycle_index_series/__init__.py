# -*- coding: utf-8 -*-
"""
Cycle index series

Reference
---------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux
  Cambridge University Press, 1998

"""
# *****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ******************************************************************************
from sage.categories.cycle_index_series import CycleIndexSeries
from sage.combinat.sf.sf import SymmetricFunctions
from sage.rings.rational_field import QQ
from sage.structure.category_object import CategoryObject
from sage.structure.dynamic_class import dynamic_class
from sage.structure.unique_representation import UniqueRepresentation


class CIS(UniqueRepresentation, CategoryObject):
    """
    Default class of cycle index series
    """

    def __init__(self):
        CategoryObject.__init__(self, category=CycleIndexSeries())

        # ##### CUSTOM ####### used to inherit the parent methods
        base = self.__class__.__base__
        self.__class__ = dynamic_class("%s_with_category" % base.__name__,
                                       (self.__class__, self.category().parent_class, ),
                                       doccls=base)
        # ####################

        self._sym_ = SymmetricFunctions(QQ)