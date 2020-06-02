#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
AUTHORS: - Bruce Westbury (2018): initial version

#*****************************************************************************
#       Copyright (C) 2018 Bruce Westbury <bruce.westbury@gmail.com>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

This implements the combinatorics of vacillating tableaux.
This module defines two classes:
    - the Parent class is VacillatingTableaux
    - the Element class is VacillatingTableau
    
The Category is PathTableaux.


"""
from six import add_metaclass

from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent

import pathtableaux

###############################################################################

@add_metaclass(InheritComparisonClasscallMetaclass)
class VacillatingTableau(AbstractTableau):


@add_metaclass(InheritComparisonClasscallMetaclass)
class VacillatingTableau(ClonableArray):

    @staticmethod
    def __classcall_private__(self, ot):
        try:
            ot = tuple([ Partition(a) for a in ot ])
        except TypeError:
            raise ValueError("%s is not a sequence of partitions." % str(ot) )
        return VacillatingTableaux()(ot)

    def _hash_(self):
        return hash(tuple(map(tuple, self)))

    
    def check(self):
        pass
    
    def _rule(x):
        raise NotImplementedError("I wish I knew!")
        
###############################################################################
        
class VacillatingTableaux(UniqueRepresentation,Parent):

    @staticmethod
    def __classcall_private__(cls):
        return super(VacillatingTableaux, cls).__classcall__(cls)

    def __init__(self):

        Parent.__init__(self, category=PathTableaux())

    def __contains__(self, ot):

        return isinstance(ot, (list, tuple, VacillatingTableau))

    def _element_constructor_(self, ot, check=True):

        if isinstance(ot, VacillatingTableaux) and ot.parent() == self:
            return ot

        return self.element_class(self, list(ot))

    Element = VacillatingTableau

