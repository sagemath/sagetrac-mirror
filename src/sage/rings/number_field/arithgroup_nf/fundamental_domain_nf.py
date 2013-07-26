r"""

Fundamental domains for Hilbert modular groups.

AUTHORS:

- Fredrik Stromberg (2013): creating initial classes.
"""

#*****************************************************************************
#       Copyright (C) 2013 Fredrik Stromberg <fredrik314@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from fundamental_domain_nf_reduction import _HilbertFundDomData
from sage.structure.sage_object import SageObject

class HilbertFundamentalDomain(SageObject):
    r"""
    Fundamental domain for Hilbert modular groups.
    """

    def __init__(self,H,**kwds):
        r"""
        Init the fundamental domain.
        
        """
        self._group = H
        self._K = H.number_field()
        if H.is_special():
            raise NotImplementedError
        if H.is_congruence():
            if H.index()<>1:
                raise NotImplementedError
        else:
            raise NotImplementedError
        if H.number_field().class_number()<>1:
            raise NotImplementedError
        if H.is_special():
            groupn = 'sl'
        else:
            groupn = 'gl'
        embs = kwds.get('embs',None)
        self._hfd = _HilbertFundDomData(self._K, embs=embs, group=groupn)
        

    def group(self):
        r"""
        Return the group for which self is a fundamental domain.

        """
        return self._group
    def __repr__(self):
        r"""
        string representation of self.
        """
        s = "Fundamental domain for {0}".format(self.group())

    def reduce(self,z,**kwds):
        r"""
        Reduce self.
        """
        
        max_rounds = kwds.get('max_rounds',0)
        bound_on_c = kwds.get('bound_on_c',None)
        return self._hfd.reduce(z, max_rounds, bound_on_c)

    def step(self,z,**kwds):
        r"""
        Apply one step in the reduction algorithm.
        """
        return self._hfd.step(z,**kwds)
        
        
    def  plot(self,**kwds):
        r"""
        Plot self.
        """
        raise NotImplementedError

