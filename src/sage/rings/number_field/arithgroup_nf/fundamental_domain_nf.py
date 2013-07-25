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
    Fundamental domain for Hilbert modular group.
    """

    def __init__(self,H,**kwds):
        r"""


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
            group = 'sl'
        else:
            group = 'gl'
        embs = kwds.get('embs',None)
        self._hfd = _HilbertFundDomData(self._K, embs=embs, group=group)
        
    


    def reduce(self,z,**kwds):
        max_round = kwds.get('max_rounds',0)
        bound_on_c = kwds.get('bound_on_c',None)
        return self._hfd.reduce(z, max_rounds, bound_on_c)

    def step(self,s,**kwds):
        return self._hfd.step(z,**kwds)
        
        
    def  plot(self,**kwds):
        raise NotImplementedError

