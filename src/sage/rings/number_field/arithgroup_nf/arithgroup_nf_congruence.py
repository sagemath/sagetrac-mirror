r"""
Congruence subgroups of `{\rm SL}_2(\mathcal{O})`
and `{\rm GL}_2(\mathcal{O})`) where `\mathcal{O}` is an (maximal) order in a number field `K`.

AUTHORS:

- Fredrik Stromberg (2013): initial version based on arithgroup_generic.py

"""

################################################################################
#
#       Copyright (C) 2009, The Sage Group -- http://www.sagemath.org/
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
#
################################################################################

from arithgroup_nf_generic import ArithmeticSubgroup_NF
from sage.misc.latex import latex

class CongruenceSubgroup_NF(ArithmeticSubgroup_NF):
    r"""
    Base class for congruence subgroups of SL(2,O)

    """

    def __init__(self,order,special=True,name='',ltx=''):
        super(CongruenceSubgroup_NF,self).__init__(order,special,name,ltx)

    def level(self):
        r"""
        Return the level (should be an ideal) of self.
        """
        if not hasattr(self,'_level'):
            return None
        return self._level
        
    def _find_cusps(self):
        r"""
        Find a set of cusp representatives.
        Note: currently only implemented for the full Hilbert modular group SL(2,O)
        """
        if self.level() <>self.base_ring().ideal(1) or not self.is_special():
            raise NotImplementedError
        K = self.base_field()
        if self.base_ring()<>self.base_field().ring_of_integers():
            # Then we have representatives given by ideal classes
            lreps = map(lambda x:x.ideal(),self._class_group.list())
            ncusps=len(lreps)
            for a in lreps:
                if self._verbose>0:
                    print "Set cusp info for a={0}".format(a)
                if a.is_trivial():
                    ca = NFCusp(K(1),K(0),lreps=lreps)
                else:
                    ag = a.gens()
                    ca = NFCusp(K,ag[0],ag[1],lreps=lreps)
                cusps.append(ca)
        else:
            raise NotImplementedError
        return cusps

class CongruenceSubgroup_Gamma_NF(CongruenceSubgroup_NF):
    r"""
    Base class for principal congruence subgroups of SL(2,O)

    """

    
    def __init__(self,order,special,level,name='',ltx=''):
        r"""

        INPUT:

        - 'ring' -- ring
        """
        assert order.ideal(level)==level
        self._level = level
        super(CongruenceSubgroup_NF,self).__init__(order,special,name,ltx)


    def __contains__(self,A):
        r"""
        Check if A is in self.

        EXAMPLES::


        
        """
        try:
            a,b,c,d = A
        except:
            return False
        if a - 1 not in self.level() or  d - 1 not in self.level():
            return False
        if b not in self.level() or  c not in self.level():
            return False
        return True

class CongruenceSubgroup_Gamma0_NF(CongruenceSubgroup_NF):
    r"""
    Base class for principal congruence subgroups of SL(2,O)

    """

    
    def __init__(self,order,level,special=True,name='',ltx=''):
        self._level = level
        super(CongruenceSubgroup_NF,self).__init__(order,special,name,ltx)

    def __contains__(self,A):
        r"""
        Check if A is in self.

        EXAMPLES::


        
        """
        try:
            a,b,c,d = A
        except:
            return False
        if a not in self.base_ring() or b not in self.base_ring() or d not in self.base_ring():
            return False
        if c not in self.level():
            return False
        return True    



def HilbertModularGroup(O,special=True,a=None):
    r"""
    Returns the Hilbert modular group `SL_{2}(O\oplus \frak{a})`
    where `O` is a maximal order in a number field and `\frak{a}`
    an ideal in `O`, consisting of matrices of the form
    `[a b // c d ]` with `a,d` in `O`, `c \in \frak{a}` and `d \in \frak{a}^{-1}`.

    INPUT:
    
    - 'O' -- Order in number field
    - 'a' -- Ideal in O. Default = None (meaning it is O itself)

    EXAMPLES::

    sage: from sage.rings.number_field.arithgroup_nf.all import *
    sage: K=QuadraticField(41); O=K.ring_of_integers()
    sage: G=HilbertModularGroup(O)
    sage: G
    Hilbert modular group `SL_{2}(O)`
    

    """
    if special:
        sn='SL'
    else:
        sn='GL'
    if a==None:
        name = 'Hilbert modular group `{sn}_{{2}}(O)`'.format(sn=sn)
        ltx = latex(name)
        return CongruenceSubgroup_Gamma0_NF(O,O,special,name,ltx)
    name = 'Hilbert modular group `{sn}_{{2}}(O+a)`'.format(sn=sn)
    ltx = latex(name)

    raise NotImplementedError
  
