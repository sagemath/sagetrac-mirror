Class for ideals of a hecke algebra.

AUTHOR:

- Preston Wake (2009-08-05)
"""
#*****************************************************************************
#       Copyright (C) 2009 Preston Wake <preston.wake@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.modular.modsym as modsym
import sage.rings.all as rings

from sage.modular.hecke.hecke_operator import is_HeckeAlgebraElement
from sage.rings.arith import lcm
from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import FreeModule

import algebra


def HeckeIdeal(hecke, gens):
    r"""
    Creates a hecke ideal from generators, or from a module.
    
    INPUT:

    - ``hecke`` - a Hecke algebra

    - ``gens`` - a list of generators of the ideal, or a module. Gens consists of elements of hecke

    OUTPUT:

    - a Hecke ideal

    EXAMPLES::

        sage: TT = ModularSymbols(11,2,1).hecke_algebra()
        sage: HeckeIdeal(TT, [TT.2 + TT.3, TT.4, TT.7])
        Ideal of Full Hecke algebra acting on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
        sage: HeckeIdeal(TT, TT.2)
        Ideal of Full Hecke algebra acting on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
    """
    try:
        B = gens.basis()
        if not is_HeckeAlgebraElement(gens[0]):
            raise TypeError, "Input must be Hecke Algebra elements."
        return HeckeIdeal_class(hecke,B)
    except:
        if not isinstance(gens, (list, tuple)):
            gens = [gens]
        if len(gens)==0: gens=(hecke(rings.ZZ(0)).matrix())
        if not is_HeckeAlgebraElement(gens[0]):
            raise TypeError, "Input must be Hecke Algebra elements."
        return HeckeIdeal_class(hecke,gens)


class HeckeIdeal_class():
    r"""
    Class for ideals of the hecke algebra. Generators are stored as coordinate vectors, not as matrices.
    """
    def __init__(self, hecke, gens):
        r"""
        INPUT:
        
        - ``hecke`` - a Hecke algebra

        - ``gens`` - a list of generators of the ideal, must be Hecke algebra elements, or coordinate vectors representing Hecke algebra elements

        NOTES:
    
        The end user should initialize Hecke ideals using the HeckeIdeal function.
        """
        self.__hecke = hecke
        try:
            gens = [ T.coordinates() for T in gens]
        except:
            pass
        self.__gens = tuple(gens)
    
    def __repr__(self):
        r"""
        String representation of self
        
        EXAMPLE::

            sage: TT = ModularSymbols(11,2,1).hecke_algebra()
            sage: HeckeIdeal(TT,TT.2).__repr__()
            'Ideal of Full Hecke algebra acting on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field'
        """
        return "Ideal of %s"%self.__hecke

    def __cmp__(self,other):
        r"""
        Uses submodule machinary.

        EXAMPLES::

            sage: TT = ModularSymbols(11,2,1).hecke_algebra()
            sage: I=HeckeIdeal(TT,TT.2)
            sage: J=HeckeIdeal(TT,-TT.2)
            sage: I==J      #indirect doctest
            True
        """
        S = self._as_submodule()
        T = other._as_submodule()
        return S.__cmp__(T)
    
    def __nonzero__(self):
        r"""
        Determines if self is the zero ideal.

        EXAMPLES::

            sage: TT = ModularSymbols(11,2,1).hecke_algebra()
            sage: HeckeIdeal(TT, TT(0)).__nonzero__()
            False
            sage: HeckeIdeal(TT, TT.2).__nonzero__()
            True
        """
        for g in self.__gens:
            if not g.is_zero():
                return True
        return False

    def _as_submodule(self):
        r"""
        Represents the ideal as a submodule of a Free module over ZZ.

        EXAMPLE::

            sage: TT = ModularSymbols(11,2,1).hecke_algebra()
            sage: HeckeIdeal(TT, TT.2)._as_submodule()
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [ 3 -1]
        """
        gens = self.__gens
        d = len(gens[0])
        MM = FreeModule(rings.ZZ,d)
        return MM.submodule(gens)

    def __contains__(self,v):
        r"""
        Checks if self contains v. Uses the submodule interpretation.

        INPUT:
    
        - ``v`` - an element of the hecke algebra, or a coordinate vector

        EXAMPLES::

            sage: TT = ModularSymbols(11,2,1).hecke_algebra()
            sage: HeckeIdeal(TT, [TT.2, TT.3]).__contains__(TT.2 + TT.3)
            True
            sage: TT.2 + TT.3 in HeckeIdeal(TT, [TT.2, TT.3])
            True
        """
        S = self._as_submodule()
        gens = self.__gens
        d = len(gens[0])
        MM = FreeModule(rings.ZZ,d)
        try:
            v=v.coordinates()
            v=MM(v)
            return S.__contains__(v)
        except:
            v=MM(v)
            return S.__contains__(v)
    
    def heckealgebra(self):
        r"""
        Returns the Hecke Algebra used to define self.

        EXAMPLES::
        
            sage: TT = ModularSymbols(11,2,1).hecke_algebra()
            sage: HeckeIdeal(TT, [TT.2, TT.3]).heckealgebra()
            Full Hecke algebra acting on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
        """
        return self.__hecke

    def module(self):
        r"""
        Returns the Hecke Module of the Hecke Algebra used to define self.

        EXAMPLES::
        
            sage: TT = ModularSymbols(11,2,1).hecke_algebra()
            sage: HeckeIdeal(TT, [TT.2, TT.3]).module()
            Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
        """
        return self.__hecke.module()

    def gens(self):
        r"""
        Returns generators of self as matrices

        EXAMPLES::

            sage: TT = ModularSymbols(11,2,1).hecke_algebra()
            sage: HeckeIdeal(TT, [TT.2+TT.7, TT.3]).gens()
            [Hecke operator on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field defined by:
            [11 -3]                                                                                                                        
            [ 0 -4], Hecke operator on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field defined by:
            [ 4 -1]                                                                                                                                
            [ 0 -1]]          
        """
        B=self.heckealgebra().basis()
        gens=[]
        for v in self.__gens:
            s=_matrix_from_coords(v,B)
            gens.append(s)                     
        return gens

    def index_in(self,other):
        r"""
        Returns the index of self in other.
    
        EXAMPLES::

        """
        S = self._as_submodule()
        return S.index_in(other._as_submodule())

    def quo(self):
        r"""
        Finds the quotient self.__hecke / self. Right now this is only implemented if [self.__hecke:self] is finite.

        OUTPUT:

        - an Abelian group

        EXAMPLES::

        """
        TT = self.__hecke
        N = self._as_submodule()
        M = HeckeIdeal(TT,TT.basis())._as_submodule()
        if not rank(N)==rank(M):
            raise NotImplementedError, "Infinite index quotients not implemented"
        G = M.coordinate_submodule(N).matrix().elementary_divisors()
        return AbelianGroup(G)

    def __add__(self,other):
        r"""
        The sum of two ideals.

        EXAMPLES::

            sage: TT = ModularSymbols(11,2,1).hecke_algebra()
            sage: I = HeckeIdeal(TT, [TT.2 + TT.3, TT.5])
            sage: J = HeckeIdeal(TT, TT.2)
            sage: I+J
            Ideal of Full Hecke algebra acting on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: (I+J).gens()
            [Hecke operator on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field defined by:
            [10 -3]
            [ 0 -5], Hecke operator on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field defined by:
            [ 9 -2]
            [ 0 -1]]
        """
        gens = [x+y for x in self.__gens for y in other.__gens]
        return HeckeIdeal_class(self.__hecke,gens)

    def __mul__(self,other):
        r"""
        The product of two ideals.

        EXAMPLES::

            sage: TT = ModularSymbols(11,2,1).hecke_algebra()
            sage: I = HeckeIdeal(TT, [TT.2 + TT.3, TT.5])
            sage: J = HeckeIdeal(TT, TT.2)
            sage: I*J
            Ideal of Full Hecke algebra acting on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: (I*J).gens()
            [Hecke operator on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field defined by:            
            [21 -3]                                                                                                                                    
            [ 0  6], Hecke operator on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field defined by:    
            [18 -4]                                                                                                                                    
            [ 0 -2]]
        """
        hecke = self.__hecke 
        gens = [hecke(x*y, check=False) for x in self.gens() for y in other.gens()]
        coords = [g.coordinates() for g in gens]
        return HeckeIdeal_class(self.__hecke,coords)
 
    def intersection(self,other):
        r"""
        The intersection of two ideals. Uses the submodule machinery.

        EXAMPLES::

            sage: TT = ModularSymbols(11,2,1).hecke_algebra()
            sage: I = HeckeIdeal(TT, [TT.2, TT.5])
            sage: J = HeckeIdeal(TT, TT.2)
            sage: I.intersection(J)
            Ideal of Full Hecke algebra acting on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: I.intersection(J).cmp(J)
            True
        """
        I = self._as_submodule()
        J = other._as_submodule()
        return HeckeIdeal_class(self.__hecke,I.intersection(J).basis())


def _matrix_from_coords(v,B):
    r"""
    Returns the hecke matrix associated to given coordinates.
    
    INPUT:
    
    - ``v`` - a vector
    
    - ``B`` - a basis of the hecke algebra

    OUTPUT:

    - a matrix representing a hecke algebra element

    EXAMPLES::
        
        sage: TT = ModularSymbols(11,2,1).hecke_algebra()
        sage: t = TT.2
        sage: sage.modular.hecke.ideal._matrix_from_coords(t.coordinates(),TT.basis())
        Hecke operator on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field defined by:             
        [ 3 -1]                                                                                                                                    
        [ 0 -2]
        sage: t.matrix()
        [ 3 -1]
        [ 0 -2]
    """
    s=sum(v[i]*B[i] for i in xrange(len(v)))
+    return s
