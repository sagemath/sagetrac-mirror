from sage.structure.element import ModuleElement
from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
from manin_map import ManinMap
import operator

from sage.categories.action import Action

M2Z = MatrixSpace_ZZ_2x2()
minusproj = M2Z([1,0,0,-1])

class ModSymAction(Action):
    def __init__(self, MSspace):
        Action.__init__(self, M2Z, MSspace, False, operator.mul)

    def _call_(self, sym, g):
        return sym.__class__(sym._map * g, sym.parent())

class ModularSymbolElement(ModuleElement):
    def __init__(self, map_data, parent):
        ModuleElement.__init__(self, parent)
        self._map = ManinMap(parent._manin_relations, map_data)

    def _add_(self, right):
        return self.__class__(self._map + right._map, self.parent())

    def _lmul_(self, right):
        return self.__class__(self._map * right, self.parent())

    def _sub_(self, right):
        return self.__class__(self._map - right._map, self.parent())

    def plus_part(self):
        r"""
        Returns the plus part of self -- i.e. self + self | [1,0,0,-1].

        Note that we haven't divided by 2.  Is this a problem?

        OUTPUT:

        self + self | [1,0,0,-1]

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: (phi.plus_part()+phi.minus_part()) == phi.scale(2)
        True
        """
        return self * minusproj + self
