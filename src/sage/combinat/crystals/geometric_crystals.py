r"""
Geometric Crystals

AUTHORS:

- Travis Scrimshaw: Initial version
"""

#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

from copy import copy, deepcopy
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.categories.classical_crystals import Crystals
from sage.matrix.constructor import matrix
from sage.matrix.matrix import is_Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.matrix_integer_sparse import Matrix_integer_sparse
from sage.rings.all import ZZ, QQ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.combinat.combinat import CombinatorialObject
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.root_system import RootSystem
from sage.combinat.root_system.weyl_group import WeylGroup

class GeometricCrystalElement(Element):
    r"""
    A geometric crystals element.

    For more information, see :class:`GeometricCrystal`.

    EXAMPLES::
    """
    def __init__(self, parent, M):
        r"""
        EXAMPLES::

            sage: C = GeometricCrystal(2, [1, 1])
            sage: mg = C.module_generators[0]
            sage: TestSuite(mg).run()
        """
        self._M = M
        Element.__init__(self, parent)

    def _repr_(self):
        r"""
        EXAMPLES::
        """
        return repr(self._M)

    def _latex_(self):
        r"""
        Generate LaTeX code for ``self``.  Requires TikZ.

        EXAMPLES::

            sage: x = GeometricCrystal(2, [1, 1]).module_generators[0]
            sage: latex(x)
        """
        return self._M._latex_()

    def __eq__(self, rhs):
        """
        Check equality.
        """
        if isinstance(rhs, GeometricCrystalElement):
            return self._M == rhs._M
        return self._M == rhs

    def __ne__(self, rhs):
        """
        Check not equals.
        """
        return not self.__eq__(rhs)

    def e(self, i, c=1):
        r"""
        Return the action of `e_i^c` on ``self``.

        EXAMPLES::
        """
        phi = self._phi(i)
        if phi == 0:
            return self
        xl = self.parent()._x(i, (c - 1) / phi)
        xr = self.parent()._x(i, (c**-1 - 1) / self._epsilon(i))
        return self.__class__(self.parent(), xl * self._M * xr)

    def f(self, i, c=1):
        r"""
        Return the action of `f_i^c` on ``self``.

        EXAMPLES::
        """
        phi = self._phi(i)
        if phi == 0:
            return self
        xl = self.parent()._x( i, (c - 1) / (c * phi) )
        xr = self.parent()._x( i, (c**-1 - 1) / (c**-1 * self._epsilon(i)) )
        return self.__class__(self.parent(), xl**-1 * self._M * xr**-1)

    def _epsilon(self, i):
        r"""
        Return `\varepsilon_i` of ``self``.
        """
        return self._M[i,i-1] / self._M[i,i] # -1 for indexing

    def _phi(self, i):
        r"""
        Return `\varphi_i` of ``self``.
        """
        return self._M[i,i-1] / self._M[i-1,i-1] # -1 for indexing

    def weight(self):
        r"""
        Return the weight of ``self`` as an element of the root lattice
        `\bigoplus_{i \in I} \ZZ \alpha_i`.

        EXAMPLES::

            sage: C = GeometricCrystal(2, [1, 1])
            sage: x = C.module_generators[0]
            sage: x.weight()
        """
        return self._M.diagonal()

    def decoration(self):
        """
        Return the decoration function `f_{G,\chi}` of ``self``.
        """
        n = self._M.nrows()
        minor = lambda J,Jprime: self._M.matrix_from_rows_and_columns(J, Jprime).determinat()
        frac = lambda i: ( minor(range(n-i,n,2), range(0,i,2))
            + minor(range(n+1-i,n,2), range(0,i+1,2)) ) \
            / minor(range(n+1-i,n,2), range(0,i,2))
        return sum(frac(i) for i in range(n-1))

    def s(self, i):
        """
        Return `s_i` of ``self``.
        """
        wt = self.weight()
        return self.e(i, wt[i-1]/wt[i])

class GeometricCrystal(Parent, UniqueRepresentation):
    r"""
    Geometric crystal.

    INPUT:

    - ``n`` -- The group `GL_n`
    - ``weight`` -- The weight

    EXAMPLES::
    """
    @staticmethod
    def __classcall_private__(cls, n, weight):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: C = GeometricCrystal(3, [3,1])
            sage: C2 = GeometricCrystal(int(3), (3,1,0))
            sage: C is C2
            True
        """
        weight = list(weight)
        if len(weight) > n:
            raise ValueError("The weight can have at most %s entries"%n)
        weight += [0]*(n-len(weight))
        return super(GeometricCrystal, cls).__classcall__(cls, n, tuple(weight))

    def __init__(self, n, weight):
        r"""
        EXAMPLES::

            sage: C = GeometricCrystal(3, [3,1,1])
            sage: TestSuite(C).run()
        """
        self._cartan_type = CartanType(['A',n-1])
        self._weight = weight
        self._MS = MatrixSpace(QQ, n, sparse=True)
        Parent.__init__(self, category=Crystals())
        self.module_generators = (self.element_class(self, self._MS.identity_matrix()),)

    def _element_constructor_(self, data):
        r"""
        Construct an element of ``self`` from ``data``.

        INPUT:

        - ``data`` -- A weight

        EXAMPLES::
        """
        return self.element_class(self, data)

    Element = GeometricCrystalElement

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: TropicalizedGeometricCrystal(2, [1, 2])
            Geometric crystal of type ['A', 2] with weight (0, 1, 2)
        """
        return "Geometric crystal of type %s with weight %s"%(self._cartan_type, self._weight)

    def _x(self, i, c):
        """
        Construct the matrix `x_i(c)`.
        """
        m = copy(self._MS.identity_matrix())
        m[i-1,i] = c # -1 for indexing
        return m

##################################################################
## Tropicalized version

class TropicalizedGeometricCrystalElement(CombinatorialObject, Element):
    r"""
    A geometric crystal element.
    """
    def __init__(self, parent, l):
        r"""
        EXAMPLES::

            sage: Y = CrystalOfMatrices(5, 2)
            sage: mg = Y.module_generators[0]
            sage: TestSuite(mg).run()
        """
        CombinatorialObject.__init__(self, l)
        Element.__init__(self, parent)

    def e(self, i):
        r"""
        Return the action of `e_i` on ``self``.

        EXAMPLES::
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        return None

    def f(self,i):
        r"""
        Return the action of `f_i` on ``self``.

        EXAMPLES::
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        return None

    def epsilon(self, i):
        r"""
        Return `\varepsilon_i` of ``self``.
        """
        M = self.parent()._theta
        pos = min(sum(s*self[j] for j,s in enumerate(exp)) for exp in M[i,i-1].exponents())
        neg = min(sum(s*self[j] for j,s in enumerate(exp)) for exp in M[i,i].exponents())
        return pos - neg

    def phi(self, i):
        r"""
        Return `\varphi_i` of ``self``.
        """
        M = self.parent()._theta
        i -= 1 # for indexing
        pos = min(sum(s*self[j] for j,s in enumerate(exp)) for exp in M[i+1,i].exponents())
        neg = min(sum(s*self[j] for j,s in enumerate(exp)) for exp in M[i,i].exponents())
        return self.parent().weight[i+1]*pos - self.parent().weight[i]*neg

class TropicalizedGeometricCrystal(Parent, UniqueRepresentation):
    r"""
    Tropicalized geometric crystal.

    INPUT:

    - ``n`` -- The group `GL_n`
    - ``weight`` -- The weight

    EXAMPLES::
    """
    @staticmethod
    def __classcall_private__(cls, n, weight):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: C = TropicalizedGeometricCrystal(3, [3,1])
            sage: C2 = TropicalizedGeometricCrystal(int(3), (3,1,0))
            sage: C is C2
            True
        """
        weight = list(weight)
        if len(weight) > n:
            raise ValueError("The weight can have at most %s entries"%n)
        weight += [0]*(n-len(weight))
        return super(TropicalizedGeometricCrystal, cls).__classcall__(cls,n,tuple(weight))

    def __init__(self, n, weight):
        r"""
        EXAMPLES::

            sage: C = TropicalizedGeometricCrystal(3, [3,1,1])
            sage: TestSuite(C).run()
        """
        self._cartan_type = CartanType(['A',n-1])
        self._weight = weight

        # Construct the matrix theta
        W = WeylGroup(self._cartan_type)
        long_elt = W.long_element()
        ell = long_elt.length()
        R = LaurentPolynomialRing(ZZ, ell, ['c%s'%(i+1) for i in range(ell)])
        #T = PolynomialRing(R, n, ['t%s'%(i+1) for i in range(n)])
        MS = MatrixSpace(R, n, sparse=True)
        self._theta = MS.identity_matrix() #MS({(i,i):t for i,t in enumerate(T.gens())})
        for i,j in enumerate(long_elt.reduced_word()):
            # j is one based
            M = copy(MS.identity_matrix())
            M[j-1,j-1] = R.gen(i)**(-1)
            M[j,j] = R.gen(i)
            M[j,j-1] = 1
            self._theta *= M
        self._theta.set_immutable()

        Parent.__init__(self, category=Crystals())
        self.module_generators = (self.element_class(self, [0]*n),)

    def _element_constructor_(self, data):
        r"""
        Construct an element of ``self`` from ``data``.

        INPUT:

        - ``data`` -- A weight

        EXAMPLES::
        """
        return self.element_class(self, data)

    Element = TropicalizedGeometricCrystalElement

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: TropicalizedGeometricCrystal(2, [1, 2])
            Tropicalized geometric crystal of type ['A', 2] with weight (0, 1, 2)
        """
        return "Tropicalized geometric crystal of type %s with weight %s"%(self._cartan_type, self._weight)

