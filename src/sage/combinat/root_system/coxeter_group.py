r"""
Coxeter groups

AUTHORS:

- Christian Stump

.. note::

    - For definitions and classification types of finite complex reflection
      groups, see :wikipedia:`Complex_reflection_group`.
    - Uses the GAP3 package *chevie*.

Version: 2011-04-26

EXAMPLES::


"""
#*****************************************************************************
#       Copyright (C) 2011 Christian Stump <christian.stump at lacim.ca>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from copy import copy
from sage.misc.all import prod
from sage.misc.cachefunc import cached_function, cached_method, cached_in_parent_method
from sage.categories.category import Category
from sage.categories.finite_permutation_groups import FinitePermutationGroups
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.combinat.root_system.cartan_type import CartanType
from sage.groups.perm_gps.permgroup import PermutationGroup_generic

from sage.rings.all import ZZ, QQ
from sage.matrix.all import Matrix, identity_matrix
from sage.matrix.matrix import is_Matrix
from sage.interfaces.gap3 import GAP3Record, gap3
from sage.interfaces.gap import gap
from sage.combinat.words.word import Word
from sage.rings.arith import gcd, lcm
from sage.combinat.root_system.complex_reflection_group import FiniteComplexReflectionGroup, IrreducibleFiniteComplexReflectionGroup, assert_chevie_available
from sage.categories.coxeter_groups import CoxeterGroups
from sage.combinat.root_system.cartan_matrix import CartanMatrix

from sage.misc.lazy_import import lazy_import
lazy_import('sage.groups.matrix_gps.coxeter_group', 'CoxeterMatrixGroup')

class FiniteCoxeterGroup(FiniteComplexReflectionGroup):
    def __init__(self, W_types, index_set=None, hyperplane_index_set=None, reflection_index_set=None):
        r"""
        TESTS::

            sage: W = CoxeterGroups().example()
            sage: TestSuite(W).run()
        """
        W_types = tuple( tuple( W_type ) if isinstance(W_type,(list,tuple)) else W_type for W_type in W_types )
        cartan_types = []
        for W_type in W_types:
            W_type = CartanType(W_type)
            assert W_type.is_finite()
            assert W_type.is_irreducible()
            cartan_types.append( W_type )
        if len(W_types) == 1:
            cls = IrreducibleFiniteComplexReflectionGroup
        else:
            cls = FiniteComplexReflectionGroup
        cls.__init__(self, W_types, index_set=index_set,
                                    hyperplane_index_set=hyperplane_index_set,
                                    reflection_index_set=reflection_index_set,
                                    is_coxeter_group = True)
        N = self.nr_reflections()
        self._is_positive_root = [None] + [ True ] * N + [False]*N

    __iter__ = CoxeterGroups.ParentMethods.__iter__.__func__

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::
        
            sage: W = CoxeterGroup(['A',3],['B',2],['I',5],['I',6]); W
            Reducible finite Coxeter group of rank 9 and type A3 x B2 x I2(5) x G2
        """
        type_str = ''
        for W_type in self._type:
            type_str += self._irrcomp_repr_(W_type)
            type_str += ' x '
        type_str = type_str[:-3]
        return 'Reducible finite Coxeter group of rank {} and type {}'.format(self._rank,type_str)

    @cached_method
    def bipartite_index_set(self):
        r"""
        Return the bipartite index set of a finite real reflection group.

        EXAMPLES::

            sage: W = CoxeterGroup(["A",5])
            sage: W.bipartite_index_set()
            [[0, 2, 4], [1, 3]]

            sage: W = CoxeterGroup(["A",5],index_set=['a','b','c','d','e'])
            sage: W.bipartite_index_set()
            [['a', 'c', 'e'], ['b', 'd']]
        """
        index_family = self._index_set
        keys = index_family.keys()
        L,R = self._gap_group.BipartiteDecomposition().sage()
        L = [ i for i in keys if index_family[i]+1 in L ]
        R = [ i for i in keys if index_family[i]+1 in R ]
        return [L,R]

    def irreducible_components(self):
        r"""
        Return a list containing the irreducible components of ``self``
        as finite reflection groups.

        EXAMPLES::

            tba
        """
        if self.nr_irreducible_components() == 1:
            irr_comps = [self]
        else:
            irr_comps = []
            for W_type in self._type:
                W_str = [ W_type["series"], W_type["rank"] ]
                if W_type["series"] == "I":
                    W_str[1] = W_type["bond"]
                irr_comps.append( CoxeterGroup(W_str) )
        return irr_comps

    def cartan_type(self):
        if len(self._type) == 1:
            ct = self._type[0]
            return CartanType([ct['series'],ct['rank']])
        else:
            return [ W.cartan_type() for W in self.irreducible_components() ]

    @cached_method
    def cartan_matrix(self):
        from sage.matrix.constructor import matrix
        return matrix(self._gap_group.CartanMat().sage())

    def simple_root(self,i):
        return self.simple_roots()[self._index_set[i]]

    def positive_roots(self):
        return self.roots()[:self.nr_reflections()]

    def almost_positive_roots(self):
        return [ -beta for beta in self.simple_roots() ] + self.positive_roots()

    def root_to_reflection(self,root):
        Phi = self.roots()
        R = self.reflections()
        i = Phi.index(root)+1
        j = Phi.index(-root)+1
        for r in R:
            if r(i) == j:
                return r
        raise ValueError("there is a bug in root_to_reflection")

    def reflection_to_positive_root(self,r):
        Phi = self.roots()
        N = len(Phi)/2
        for i in range(1,N+1):
            if r(i) == i+N:
                return Phi[i-1]
        raise ValueError("there is a bug in reflection_to_positive_root")

    @cached_method
    def fundamental_weights(self):
        m = self.cartan_matrix().transpose().inverse()
        S = self.simple_roots()
        zero = S[0] - S[0]
        weights = [ sum( [ m[i,j] * S[j] for j in range(len(S)) ], zero )
                    for i in range(len(S)) ]
        for weight in weights:
            weight.set_immutable()
        return weights

    def fundamental_weight(self,i):
        return self.fundamental_weights()[self._index_set[i]]

    def permutahedron(self,coefficients=None):
        n = self.rank()
        weights = self.fundamental_weights()
        if coefficients is None:
            coefficients = [1]*n
        v = sum( coefficients[i] * weights[i] for i in range(n) )
        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron( vertices=[ v*(~w).as_matrix() for w in self] )

    class Element(FiniteComplexReflectionGroup.Element):

        has_descent = CoxeterGroups.ElementMethods.has_descent.__func__
        reduced_word = cached_in_parent_method(CoxeterGroups.ElementMethods.reduced_word.__func__)

        def has_left_descent(self, i):
            r"""
            EXAMPLES::

                sage: W = CoxeterGroup(["A",3])
                sage: s = W.simple_reflections()
                sage: (s[1]*s[2]).has_left_descent(1)
                True
                sage: (s[1]*s[2]).has_left_descent(2)
                False
            """
            W = self.parent()
            assert i in W.index_set()
            return not W._is_positive_root[self(W._index_set[i]+1)]

        def act_on_root(self,root):
            Phi = self.parent().roots()
            return Phi[ (~self)(Phi.index(root)+1)-1 ]

        def inversion_set(self):
            Phi_plus = set(self.parent().positive_roots())
            return [ root for root in Phi_plus if self.act_on_root(root) not in Phi_plus ]

class IrreducibleFiniteCoxeterGroup(FiniteCoxeterGroup, IrreducibleFiniteComplexReflectionGroup):
    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: for i in [2..7]: print CoxeterGroup(["I",i])
            Reducible finite Coxeter group of rank 2 and type A1 x A1
            Irreducible finite Coxeter group of rank 2 and type A2
            Irreducible finite Coxeter group of rank 2 and type B2
            Irreducible finite Coxeter group of rank 2 and type I2(5)
            Irreducible finite Coxeter group of rank 2 and type G2
            Irreducible finite Coxeter group of rank 2 and type I2(7)
        """
        type_str = self._irrcomp_repr_(self._type[0])
        return 'Irreducible finite Coxeter group of rank {} and type {}'.format(self._rank,type_str)

    class Element(FiniteCoxeterGroup.Element,IrreducibleFiniteComplexReflectionGroup.Element):
        pass

def CoxeterGroup(data, implementation="reflection", base_ring=None, index_set=None):
    """
    Return an implementation of the Coxeter group given by ``data``.

    INPUT:

    - ``data`` -- a Cartan type (or coercible into; see :class:`CartanType`)
      or a Coxeter matrix or graph

    - ``implementation`` -- (default: ``'reflection'``) can be one of
      the following:

      * ``'permutation'`` - as a permutation representation
      * ``'matrix'`` - as a Weyl group (as a matrix group acting on the
        root space); if this is not implemented, this uses the "reflection"
        implementation
      * ``'coxeter3'`` - using the coxeter3 package
      * ``'reflection'`` - as elements in the reflection representation; see
        :class:`~sage.groups.matrix_gps.coxeter_groups.CoxeterMatrixGroup`

    - ``base_ring`` -- (optional) the base ring for the ``'reflection'``
      implementation

    - ``index_set`` -- (optional) the index set for the ``'reflection'``
      implementation

    EXAMPLES:

    Now assume that ``data`` represents a Cartan type. If
    ``implementation`` is not specified, the reflection representation
    is returned::

        sage: W = CoxeterGroup(["A",2])
        sage: W
        Finite Coxeter group over Universal Cyclotomic Field with Coxeter matrix:
        [1 3]
        [3 1]

        sage: W = CoxeterGroup(["A",3,1]); W
        Coxeter group over Universal Cyclotomic Field with Coxeter matrix:
        [1 3 2 3]
        [3 1 3 2]
        [2 3 1 3]
        [3 2 3 1]

        sage: W = CoxeterGroup(['H',3]); W
        Finite Coxeter group over Universal Cyclotomic Field with Coxeter matrix:
        [1 3 2]
        [3 1 5]
        [2 5 1]

    We now use the ``implementation`` option::

        sage: W = CoxeterGroup(["A",2], implementation = "permutation") # optional - chevie
        sage: W                                                         # optional - chevie
        Permutation Group with generators [(1,3)(2,5)(4,6), (1,4)(2,3)(5,6)]
        sage: W.category()                       # optional - chevie
        Join of Category of finite permutation groups and Category of finite coxeter groups

        sage: W = CoxeterGroup(["A",2], implementation="matrix")
        sage: W
        Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)

        sage: W = CoxeterGroup(["H",3], implementation="matrix")
        sage: W
        Finite Coxeter group over Universal Cyclotomic Field with Coxeter matrix:
        [1 3 2]
        [3 1 5]
        [2 5 1]

        sage: W = CoxeterGroup(["H",3], implementation="reflection")
        sage: W
        Finite Coxeter group over Universal Cyclotomic Field with Coxeter matrix:
        [1 3 2]
        [3 1 5]
        [2 5 1]

        sage: W = CoxeterGroup(["A",4,1], implementation="permutation")
        Traceback (most recent call last):
        ...
        NotImplementedError: Coxeter group of type ['A', 4, 1] as permutation group not implemented

    We use the different options for the "reflection" implementation::

        sage: W = CoxeterGroup(["H",3], implementation="reflection", base_ring=RR)
        sage: W
        Finite Coxeter group over Real Field with 53 bits of precision with Coxeter matrix:
        [1 3 2]
        [3 1 5]
        [2 5 1]
        sage: W = CoxeterGroup([[1,10],[10,1]], implementation="reflection", index_set=['a','b'], base_ring=SR)
        sage: W
        Finite Coxeter group over Symbolic Ring with Coxeter matrix:
        [ 1 10]
        [10  1]

    TESTS::

        sage: W = groups.misc.CoxeterGroup(["H",3])
    """
    if implementation not in ["permutation", "matrix", "coxeter3",
                              "reflection", None]:
        raise ValueError("invalid type implementation")

    try:
        cartan_type = CartanType(data)
    except (TypeError, ValueError): # If it is not a Cartan type, try to see if we can represent it as a matrix group
        return CoxeterMatrixGroup(data, base_ring, index_set)

    if implementation is None:
        implementation = "matrix"

    if implementation == "reflection":
        return CoxeterMatrixGroup(cartan_type, base_ring, index_set)
    if implementation == "coxeter3":
        try:
            from sage.libs.coxeter3.coxeter_group import CoxeterGroup
        except ImportError:
            raise RuntimeError("coxeter3 must be installed")
        else:
            return CoxeterGroup(cartan_type)
    if implementation == "permutation" and assert_chevie_available() and \
       cartan_type.is_finite() and cartan_type.is_irreducible():
        return CoxeterGroupAsPermutationGroup(cartan_type)
    elif implementation == "matrix":
        if cartan_type.is_crystallographic():
            return WeylGroup(cartan_type)
        return CoxeterMatrixGroup(cartan_type, base_ring, index_set)

    raise NotImplementedError("Coxeter group of type {} as {} group not implemented".format(cartan_type, implementation))

