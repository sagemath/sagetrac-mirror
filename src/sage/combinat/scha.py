r"""
The Hopf algebra of super characters of upper triangular unimodular matrices
over a finite field.

Bases
-----

This module implements the supercharacter, superclass, powersum, homogeneous,
and elementary bases of this algebra. In addition, the following change of
bases are implemented:

- between the supercharacter basis and the superclass basis

- between the superclass basis and the powersum basis

- between the powersum basis and the elementary basis

- from the powersum basis to the homogeneous basis

- from the homogeneous basis to the superclass basis

Here is an illustration of these change of bases:

                            |<--------------------------->|
                            |                             |
    supercharacter <--> superclass <-- homogeneous <-- powersum <--> elementary


Products and coproducts
-----------------------

Products are explicitly implemented for supercharacters and powersums, and
coproducts are explicitly implemented for supercharacters and superclasses.
Products and coproducts in other bases are computed using a change of basis to
a basis in which it is implemented.

TODO:

- decide whether it is more efficient to convert to the powersum basis for
  products and to superclasses for coproducts

"""
#*****************************************************************************
#       Copyright (C) 2010 Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.all import Rings, Realizations
from sage.categories.all import tensor
from sage.categories.category_types import Category_realization
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.combinat import CombinatorialObject, CombinatorialClass
from sage.combinat.subset import Subsets
from sage.misc.cachefunc import cached_method
from sage.rings.finite_rings.constructor import FiniteField
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.rings.rational_field import QQ
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.finite_rings.constructor import GF
from sage.combinat.set_partition import SetPartitions
from sage.sets.set import Set, Set_object_enumerated
from sage.combinat.set_partition_ordered import OrderedSetPartitions
from sage.misc.misc_c import prod
from sage.functions.other import factorial
from sage.combinat.posets.lattices import FiniteLatticePoset, LatticePoset
from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject
from sage.sets.family import Family

##### Labelled set partitions
def LabelledSetPartitions(arg1, arg2=None):
    r"""
    Function that constructs families of labelled set partitions.

    EXAMPLES:

    With two arguments, the labelled set partitions of a given size are
    returned::

        sage: LabelledSetPartitions(3,2)
        Set partitions of [3] with arcs labelled by elements of Finite Field of size 2

    With one argument, the family of all set partitions is returned::

        sage: LabelledSetPartitions(3)
        Set partitions with arcs labelled by elements of Finite Field of size 3
    """
    if arg2 is None:
        return LabelledSetPartitions_all(q=arg1)
    else:
        return LabelledSetPartitions_n(n=arg1,q=arg2)

class LabelledSetPartition(CombinatorialObject):
    r"""
    Class modelling a set partition with labelled arcs. These objects are
    essentially modelled as a list of arcs.
    
    EXAMPLES:
    
    The set partition `\{ \{1,3,5\}, \{2,4\} \}` has arcs `(1,3)`, `(3,5)` and
    `(2,4)`. We label the arcs by 1, 3, and 4, respectively::

        sage: LabelledSetPartition(5, [(1,3,1), (3,5,3), (2,4,4)])
        [5, [(1, 3, 1), (2, 4, 4), (3, 5, 3)]]
    """
    def __init__(self, n, arcs):
        r"""
        EXAMPLES:
        
        The set partition `\{ \{1,3,5\}, \{2,4\} \}` has arcs `(1,3)`, `(3,5)` and
        `(2,4)`. We label the arcs by 1, 3, and 4, respectively::

            sage: LabelledSetPartition(5, [(1,3,1), (3,5,3), (2,4,4)])
            [5, [(1, 3, 1), (2, 4, 4), (3, 5, 3)]]
        """
        self._n = n
        _arcs = []
        for arc in arcs:
            if len(arc) == 3:
                _arcs.append(tuple(arc))
            elif len(arc) == 2:
                _arcs.append((arc[0],arc[1],1))
            else:
                raise ValueError, "arcs must be tuples of length 2 or 3"
        super(LabelledSetPartition, self).__init__([Integer(n),sorted(_arcs)])

    def __repr__(self):
        r"""
        EXAMPLES::
        
            sage: LabelledSetPartition(5, [(1,3,1), (3,5,3), (2,4,4)])
            [5, [(1, 3, 1), (2, 4, 4), (3, 5, 3)]]
        """
        return "[%s, %s]" % (self[0],self[1])

    def size(self):
        r"""
        Return the size of the underlying set (the size of the set partitioned
        by this set partition).

        EXAMPLES::

            sage: phi = LabelledSetPartition(5, [(1,3,1), (3,5,3), (2,4,4)])
            sage: phi.size()
            5
        """
        return self._n

    def arcs(self):
        r"""
        The labelled arcs, as a list.

        EXAMPLES::

            sage: phi = LabelledSetPartition(5, [(1,3,1), (3,5,3), (2,4,4)])
            sage: phi.arcs()
            [(1, 3, 1), (2, 4, 4), (3, 5, 3)]
        """
        return self[1][:]

    def arcs_dict(self):
        r"""
        The labelled arcs, as a dictionary. They keys are the arcs (as tuples)
        and the values are the labels.

        EXAMPLES::

            sage: phi = LabelledSetPartition(5, [(1,3,1), (3,5,3), (2,4,4)])
            sage: phi.arcs_dict()
            {(1, 3): 1, (2, 4): 4, (3, 5): 3}
        """
        return dict(((i,j),a) for (i,j,a) in self.arcs())

    def to_set_partition(self):
        r"""
        Return the corresponding unlabelled set partition.

        EXAMPLES::

            sage: phi = LabelledSetPartition(5, [(1,3,1), (3,5,3), (2,4,4)])
            sage: phi.to_set_partition()
            {{2, 4}, {1, 3, 5}}
        """
        from sage.graphs.graph import Graph
        G = Graph(dict((i,[]) for i in range(1,self._n+1)))
        for arc in self.arcs():
            G.add_edge(*arc)
        partition = G.connected_components()
        return Set(map(Set,partition))

    def __cmp__(self, other):
        r"""
        Total order on labelled set partitions.
        
        First, the size of the underlying sets are compared, then the number of
        arcs, and finally the (sorted) list of arcs.

        EXAMPLES::

            sage: phi = LabelledSetPartition(5, [(1,3,1), (3,5,3), (2,4,4)])
            sage: psi = LabelledSetPartition(5, [(1,3,1), (3,5,3)])
            sage: chi = LabelledSetPartition(7, [(1,3,1), (3,5,3), (2,4,4)])
            sage: cmp(phi,psi)
            1
            sage: cmp(psi,phi)
            -1
            sage: cmp(chi,phi)
            1
            sage: cmp(phi,chi)
            -1
        """
        s = cmp(self.size(),other.size())
        if s != 0:
            return s
        s = cmp(len(self.arcs()),len(other.arcs()))
        if s != 0:
            return s
        return -cmp(self.arcs(),other.arcs())

    def _latex_(self):
        from sage.misc.latex import latex
        return latex(self.to_set_partition())

class LabelledSetPartitions_n(UniqueRepresentation,CombinatorialClass):
    r"""
    Set partitions of [n] with arcs labelled by elements of the Finite Field of
    size q.

    EXAMPLES::

        sage: LabelledSetPartitions(3,2)
        Set partitions of [3] with arcs labelled by elements of Finite Field of size 2
    """
    def __init__(self, n, q):
        r"""
        Set partitions of [n] with arcs labelled by elements of the Finite Field of size q

        EXAMPLES::

            sage: LabelledSetPartitions(3,2)
            Set partitions of [3] with arcs labelled by elements of Finite Field of size 2
        """
        self._n = n
        self._q = q
        self._field = GF(q, 'a')
    def q(self):
        r"""
        Return the size of the finite field used for labels.

        EXAMPLES::

            sage: LabelledSetPartitions(3,2).q()
            2
        """
        return self._q
    def __repr__(self):
        r"""
        EXAMPLES::

            sage: LabelledSetPartitions(3,2)
            Set partitions of [3] with arcs labelled by elements of Finite Field of size 2
        """
        return "Set partitions of [%s] with arcs labelled by elements of %s" % (self._n, self._field)
    def __iter__(self):
        r"""
        Iterate through all the labelled set partitions in this family.

        EXAMPLES::

            sage: for phi in LabelledSetPartitions(3,2): print phi
            [3, [(1, 2, 1), (2, 3, 1)]]
            [3, [(2, 3, 1)]]
            [3, [(1, 3, 1)]]
            [3, [(1, 2, 1)]]
            [3, []]
            sage: for phi in LabelledSetPartitions(2,4): print phi
            [2, [(1, 2, a)]]
            [2, [(1, 2, a + 1)]]
            [2, [(1, 2, 1)]]
            [2, []]
        """
        for partition in SetPartitions(self._n):
            # convert set partition to arcs
            arcs = []
            for part in partition:
                sorted_part = sorted(part)
                for i in range(len(sorted_part)-1):
                    arcs.append(sorted_part[i:i+2])
            # label arcs with field elements
            for labels in self._field**len(arcs):
                if all(labels):
                    yield LabelledSetPartition(self._n,[(i,j,a) for ((i,j),a) in zip(arcs,labels)])

    def __contains__(self, phi):
        r"""
        Test whether ``phi`` belongs to this family.

        EXAMPLES::

            sage: phi = LabelledSetPartition(3, [(1, 3, 1)])
            sage: phi in LabelledSetPartitions(3, 2)
            True
            sage: phi in LabelledSetPartitions(4, 2)
            False
            sage: phi in LabelledSetPartitions(3, 3)
            True

        ::

            sage: F.<a> = GF(4,'a')
            sage: phi = LabelledSetPartition(3, [(1, 2, a + 1)])
            sage: phi in LabelledSetPartitions(3,4)
            True
            sage: phi in LabelledSetPartitions(3,3)
            False
            sage: phi in LabelledSetPartitions(3,2)
            False
        """
        return isinstance(phi, LabelledSetPartition) and phi.size() == self._n and all(arc[2] in self._field for arc in phi.arcs())

    @cached_method
    def rank(self, phi):
        r"""
        Return the position in which ``phi`` occurs in the iterator.


        EXAMPLES::

            sage: phi = LabelledSetPartition(3, [(1, 3, 1)])
            sage: LabelledSetPartitions(3,2).rank(phi)
            2

        ::

            sage: F.<a> = GF(4,'a')
            sage: phi = LabelledSetPartition(3, [(1, 2, a + 1)])
            sage: LabelledSetPartitions(3,4).rank(phi)
            16
        """
        return super(LabelledSetPartitions_n, self).rank(phi)

class LabelledSetPartitions_all(UniqueRepresentation,CombinatorialClass):
    r"""
    All set partitions with arcs labelled by elements of the Finite Field of
    size q.

    EXAMPLES::

        sage: LabelledSetPartitions(3)
        Set partitions with arcs labelled by elements of Finite Field of size 3
    """
    def __init__(self, q):
        r"""
        All set partitions with arcs labelled by elements of the Finite Field
        of size q.

        EXAMPLES::

            sage: LabelledSetPartitions(3)
            Set partitions with arcs labelled by elements of Finite Field of size 3
        """
        self._q = q
        self._field = GF(q, 'a')

    def q(self):
        r"""
        Return the size of the finite field used for labels.

        EXAMPLES::

            sage: LabelledSetPartitions(3).q()
            3
        """
        return self._q

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: LabelledSetPartitions(3)
            Set partitions with arcs labelled by elements of Finite Field of size 3
        """
        return "Set partitions with arcs labelled by elements of %s" % (self._field,)

    def __iter__(self):
        r"""
        Iterate through all the labelled set partitions.

        EXAMPLES::

            sage: it = iter(LabelledSetPartitions(3))
            sage: for i in range(7): print it.next()
            [0, []]
            [1, []]
            [2, [(1, 2, 1)]]
            [2, [(1, 2, 2)]]
            [2, []]
            [3, [(1, 2, 1), (2, 3, 1)]]
            [3, [(1, 2, 2), (2, 3, 1)]]

        ::

            sage: it = iter(LabelledSetPartitions(4))
            sage: for i in range(7): print it.next()
            [0, []]
            [1, []]
            [2, [(1, 2, a)]]
            [2, [(1, 2, a + 1)]]
            [2, [(1, 2, 1)]]
            [2, []]
            [3, [(1, 2, a), (2, 3, a)]]
        """
        from sage.sets.non_negative_integers import NonNegativeIntegers
        for n in NonNegativeIntegers():
            for lsp in LabelledSetPartitions(n,self.q()):
                yield lsp

    def __contains__(self, phi):
        r"""
        Test whether ``phi`` belongs to this family.

        EXAMPLES::

            sage: phi = LabelledSetPartition(3, [(1, 3, 1)])
            sage: phi in LabelledSetPartitions(4)
            True

        ::

            sage: F.<a> = GF(4,'a')
            sage: phi = LabelledSetPartition(3, [(1, 2, a + 1)])
            sage: phi in LabelledSetPartitions(4)
            True
            sage: phi in LabelledSetPartitions(3)
            False
        """
        return isinstance(phi, LabelledSetPartition) and all(arc[2] in self._field for arc in phi.arcs())

##### The Hopf algebra of supercharacters
class GenericCodeCollector(UniqueRepresentation, Parent):
    r"""
    This class collects some generic code that should be in other places in
    Sage.

    TODO: move these to the correct location in Sage
    """
    def __contains__(self, x):
        return any(x in realization for realization in self.realizations())

    def _an_element_(self):
        return self.supercharacter_basis().an_element()

    def one(self):
        return self.supercharacter_basis().one()

    def zero(self):
        return self.supercharacter_basis().one()

    def base_ring(self):
        return self._base

def SupercharacterHopfAlgebra(base_ring=None, q=None):
    r"""
    The Hopf algebra of supercharacters.

    EXAMPLES::

        sage: SupercharacterHopfAlgebra(2)
        Hopf algebra of supercharacters at q=2 over Rational Field
        sage: SupercharacterHopfAlgebra(3)
        Hopf algebra of supercharacters at q=3 over Cyclotomic Field of order 3 and degree 2
        sage: SupercharacterHopfAlgebra(SR, 4)
        Hopf algebra of supercharacters at q=4 over Symbolic Ring
    """
    if q is None:
        q = base_ring
        base_ring = None
    if q == 2:
        return SupercharacterHopfAlgebra_q2(q, base_ring)
    else:
        return SupercharacterHopfAlgebra_generic_q(q, base_ring)

class SupercharacterHopfAlgebra_generic_q(GenericCodeCollector):
    r"""
    The Hopf algebra of supercharacters for an arbitrary q.

    This is the entry point for the Hopf algebra of supercharacters: all bases,
    etc., are constructed from here. This class also deals with the change of
    bases routines (see the ``__init_extra__`` method).

    EXAMPLES::

        sage: SupercharacterHopfAlgebra(3)
        Hopf algebra of supercharacters at q=3 over Cyclotomic Field of order 3 and degree 2
        sage: SupercharacterHopfAlgebra(SR, 4)
        Hopf algebra of supercharacters at q=4 over Symbolic Ring
    """
    def __init__(self, q, base_ring = None):
        r"""
        The Hopf algebra of supercharacters for the parameter ``q``.

        EXAMPLES::

            sage: SupercharacterHopfAlgebra(3)
            Hopf algebra of supercharacters at q=3 over Cyclotomic Field of order 3 and degree 2
            sage: SupercharacterHopfAlgebra(SR, 4)
            Hopf algebra of supercharacters at q=4 over Symbolic Ring

        TESTS::

            sage: X = SupercharacterHopfAlgebra(q=3).X()
            sage: TestSuite(X).run(skip=["_test_len","_test_prod","_test_coproduct","_test_associativity"]) # long time
            sage: X = SupercharacterHopfAlgebra(q=4).X()
            sage: TestSuite(X).run(skip=["_test_len","_test_prod","_test_coproduct","_test_associativity"]) # long time
        """
        q = Integer(q)
        if base_ring is None:
            if q == 2:
                base_ring = QQ
            else:
                base_ring = CyclotomicField(q)
        assert(base_ring in Rings())
        self._q = q
        self._field = GF(q, 'a')
        # TODO: the following line won't be needed when CategoryObject won't override base_ring
        self._base = base_ring 
        # TODO: Why doesn't GradedHopfAlgebrasWithBasis work here?
        #Parent.__init__(self, category = GradedHopfAlgebrasWithBasis(base_ring).WithRealizations())
        Parent.__init__(self, category = GradedHopfAlgebras(base_ring).WithRealizations())

    def q(self):
        r"""
        Return the parameter ``q``.

        EXAMPLES::

            sage: SupercharacterHopfAlgebra(3).q()
            3
        """
        return self._q

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: SupercharacterHopfAlgebra(3)._repr_()
            'Hopf algebra of supercharacters at q=3 over Cyclotomic Field of order 3 and degree 2'
        """
        return "Hopf algebra of supercharacters at q=%s over %s"%(self.q(), self.base_ring())

    def supercharacter_basis(self): 
        r"""
        The supercharacter basis of the Hopf algebra of supercharacters.

        EXAMPLES:

        There are various synonyms for this basis::

            sage: scha = SupercharacterHopfAlgebra(2)
            sage: scha.supercharacter_basis()
            Hopf algebra of supercharacters at q=2 over Rational Field on the supercharacter basis
            sage: scha.X()
            Hopf algebra of supercharacters at q=2 over Rational Field on the supercharacter basis

        This basis is also defined for `q > 2`::

            sage: scha = SupercharacterHopfAlgebra(3)
            sage: scha.supercharacter_basis()
            Hopf algebra of supercharacters at q=3 over Cyclotomic Field of order 3 and degree 2 on the supercharacter basis
            sage: scha.X()
            Hopf algebra of supercharacters at q=3 over Cyclotomic Field of order 3 and degree 2 on the supercharacter basis
        """
        return SupercharacterBasis(self)

    X = supercharacter_basis

    def superclass_basis(self):
        r"""
        The superclass basis of the Hopf algebra of supercharacters.

        EXAMPLES:

        There are various synonyms for this basis::

            sage: scha = SupercharacterHopfAlgebra(2)
            sage: scha.superclass_basis()
            Hopf algebra of supercharacters at q=2 over Rational Field on the superclass basis
            sage: scha.K()
            Hopf algebra of supercharacters at q=2 over Rational Field on the superclass basis

        This basis is also defined for `q > 2`::

            sage: scha = SupercharacterHopfAlgebra(3)
            sage: scha.superclass_basis()
            Hopf algebra of supercharacters at q=3 over Cyclotomic Field of order 3 and degree 2 on the superclass basis
            sage: scha.kappa()
            Hopf algebra of supercharacters at q=3 over Cyclotomic Field of order 3 and degree 2 on the superclass basis
            sage: scha.K()
            Hopf algebra of supercharacters at q=3 over Cyclotomic Field of order 3 and degree 2 on the superclass basis
        """
        return SuperclassBasis(self)

    K = kappa = superclass_basis

    def __init_extra__(self):
        r"""
        Register coercions between the supercharacter and superclass bases.
        """
        category = Bases(self)
        X = self.supercharacter_basis()
        K = self.superclass_basis()

        # supercharacter to superclass basis
        X.module_morphism(K._supercharacter_to_superclass_on_basis, codomain=K, category=category).register_as_coercion()

        # superclass to supercharacter basis
        K.module_morphism(K._superclass_to_supercharcter_on_basis, codomain=X, category=category).register_as_coercion()

    _shorthands = set(['X', 'K'])

    def inject_shorthands(self, shorthands=_shorthands):
        r"""
        Imports shorthands into the global namespace.

        INPUT:

         - shorthands - a list (or iterable) of strings (default: ['X', 'K'])

        EXAMPLES::

            sage: scha = SupercharacterHopfAlgebra(3)
            sage: scha.inject_shorthands()

            sage: X.one()
            X0[]
            sage: K.one()
            K0[]

        """
        from sage.misc.misc import inject_variable
        for shorthand in shorthands:
            assert shorthand in self._shorthands
            inject_variable(shorthand, getattr(self, shorthand)())

    ##### supercharacter table
    def supercharacter_table(self, n):
        r"""
        The supercharacter table for degree ``n`` (as a matrix).

        EXAMPLES::

            sage: scha = SupercharacterHopfAlgebra(q=2)
            sage: scha.supercharacter_table(3)
            [ 1 -1  1 -1  1]
            [-1 -1  1  1  1]
            [ 0  0 -2  0  2]
            [-1  1  1 -1  1]
            [ 1  1  1  1  1]

        ::

            sage: scha = SupercharacterHopfAlgebra(q=3)
            sage: scha.supercharacter_table(2)
            [     zeta3 -zeta3 - 1          1]
            [-zeta3 - 1      zeta3          1]
            [         1          1          1]
        """
        return SupercharacterTable(self.supercharacter_basis()).table(n)

class SupercharacterHopfAlgebra_q2(SupercharacterHopfAlgebra_generic_q):
    r"""
    The Hopf algebra of supercharacters for `q=2`. For `q=2`, we have two
    additional bases: the powersum basis and that homogeneous basis.

    EXAMPLES::

        sage: SupercharacterHopfAlgebra(2)
        Hopf algebra of supercharacters at q=2 over Rational Field
        sage: SupercharacterHopfAlgebra(SR, 3)
        Hopf algebra of supercharacters at q=3 over Symbolic Ring
    """
    def powersum_basis(self):
        r"""
        The powersum basis of the Hopf algebra of supercharacters.

        EXAMPLES:

        There are various synonyms for this basis::

            sage: scha = SupercharacterHopfAlgebra(2)
            sage: scha.powersum_basis()
            Hopf algebra of supercharacters at q=2 over Rational Field on the powersum basis
            sage: scha.P()
            Hopf algebra of supercharacters at q=2 over Rational Field on the powersum basis

        Note that this basis is note defined for `q > 2`.
        """
        return PowersumBasis(self)

    P = powersum_basis

    def homogeneous_basis(self):
        r"""
        The homogeneous basis of the Hopf algebra of supercharacters.

        EXAMPLES:

        There are various synonyms for this basis::

            sage: scha = SupercharacterHopfAlgebra(2)
            sage: scha.homogeneous_basis()
            Hopf algebra of supercharacters at q=2 over Rational Field on the homogeneous basis
            sage: scha.H()
            Hopf algebra of supercharacters at q=2 over Rational Field on the homogeneous basis

        Note that this basis is note defined for `q > 2`.
        """
        return HomogeneousBasis(self)

    H = homogeneous_basis

    def elementary_basis(self):
        r"""
        The elementary basis of the Hopf algebra of supercharacters.

        EXAMPLES:

        There are various synonyms for this basis::

            sage: SupercharacterHopfAlgebra(2).E()
            Hopf algebra of supercharacters at q=2 over Rational Field on the elementary basis

        Note that this basis is note defined for `q > 2`.
        """
        return ElementaryBasis(self)

    E = elementary_basis

    def __init_extra__(self):
        r"""
        Register coercions: 

        -  between the powersum basis and the superclass basis
        
        -  from the homogeneous basis to the superclass basis

        -  from the powersum basis to the homogeneous basis

        -  between the powersum basis and the elementary basis

        TESTS:

        This example is taken from the paper of Rosas and Sagan::

            sage: scha = SupercharacterHopfAlgebra(2)
            sage: K = scha.K() 
            sage: H = scha.H()
            sage: phi = LabelledSetPartition(4,[(1,3,1),(2,4,1)])
            sage: K(H[phi])
            K[1|2|3|4] + K[12|3|4] + 2*K[123|4] + 4*K[1234] + 2*K[124|3] + K[12|34] + 2*K[13|2|4] + 4*K[13|24] + 2*K[134|2] + K[14|2|3] + K[14|23] + K[1|23|4] + 2*K[1|234] + 2*K[1|24|3] + K[1|2|34]
        """
        super(SupercharacterHopfAlgebra_q2,self).__init_extra__()
        K = self.superclass_basis()
        P = self.powersum_basis()
        H = self.homogeneous_basis()
        E = self.elementary_basis()

        category = Bases(self)

        # powersum to superclass
        P2K = P.module_morphism(P._powersum_to_superclass_on_basis, codomain=K, category=category, triangular='lower', unitriangular=True, cmp=cmp)
        P2K.register_as_coercion()

        # superclass to powersum
        (~P2K).register_as_coercion()

        # homogeneous basis to superclass basis
        H.module_morphism(H._homogeneous_to_superclass_on_basis, codomain=K, category=category).register_as_coercion()
        
        # powersum to homogeneous
        P.module_morphism(H._powersum_to_homogeneous_on_basis, codomain=H, category=category).register_as_coercion()

        # elementary to powersum
        E.module_morphism(E._elementary_to_powersum_on_basis, codomain=P, category=category).register_as_coercion()

        # powersum to elementary
        P.module_morphism(E._powersum_to_elementary_on_basis, codomain=E, category=category).register_as_coercion()

    _shorthands = set(['X','K','P','H','E'])
    def inject_shorthands(self, shorthands=_shorthands):
        r"""
        Imports shorthands into the global namespace.

        INPUT:

         - shorthands - a list (or iterable) of strings (default: ['X', 'K', 'P', 'H', 'E'])

        EXAMPLES::

            sage: scha = SupercharacterHopfAlgebra(2)
            sage: scha.inject_shorthands()
            sage: X.one()
            X[]
            sage: K.one()
            K[]
            sage: P.one()
            P[]
            sage: H.one()
            H[]
        """
        from sage.misc.misc import inject_variable
        for shorthand in shorthands:
            assert shorthand in self._shorthands
            inject_variable(shorthand, getattr(self, shorthand)())


##### category of bases of the hopf algebra of supercharacters
class Bases(Category_realization):
    r"""
    The category of bases of the Hopf algebra of supercharacters.

    This class collects code common to all the various bases that get
    implemented (these methods can be found in the attribute
    :method:`ParentMethods`), in addition to methods for elements of the
    algebra (see :method:`ElementMethods`). In most cases, these are just
    default implementations that will get specialized in a basis.

    EXAMPLES::

        sage: from sage.combinat.scha import Bases
        sage: Bases(SupercharacterHopfAlgebra(2))
        The category of bases of Hopf algebra of supercharacters at q=2 over Rational Field

    ::

        sage: Bases(SupercharacterHopfAlgebra(3))
        The category of bases of Hopf algebra of supercharacters at q=3 over Cyclotomic Field of order 3 and degree 2

    """
    def super_categories(self):
        return [Realizations(self.base()), GradedHopfAlgebrasWithBasis(self.base().base_ring())]

    class ParentMethods:
        r"""
        This class collects code common to all the various bases. In most
        cases, these are just default implementations that will get
        specialized in a basis.
        """
        def _repr_(self):
            return "%s %s"%(self.realization_of(), "on the %s basis"%(self._basis_name))

        def q(self):
            return self.realization_of().q()

        @cached_method
        def algebra_generators(self):
            return self.basis()

        @cached_method
        def counit_on_basis(self, phi):
            if phi.size() > 0:
                return self.base_ring().zero()
            else:
                return self.base_ring().one()

        @cached_method
        def one_basis(self):
            return LabelledSetPartition(0,[])

        @cached_method
        def some_elements(self):
            return map(self, [LabelledSetPartition(0,[]),
                LabelledSetPartition(1,[]),
                LabelledSetPartition(2,[(1,2)]),
                LabelledSetPartition(2,[]),
                ]) + [self.an_element()]

        @cached_method
        def product_on_basis(self, phi, psi):
            r"""
            Default implementation: compute in the powersum basis for
            ``q==2`` and in the supercharacter basis otherwise.

            TESTS::

                sage: X = SupercharacterHopfAlgebra(q=2).X()
                sage: phi = LabelledSetPartition(3, [(1,3,1)])
                sage: psi = LabelledSetPartition(2, [])
                sage: X.product_on_basis(phi, psi)
                X[13|2|4|5]

            ::

                sage: K = SupercharacterHopfAlgebra(q=3).K()
                sage: psi = LabelledSetPartition(1, [])
                sage: K.product_on_basis(psi, psi)
                K2[] + K2[(1, 2, 1)] + K2[(1, 2, 2)]
            """
            if self.q() == 2:
                B = self.realization_of().powersum_basis()
            else:
                B = self.realization_of().supercharacter_basis()
            return self(B.product(B(self[phi]), B(self[psi])))

        @cached_method
        def coproduct_on_basis(self, phi):
            r"""
            Default implementation: compute in the superclass basis.

            TESTS::

                sage: X = SupercharacterHopfAlgebra(q=2).X()
                sage: psi = LabelledSetPartition(2, [])
                sage: X.coproduct_on_basis(psi)
                X[] # X[1|2] + 2*X[1] # X[1] + X[1|2] # X[]

            ::

                sage: X = SupercharacterHopfAlgebra(q=3).X()
                sage: psi = LabelledSetPartition(2, [])
                sage: X.coproduct_on_basis(psi)
                X0[] # X2[] + 2*X1[] # X1[] + X2[] # X0[]
            """
            K = self.realization_of().superclass_basis()
            coprod = K(self[phi]).coproduct()
            morphism = K.tensor_square().module_morphism(lambda (phi,psi): tensor([self(K[phi]),self(K[psi])]), codomain=self.tensor_square())
            return morphism(coprod)

        @cached_method
        def antipode_on_basis(self, phi):
            r"""
            Compute the antipode of the basis element indexed by ``phi``.

            This is a default implementation that uses the recursive
            construction of the antipode of a graded connected Hopf algebra.

            EXAMPLES::

                sage: X = SupercharacterHopfAlgebra(2).X()
                sage: for n in range(4):
                ...       for phi in X.basis(n).keys():
                ...           print "S(%s) = %s" % (X.basis()[phi],X.antipode_on_basis(phi))
                S(X[]) = X[]
                S(X[1]) = -X[1]
                S(X[12]) = 2*X[1|2] - X[12]
                S(X[1|2]) = X[1|2]
                S(X[123]) = -4*X[1|2|3] + 2*X[12|3] - X[123] + 2*X[1|23]
                S(X[1|23]) = -2*X[1|2|3] + X[12|3]
                S(X[13|2]) = -8*X[1|2|3] + 4*X[12|3] - X[13|2] + 4*X[1|23]
                S(X[12|3]) = -2*X[1|2|3] + X[1|23]
                S(X[1|2|3]) = -X[1|2|3]

            Here is the matrix of antipode in degree 4 for `q=2`::

                sage: X = SupercharacterHopfAlgebra(2).X()
                sage: m = []
                sage: for phi in X.basis(4).keys():
                ...       Sphi = X.antipode_on_basis(phi)
                ...       m.append([Sphi.coefficient(psi) for psi in X.basis(4).keys()])
                sage: matrix(m)
                [ -1   2   0   0   2   2   0   0  -4   0  -4   0   0  -4   8]
                [  0   0   0   0   1   0   0   0   0   0  -2   0   0  -2   4]
                [  0   3  -1   0   3   2   0   0  -7   1 -10   0   1  -7  16]
                [  0   3   0  -1   3   2   0   0  -7   1 -10   0   1  -7  16]
                [  0   1   0   0   0   0   0   0  -2   0  -2   0   0   0   4]
                [  0   0   0   0   0   1   0   0  -2   0   0   0   0  -2   4]
                [  0   4   0   0   4  10  -1   0 -22   4 -24   0   4 -22  42]
                [  0   2   0   0   2  10   0  -1 -24   6 -28   0   6 -24  46]
                [  0   0   0   0   0   0   0   0   0   0   0   0   0  -1   2]
                [  0   0   0   0   0   0   0   0   0   0  -4   0   1  -4   8]
                [  0   0   0   0   0   0   0   0   0   0  -1   0   0   0   2]
                [  0   0   0   0   0   6   0   0 -18   6 -24  -1   6 -18  38]
                [  0   0   0   0   0   0   0   0  -4   1  -4   0   0   0   8]
                [  0   0   0   0   0   0   0   0  -1   0   0   0   0   0   2]
                [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   1]
                sage: matrix(m)**2 == 1
                True

            We compute a few antipodes for `q=3`:
            
                sage: X = SupercharacterHopfAlgebra(3).X()
                sage: X.antipode_on_basis(LabelledSetPartition(0,[]))
                X0[]
                sage: X.antipode_on_basis(LabelledSetPartition(1,[]))
                (-1)*X1[]
                sage: X.antipode_on_basis(LabelledSetPartition(2,[]))
                X2[]
                sage: X.antipode_on_basis(LabelledSetPartition(2,[(1,2,1)]))
                2*X2[] + (-1)*X2[(1, 2, 1)]
            """
            if phi.size() == 0:
                return self.basis()[phi]
            near_delta_bar = self.coproduct_on_basis(phi) - tensor([self.basis()[phi],self.one()])
            res = self.zero()
            for (mono, coeff) in near_delta_bar.monomial_coefficients().iteritems():
                res += coeff*self.product(self.antipode_on_basis(mono[0]),self.basis()[mono[1]])
            return -res

        def omega(self, x):
            r"""
            The image of ``x`` under `\omega`.

            `\omega` is the endomorphism of the Hopf algebra of
            supercharacters defined by mapping the elementary basis element
            indexed by the set partition ``phi`` to the homogeneous basis
            element indexed by ``phi``.

            .. note::

                The map :method:`omega1` is another implementation of this
                endomorphism.

            EXAMPLES::

                sage: E = SupercharacterHopfAlgebra(2).elementary_basis()
                sage: H = SupercharacterHopfAlgebra(2).homogeneous_basis()
                sage: x = E.an_element(); x
                5*E[1|2] - 3*E[13|2]
                sage: H(E.omega(x))
                5*H[1|2] - 3*H[13|2]

            We confirm that the powersum basis elements are eigenvectors of
            `\omega` with eigenvalues `\pm1`::

                sage: P = SupercharacterHopfAlgebra(2).powersum_basis()
                sage: for a in P.basis(4):
                ...       P.omega(a)
                -P[1234]
                P[1|234]
                P[134|2]
                P[124|3]
                P[123|4]
                P[12|34]
                P[13|24]
                P[14|23]
                -P[1|2|34]
                -P[1|24|3]
                -P[1|23|4]
                -P[14|2|3]
                -P[13|2|4]
                -P[12|3|4]
                P[1|2|3|4]
            """
            E = self.realization_of().elementary_basis()
            return self(E.omega(E(x)))

        def omega1(self, x):
            r"""
            The image of ``x`` under `\omega_1`.

            `\omega_1` is the endomorphism of the Hopf algebra of
            supercharacters defined on the powersum basis by
            `\omega_1(P_\phi) = (-1)^{|\phi|-\ell(phi)} P_\phi`.

            .. note::

                The map :method:`omega` is another implementation of this
                same endomorphism.

            EXAMPLES::

                sage: X = SupercharacterHopfAlgebra(q=2).X()
                sage: x = X[4, [(3,4)]]
                sage: X.omega1(x)
                2*X[1|2|3|4] - X[1|2|34]
            """
            P = self.realization_of().powersum_basis()
            return self(P.omega1(P(x)))

        @lazy_attribute
        def omega2(self):
            r"""
            The endomorphism of the Hopf algebra of supercharacters
            obtained by linearly extending the `\omega_2` map.

            `\omega_2` is the endomorphism of the Hopf algebra of
            supercharacters defined by taking the inner tensor product
            with `\sum_n X_{\{[n]\}}`.

            EXAMPLES::

                sage: X = SupercharacterHopfAlgebra(q=2).X()
                sage: x = X[4, [(3,4)]]
                sage: y = X[4, [(1,3),(2,4)]]
                sage: X.omega2(x - 3*y)
                X[123|4] - 3*X[13|24]
            """
            return self._module_morphism(self.omega2_on_basis, codomain=self)

        @cached_method
        def omega2_on_basis(self, phi):
            r"""
            The endomorphism of ``self`` defined by taking the inner
            tensor product with `\sum_n X_{\{[n]\}}`.

            EXAMPLES::

                sage: X = SupercharacterHopfAlgebra(q=2).X()
                sage: for m in range(4):
                ...       for phi in X.basis(m).keys():
                ...           print "%s --> %s" % (X[phi], X.omega2_on_basis(phi))
                X[] --> X[]
                X[1] --> X[1]
                X[12] --> X[1|2]
                X[1|2] --> X[12]
                X[123] --> X[1|2|3]
                X[1|23] --> X[12|3]
                X[13|2] --> X[13|2]
                X[12|3] --> X[1|23]
                X[1|2|3] --> X[123]
            """
            X = self.realization_of().supercharacter_basis()
            x = X[LabelledSetPartitions(phi.size(),2).first()]
            return self.inner_tensor_product(self[phi], self(x))

        def omega3(self, x):
            r"""
            The image of ``x`` under `\omega_3`.

            `\omega_3` is the endomorphism of the Hopf algebra of
            supercharacters defined on the superclass basis by
            `\omega_3(K_\phi) = (-1)^{|\phi|-\ell(phi)} K_\phi`.

            EXAMPLES::

                sage: X = SupercharacterHopfAlgebra(q=2).X()
                sage: x = X[3, [(1,3)]]
                sage: X.omega3(x)
                1/2*X[1|2|3] + 1/2*X[12|3] + 1/2*X[123] + 1/2*X[1|23]
            """
            K = self.realization_of().superclass_basis()
            return self(K.omega3(K(x)))

        def inner_tensor_product(self, x, y):
            r"""
            Default implementation: compute in the supercharacter basis.
            """
            X = self.realization_of().supercharacter_basis()
            return self(X.inner_tensor_product(X(x),X(y)))

        def _test_antipode(self,**options):
            r"""
            Tests that the antipode satisfies its defining equations:

                product * (S # 1) * coproduct == unit * counit

                product * (1 # S) * coproduct == unit * counit

            and that it is an involution.

            TESTS::

                sage: X = SupercharacterHopfAlgebra(2).X()
                sage: X._test_antipode()
            """
            tester = self._tester(**options)
            T = self.tensor(self)
            coprod = self.coproduct
            S_Id = T.module_morphism(lambda (phi,psi): tensor([self.antipode_on_basis(phi),self.basis()[psi]]), codomain=T)
            Id_S = T.module_morphism(lambda (phi,psi): tensor([self.basis()[phi],self.antipode_on_basis(psi)]), codomain=T)
            prod = T.module_morphism(lambda (phi,psi): self.product_on_basis(phi,psi), codomain=self)
            for x in tester.some_elements():
                tester.assert_(prod(S_Id(coprod(x)))==self.counit(x)*self.one())
                tester.assert_(prod(Id_S(coprod(x)))==self.counit(x)*self.one())
                tester.assert_(self.antipode(self.antipode(x)) == x)

        def _test_coproduct(self, **options):
            r"""
            Run some tests on the coproduct.

            TESTS::

                sage: X = SupercharacterHopfAlgebra(2).X()
                sage: X._test_coproduct()

            ::

                sage: X = SupercharacterHopfAlgebra(3).X()
                sage: X._test_coproduct()
            """
            tester = self._tester(**options)
            coprod = self.coproduct
            # test the zero
            zero = self.zero()
            tester.assert_(coprod(zero)==tensor([zero,zero]))
            # test coproduct is an algebra morphism
            for x in tester.some_elements():
                for y in tester.some_elements():
                    tester.assert_(coprod(x*y) == coprod(x)*coprod(y))

        @cached_method
        def rho_on_basis(self, phi):
            r"""
            The projection of the basis element indexed by ``phi`` into the
            Hopf algebra of symmetric functions.

            EXAMPLES::

                sage: SupercharacterHopfAlgebra(q=2).inject_shorthands()
                sage: phi = LabelledSetPartition(5,[(1,3),(2,4)])
                sage: K.rho_on_basis(phi)
                2*m[2, 2, 1]
                sage: P.rho_on_basis(phi)
                p[2, 2, 1]
                sage: E.rho_on_basis(phi)
                4*e[2, 2, 1]
                sage: H.rho_on_basis(phi)
                4*h[2, 2, 1]
                sage: X.rho_on_basis(phi)
                480*m[1, 1, 1, 1, 1] + 72*m[2, 1, 1, 1] - 8*m[2, 2, 1] - 8*m[3, 1, 1] + 8*m[3, 2]
            """
            if self.q() != 2:
                raise NotImplementedError

            from sage.combinat.sf import sfa
            from sage.combinat.partition import Partition
            pi = phi.to_set_partition()
            la = Partition(sorted(map(len,pi),reverse=True))
            shca = self.realization_of()
            R = self.base_ring()

            if self is shca.superclass_basis():
                coeff = prod(factorial(i) for i in la.to_exp())
                return coeff * sfa.SFAMonomial(R)[la]
            elif self is shca.powersum_basis():
                return sfa.SFAPower(R)[la]
            elif self is shca.elementary_basis():
                coeff = prod(factorial(i) for i in la)
                return coeff * sfa.SFAElementary(R)[la]
            elif self is shca.homogeneous_basis():
                coeff = prod(factorial(i) for i in la)
                return coeff * sfa.SFAHomogeneous(R)[la]
            else:
                return shca.superclass_basis()(self[phi]).rho()

        to_symmetric_function_on_basis = rho_on_basis

        @lazy_attribute
        def rho(self):
            r"""
            The projection onto the Hopf algebra of symmetric functions.

            EXAMPLES::

                sage: SupercharacterHopfAlgebra(q=2).inject_shorthands()
                sage: phi = LabelledSetPartition(5,[(1,3),(2,4)])
                sage: K.rho(K[phi])
                2*m[2, 2, 1]
                sage: P.rho(P[phi])
                p[2, 2, 1]
                sage: E.rho(E[phi])
                4*e[2, 2, 1]
                sage: H.rho(H[phi])
                4*h[2, 2, 1]
                sage: X.rho(X[phi])
                480*m[1, 1, 1, 1, 1] + 72*m[2, 1, 1, 1] - 8*m[2, 2, 1] - 8*m[3, 1, 1] + 8*m[3, 2]
            """
            if self.q() != 2:
                raise NotImplementedError
            from sage.combinat.sf import sfa
            shca = self.realization_of()
            R = self.base_ring()

            codomain = None
            if self is shca.superclass_basis():
                codomain = sfa.SFAMonomial(R)
            elif self is shca.powersum_basis():
                codomain = sfa.SFAPower(R)
            elif self is shca.elementary_basis():
                codomain = sfa.SFAElementary(R)
            elif self is shca.homogeneous_basis():
                codomain = sfa.SFAHomogeneous(R)

            if codomain is not None:
                return self.module_morphism(self.rho_on_basis, codomain=codomain)
            else:
                codomain = sfa.SFAMonomial(R)
                K = shca.superclass_basis()
                morph1 = K.module_morphism(K.rho_on_basis, codomain=codomain)
                morph2 = self.module_morphism(K.coerce_map_from(self).on_basis(), codomain=K)
                return morph1 * morph2

        to_symmetric_function = rho

    class ElementMethods:
        r"""
        Methods for elements of the algebra.
        """
        def restriction(self, A):
            r"""
            Restriction of the corresponding supercharacter to the set
            ``A``.

            This is a default implementation: the restriction is computed
            in the supercharacter basis.

            EXAMPLES::

                sage: X = SupercharacterHopfAlgebra(q=4).X()
                sage: a = X[7, [(2,5,1)]] + 3 * X[7, [(5,7,1)]]
                sage: a.restriction(Set([2,3,4,5]))
                12*X4[] + X4[(1, 4, 1)]
            """
            return self.parent().restriction(self,A)

        def inner_tensor_product(self, x):
            r"""
            The inner tensor product of this element with the element
            ``x``.

            This is a default implementation: the restriction is computed
            in the supercharacter basis.

            EXAMPLES::

                sage: X = SupercharacterHopfAlgebra(5).X()
                sage: y = X[5, [(2,3,2)]]
                sage: x = X[5, [(1,4,1)]]
                sage: y.inner_tensor_product(x)
                X5[(1, 4, 1), (2, 3, 2)]
                sage: z = X[5, [(1,4,2)]]
                sage: z.inner_tensor_product(x)
                9*X5[(1, 4, 3)] + 4*X5[(1, 4, 3), (2, 3, 1)] + 4*X5[(1, 4, 3), (2, 3, 2)] + 4*X5[(1, 4, 3), (2, 3, 3)] + 4*X5[(1, 4, 3), (2, 3, 4)]
            """
            return self.parent().inner_tensor_product(self, x)

        def omega(self):
            r"""
            The image of this element under `\omega`.

            `\omega` is the endomorphism of the Hopf algebra of
            supercharacters defined by mapping the elementary basis element
            indexed by the set partition ``phi`` to the homogeneous basis
            element indexed by ``phi``.

            EXAMPLES::

                sage: E = SupercharacterHopfAlgebra(2).elementary_basis()
                sage: H = SupercharacterHopfAlgebra(2).homogeneous_basis()
                sage: a = E.an_element()
                sage: a
                5*E[1|2] - 3*E[13|2]
                sage: H(a.omega())
                5*H[1|2] - 3*H[13|2]
            """
            return self.parent().omega(self)

        def omega1(self):
            r"""
            The image of this element under `\omega_1`.

            `\omega_1` is the endomorphism of the Hopf algebra of
            supercharacters defined on the powersum basis by
            `\omega_1(P_\phi) = (-1)^{|\phi|-\ell(phi)} P_\phi`.

            EXAMPLES::

                sage: P = SupercharacterHopfAlgebra(q=2).powersum_basis()
                sage: a = P[4, [(3,4)]]
                sage: a
                P[1|2|34]
                sage: a.omega1()
                -P[1|2|34]
            """
            return self.parent().omega1(self)

        def omega2(self):
            r"""
            The image of this element under `\omega_2`.

            `\omega_2` is the endomorphism of the Hopf algebra of
            supercharacters defined by taking the inner tensor product
            with `\sum_n X_{\{[n]\}}`.

            EXAMPLES::

                sage: X = SupercharacterHopfAlgebra(q=2).X()
                sage: x = X[4, [(3,4)]]
                sage: y = X[4, [(1,3),(2,4)]]
                sage: (x - 3*y).omega2()
                X[123|4] - 3*X[13|24]
            """
            return self.parent().omega2(self)

        def omega3(self):
            r"""
            The image of this element under `\omega_3`.

            `\omega_3` is the endomorphism of the Hopf algebra of
            supercharacters defined on the superclass basis by
            `\omega_3(K_\phi) = (-1)^{|\phi|-\ell(phi)} K_\phi`.

            EXAMPLES::

                sage: K = SupercharacterHopfAlgebra(q=2).K()
                sage: x = K[4, [(3,4,1)]]
                sage: y = K[4, [(1,3,1),(2,4,1)]]
                sage: (x + 3*y).omega3()
                3*K[13|24] - K[1|2|34]
            """
            return self.parent().omega3(self)

        def rho(self):
            r"""
            The projection of this element into the Hopf algebra of
            symmetric functions.

            EXAMPLES::

                sage: SupercharacterHopfAlgebra(q=2).inject_shorthands()
                sage: phi = LabelledSetPartition(5,[(1,3),(2,4)])
                sage: K[phi].rho()
                2*m[2, 2, 1]
                sage: P[phi].rho()
                p[2, 2, 1]
                sage: E[phi].rho()
                4*e[2, 2, 1]
                sage: H[phi].rho()
                4*h[2, 2, 1]
                sage: X[phi].rho()
                480*m[1, 1, 1, 1, 1] + 72*m[2, 1, 1, 1] - 8*m[2, 2, 1] - 8*m[3, 1, 1] + 8*m[3, 2]
            """
            return self.parent().rho(self)

        to_symmetric_function = rho

class SupercharacterHopfAlgebraBasis(CombinatorialFreeModule):
    r"""
    Code shared by all bases of the Hopf algebra of supercharacters.
    """
    def __init__(self, SCHA, category=None):
        if category is None:
            category = Bases(SCHA)
        self._field = SCHA._field
        self.abstract_algebra = SCHA
        CombinatorialFreeModule.__init__(self, SCHA.base_ring(),
                LabelledSetPartitions(SCHA.q()), category=category,
                prefix=self._prefix)

    def __getitem__(self, key):
        r"""
        Override the ``__getitem__`` method to allow passing of arguments to
        LabelledSetPartition.

        EXAMPLES:

        Define a basis element by specifying the degree and a list of
        arcs::

            sage: X = SupercharacterHopfAlgebra(2).X()
            sage: X[0,[]]
            X[]
            sage: X[4,[]]
            X[1|2|3|4]

        Alternatively, define a basis element by passing data that
        describes a set partition (the labels are taken to be 1)::

            sage: X = SupercharacterHopfAlgebra(2).X()
            sage: X[[1,2,3,4]]
            X[1234]
            sage: X[[1,3],[2]]
            X[13|2]
            sage: X[[1],[2],[3],[4]]
            X[1|2|3|4]
            sage: x = Set([Set([1,3]), Set([2])]); x
            {{1, 3}, {2}}
            sage: X[x]
            X[13|2]

        ::

            sage: X = SupercharacterHopfAlgebra(3).X()
            sage: X[0,[]]
            X0[]
            sage: X[4,[]]
            X4[]
            sage: X[4,[(1,2,1)]]
            X4[(1, 2, 1)]
            sage: X[[1,3],[2]]
            X3[(1, 3, 1)]
            sage: x = Set([Set([1,3]), Set([2])]); x
            {{1, 3}, {2}}
            sage: X[x]
            X3[(1, 3, 1)]
        """
        if isinstance(key, LabelledSetPartition):
            phi = key
        elif isinstance(key, list):
            phi = label_set_partition([key])
        elif isinstance(key, Set_object_enumerated):
            phi = label_set_partition(key)
        elif isinstance(key[0], (int,Integer)):
            phi = LabelledSetPartition(*key)
        else:
            phi = label_set_partition(key)

        return self.monomial(phi)

    ##### TODO: why can't these be put in ParentMethods of Bases?
    def _repr_term(self, phi):
        r"""
        The string representation of a basis element.

        EXAMPLES:

        For ``q==2``, we use a compatified notation of the set partition::

            sage: X = SupercharacterHopfAlgebra(2).X()
            sage: X[0,[]]
            X[]
            sage: X[4,[]]
            X[1|2|3|4]
            sage: X[5,[(1,3),(3,5),(2,4)]]
            X[135|24]

        For ``q>2``, the representation includes the size of the underlying set
        and lists the labelled arcs::

            sage: X = SupercharacterHopfAlgebra(3).X()
            sage: X[0,[]]
            X0[]
            sage: X[4,[]]
            X4[]
            sage: X[4,[(1,2,1)]]
            X4[(1, 2, 1)]
            sage: X[5,[(1,3,1),(3,5,2),(2,4,1)]]
            X5[(1, 3, 1), (2, 4, 1), (3, 5, 2)]
        """
        if self.q() == 2:
            return "%s[%s]" % (self._prefix, "|".join("".join(map(str,block)) for block in sorted(phi.to_set_partition(),key=min)))
        else:
            return "%s%s%s" % (self._prefix, phi.size(), phi.arcs())

    def _an_element_(self):
        r"""
        An element of the algebra.

        EXAMPLES:

            sage: X = SupercharacterHopfAlgebra(2).X()
            sage: X.an_element()
            5*X[1|2] - 3*X[13|2]

        ::

            sage: P = SupercharacterHopfAlgebra(2).P()
            sage: P.an_element()
            5*P[1|2] - 3*P[13|2]

        ::

            sage: K = SupercharacterHopfAlgebra(4).K()
            sage: K.an_element()
            5*K2[] + (-3)*K3[(1, 3, a)]

        """
        R = self.base_ring()
        a = self._field.gen()
        phi = LabelledSetPartition(3, [(1,3,a)])
        psi = LabelledSetPartition(2, [])
        return self.sum_of_terms([(phi,R(-3)),(psi,R(5))])

    def basis(self, d=None):
        r"""
        Return the basis of the homogeneous component of degree ``d``.
        If ``d`` is None, then a basis of the algebra is returned.

        EXAMPLES::

            sage: X = SupercharacterHopfAlgebra(QQ, 2).supercharacter_basis()
            sage: X.basis()
            Lazy family (Term map from Set partitions with arcs labelled by elements of Finite Field of size 2 to Hopf algebra of supercharacters at q=2 over Rational Field on the supercharacter basis(i))_{i in Set partitions with arcs labelled by elements of Finite Field of size 2}
            sage: X.basis(3).list()
            [X[123], X[1|23], X[13|2], X[12|3], X[1|2|3]]

        ::

            sage: P = SupercharacterHopfAlgebra(QQ, 2).powersum_basis()
            sage: P.basis(3).list()
            [P[123], P[1|23], P[13|2], P[12|3], P[1|2|3]]

        ::

            sage: H = SupercharacterHopfAlgebra(QQ, 2).homogeneous_basis()
            sage: H.basis(3).list()
            [H[123], H[1|23], H[13|2], H[12|3], H[1|2|3]]

        ::

            sage: K = SupercharacterHopfAlgebra(QQ, 2).superclass_basis()
            sage: K.basis(3).list()
            [K[123], K[1|23], K[13|2], K[12|3], K[1|2|3]]

        ::

            sage: E = SupercharacterHopfAlgebra(QQ, 2).elementary_basis()
            sage: E.basis(3).list()
            [E[123], E[1|23], E[13|2], E[12|3], E[1|2|3]]

        """
        if d is None:
            # BUG : the following replaces my basis with the one from CombinatorialFreeModule
            #return super(SupercharacterHopfAlgebraBasis,self).basis()
            return Family(self._basis_keys, self.monomial)
        else:
            return Family(LabelledSetPartitions(d,self.q()), self.monomial)

##### various basis
class SupercharacterBasis(SupercharacterHopfAlgebraBasis):
    r"""
    The Hopf algebra of supercharacters on the supercharacter basis.

    EXAMPLES::

        sage: scha = SupercharacterHopfAlgebra(2)
        sage: scha
        Hopf algebra of supercharacters at q=2 over Rational Field
        sage: X = scha.X() # or X = scha.supercharacter_basis()
        sage: X
        Hopf algebra of supercharacters at q=2 over Rational Field on the supercharacter basis

    Elements of the algebra look like::

        sage: X.an_element()
        5*X[1|2] - 3*X[13|2]
        sage: X.some_elements()
        [X[], X[1], X[12], X[1|2], 5*X[1|2] - 3*X[13|2]]

    The algebra is indexed by labelled set partitions, so to create a basis
    element of the algebra, we must first create a labelled set partition (see
    below for a shorthand). For example, the set partition `\{ \{1,3,5\},
    \{2,4\} \}' with arcs `(1,3)`, `(3,5)` and `(2,4)` labelled by `1`, `5`
    and `4`, respectively is created using the command::

        sage: phi = LabelledSetPartition(5, [(1,3,1), (3,5,3), (2,4,4)])
        sage: phi
        [5, [(1, 3, 1), (2, 4, 4), (3, 5, 3)]]

    These objects are used to access the corresponding basis element of the
    algebra::

        sage: phi = LabelledSetPartition(5, [(1,2,1),(2,5,1)])
        sage: X[phi]
        X[125|3|4]

    For simplicity, we can pass to ``X`` the arguments needed to construct a
    labelled set partition (if no label is given, then it is taken to be 1)::

        sage: X[5, [(1,2),(2,5)]]
        X[125|3|4]
        sage: X[5,[]]
        X[1|2|3|4|5]

    Here is an example ``q==3``:: 

        sage: scha = SupercharacterHopfAlgebra(3)
        sage: scha
        Hopf algebra of supercharacters at q=3 over Cyclotomic Field of order 3 and degree 2
        sage: X = scha.X() # or X = scha.supercharacter_basis()
        sage: X
        Hopf algebra of supercharacters at q=3 over Cyclotomic Field of order 3 and degree 2 on the supercharacter basis

    A typical basis element looks like this::

        sage: X[5, [(1,3,1), (3,5,3), (2,4,4)]]
        X5[(1, 3, 1), (2, 4, 4), (3, 5, 3)]

    Here, the number following the ``X`` denotes the size of the set partition
    and the list of elements are the arcs.

    TESTS::

        sage: X = SupercharacterHopfAlgebra(q=2).X()
        sage: TestSuite(X).run(skip=["_test_len"])

    """
    _basis_name = "supercharacter"
    _prefix = "X"

    ##### algebra product
    @cached_method
    def product_on_basis(self, phi, psi):
        r"""
        The product of the supercharacters indexed by ``phi`` and ``psi``.

        This is defined as the shifted concatenation of the labelled set
        partitions.

        EXAMPLES::

            sage: X = SupercharacterHopfAlgebra(2).supercharacter_basis()
            sage: X[3,[]] * X[4,[]]
            X[1|2|3|4|5|6|7]
            sage: X[3,[(1,3)]] * X[4,[(1,3),(2,4)]]
            X[13|2|46|57]

        ::

            sage: X = SupercharacterHopfAlgebra(3).supercharacter_basis()
            sage: X[3,[]] * X[4,[]]
            X7[]
            sage: X[3,[(1,3)]] * X[4,[(1,3),(2,4)]]
            X7[(1, 3, 1), (4, 6, 1), (5, 7, 1)]

        """
        shift = phi.size()
        lsp = phi.arcs() + [(i+shift,j+shift,l) for (i,j,l) in psi.arcs()]
        return self.basis()[LabelledSetPartition(phi.size()+psi.size(), lsp)]

    ##### inner tensor product of characters
    @lazy_attribute
    def inner_tensor_product(self):
        r"""
        The inner tensor product of supercharacters.

        See [1] for details.

        EXAMPLES::

            sage: X = SupercharacterHopfAlgebra(5).X()
            sage: x = X[LabelledSetPartition(5, [(1,4,1)])]
            sage: y = X[LabelledSetPartition(5, [(2,3,2)])]
            sage: X.inner_tensor_product(x,y)
            X5[(1, 4, 1), (2, 3, 2)]
            sage: z = X[LabelledSetPartition(5, [(1,4,2)])]
            sage: X.inner_tensor_product(x,z)
            9*X5[(1, 4, 3)] + 4*X5[(1, 4, 3), (2, 3, 1)] + 4*X5[(1, 4, 3), (2, 3, 2)] + 4*X5[(1, 4, 3), (2, 3, 3)] + 4*X5[(1, 4, 3), (2, 3, 4)]

        REFERENCES:

        -   [1] Nathaniel Thiem, Branching rules in the ring of superclass
            functions of unipotent upper-triangular matrices, J Algebr Comb
            (2010) 31: 267--298
        """
        return self._module_morphism(self._module_morphism(self.inner_tensor_on_basis, position = 0, codomain=self), position = 1)

    @cached_method
    def inner_tensor_on_basis(self, phi, psi):
        r"""
        The inner tensor product of the supercharacters indexed by ``phi`` and
        ``psi``.

        EXAMPLES:

        These exxamples were taken from [1]::

            sage: q = 5
            sage: F = GF(q)
            sage: a,b,c,d = map(F, range(1,5)) 
            sage: X = SupercharacterHopfAlgebra(q).X()
            sage: mu = LabelledSetPartition(6, [(1,6,a)])
            sage: nu = LabelledSetPartition(6, [(2,4,b)])
            sage: X.inner_tensor_on_basis(mu,nu)
            X6[(1, 6, 1), (2, 4, 2)]

        ::

            sage: mu = LabelledSetPartition(6, [(1,4,c)])
            sage: nu = LabelledSetPartition(6, [(2,5,d)])
            sage: X.inner_tensor_on_basis(mu,nu)
            X6[(1, 4, 3), (2, 5, 4)]

        ::

            sage: a = LabelledSetPartition(6, [(1,6,F(1))])
            sage: b = LabelledSetPartition(6, [(2,5,F(2))])
            sage: c = LabelledSetPartition(6, [(1,4,F(3))])
            sage: d = LabelledSetPartition(6, [(2,5,F(4))])
            sage: X.inner_tensor_on_basis(b,d)
            9*X6[(2, 5, 1)] + 4*X6[(2, 5, 1), (3, 4, 1)] + 4*X6[(2, 5, 1), (3, 4, 2)] + 4*X6[(2, 5, 1), (3, 4, 3)] + 4*X6[(2, 5, 1), (3, 4, 4)]
            sage: X.inner_tensor_on_basis(a,c)
            X6[(1, 6, 1)] + X6[(1, 6, 1), (2, 4, 1)] + X6[(1, 6, 1), (2, 4, 2)] + X6[(1, 6, 1), (2, 4, 3)] + X6[(1, 6, 1), (2, 4, 4)] + X6[(1, 6, 1), (3, 4, 1)] + X6[(1, 6, 1), (3, 4, 2)] + X6[(1, 6, 1), (3, 4, 3)] + X6[(1, 6, 1), (3, 4, 4)]

        ::

            sage: mu = LabelledSetPartition(6, [(1,6,F(1)),(2,5,F(2))])
            sage: nu = LabelledSetPartition(6, [(1,4,F(3)),(2,5,F(4))])
            sage: X.inner_tensor_on_basis(mu,nu)
            125*X6[(1, 6, 1), (2, 5, 1)] + 125*X6[(1, 6, 1), (2, 5, 1), (3, 4, 1)] + 125*X6[(1, 6, 1), (2, 5, 1), (3, 4, 2)] + 125*X6[(1, 6, 1), (2, 5, 1), (3, 4, 3)] + 125*X6[(1, 6, 1), (2, 5, 1), (3, 4, 4)]

        REFERENCES:

        -   [1] Nathaniel Thiem, Branching rules in the ring of superclass
            functions of unipotent upper-triangular matrices, J Algebr Comb
            (2010) 31: 267--298
        """
        if phi.size() != psi.size():
            return self.zero()
        if phi.arcs() == []:
            return self.basis()[psi]
        arcs_psi = psi.arcs()
        if arcs_psi == []:
            return self.basis()[phi]
        return self.inner_tensor_product(self._inner_tensor_on_basis_by_arc(phi,arcs_psi.pop()),
                                     self.basis()[LabelledSetPartition(psi.size(), arcs_psi)])

    def _inner_tensor_on_basis_by_arc(self, phi, (i,j,a)):
        r"""
        The inner tensor product of the supercharacter indexed by ``phi`` by
        the arc from ``i`` to ``j`` with label ``a``.

        EXAMPLES::

            sage: q = 5
            sage: F = GF(q)
            sage: a,b,c,d = map(F, range(1,5)) 
            sage: X = SupercharacterHopfAlgebra(q).X()
            sage: mu = LabelledSetPartition(6, [(1,6,a)])
            sage: X._inner_tensor_on_basis_by_arc(mu,(2,4,b))
            X6[(1, 6, 1), (2, 4, 2)]
            sage: mu = LabelledSetPartition(6, [(1,4,c)])
            sage: X._inner_tensor_on_basis_by_arc(mu,(2,5,d))
            X6[(1, 4, 3), (2, 5, 4)]
        """
        arcs_phi = phi.arcs()
        for (k,arc) in enumerate(arcs_phi):
            if arc[0]==i or arc[1]==j:
                del arcs_phi[k]
                return self.inner_tensor_product(
                        self.basis()[LabelledSetPartition(phi.size(),arcs_phi)],
                        self._inner_tensor_of_two_arcs(arc,(i,j,a),phi.size())
                        )
        else:
            return self.basis()[LabelledSetPartition(phi.size(), arcs_phi+[(i,j,a)])]

    def _inner_tensor_of_two_arcs(self, arc1, arc2, n):
        r"""
        The inner tensor product of the supercharacters indexed by two labelled
        set partitions of ``[n]`` containing the labelled arcs ``arc1`` and
        ``arc2``, respectively.

        This method implements Corollary 4.9 of [1].

        EXAMPLES::

            sage: X = SupercharacterHopfAlgebra(q=2).X()
            sage: X._inner_tensor_of_two_arcs( (1,2,1), (3,4,1), 5 )
            X[12|34|5]
            sage: X._inner_tensor_of_two_arcs( (1,2,1), (2,4,1), 5 )
            X[124|3|5]
            sage: X._inner_tensor_of_two_arcs( (1,2,1), (1,4,1), 5 )
            X[14|2|3|5]
            sage: X._inner_tensor_of_two_arcs( (1,3,1), (1,4,1), 5 )
            X[14|2|3|5] + X[14|23|5]
            sage: X._inner_tensor_of_two_arcs( (1,4,1), (2,4,1), 5 )
            X[14|2|3|5] + X[14|23|5]
            sage: X._inner_tensor_of_two_arcs( (1,4,1), (1,4,1), 5 )
            3*X[14|2|3|5] + X[14|23|5]

        ::

            sage: X = SupercharacterHopfAlgebra(q=3).X()
            sage: X._inner_tensor_of_two_arcs( (1,2,1), (3,4,2), 5 )
            X5[(1, 2, 1), (3, 4, 2)]
            sage: X._inner_tensor_of_two_arcs( (1,2,1), (2,4,2), 5 )
            X5[(1, 2, 1), (2, 4, 2)]
            sage: X._inner_tensor_of_two_arcs( (1,2,1), (1,4,2), 5 )
            X5[(1, 4, 2)]
            sage: X._inner_tensor_of_two_arcs( (1,3,1), (1,4,2), 5 )
            X5[(1, 4, 2)] + X5[(1, 4, 2), (2, 3, 1)] + X5[(1, 4, 2), (2, 3, 2)]
            sage: X._inner_tensor_of_two_arcs( (1,4,1), (2,4,2), 5 )
            X5[(1, 4, 1)] + X5[(1, 4, 1), (2, 3, 1)] + X5[(1, 4, 1), (2, 3, 2)]
            sage: X._inner_tensor_of_two_arcs( (1,4,1), (1,4,2), 5 )
            5*X5[(1, 4, 3)] + 2*X5[(1, 4, 3), (2, 3, 1)] + 2*X5[(1, 4, 3), (2, 3, 2)]

        TESTS::

            sage: X = SupercharacterHopfAlgebra(q=5).X()
            sage: F = GF(5)
            sage: X._inner_tensor_of_two_arcs( (2,5,F(2)), (2,5,F(4)), 6 )
            9*X6[(2, 5, 1)] + 4*X6[(2, 5, 1), (3, 4, 1)] + 4*X6[(2, 5, 1), (3, 4, 2)] + 4*X6[(2, 5, 1), (3, 4, 3)] + 4*X6[(2, 5, 1), (3, 4, 4)]

        REFERENCES:

        -   [1] Nathaniel Thiem, Branching rules in the ring of superclass
            functions of unipotent upper-triangular matrices, J Algebr Comb
            (2010) 31: 267--298
        """
        (i,k,a), (j,l,b) = sorted([arc1,arc2])
        assert(i<k and j<l)
        q = self.q()
        Bmap = lambda arcs: self[LabelledSetPartition(n,arcs)]
        if (i,k) != (j,l):
            if set([i,k]).intersection(set([j,l])) == set([]):
                return Bmap([(i,k,a),(j,l,b)])
            elif i < j == k < l:
                return Bmap([(i,j,a),(j,l,b)])
            elif i == j < k < l:
                return Bmap([(i,l,b)]) + self.sum(Bmap([(jj,k,c),(i,l,b)])
                        for jj in range(i+1,k)
                        for c in self._field if c != 0)
            elif i < j < k == l:
                return Bmap([(i,l,a)]) + self.sum(Bmap([(i,l,a),(j,kk,c)])
                        for kk in range(j+1,l)
                        for c in self._field if c != 0)
        else:
            if b == -a:
                return Bmap([]) \
                    + self.sum(Bmap([(i,jj,c)])
                            for jj in range(i+1,l)
                            for c in self._field if c != 0) \
                    + self.sum(Bmap([(kk,l,c)]) 
                            for kk in range(i+1,l)
                            for c in self._field if c != 0) \
                    + self.sum(Bmap([(i,jj,c),(kk,l,d)])
                            for jj in range(i+1,l)
                            for kk in range(i+1,l)
                            for c in self._field if c != 0
                            for d in self._field if d != 0)
            else:
                return ((l-i-1)*(q-1)+1) * Bmap([(i,l,a+b)]) \
                        + (q-1) * self.sum(Bmap([(jj,kk,c),(i,l,a+b)])
                                for jj in range(i+1,l-1)
                                for kk in range(jj+1,l)
                                for c in self._field if c != 0)

    ##### restriction
    def _arc_to_partition(self, arc, A):
        r"""
        The standardization of the set partition of the set ``A`` defined by
        the ``arc``.

        EXAMPLES::

            sage: X = SupercharacterHopfAlgebra(q=2).X()
            sage: X._arc_to_partition( (3,7,2), [3,5,7])
            [3, [(1, 3, 2)]]

        """
        if arc == ():
            return LabelledSetPartition(len(A), ())
        # standardize
        std = dict((j,i+1) for (i,j) in enumerate(sorted(A)))
        return LabelledSetPartition(len(A),[(std[arc[0]], std[arc[1]], arc[2])])

    def _restriction_on_arc(self, arc, A):
        r"""
        Compute the standardized restriction of the supercharacter indexed by
        arc to `X_A` for a subset `A` of `[n]`.

        This method is based on Theorem 4.6 of [1].

        EXAMPLES:

        These examples appear on pages 285--286 of [1]::

            sage: F = GF(4,'a')
            sage: a = F.gen()
            sage: X = SupercharacterHopfAlgebra(q=4).X()
            sage: A = Set([2,3,4,5])
            sage: X._restriction_on_arc((2,5,a),A)
            X4[(1, 4, a)]
            sage: X._restriction_on_arc([1,5,F(1)],A)
            X4[] + X4[(1, 4, 1)] + X4[(1, 4, a)] + X4[(1, 4, a + 1)] + X4[(2, 4, 1)] + X4[(2, 4, a)] + X4[(2, 4, a + 1)] + X4[(3, 4, 1)] + X4[(3, 4, a)] + X4[(3, 4, a + 1)]
            sage: X._restriction_on_arc([1,7,F(1)],A)
            52*X4[] + 12*X4[(1, 2, 1)] + 12*X4[(1, 2, a)] + 12*X4[(1, 2, a + 1)] + 12*X4[(1, 3, 1)] + 12*X4[(1, 3, a)] + 12*X4[(1, 3, a + 1)] + 12*X4[(1, 4, 1)] + 12*X4[(1, 4, a)] + 12*X4[(1, 4, a + 1)] + 12*X4[(2, 3, 1)] + 12*X4[(2, 3, a)] + 12*X4[(2, 3, a + 1)] + 12*X4[(2, 4, 1)] + 12*X4[(2, 4, a)] + 12*X4[(2, 4, a + 1)] + 12*X4[(3, 4, 1)] + 12*X4[(3, 4, a)] + 12*X4[(3, 4, a + 1)]
            sage: X._restriction_on_arc([5,7,F(1)],A)
            4*X4[]

        REFERENCES:

        -   [1] Nathaniel Thiem, Branching rules in the ring of superclass
            functions of unipotent upper-triangular matrices, J Algebr Comb
            (2010) 31: 267--298
        """
        q = self.q()
        if len(arc) == 0:
            raise NotImplementedError
        (i,l,a) = arc
        Bmap = lambda arc: self[self._arc_to_partition(arc,A)]
        if i in A:
            if l in A:
                new_term = Bmap((i,l,a))
            else:
                new_term = Bmap(()) + \
                    self.sum(Bmap((i,k,b))
                            for k in range(i+1,l) if k in A
                            for b in self._field if b != 0)
        else:
            if l in A:
                new_term = Bmap(()) + \
                    self.sum(Bmap((j,l,b)) for j in range(i+1,l) if j in A
                                    for b in self._field if b != 0)
            else:
                new_term = (sum(1 for j in range(i,l) if j in A)*(q-1) + 1) * Bmap(())
                for j in A:
                    if i < j < l:
                        for k in A:
                            if j < k < l:
                                new_term += (q-1) * self.sum(Bmap((j,k,c)) for c in self._field if c!= 0)
        return q**sum(1 for k in range(i+1,l) if k not in A) * new_term

    @cached_method
    def restriction_on_basis(self, phi, A):
        r"""
        Compute the restriction of the supercharacter indexed by ``phi`` to the
        subset ``A``.


        EXAMPLES::

            sage: X = SupercharacterHopfAlgebra(2).X()
            sage: arc = LabelledSetPartition(3,[(1,3,GF(2)(1))])
            sage: for A in Subsets(3):
            ...       print "Res_%s(%s) = %s" % (A,X[arc],X.restriction_on_basis(arc,A))
            Res_{}(X[13|2]) = 2*X[]
            Res_{1}(X[13|2]) = 2*X[1]
            Res_{2}(X[13|2]) = 2*X[1]
            Res_{3}(X[13|2]) = 2*X[1]
            Res_{1, 2}(X[13|2]) = X[1|2] + X[12]
            Res_{1, 3}(X[13|2]) = 2*X[12]
            Res_{2, 3}(X[13|2]) = X[1|2] + X[12]
            Res_{1, 2, 3}(X[13|2]) = X[13|2]

        These examples appear on pages 285--286 of [1]::

            sage: F = GF(4,'a')
            sage: a = F.gen()
            sage: X = SupercharacterHopfAlgebra(q=4).X()
            sage: A = Set([2,3,4,5])
            sage: X.restriction_on_basis(LabelledSetPartition(5,[(2,5,a)]),A)
            X4[(1, 4, a)]
            sage: X.restriction_on_basis(LabelledSetPartition(5,[(1,5,F(1))]),A)
            X4[] + X4[(1, 4, 1)] + X4[(1, 4, a)] + X4[(1, 4, a + 1)] + X4[(2, 4, 1)] + X4[(2, 4, a)] + X4[(2, 4, a + 1)] + X4[(3, 4, 1)] + X4[(3, 4, a)] + X4[(3, 4, a + 1)]
            sage: X.restriction_on_basis(LabelledSetPartition(5,[(1,7,F(1))]),A)
            52*X4[] + 12*X4[(1, 2, 1)] + 12*X4[(1, 2, a)] + 12*X4[(1, 2, a + 1)] + 12*X4[(1, 3, 1)] + 12*X4[(1, 3, a)] + 12*X4[(1, 3, a + 1)] + 12*X4[(1, 4, 1)] + 12*X4[(1, 4, a)] + 12*X4[(1, 4, a + 1)] + 12*X4[(2, 3, 1)] + 12*X4[(2, 3, a)] + 12*X4[(2, 3, a + 1)] + 12*X4[(2, 4, 1)] + 12*X4[(2, 4, a)] + 12*X4[(2, 4, a + 1)] + 12*X4[(3, 4, 1)] + 12*X4[(3, 4, a)] + 12*X4[(3, 4, a + 1)]
            sage: X.restriction_on_basis(LabelledSetPartition(5,[(5,7,F(1))]),A)
            4*X4[]

        REFERENCES:

        -   [1] Nathaniel Thiem, Branching rules in the ring of superclass
            functions of unipotent upper-triangular matrices, J Algebr Comb
            (2010) 31: 267--298
        """
        if phi.arcs() == []:
            return self.basis()[LabelledSetPartition(len(A),[])]
        res_on_arcs = self._restriction_on_arc
        return reduce(self.inner_tensor_product,[res_on_arcs(arc,A) for arc in phi.arcs()])


    @lazy_attribute
    def restriction(self):
        r"""
        The endomorphism of the Hopf algebra of supercharacters obtained by
        linearly extending the restriction map.

        EXAMPLES::

            sage: X = SupercharacterHopfAlgebra(q=4).X()
            sage: x = X[7, [(2,5,1)]]
            sage: y = X[7, [(5,7,1)]]
            sage: X.restriction(x + 3*y, Set([2,3,4,5]))
            12*X4[] + X4[(1, 4, 1)]
        """
        return self._module_morphism(self.restriction_on_basis, position=0, codomain=self)

    ##### coproduct
    @cached_method
    def coproduct_on_basis(self, phi):
        r"""
        The coproduct of the supercharacter indexed by ``phi``.

        EXAMPLES::

            sage: q = 2
            sage: X = SupercharacterHopfAlgebra(q).X()
            sage: mu = LabelledSetPartition(3, [(1,3)])
            sage: x = X.basis()[mu]
            sage: X.coproduct(x)
            X[] # X[13|2] + 2*X[1] # X[1|2] + 4*X[1] # X[12] + 2*X[1|2] # X[1] + 4*X[12] # X[1] + X[13|2] # X[]

        TESTS::

            sage: X.coproduct(X.one())
            X[] # X[]
            sage: X.coproduct(X.zero())
            0

        Let's test that this behaves correctly on some labelled set partitions
        with no arcs (and hence no labels)::

            sage: arcless = lambda n : LabelledSetPartition(n,[])
            sage: X.coproduct_on_basis(arcless(0))
            X[] # X[]
            sage: X.coproduct_on_basis(arcless(1))
            X[] # X[1] + X[1] # X[]
            sage: X.coproduct_on_basis(arcless(2))
            X[] # X[1|2] + 2*X[1] # X[1] + X[1|2] # X[]
            sage: X.coproduct_on_basis(arcless(3))
            X[] # X[1|2|3] + 3*X[1] # X[1|2] + 3*X[1|2] # X[1] + X[1|2|3] # X[]
        """
        if phi.size() == 0:
            return tensor([self.one(),self.one()])
        superrestriction = self.restriction_on_basis
        res = tensor([self.one(),self.basis()[phi]]) \
                + tensor([self.basis()[phi],self.one()])
        q = self.q()
        for (J,Jcomp) in OrderedSetPartitions(phi.size(), 2):
            res += q**(-sum(l-i-1 for (i,l,a) in phi.arcs())) \
                    * tensor([superrestriction(phi,J), superrestriction(phi,Jcomp)])
        return res

class SuperclassBasis(SupercharacterHopfAlgebraBasis):
    r"""
    The Hopf algebra of supercharacters on the superclass basis.

    EXAMPLES::

        sage: scha = SupercharacterHopfAlgebra(2)
        sage: scha
        Hopf algebra of supercharacters at q=2 over Rational Field
        sage: K = scha.K() # or K = scha.superclass_basis()
        sage: K
        Hopf algebra of supercharacters at q=2 over Rational Field on the superclass basis

    Elements of the algebra look like::

        sage: K.an_element()
        5*K[1|2] - 3*K[13|2]
        sage: K.some_elements()
        [K[], K[1], K[12], K[1|2], 5*K[1|2] - 3*K[13|2]]

    There are change of bases implemented between the superclass and
    supercharacter bases::

        sage: X = scha.supercharacter_basis()
        sage: a = X[3,[]]; a
        X[1|2|3]
        sage: K(a)
        K[1|2|3] + K[12|3] + K[123] + K[13|2] + K[1|23]
        sage: X(K(a))
        X[1|2|3]

    TESTS:

    We test that converting between the supercharacter and superclass basis
    works properly::

        sage: q = 2
        sage: scha = SupercharacterHopfAlgebra(q)
        sage: X = scha.supercharacter_basis()
        sage: K = scha.superclass_basis()
        sage: for n in range(5):
        ...       for phi in X.basis(n).keys():
        ...           assert(K(X(K[phi]))==K[phi])
        ...           assert(X(K(X[phi]))==X[phi])

    ::

        sage: q = 4
        sage: scha = SupercharacterHopfAlgebra(q)
        sage: X = scha.supercharacter_basis()
        sage: K = scha.superclass_basis()
        sage: for n in range(4):
        ...       for phi in X.basis(n).keys():
        ...           assert(K(X(K[phi]))==K[phi])
        ...           assert(X(K(X[phi]))==X[phi])

    We run the test suites::

        sage: K = SupercharacterHopfAlgebra(q=2).K()
        sage: TestSuite(K).run(skip=["_test_len"]) # long time
        sage: K = SupercharacterHopfAlgebra(q=3).K()
        sage: TestSuite(K).run(skip=["_test_len","_test_prod","_test_coproduct","_test_associativity"]) # long time
    """
    _basis_name = "superclass"
    _prefix = "K"

    @cached_method
    def _supercharacter_to_superclass_on_basis(self, phi):
        r"""
        Expand the supercharacter indexed by ``phi`` in the superclass basis.

        EXAMPLES::

            sage: K = SupercharacterHopfAlgebra(2).superclass_basis()
            sage: phi = LabelledSetPartition(3,[(1,3,1)])
            sage: a = K._supercharacter_to_superclass_on_basis(phi); a
            2*K[1|2|3] - 2*K[13|2]

        TESTS::

            sage: X = SupercharacterHopfAlgebra(2).supercharacter_basis()
            sage: X(a)
            X[13|2]
        """
        scf = SupercharacterTable(self).supercharacter_formula
        return self.sum_of_terms((psi,scf(phi,psi)) for psi in self.basis(phi.size()).keys())

    @cached_method
    def _superclass_to_supercharcter_on_basis(self, phi):
        r"""
        Expand the superclass indexed by ``phi`` in the supercharacter basis.

        EXAMPLES::

            sage: K = SupercharacterHopfAlgebra(2).superclass_basis()
            sage: phi = LabelledSetPartition(3,[])
            sage: psi = LabelledSetPartition(3,[(1,3,1)])
            sage: K._superclass_to_supercharcter_on_basis(phi)
            1/8*X[1|2|3] + 1/8*X[12|3] + 1/8*X[123] + 1/4*X[13|2] + 1/8*X[1|23]
            sage: K._superclass_to_supercharcter_on_basis(psi)
            1/8*X[1|2|3] + 1/8*X[12|3] + 1/8*X[123] - 1/4*X[13|2] + 1/8*X[1|23]
        """
        n, q = phi.size(), self.q()
        X = self.realization_of().supercharacter_basis()
        trans_mat = SupercharacterTable(self).table_inverse(n)
        vec = trans_mat.row(self.basis(n).keys().rank(phi))
        return X.sum_of_terms(zip(self.basis(n).keys(),vec))

    @cached_method
    def coproduct_on_basis(self, phi):
        r"""
        The coproduct of the superclass indexed by ``phi``.

        .. note::

            This only works for `q=2`.

        EXAMPLES::

            sage: K = SupercharacterHopfAlgebra(2).superclass_basis()
            sage: K.coproduct_on_basis(LabelledSetPartition(3,[]))
            K[] # K[1|2|3] + 3*K[1] # K[1|2] + 3*K[1|2] # K[1] + K[1|2|3] # K[]

        TESTS:

        Test that coproduct agrees with the coproduct on supercharacters::

            sage: K = SupercharacterHopfAlgebra(2).superclass_basis()
            sage: X = SupercharacterHopfAlgebra(2).supercharacter_basis()
            sage: KK = K.tensor_square()
            sage: XX = X.tensor_square()
            sage: XX.module_morphism(lambda (phi,psi): tensor([K(X[phi]),K(X[psi])]), codomain=KK).register_as_coercion()
            sage: for n in range(5):
            ...       for a in K.basis(n):
            ...           assert(a.coproduct() == KK(X(a).coproduct()))
        """
        if self.q() == 2:
            A = list(phi.to_set_partition())
            T = self.tensor_square()
            if phi.size() == 0:
                return tensor([self.one(),self.one()])
            res = tensor([self.one(),self[phi]]) + tensor([self[phi],self.one()])
            for (left, right) in OrderedSetPartitions(A,2):
                left = label_set_partition(standardize_set_partition(left))
                right = label_set_partition(standardize_set_partition(right))
                res = res + T[left, right]
            return res
        else:
            return super(SuperclassBasis, self).coproduct_on_basis(phi)

    @lazy_attribute
    def omega3(self):
        r"""
        The endomorphism of the Hopf algebra of supercharacters obtained by
        linearly extending the `\omega_3` map.

        EXAMPLES::

            sage: K = SupercharacterHopfAlgebra(q=2).K()
            sage: x = K[4, [(3,4,1)]]
            sage: y = K[4, [(1,3,1),(2,4,1)]]
            sage: K.omega3(x + 3*y)
            3*K[13|24] - K[1|2|34]
        """
        return self._module_morphism(self.omega3_on_basis, codomain=self)

    @cached_method
    def omega3_on_basis(self, phi):
        r"""
        The endomorphism of ``self`` defined on the superclass basis
        by `\omega_3(K_\phi) = (-1)^{|\phi|-\ell(phi)} K_\phi`.

        EXAMPLES::

            sage: K = SupercharacterHopfAlgebra(q=2).K()
            sage: K.omega3_on_basis(LabelledSetPartition(4, [(3,4,1)]))
            -K[1|2|34]
        """
        return self((-1)**(phi.size()-len(phi.to_set_partition())) * self.basis()[phi])

class PowersumBasis(SupercharacterHopfAlgebraBasis):
    r"""
    The Hopf algebra of supercharacters on the powersum basis.

    The definition of this basis is based on Theorem 3.1(i) of [NCSym].

    EXAMPLES::

        sage: scha = SupercharacterHopfAlgebra(2)
        sage: scha
        Hopf algebra of supercharacters at q=2 over Rational Field
        sage: P = scha.P() # or P = scha.powersum_basis()
        sage: P
        Hopf algebra of supercharacters at q=2 over Rational Field on the powersum basis

    Elements of the algebra look like::

        sage: P.an_element()
        5*P[1|2] - 3*P[13|2]
        sage: P.some_elements()
        [P[], P[1], P[12], P[1|2], 5*P[1|2] - 3*P[13|2]]

    There are change of bases routines implemented between the powersum and
    superclass bases::

        sage: a = P[4,[(1,3),(2,4)]]; a
        P[13|24]
        sage: K = scha.superclass_basis()
        sage: K(a)
        K[1234] + K[13|24]
        sage: P(K(a))
        P[13|24]

    As a result, we can convert between the powersum basis and the
    supercharacter basis (an implicit conversion to the superclass basis is
    performed in the background).

    .. note::
    
        The powersum basis is only defined for `q=2`.

    TESTS:

    We test that converting between the supercharacter and powersum basis works
    properly::

        sage: q = 2
        sage: scha = SupercharacterHopfAlgebra(q)
        sage: X = scha.supercharacter_basis()
        sage: P = scha.powersum_basis()
        sage: for n in range(5):
        ...       for phi in X.basis(n).keys():
        ...           assert(P(X(P[phi]))==P[phi])
        ...           assert(X(P(X[phi]))==X[phi])

    We run the test suites::

        sage: P = SupercharacterHopfAlgebra(q=2).P()
        sage: TestSuite(P).run(skip=["_test_len"]) # long time

    REFERENCES:

    -   [NCSym] Mercedes H. Rosas and Bruce E. Sagan, Symmetric functions in
        noncommuting variables, Transactions of the American Mathematical
        Society, Volume 358, Number 1, Pages 215--232
    """
    _basis_name = "powersum"
    _prefix = "P"

    @cached_method
    def _powersum_to_superclass_on_basis(self, phi):
        r"""
        Expand the powersum indexed by ``phi`` in the superclass basis.

        The definition of this basis is based on Theorem 3.1(i) of [NCSym].

        EXAMPLES::

            sage: P = SupercharacterHopfAlgebra(2).powersum_basis()
            sage: phi = LabelledSetPartition(3,[])
            sage: psi = LabelledSetPartition(3,[(1,3,1)])
            sage: P._powersum_to_superclass_on_basis(phi)
            K[1|2|3] + K[12|3] + K[123] + K[13|2] + K[1|23]
            sage: P._powersum_to_superclass_on_basis(psi)
            K[123] + K[13|2]

        REFERENCES:

        -   [NCSym] Mercedes H. Rosas and Bruce E. Sagan, Symmetric functions in
            noncommuting variables, Transactions of the American Mathematical
            Society, Volume 358, Number 1, Pages 215--232
        """
        K = self.realization_of().superclass_basis()
        L = LatticeOfSetPartitions(phi.size())
        return K.sum_of_monomials([label_set_partition(x.element) for x in L.order_filter([L(phi.to_set_partition())])])

    @cached_method
    def product_on_basis(self, phi, psi):
        r"""
        The product of the powersums indexed by ``phi`` and ``psi``.

        This is defined as the shifted concatenation of the set partitions.

        EXAMPLES::

            sage: P = SupercharacterHopfAlgebra(2).powersum_basis()
            sage: P[3,[]] * P[4,[]]
            P[1|2|3|4|5|6|7]
            sage: P[3,[(1,3)]] * P[4,[(1,3),(2,4)]]
            P[13|2|46|57]

        Test that the product agrees with the product on supercharacters::

            sage: P = SupercharacterHopfAlgebra(2).powersum_basis()
            sage: X = SupercharacterHopfAlgebra(2).supercharacter_basis()
            sage: for a in P.basis(3):
            ...       for b in P.basis(2):
            ...           assert(X(a*b) == X(a)*X(b))
        """
        shift = phi.size()
        lsp = phi.arcs() + [(i+shift,j+shift,l) for (i,j,l) in psi.arcs()]
        return self.basis()[LabelledSetPartition(phi.size()+psi.size(), lsp)]


    @lazy_attribute
    def omega1(self):
        r"""
        The endomorphism of the Hopf algebra of supercharacters obtained by
        linearly extending the `\omega_1` map.

        .. note::

            The map :method:`omega` is another implementation of this
            same endomorphism.

        EXAMPLES::

            sage: P = SupercharacterHopfAlgebra(q=2).P()
            sage: a = 3*P[4,[(1,3),(2,4)]] + P[4,[(3,4)]]
            sage: a
            3*P[13|24] + P[1|2|34]
            sage: P.omega1(a)
            3*P[13|24] - P[1|2|34]
        """
        return self._module_morphism(self.omega1_on_basis, codomain=self)

    @cached_method
    def omega1_on_basis(self, phi):
        r"""
        The endomorphism of ``self`` defined on the powersum basis
        by `\omega_1(P_\phi) = (-1)^{|\phi|-\ell(phi)} P_\phi`.

        .. note::

            The map :method:`omega` is another implementation of this
            same endomorphism.

        EXAMPLES::

            sage: P = SupercharacterHopfAlgebra(q=2).P()
            sage: P.omega1_on_basis(LabelledSetPartition(4, [(3,4,1)]))
            -P[1|2|34]
        """
        return self((-1)**(phi.size()-len(phi.to_set_partition())) * self.basis()[phi])

class HomogeneousBasis(SupercharacterHopfAlgebraBasis):
    r"""
    The Hopf algebra of supercharacters on the homogeneous basis.

    EXAMPLES::

        sage: scha = SupercharacterHopfAlgebra(2)
        sage: scha
        Hopf algebra of supercharacters at q=2 over Rational Field
        sage: H = scha.H() # or H = scha.homogeneous_basis()
        sage: H
        Hopf algebra of supercharacters at q=2 over Rational Field on the homogeneous basis

    Elements of the algebra look like::

        sage: H.an_element()
        5*H[1|2] - 3*H[13|2]
        sage: H.some_elements()
        [H[], H[1], H[12], H[1|2], 5*H[1|2] - 3*H[13|2]]

    There are routines for changing bases from the homogeneous basis to the
    superclass basis and from the powersum basis to the homogeneous basis::
    
        sage: a = H[4, [(1,3),(2,4)]]; a
        H[13|24]
        sage: K = scha.superclass_basis()
        sage: K(a)
        K[1|2|3|4] + K[12|3|4] + 2*K[123|4] + 4*K[1234] + 2*K[124|3] + K[12|34] + 2*K[13|2|4] + 4*K[13|24] + 2*K[134|2] + K[14|2|3] + K[14|23] + K[1|23|4] + 2*K[1|234] + 2*K[1|24|3] + K[1|2|34]

    As a result, we can convert between any of the homogeneous basis and the
    supercharacter basis (implicit conversions to intermediate bases are
    performed in the background).

    Note: the homogeneous basis is only defined for `q=2`.

    TESTS:

    We test that converting between the supercharacter and homogeneous basis
    works properly::

        sage: q = 2
        sage: scha = SupercharacterHopfAlgebra(q)
        sage: X = scha.supercharacter_basis()
        sage: H = scha.homogeneous_basis()
        sage: for n in range(4):
        ...       for phi in X.basis(n).keys():
        ...           assert(H(X(H[phi]))==H[phi])
        ...           assert(X(H(X[phi]))==X[phi])

    We run the test suites::

        sage: H = SupercharacterHopfAlgebra(q=2).H()
        sage: TestSuite(H).run(skip=["_test_len","_test_prod","_test_coproduct","_test_associativity"])
    """
    _basis_name = "homogeneous"
    _prefix = "H"

    @cached_method
    def _homogeneous_to_superclass_on_basis(self, phi):
        r"""
        Expand the homogeneous basis element indexed by ``phi`` in the
        superclass basis.

        The definition of this basis is based on Theorem 3.1(iii) of [NCSym].

        EXAMPLES::

            sage: H = SupercharacterHopfAlgebra(2).homogeneous_basis()
            sage: phi = LabelledSetPartition(3,[])
            sage: psi = LabelledSetPartition(3,[(1,3,1)])
            sage: H._homogeneous_to_superclass_on_basis(phi)
            K[1|2|3] + K[12|3] + K[123] + K[13|2] + K[1|23]
            sage: H._homogeneous_to_superclass_on_basis(psi)
            K[1|2|3] + K[12|3] + 2*K[123] + 2*K[13|2] + K[1|23]

        REFERENCES:

        -   [NCSym] Mercedes H. Rosas and Bruce E. Sagan, Symmetric functions in
            noncommuting variables, Transactions of the American Mathematical
            Society, Volume 358, Number 1, Pages 215--232
        """
        K = self.realization_of().superclass_basis()
        return K.sum(prod(factorial(len(a.intersection(b)))
                    for a in phi.to_set_partition()
                    for b in psi.to_set_partition()) * K[psi]
                    for psi in K.basis(phi.size()).keys())

    @cached_method
    def _powersum_to_homogeneous_on_basis(self, phi):
        r"""
        Expand the powersum indexed by ``phi`` in the homogeneous basis.

        This change of basis formula is Theorem 3.4 of [NCSym].

        EXAMPLES::

            sage: H = SupercharacterHopfAlgebra(2).homogeneous_basis()
            sage: phi = LabelledSetPartition(3,[])
            sage: psi = LabelledSetPartition(3,[(1,3,1)])
            sage: H._powersum_to_homogeneous_on_basis(phi)
            H[1|2|3]
            sage: H._powersum_to_homogeneous_on_basis(psi)
            -H[1|2|3] + H[13|2]

        REFERENCES:

        -   [NCSym] Mercedes H. Rosas and Bruce E. Sagan, Symmetric functions in
            noncommuting variables, Transactions of the American Mathematical
            Society, Volume 358, Number 1, Pages 215--232
        """
        H = self.realization_of().homogeneous_basis()
        L = LatticeOfSetPartitions(phi.size())
        phi = L(phi.to_set_partition())
        mu = L.mobius_function
        return 1/abs(mu(L.bottom(),phi)) * H.sum(mu(sigma,phi)*H[label_set_partition(sigma.element)]
                for sigma in L.order_ideal([phi]))

class ElementaryBasis(SupercharacterHopfAlgebraBasis):
    r"""
    The Hopf algebra of supercharacters on the elementary basis.

    The definition of this basis is based on Theorem 3.4 of [NCSym].

    EXAMPLES::

        sage: E = SupercharacterHopfAlgebra(2).E()
        sage: E
        Hopf algebra of supercharacters at q=2 over Rational Field on the elementary basis

    Elements of the algebra look like::

        sage: E.an_element()
        5*E[1|2] - 3*E[13|2]
        sage: E.some_elements()
        [E[], E[1], E[12], E[1|2], 5*E[1|2] - 3*E[13|2]]

    There are change of bases routines implemented between the powersum and
    superclass bases::

        sage: a = E[4,[(1,3),(2,4)]]; a
        E[13|24]
        sage: P = SupercharacterHopfAlgebra(2).powersum_basis()
        sage: P(a)
        P[1|2|3|4] - P[13|2|4] + P[13|24] - P[1|24|3]
        sage: E(P(a))
        E[13|24]

    As a result, we can convert between the powersum basis and the
    supercharacter basis (an implicit conversion to the superclass basis is
    performed in the background)::

        sage: a = E[4,[(1,3),(2,4)]]; a
        E[13|24]
        sage: K = SupercharacterHopfAlgebra(2).superclass_basis()
        sage: K(a)
        K[1|2|3|4] + K[12|3|4] + K[12|34] + K[14|2|3] + K[14|23] + K[1|23|4] + K[1|2|34]

    .. note::

        The elementary basis is only defined for `q=2`.

    TESTS:

    We test that converting between the supercharacter and powersum basis works
    properly::

        sage: X = SupercharacterHopfAlgebra(2).supercharacter_basis()
        sage: E = SupercharacterHopfAlgebra(2).elementary_basis()
        sage: for n in range(5):
        ...       for phi in X.basis(n).keys():
        ...           assert(E(X(E[phi]))==E[phi])
        ...           assert(X(E(X[phi]))==X[phi])

    We run the test suites::

        sage: E = SupercharacterHopfAlgebra(q=2).E()
        sage: TestSuite(E).run(skip=["_test_len","_test_prod","_test_coproduct","_test_associativity"])

    REFERENCES:

    -   [NCSym] Mercedes H. Rosas and Bruce E. Sagan, Symmetric functions in
        noncommuting variables, Transactions of the American Mathematical
        Society, Volume 358, Number 1, Pages 215--232
    """
    _basis_name = "elementary"
    _prefix = "E"

    @cached_method
    def _elementary_to_powersum_on_basis(self, phi):
        r"""
        Expand the elementary basis element indexed by ``phi`` in the
        powersum basis.

        The definition of this basis is based on Theorem 3.4 of [NCSym].

        EXAMPLES::

            sage: SupercharacterHopfAlgebra(2).inject_shorthands()
            sage: phi = LabelledSetPartition(3,[])
            sage: psi = LabelledSetPartition(3,[(1,3,1)])
            sage: zeta = LabelledSetPartition(0,[])
            sage: E._elementary_to_powersum_on_basis(phi)
            P[1|2|3]
            sage: E._elementary_to_powersum_on_basis(psi)
            P[1|2|3] - P[13|2]
            sage: E._elementary_to_powersum_on_basis(zeta)
            P[]

        TESTS::

            sage: X.convert_map_from(E)
            sage: X(E.one())

        REFERENCES:

        -   [NCSym] Mercedes H. Rosas and Bruce E. Sagan, Symmetric functions in
            noncommuting variables, Transactions of the American Mathematical
            Society, Volume 358, Number 1, Pages 215--232
        """
        P = self.realization_of().powersum_basis()
        L = LatticeOfSetPartitions(phi.size())
        mu = L.mobius_function
        pi = L(phi.to_set_partition())
        hat0 = L.bottom()
        return P.sum(mu(hat0,sigma)*P[label_set_partition(sigma.element)]
                    for sigma in L.order_ideal([pi]))

    @cached_method
    def _powersum_to_elementary_on_basis(self, phi):
        r"""
        Expand the powersum indexed by ``phi`` in the elementary basis.

        The definition of this basis is based on Theorem 3.4 of [NCSym].

        EXAMPLES::

            sage: E = SupercharacterHopfAlgebra(2).elementary_basis()
            sage: phi = LabelledSetPartition(3,[])
            sage: psi = LabelledSetPartition(3,[(1,3,1)])
            sage: E._powersum_to_elementary_on_basis(phi)
            E[1|2|3]
            sage: E._powersum_to_elementary_on_basis(psi)
            E[1|2|3] - E[13|2]

        REFERENCES:

        -   [NCSym] Mercedes H. Rosas and Bruce E. Sagan, Symmetric functions in
            noncommuting variables, Transactions of the American Mathematical
            Society, Volume 358, Number 1, Pages 215--232
        """
        E = self.realization_of().elementary_basis()
        L = LatticeOfSetPartitions(phi.size())
        pi = L(phi.to_set_partition())
        mu = L.mobius_function
        return 1/mu(L.bottom(),pi) * E.sum(mu(sigma,pi)*E[label_set_partition(sigma.element)]
                for sigma in L.order_ideal([pi]))

    @lazy_attribute
    def omega(self):
        r"""
        The endomorphism of the Hopf algebra of supercharacters obtained by
        linearly extending the `\omega` map.

        `\omega` is the endomorphism of the Hopf algebra of supercharacters
        defined by mapping the elementary basis element indexed by the set
        partition ``phi`` to the homogeneous basis element indexed by
        ``phi``.

        .. note::

            The map :method:`omega1` is another implementation of this
            same endomorphism.

        EXAMPLES::

            sage: E = SupercharacterHopfAlgebra(q=2).E()
            sage: a = 3*E[4,[(1,3),(2,4)]] + E[4,[(3,4)]]
            sage: a
            3*E[13|24] + E[1|2|34]
            sage: E.omega(a)
            3*H[13|24] + H[1|2|34]
        """
        H = self.realization_of().homogeneous_basis()
        return self._module_morphism(H, codomain=H)

    @cached_method
    def omega1_on_basis(self, phi):
        r"""
        The endomorphism of ``self`` defined on the powersum basis
        by `\omega_1(P_\phi) = (-1)^{|\phi|-\ell(phi)} P_\phi`.

        .. note::

            The map :method:`omega` is another implementation of this
            same endomorphism.

        EXAMPLES::

            sage: P = SupercharacterHopfAlgebra(q=2).P()
            sage: P.omega1_on_basis(LabelledSetPartition(4, [(3,4,1)]))
            -P[1|2|34]
        """
        return self((-1)**(phi.size()-len(phi.to_set_partition())) * self.basis()[phi])

##### miscellaneous code
class SupercharacterTable(UniqueRepresentation, SageObject):
    r"""
    The supercharacter table. This is a lazy implementation, in that no values
    are computed until requested.

    EXAMPLES::

        sage: from sage.combinat.scha import SupercharacterTable
        sage: SupercharacterTable(3)
        Supercharacter table of supercharacters of unipotent upper-triangular matrices over Finite Field of size 3

    To get the table itself, as a matrix, use the :method:`table`::

        sage: SupercharacterTable(q=2).table(3)
        [ 1 -1  1 -1  1]
        [-1 -1  1  1  1]
        [ 0  0 -2  0  2]
        [-1  1  1 -1  1]
        [ 1  1  1  1  1]

    ::

        sage: SupercharacterTable(q=3).table(2)
        [     zeta3 -zeta3 - 1          1]
        [-zeta3 - 1      zeta3          1]
        [         1          1          1]

    One can also pass the Hopf algebra of supercharacters to get the
    corresponding supercharacter table::

        sage: X = SupercharacterHopfAlgebra(2).supercharacter_basis()
        sage: SupercharacterTable(X).table(3)
        [ 1 -1  1 -1  1]
        [-1 -1  1  1  1]
        [ 0  0 -2  0  2]
        [-1  1  1 -1  1]
        [ 1  1  1  1  1]
    """
    def __init__(self, q):
        r"""
        The supercharacter table. This is a lazy implementation, in that no
        values are computed until requested.

        EXAMPLES::

            sage: from sage.combinat.scha import SupercharacterTable
            sage: SupercharacterTable(3)
            Supercharacter table of supercharacters of unipotent upper-triangular matrices over Finite Field of size 3
        """
        if isinstance(q, SupercharacterHopfAlgebraBasis):
            scha = q
            self._q = scha.q()
            self._field = scha._field
            self._base_ring = scha.base_ring()
            if self._q == 2:
                self._zeta = self._base_ring(-1)
                self.theta = self._zeta.__pow__
            else:
                self._zeta = CyclotomicField(self._q.prime_factors()[0]).gen()
        else:
            self._q = q
            self._field = GF(q, 'a')
            if q == 2:
                self._base_ring = QQ
                self._zeta = QQ(-1)
                self.theta = self._zeta.__pow__
            else:
                self._base_ring = CyclotomicField(q)
                self._zeta = CyclotomicField(q.prime_factors()[0]).gen()

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.scha import SupercharacterTable
            sage: SupercharacterTable(2)._repr_()
            'Supercharacter table of supercharacters of unipotent upper-triangular matrices over Finite Field of size 2'
            sage: SupercharacterTable(4)._repr_()
            'Supercharacter table of supercharacters of unipotent upper-triangular matrices over Finite Field in a of size 2^2'
        """
        return "Supercharacter table of supercharacters of unipotent upper-triangular matrices over %s" % (self._field,)

    @cached_method
    def theta(self, x):
        r"""
        This method implements a nontrivial homomorphism from the
        multiplicative group of a finite field into the complex numbers.

        This implementation is defined for `q=p^s` with `p` prime by mapping
        elements of the finite field into the prime subfield (by taking the
        trace), and returns the corresponding power of a primitive `p`-th root
        of unity.

        EXAMPLES::
        
            sage: from sage.combinat.scha import SupercharacterTable
            sage: theta = SupercharacterTable(2).theta
            sage: map(theta, GF(2))
            [1, -1]

        ::

            sage: theta = SupercharacterTable(9).theta
            sage: map(theta, GF(9,'a'))
            [1, -zeta3 - 1, 1, -zeta3 - 1, zeta3, zeta3, 1, zeta3, -zeta3 - 1]
        """
        return self._zeta**(x.trace())

    @cached_method
    def supercharacter_formula(self, phi, psi):
        r"""
        Evaluate the supercharacter indexed by ``phi`` at the superclass
        indexed by ``psi``.

        This implementation is based on Formula (2.1) of [1].

        EXAMPLES::

            sage: from sage.combinat.scha import SupercharacterTable
            sage: F = GF(3)
            sage: sct = SupercharacterTable(3)
            sage: phi = LabelledSetPartition(2, [(1,2,F(1))])
            sage: psi = LabelledSetPartition(2, [(1,2,F(2))])
            sage: sct(phi,phi)
            zeta3
            sage: sct(phi,psi)
            -zeta3 - 1

        REFERENCES:

        -   [1] Nathaniel Thiem, Branching rules in the ring of superclass
            functions of unipotent upper-triangular matrices, J Algebr Comb
            (2010) 31: 267--298
        """
        if phi.size() != psi.size():
            return self._base_ring.zero()
        phi_arcs = phi.arcs_dict()
        psi_arcs = psi.arcs_dict()
        res = self._base_ring.one()
        field_zero = self._field.zero()
        q = self._q
        for (i,l) in phi_arcs:
            if any((i,j) in psi_arcs for j in range(i+1,l)) or \
               any((j,l) in psi_arcs for j in range(i+1,l)):
                   res = self._base_ring.zero()
                   break
            else:
                res = res * q**(l-i-1) * \
                        self.theta(phi_arcs.get((i,l),field_zero)*psi_arcs.get((i,l),field_zero)) \
                        / q**sum(1 for (j,k) in psi_arcs if i<j<l and j<k<l)
        return res

    __call__ = supercharacter_formula

    @cached_method
    def table(self, n):
        r"""
        The supercharacter table for degree ``n`` as a matrix.

        .. note:

            The matrix is computed only once and then cached.
            
        EXAMPLES::

            sage: from sage.combinat.scha import SupercharacterTable
            sage: SupercharacterTable(q=2).table(3)
            [ 1 -1  1 -1  1]
            [-1 -1  1  1  1]
            [ 0  0 -2  0  2]
            [-1  1  1 -1  1]
            [ 1  1  1  1  1]

        ::

            sage: SupercharacterTable(q=3).table(2)
            [     zeta3 -zeta3 - 1          1]
            [-zeta3 - 1      zeta3          1]
            [         1          1          1]

        ::

            sage: SupercharacterTable(q=2).table(4)
            [-1  1 -1 -1  1  1  1 -1 -1  1 -1  1  1 -1  1]
            [ 1  1 -1  1 -1 -1  1 -1 -1  1 -1  1  1  1  1]
            [ 0  0  2  0  0  0 -2  0 -2  2  0  2 -2  0  2]
            [ 0  0  0  2  0  0 -2  0  0 -2  0  2  2 -2  2]
            [ 1 -1  1 -1  1 -1  1 -1  1  1 -1  1  1 -1  1]
            [ 1 -1 -1 -1 -1  1  1  1 -1  1  1  1  1 -1  1]
            [ 0  0  0  0  0  0  4  0  0 -4  0  4 -4  0  4]
            [ 0  0  0  0  0  0  0  2  0  0 -2 -4  0  0  4]
            [-1 -1 -1  1  1 -1  1  1 -1  1  1  1  1  1  1]
            [ 0  0  0 -2  0  0 -2  0  0 -2  0  2  2  2  2]
            [-1 -1  1  1 -1  1  1 -1  1  1 -1  1  1  1  1]
            [ 0  0  0  0  0  0  0 -2  0  0  2 -4  0  0  4]
            [ 0  0 -2  0  0  0 -2  0  2  2  0  2 -2  0  2]
            [-1  1  1 -1 -1 -1  1  1  1  1  1  1  1 -1  1]
            [ 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1]
        """
        from sage.matrix.constructor import matrix
        supercharacter_formula = self.supercharacter_formula
        lsp = LabelledSetPartitions(n, self._q)
        M = {}
        for (i,phi) in enumerate(lsp):
            for (j,psi) in enumerate(lsp):
                M[i,j] = supercharacter_formula(phi,psi)
        return matrix(self._base_ring, M)

    @cached_method
    def table_inverse(self, n):
        r"""
        The inverse of supercharacter table for degree ``n``.

        .. note::

            The matrix is computed only once and then cached.
            
        EXAMPLES::

            sage: from sage.combinat.scha import SupercharacterTable
            sage: SupercharacterTable(2).table_inverse(3)
            [ 1/4 -1/4    0 -1/4  1/4]
            [-1/4 -1/4    0  1/4  1/4]
            [ 1/8  1/8 -1/4  1/8  1/8]
            [-1/4  1/4    0 -1/4  1/4]
            [ 1/8  1/8  1/4  1/8  1/8]

        ::

            sage: SupercharacterTable(3).table_inverse(2)
            [-1/3*zeta3 - 1/3        1/3*zeta3              1/3]
            [       1/3*zeta3 -1/3*zeta3 - 1/3              1/3]
            [             1/3              1/3              1/3]
        """
        return self.table(n).inverse()

##### Lattice of Set Partitions
def LatticeOfSetPartitions(n):
    r"""
    The lattice of set partitions.

    .. note::

        This is need expand the powersum basis in the homogeneous basis.

    TODO: 

    -  this should be improved and eventually included in Sage someplace

    -  Example 3.10.4 of Stanley's Enumerative Combinatorics I gives an
       explicit description of the Mobius function of the lattice of set
       partitions. Use it to implement :meth:`mobius_function`.

    EXAMPLES::

        sage: from sage.combinat.scha import LatticeOfSetPartitions
        sage: LatticeOfSetPartitions(3)
        Finite lattice containing 5 elements

        sage: from sage.combinat.scha import LatticeOfSetPartitions
        sage: L = LatticeOfSetPartitions(2)
        sage: for T in L:
        ...      for S in L:
        ...         print (T,S), L.is_lequal(T,S)
        ({{2}, {1}}, {{2}, {1}}) True
        ({{2}, {1}}, {{1, 2}}) True
        ({{1, 2}}, {{2}, {1}}) False
        ({{1, 2}}, {{1, 2}}) True
    """
    elements = SetPartitions(n)

    def upper_covers_iter(x):
        l = list(x.element) if hasattr(x,'element') else x
        for (s,t) in Subsets(l,2):
            m = l[:]
            m.remove(s)
            m.remove(t)
            m.append(s.union(t))
            yield Set(m)

    relns = [ (S,T) for S in elements for T in upper_covers_iter(S) ]
    L = LatticePoset((elements,relns), cover_relations=False)
    # pre-compute the mobius function for speedup
    L.mobius_function_matrix()

    def partial_order(T,S):
        T, S = T.element, S.element
        return all(any(set(t).issubset(set(s)) for s in S) for t in T)

    L.is_lequal = partial_order

    return L

##### Miscellaneous
def set_partition_to_arcs(partition):
    r"""
    Construct a list of arcs labelled by 1 from a set partition.

    EXAMPLES::

        sage: from sage.combinat.scha import set_partition_to_arcs
        sage: x = Set([Set([1,3]), Set([2])]); x
        {{1, 3}, {2}}
        sage: set_partition_to_arcs(x)
        [(1, 3, 1)]
    """
    one = Integer(1)
    arcs = []
    for part in partition:
        sorted_part = sorted(part)
        for i in range(len(sorted_part)-1):
            arcs.append(tuple(sorted_part[i:i+2]+[one]))
    return arcs

def label_set_partition(partition):
    r"""
    Create a labelled set partition from a set partition; the labels of the
    arcs are taken to be 1.

    .. note::

        This is need to expand the powersum basis in the homogeneous basis.

    EXAMPLES::

        sage: from sage.combinat.scha import label_set_partition
        sage: x = Set([Set([1,3]), Set([2])]); x
        {{1, 3}, {2}}
        sage: label_set_partition(x)
        [3, [(1, 3, 1)]]
    """
    arcs = set_partition_to_arcs(partition)
    return LabelledSetPartition(sum(map(len,partition)), arcs)

def standardize_set_partition(S):
    underlying_set = sorted(sum(map(list,S),[]))
    std = dict((j,i+1) for (i,j) in enumerate(underlying_set))
    return Set([Set([std[i] for i in a]) for a in S])
