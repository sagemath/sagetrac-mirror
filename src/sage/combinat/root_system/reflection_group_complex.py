r"""
Finite complex reflection groups
----------------------------------

Let `V` be a finite-dimensional complex vector space. A reflection of
`V` is an operator `r \in \operatorname{GL}(V)` that has finite order
and fixes pointwise a hyperplane in `V`.

For more definitions and classification types of finite complex
reflection groups, see :wikipedia:`Complex_reflection_group`.

The point of entry to work with reflection groups is :func:`~sage.combinat.root_system.reflection_group_real.ReflectionGroup`
which can be used with finite Cartan-Killing types::

    sage: ReflectionGroup(['A',2])
    Irreducible real reflection group of rank 2 and type A2
    sage: ReflectionGroup(['F',4])
    Irreducible real reflection group of rank 4 and type F4
    sage: ReflectionGroup(['H',3])
    Irreducible real reflection group of rank 3 and type H3

or with Shephard-Todd types::

    sage: ReflectionGroup((1,1,3))
    Irreducible real reflection group of rank 2 and type A2
    sage: ReflectionGroup((2,1,3))
    Irreducible real reflection group of rank 3 and type B3
    sage: ReflectionGroup((3,1,3))
    Irreducible complex reflection group of rank 3 and type G(3,1,3)
    sage: ReflectionGroup((4,2,3))
    Irreducible complex reflection group of rank 3 and type G(4,2,3)
    sage: ReflectionGroup(4)
    Irreducible complex reflection group of rank 2 and type ST4
    sage: ReflectionGroup(31)
    Irreducible complex reflection group of rank 4 and type ST31

Also reducible types are allowed using concatenation::

    sage: ReflectionGroup(['A',3],(4,2,3))
    Reducible complex reflection group of rank 6 and type A3 x G(4,2,3)

Some special cases also occur, among them are::

    sage: W = ReflectionGroup((2,2,2)); W
    Reducible real reflection group of rank 2 and type A1 x A1
    sage: W = ReflectionGroup((2,2,3)); W
    Irreducible real reflection group of rank 3 and type A3

.. WARNING:: Uses the GAP3 package *Chevie* which is available as an
             experimental package (installed by ``sage -i gap3``) or to
             download by hand at `Jean Michel's website <http://webusers.imj-prg.fr/~jean.michel/gap3/>`_.

A guided tour
-------------

We start with the example type `B_2`::

    sage: W = ReflectionGroup(['B',2]); W
    Irreducible real reflection group of rank 2 and type B2

Most importantly, observe that the group elements are usually represented
by permutations of the roots::

    sage: for w in W: print w
    ()
    (1,3)(2,6)(5,7)
    (1,5)(2,4)(6,8)
    (1,7,5,3)(2,4,6,8)
    (1,3,5,7)(2,8,6,4)
    (2,8)(3,7)(4,6)
    (1,7)(3,5)(4,8)
    (1,5)(2,6)(3,7)(4,8)

This has the drawback that one can hardly see anything. Usually, one
would look at elements with either of the following methods::

    sage: for w in W: w.reduced_word()
    []
    [2]
    [1]
    [1, 2]
    [2, 1]
    [2, 1, 2]
    [1, 2, 1]
    [2, 1, 2, 1]

    sage: for w in W: w.reduced_word_in_reflections()
    []
    [2]
    [1]
    [1, 2]
    [1, 4]
    [3]
    [4]
    [1, 3]

    sage: for w in W: w.reduced_word(); w.to_matrix(); print("")
    []
    [1 0]
    [0 1]
    <BLANKLINE>
    [2]
    [ 1  1]
    [ 0 -1]
    <BLANKLINE>
    [1]
    [-1  0]
    [ 2  1]
    <BLANKLINE>
    [1, 2]
    [-1 -1]
    [ 2  1]
    <BLANKLINE>
    [2, 1]
    [ 1  1]
    [-2 -1]
    <BLANKLINE>
    [2, 1, 2]
    [ 1  0]
    [-2 -1]
    <BLANKLINE>
    [1, 2, 1]
    [-1 -1]
    [ 0  1]
    <BLANKLINE>
    [2, 1, 2, 1]
    [-1  0]
    [ 0 -1]
    <BLANKLINE>

The standard references for actions of complex reflection groups have
the matrices acting on the left, so::

    sage: W.simple_reflection(1).to_matrix()
    [-1  0]
    [ 2  1]

sends the simple root `\alpha_0` to its negative, while sending `\alpha_1`
to `2\alpha_0+\alpha_1`.

.. TODO::

    - properly provide root systems for real reflection groups
    - element class should be unique to be able to work with large groups
      without creating elements multiple times.
    - is_shephard_group, is_generalized_coxeter_group
    - exponents & coexponents
    - coinvariant ring:

      * fake degrees from Torsten Hoge
      * operation of linear characters on all characters
      * harmonic polynomials

    - linear forms for hyperplanes
    - field of definition
    - intersection lattice and characteristic polynomial::

        X = [ alpha(t) for t in W.distinguished_reflections() ]
        X = Matrix(CF,X).transpose()
        Y = Matroid(X)

    - linear characters
    - permutation pi on irreducibles
    - hyperplane orbits (76.13 in Gap Manual)
    - improve invariant_form with a code similar to the one in reflection_group_real.py
    - add a method reflection_to_root or distinguished_reflection_to_positive_root.
    - diagrams in ASCII-art (76.15)
    - standard (BMR) presentations
    - character table directly from Chevie
    - GenericOrder (76.20), TorusOrder (76.21)
    - correct fundamental invariants for G34, check the others
    - copy hardcoded data (degrees, invariants, braid relations...) into sage
    - add other hardcoded data from the tables in chevie (location is SAGEDIR/local/gap3/gap-jm5-2015-02-01/gap3/pkg/chevie/tbl): basic derivations, discriminant, ...
    - transfer code for reduced_word_in_reflections into Gap4 or Sage
    - list of reduced words for an element
    - list of reduced words in reflections for an element
    - Hurwitz action?
    - is_crystallographic should be hardcoded

AUTHORS:

- Christian Stump (initial version 2011--2015)
"""
#*****************************************************************************
#       Copyright (C) 2011-2016 Christian Stump <christian.stump at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method, cached_in_parent_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.flatten import flatten
from sage.misc.misc_c import prod
from sage.categories.category import Category
from sage.categories.permutation_groups import PermutationGroups
from sage.categories.complex_reflection_groups import ComplexReflectionGroups
from sage.categories.coxeter_groups import CoxeterGroups
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.sets.family import Family
from sage.structure.unique_representation import UniqueRepresentation
from sage.groups.perm_gps.permgroup import PermutationGroup_generic
from sage.rings.all import ZZ, QQ
from sage.matrix.all import Matrix, identity_matrix
from sage.matrix.matrix import is_Matrix
from sage.interfaces.gap3 import gap3
from sage.combinat.words.word import Word
from sage.rings.universal_cyclotomic_field import E
from sage.arith.misc import lcm
from sage.modules.free_module_element import vector
from sage.combinat.root_system.cartan_matrix import CartanMatrix

from sage.misc.sage_eval import sage_eval

class ComplexReflectionGroup(UniqueRepresentation, PermutationGroup_generic):
    """
    A complex reflection group given as a permutation group.

    .. SEEALSO::

        :func:`ReflectionGroup`
    """
    def __init__(self, W_types, index_set=None, hyperplane_index_set=None, reflection_index_set=None):
        r"""
        TESTS::

            sage: from sage.categories.complex_reflection_groups import ComplexReflectionGroups
            sage: W = ComplexReflectionGroups().example()
            sage: TestSuite(W).run()
        """
        W_components = []
        reflection_type = []
        for W_type in W_types:
            if W_type == (1,1,1):
                raise ValueError("the one element group is not considered a reflection group")
            elif W_type in ZZ:
                call_str = 'ComplexReflectionGroup(%s)'%W_type
            elif isinstance(W_type,CartanMatrix):
                call_str = 'PermRootGroup(IdentityMat(%s),%s)'%(W_type._rank,str(W_type._M._gap_()))
            elif is_Matrix(W_type):
                call_str = 'PermRootGroup(IdentityMat(%s),%s)'%(W_type._rank,str(W_type._gap_()))
            elif W_type in ZZ or ( isinstance(W_type, tuple) and len(W_type) == 3 ):
                call_str = 'ComplexReflectionGroup%s'%str(W_type)
            else:
                if W_type[0] == "I":
                    call_str = 'CoxeterGroup("I",2,%s)'%W_type[1]
                else:
                    call_str = 'CoxeterGroup("%s",%s)'%W_type

            W_components.append(gap3(call_str))
            X = list(W_components[-1].ReflectionType())
            if len(X) > 1:
                raise ValueError("input data %s is invalid"%W_type)
            X = X[0]
            type_dict = dict()
            type_dict["series"] = X.series.sage()
            type_dict["rank"] = X.rank.sage()
            type_dict["indices"] = X.indices.sage()
            if hasattr(X.ST,"sage"):
                type_dict["ST"] = X.ST.sage()
            elif hasattr(X.p,"sage") and hasattr(X.q,"sage"):
                type_dict["ST"] = ( X.p.sage(), X.q.sage(), X.rank.sage() )
            elif hasattr(X.bond,"sage"):
                type_dict["bond"] = X.bond.sage()
            if type_dict["series"] == "B" and (X.cartanType.sage() == 1 or X.indices.sage() == [2,1]):
                type_dict["series"] = "C"
            reflection_type.append( type_dict )

        self._type = reflection_type
        self._gap_group = prod(W_components)
        generators = [str(x) for x in self._gap_group.generators]
        self._index_set = index_set
        self._hyperplane_index_set = hyperplane_index_set
        self._reflection_index_set = reflection_index_set

        self._conjugacy_classes = {}
        self._conjugacy_classes_representatives = None
        self._reflection_representation = None

        self._rank = self._gap_group.rank.sage()
        if len(generators) == self._rank:
            category = ComplexReflectionGroups().Finite().WellGenerated()
            if all(str(W_comp).find('CoxeterGroup') >= 0 for W_comp in W_components):
                category = Category.join([category,CoxeterGroups()])
        else:
            category = ComplexReflectionGroups().Finite()
        if len(self._type) == 1:
            category = category.Irreducible()

        category = Category.join([category,PermutationGroups()]).Finite()

        PermutationGroup_generic.__init__(self, gens=generators,
                                          canonicalize=False,
                                          category=category)

        l_set = range(1, len(self.gens())+1)
        if self._index_set is None:
            self._index_set = tuple(l_set)
        else:
            if len(self._index_set) != len(l_set):
                raise ValueError("the given index set (= %s) does not have the right size"%self._index_set.values())
        self._index_set_inverse = {i: ii for ii,i in enumerate(self._index_set)}
        Nstar_set = range(1,self.number_of_reflecting_hyperplanes()+1)
        if self._hyperplane_index_set is None:
            self._hyperplane_index_set = tuple(Nstar_set)
        else:
            if len(self._hyperplane_index_set) != len(Nstar_set):
                raise ValueError("the given hyperplane index set (= %s) does not have the right size"%self._index_set.values())
        self._hyperplane_index_set_inverse = {i: ii for ii,i in enumerate(self._hyperplane_index_set)}
        N_set = range(1, self.number_of_reflections()+1)
        if self._reflection_index_set is None:
            self._reflection_index_set = tuple(N_set)
        else:
            if len(self._reflection_index_set) != len(N_set):
                raise ValueError("the given reflection index set (= %s) does not have the right size"%self._index_set.values())
        self._reflection_index_set_inverse = {i: ii for ii,i in enumerate(self._reflection_index_set)}

    def _irrcomp_repr_(self,W_type):
        r"""
        Return the string representation of an irreducible component
        of ``self``.

        TESTS::

            sage: W = ReflectionGroup(25,[4,1,4],[1,1,4],[5,5,2]); W
            Reducible complex reflection group of rank 12 and type ST25 x G(4,1,4) x A3 x I2(5)
            sage: for W_type in W._type: print(W._irrcomp_repr_(W_type))
            ST25
             G(4,1,4)
             A3
             I2(5)
        """
        type_str = ''
        if "ST" in W_type:
            if W_type["ST"] in ZZ:
                type_str += "ST" + str(W_type["ST"])
            else:
                type_str += 'G' + str(W_type["ST"]).replace(' ','')
        else:
            type_str += str(W_type["series"])
            if W_type["series"] == "I":
                type_str += '2(' + str(W_type["bond"]) + ')'
            else:
                type_str += str(W_type["rank"])
        return type_str

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(25, [4,1,4],[1,1,4],[5,5,2]); W
            Reducible complex reflection group of rank 12 and type ST25 x G(4,1,4) x A3 x I2(5)
        """
        type_str = ''
        for W_type in self._type:
            type_str += self._irrcomp_repr_(W_type)
            type_str += ' x '
        type_str = type_str[:-3]
        return 'Reducible complex reflection group of rank %s and type %s'%(self._rank,type_str)

    def __iter__(self):
        r"""
        Return an iterator going through all elements in ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: for w in W: w
            ()
            (1,3)(2,5)(4,6)
            (1,4)(2,3)(5,6)
            (1,6,2)(3,5,4)
            (1,2,6)(3,4,5)
            (1,5)(2,4)(3,6)
        """
        for w,word in self._iterator_tracking_words():
            if w._reduced_word is None:
                w._reduced_word = word
            yield w

    # This is the default implementation for any group with generators
    # I leave it here in case it is needed at some point.
    #def _iterator_tracking_words(self):
        #I = self.gens()
        #index_list = range(len(I))
        #elements = [ (self.one(),tuple()) ]
        #elements_set = set( x[0] for x in elements )
        #while elements:
            #x,word = elements.pop()
            #yield x,word
            #for i in index_list:
                #y = x._mul_(I[i])
                #if y not in elements_set:
                    #elements.append((y,word+tuple([i])))
                    #elements_set.add(y)

    def _iterator_tracking_words(self):
        r"""
        Return an iterator through the elements of ``self`` together
        with the words in the simple generators.

        The iterator is a breadth first search through the graph of the
        elements of the group with generators.

        EXAMPLES::

            sage: W = ReflectionGroup(4)
            sage: for w in W._iterator_tracking_words(): w
            ((), [])
            ((1,3,9)(2,4,7)(5,10,18)(6,11,16)(8,12,19)(13,15,20)(14,17,21)(22,23,24), [0])
            ((1,5,13)(2,6,10)(3,7,14)(4,8,15)(9,16,22)(11,12,17)(18,19,23)(20,21,24), [1])
            ((1,9,3)(2,7,4)(5,18,10)(6,16,11)(8,19,12)(13,20,15)(14,21,17)(22,24,23), [0, 0])
            ((1,7,6,12,23,20)(2,8,17,24,9,5)(3,16,10,19,15,21)(4,14,11,22,18,13), [0, 1])
            ((1,10,4,12,21,22)(2,11,19,24,13,3)(5,15,7,17,16,23)(6,18,8,20,14,9), [1, 0])
            ((1,13,5)(2,10,6)(3,14,7)(4,15,8)(9,22,16)(11,17,12)(18,23,19)(20,24,21), [1, 1])
            ((1,16,12,15)(2,14,24,18)(3,5,19,17)(4,6,22,20)(7,8,23,9)(10,13,21,11), [0, 0, 1])
            ((1,2,12,24)(3,6,19,20)(4,17,22,5)(7,11,23,13)(8,21,9,10)(14,16,18,15), [0, 1, 0])
            ((1,14,12,18)(2,15,24,16)(3,22,19,4)(5,6,17,20)(7,10,23,21)(8,11,9,13), [0, 1, 1])
            ((1,18,12,14)(2,16,24,15)(3,4,19,22)(5,20,17,6)(7,21,23,10)(8,13,9,11), [1, 0, 0])
            ((1,15,12,16)(2,18,24,14)(3,17,19,5)(4,20,22,6)(7,9,23,8)(10,11,21,13), [1, 1, 0])
            ((1,6,23)(2,17,9)(3,10,15)(4,11,18)(5,8,24)(7,12,20)(13,14,22)(16,19,21), [0, 0, 1, 0])
            ((1,22,21,12,4,10)(2,3,13,24,19,11)(5,23,16,17,7,15)(6,9,14,20,8,18), [0, 0, 1, 1])
            ((1,4,21)(2,19,13)(3,11,24)(5,7,16)(6,8,14)(9,18,20)(10,12,22)(15,17,23), [0, 1, 0, 0])
            ((1,17,13,12,5,11)(2,20,10,24,6,21)(3,23,14,19,7,18)(4,9,15,22,8,16), [0, 1, 1, 0])
            ((1,19,9,12,3,8)(2,22,7,24,4,23)(5,21,18,17,10,14)(6,13,16,20,11,15), [1, 0, 0, 1])
            ((1,20,23,12,6,7)(2,5,9,24,17,8)(3,21,15,19,10,16)(4,13,18,22,11,14), [1, 1, 0, 0])
            ((1,11,5,12,13,17)(2,21,6,24,10,20)(3,18,7,19,14,23)(4,16,8,22,15,9), [0, 0, 1, 0, 0])
            ((1,23,6)(2,9,17)(3,15,10)(4,18,11)(5,24,8)(7,20,12)(13,22,14)(16,21,19), [0, 0, 1, 1, 0])
            ((1,8,3,12,9,19)(2,23,4,24,7,22)(5,14,10,17,18,21)(6,15,11,20,16,13), [0, 1, 0, 0, 1])
            ((1,21,4)(2,13,19)(3,24,11)(5,16,7)(6,14,8)(9,20,18)(10,22,12)(15,23,17), [0, 1, 1, 0, 0])
            ((1,12)(2,24)(3,19)(4,22)(5,17)(6,20)(7,23)(8,9)(10,21)(11,13)(14,18)(15,16), [0, 0, 1, 0, 0, 1])
            ((1,24,12,2)(3,20,19,6)(4,5,22,17)(7,13,23,11)(8,10,9,21)(14,15,18,16), [0, 0, 1, 1, 0, 0])
        """
        I = tuple(self.simple_reflections())
        index_list = range(len(I))

        level_set_cur = [(self.one(), list())]
        level_set_old = set([ self.one() ])
        while level_set_cur:
            level_set_new = []
            for x, word in level_set_cur:
                yield x, word
                for i in index_list:
                    y = x._mul_(I[i])
                    if y not in level_set_old:
                        level_set_old.add(y)
                        level_set_new.append((y,word+[i]))
            level_set_cur = level_set_new

    __len__ = ComplexReflectionGroups.Finite.ParentMethods.cardinality.__func__

    @cached_method
    def index_set(self):
        r"""
        Return the index set of the simple reflections of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,4))
            sage: W.index_set()
            (1, 2, 3)
            sage: W = ReflectionGroup((1,1,4), index_set=[1,3,'asdf'])
            sage: W.index_set()
            (1, 3, 'asdf')
            sage: W = ReflectionGroup((1,1,4), index_set=('a', 'b', 'c'))
            sage: W.index_set()
            ('a', 'b', 'c')
        """
        return self._index_set

    def series(self):
        r"""
        Return the series of the classification type to which ``self``
        belongs.

        For real reflection groups, these are the Cartan-Killing
        classification types "A","B","C","D","E","F","G","H","I", and
        for complx non-real reflection groups these are the
        Shephard-Todd classification type "ST".

        EXAMPLES::

            sage: ReflectionGroup((1,1,3)).series()
            ['A']
            sage: ReflectionGroup((3,1,3)).series()
            ['ST']
        """
        return [self._type[i]['series'] for i in range(len(self._type))]

    @cached_method
    def hyperplane_index_set(self):
        r"""
        Return the index set of the hyperplanes of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,4))
            sage: W.hyperplane_index_set()
            (1, 2, 3, 4, 5, 6)
            sage: W = ReflectionGroup((1,1,4), hyperplane_index_set=[1,3,'asdf',7,9,11])
            sage: W.hyperplane_index_set()
            (1, 3, 'asdf', 7, 9, 11)
            sage: W = ReflectionGroup((1,1,4),hyperplane_index_set=('a','b','c','d','e','f'))
            sage: W.hyperplane_index_set()
            ('a', 'b', 'c', 'd', 'e', 'f')
        """
        return self._hyperplane_index_set

    @cached_method
    def distinguished_reflections(self):
        r"""
        Return a finite family containing the distinguished reflections
        of ``self`` indexed by :meth:`hyperplane_index_set`.

        These are the reflections in ``self`` acting on the complement
        of the fixed hyperplane `H` as `\operatorname{exp}(2 \pi i / n)`,
        where `n` is the order of the reflection subgroup fixing `H`.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.distinguished_reflections()
            Finite family {1: (1,4)(2,3)(5,6), 2: (1,3)(2,5)(4,6), 3: (1,5)(2,4)(3,6)}

            sage: W = ReflectionGroup((1,1,3),hyperplane_index_set=['a','b','c'])
            sage: W.distinguished_reflections()
            Finite family {'a': (1,4)(2,3)(5,6), 'c': (1,5)(2,4)(3,6), 'b': (1,3)(2,5)(4,6)}

            sage: W = ReflectionGroup((3,1,1))
            sage: W.distinguished_reflections()
            Finite family {1: (1,2,3)}

            sage: W = ReflectionGroup((1,1,3),(3,1,2))
            sage: W.distinguished_reflections()
            Finite family {1: (1,6)(2,5)(7,8), 2: (1,5)(2,7)(6,8),
             3: (3,9,15)(4,10,16)(12,17,23)(14,18,24)(20,25,29)(21,22,26)(27,28,30),
             4: (3,11)(4,12)(9,13)(10,14)(15,19)(16,20)(17,21)(18,22)(23,27)(24,28)(25,26)(29,30),
             5: (1,7)(2,6)(5,8),
             6: (3,19)(4,25)(9,11)(10,17)(12,28)(13,15)(14,30)(16,18)(20,27)(21,29)(22,23)(24,26),
             7: (4,21,27)(10,22,28)(11,13,19)(12,14,20)(16,26,30)(17,18,25)(23,24,29),
             8: (3,13)(4,24)(9,19)(10,29)(11,15)(12,26)(14,21)(16,23)(17,30)(18,27)(20,22)(25,28)}
        """
        # makes sure that the simple reflections come first
        gens = self.gens()
        R = [t for t in gens]
        # Then import all distinguished reflections from gap,
        #   the Set is used as every such appears multiple times.
        for r in self._gap_group.Reflections():
            t = self(str(r))
            if t not in R:
                R.append(t)
        return Family(self._hyperplane_index_set,
                      lambda i: R[self._hyperplane_index_set_inverse[i]])

    def distinguished_reflection(self, i):
        r"""
        Return the ``i``-th distinguished reflection of ``self``.

        These are the reflections in ``self`` acting on the complement
        of the fixed hyperplane `H` as `\operatorname{exp}(2 \pi i / n)`,
        where `n` is the order of the reflection subgroup fixing `H`.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.distinguished_reflection(1)
            (1,4)(2,3)(5,6)
            sage: W.distinguished_reflection(2)
            (1,3)(2,5)(4,6)
            sage: W.distinguished_reflection(3)
            (1,5)(2,4)(3,6)

            sage: W = ReflectionGroup((3,1,1),hyperplane_index_set=['a'])
            sage: W.distinguished_reflection('a')
            (1,2,3)

            sage: W = ReflectionGroup((1,1,3),(3,1,2))
            sage: for i in range(W.number_of_reflecting_hyperplanes()):
            ....:     W.distinguished_reflection(i+1)
            (1,6)(2,5)(7,8)
            (1,5)(2,7)(6,8)
            (3,9,15)(4,10,16)(12,17,23)(14,18,24)(20,25,29)(21,22,26)(27,28,30)
            (3,11)(4,12)(9,13)(10,14)(15,19)(16,20)(17,21)(18,22)(23,27)(24,28)(25,26)(29,30)
            (1,7)(2,6)(5,8)
            (3,19)(4,25)(9,11)(10,17)(12,28)(13,15)(14,30)(16,18)(20,27)(21,29)(22,23)(24,26)
            (4,21,27)(10,22,28)(11,13,19)(12,14,20)(16,26,30)(17,18,25)(23,24,29)
            (3,13)(4,24)(9,19)(10,29)(11,15)(12,26)(14,21)(16,23)(17,30)(18,27)(20,22)(25,28)
        """
        return self.distinguished_reflections()[i]

    @cached_method
    def reflecting_hyperplanes(self, as_linear_functionals=False):
        r"""
        Return the list of all reflecting hyperplanes of ``self``,
        either as a codimension 1 space, or as its linear functional.

        INPUT:

        - ``as_linear_functionals`` -- (default:``False``) flag whether
          to return the hyperplane or its linear functional in the basis
          dual to the given root basis

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: for H in W.reflecting_hyperplanes(): H
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 2]
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [  1 1/2]
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [ 1 -1]

            sage: for H in W.reflecting_hyperplanes(as_linear_functionals=True): H
            (1, -1/2)
            (1, -2)
            (1, 1)


            sage: W = ReflectionGroup((2,1,2))
            sage: for H in W.reflecting_hyperplanes(): H
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 1]
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [  1 1/2]
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 0]
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [0 1]

            sage: for H in W.reflecting_hyperplanes(as_linear_functionals=True): H
            (1, -1)
            (1, -2)
            (0, 1)
            (1, 0)
        """
        Hs = []
        for r in self.distinguished_reflections():
            mat = (r.to_matrix().transpose() - identity_matrix(self.rank()))
            if as_linear_functionals:
                Hs.append( mat.row_space().gen() )
            else:
                Hs.append( mat.right_kernel() )
        return Family(self._hyperplane_index_set,
                      lambda i: Hs[self._hyperplane_index_set_inverse[i]])

    def reflecting_hyperplane(self, i, as_linear_functional=False):
        r"""
        Return the ``i``-th reflecting hyperplane of ``self``.

        The ``i``-th reflecting hyperplane corresponds to the ``i``
        distinguished reflection.

        INPUT:

        - ``i`` -- an index in the index set
        - ``as_linear_functionals`` -- (default:``False``) flag whether
          to return the hyperplane or its linear functional in the basis
          dual to the given root basis

        EXAMPLES::

            sage: W = ReflectionGroup((2,1,2))
            sage: W.reflecting_hyperplane(3)
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 0]

        One can ask for the result as a linear form::

            sage: W.reflecting_hyperplane(3, True)
            (0, 1)
        """
        return self.reflecting_hyperplanes(as_linear_functionals=as_linear_functional)[i]

    @cached_method
    def reflection_index_set(self):
        r"""
        Return the index set of the reflections of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,4))
            sage: W.reflection_index_set()
            (1, 2, 3, 4, 5, 6)
            sage: W = ReflectionGroup((1,1,4), reflection_index_set=[1,3,'asdf',7,9,11])
            sage: W.reflection_index_set()
            (1, 3, 'asdf', 7, 9, 11)
            sage: W = ReflectionGroup((1,1,4), reflection_index_set=('a','b','c','d','e','f'))
            sage: W.reflection_index_set()
            ('a', 'b', 'c', 'd', 'e', 'f')
        """
        return self._reflection_index_set

    @cached_method
    def reflections(self):
        r"""
        Return a finite family containing the reflections of ``self``,
        indexed by :meth:`self.reflection_index_set`.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.reflections()
            Finite family {1: (1,4)(2,3)(5,6), 2: (1,3)(2,5)(4,6), 3: (1,5)(2,4)(3,6)}

            sage: W = ReflectionGroup((1,1,3),reflection_index_set=['a','b','c'])
            sage: W.reflections()
            Finite family {'a': (1,4)(2,3)(5,6), 'c': (1,5)(2,4)(3,6), 'b': (1,3)(2,5)(4,6)}

            sage: W = ReflectionGroup((3,1,1))
            sage: W.reflections()
            Finite family {1: (1,2,3), 2: (1,3,2)}

            sage: W = ReflectionGroup((1,1,3),(3,1,2))
            sage: W.reflections()
            Finite family {1: (1,6)(2,5)(7,8), 2: (1,5)(2,7)(6,8),
                           3: (3,9,15)(4,10,16)(12,17,23)(14,18,24)(20,25,29)(21,22,26)(27,28,30),
                           4: (3,11)(4,12)(9,13)(10,14)(15,19)(16,20)(17,21)(18,22)(23,27)(24,28)(25,26)(29,30),
                           5: (1,7)(2,6)(5,8),
                           6: (3,19)(4,25)(9,11)(10,17)(12,28)(13,15)(14,30)(16,18)(20,27)(21,29)(22,23)(24,26),
                           7: (4,21,27)(10,22,28)(11,13,19)(12,14,20)(16,26,30)(17,18,25)(23,24,29),
                           8: (3,13)(4,24)(9,19)(10,29)(11,15)(12,26)(14,21)(16,23)(17,30)(18,27)(20,22)(25,28),
                           9: (3,15,9)(4,16,10)(12,23,17)(14,24,18)(20,29,25)(21,26,22)(27,30,28),
                           10: (4,27,21)(10,28,22)(11,19,13)(12,20,14)(16,30,26)(17,25,18)(23,29,24)}
        """
        T = self.distinguished_reflections().values()
        for i in range(self.number_of_reflecting_hyperplanes()):
            for j in range(2, T[i].order()):
                T.append(T[i]**j)
        return Family(self._reflection_index_set,
                      lambda i: T[self._reflection_index_set_inverse[i]])

    def reflection(self,i):
        r"""
        Return the ``i``-th reflection of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.reflection(1)
            (1,4)(2,3)(5,6)
            sage: W.reflection(2)
            (1,3)(2,5)(4,6)
            sage: W.reflection(3)
            (1,5)(2,4)(3,6)

            sage: W = ReflectionGroup((3,1,1),reflection_index_set=['a','b'])
            sage: W.reflection('a')
            (1,2,3)
            sage: W.reflection('b')
            (1,3,2)
        """
        return self.reflections()[i]

    def reflection_character(self):
        r"""
        Return the reflection characters of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.reflection_character()
            [2, 0, -1]
        """
        return self._gap_group.ReflectionCharacter().sage()

    def is_crystallographic(self):
        r"""
        Return ``True`` if self is crystallographic.

        This is, if the field of definition is the rational field.

        .. TODO::

            Make this more robust and do not use the matrix
            representation of the simple reflections.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3)); W
            Irreducible real reflection group of rank 2 and type A2
            sage: W.is_crystallographic()
            True

            sage: W = ReflectionGroup((2,1,3)); W
            Irreducible real reflection group of rank 3 and type B3
            sage: W.is_crystallographic()
            True

            sage: W = ReflectionGroup(23); W
            Irreducible real reflection group of rank 3 and type H3
            sage: W.is_crystallographic()
            False

            sage: W = ReflectionGroup((3,1,3)); W
            Irreducible complex reflection group of rank 3 and type G(3,1,3)
            sage: W.is_crystallographic()
            False

            sage: W = ReflectionGroup((4,2,2)); W
            Irreducible complex reflection group of rank 2 and type G(4,2,2)
            sage: W.is_crystallographic()
            False
        """
        return self.is_real() and all(t.to_matrix().base_ring() is QQ for t in self.simple_reflections())

    def _element_class(self):
        r"""
        A temporary workaround for compatibility with Sage's
        permutation groups.

        TESTS::

            sage: W = ReflectionGroup(23)                   # optional - chevie
            sage: W._element_class() is W.element_class     # optional - chevie
            True
        """
        return self.element_class

    def number_of_irreducible_components(self):
        r"""
        Return the number of irreducible components of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.number_of_irreducible_components()
            1

            sage: W = ReflectionGroup((1,1,3),(2,1,3))
            sage: W.number_of_irreducible_components()
            2
        """
        return len(self._type)

    def irreducible_components(self):
        r"""
        Return a list containing the irreducible components of ``self``
        as finite reflection groups.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.irreducible_components()
            [Irreducible real reflection group of rank 2 and type A2]

            sage: W = ReflectionGroup((1,1,3),(2,1,3))
            sage: W.irreducible_components()
            [Irreducible real reflection group of rank 2 and type A2,
            Irreducible real reflection group of rank 3 and type B3]
        """
        from sage.combinat.root_system.reflection_group_real import ReflectionGroup
        irr_comps = []
        for W_type in self._type:
            if W_type["series"] in ["A","B","D","E","F","G","H","I"]:
                W_str = (W_type["series"],W_type["rank"])
            elif "ST" in W_type:
                W_str = W_type["ST"]
            irr_comps.append(ReflectionGroup(W_str))
        return irr_comps

    @cached_method
    def conjugacy_classes_representatives(self):
        r"""
        Return the shortest representatives of the conjugacy classes of
        ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: [w.reduced_word() for w in W.conjugacy_classes_representatives()]
            [[], [1], [1, 2]]

            sage: W = ReflectionGroup((1,1,4))
            sage: [w.reduced_word() for w in W.conjugacy_classes_representatives()]
            [[], [1], [3, 1], [1, 2], [3, 1, 2]]

            sage: W = ReflectionGroup((3,1,2))
            sage: [w.reduced_word() for w in W.conjugacy_classes_representatives()]
            [[], [1], [1, 1], [2, 1, 2, 1], [2, 1, 2, 1, 1],
             [2, 1, 1, 2, 1, 1], [2], [1, 2], [1, 1, 2]]

            sage: W = ReflectionGroup(23)
            sage: [w.reduced_word() for w in W.conjugacy_classes_representatives()]
            [[], [1], [1, 2], [3, 1], [2, 3], [1, 2, 3], [1, 2, 1, 2],
             [1, 2, 1, 2, 3], [1, 2, 3, 1, 2, 1, 2, 3, 2],
             [3, 2, 1, 2, 3, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1]]
        """
        if self._conjugacy_classes_representatives is None:
            S = str(gap3('List(ConjugacyClasses(%s),Representative)'%self._gap_group._name))
            exec('self._conjugacy_classes_representatives=' + _gap_return(S))
        return self._conjugacy_classes_representatives

    @cached_method
    def conjugacy_classes(self):
        r"""
        Return the conjugacy classes of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: for C in W.conjugacy_classes(): sorted(C)
            [()]
            [(1,3)(2,5)(4,6), (1,4)(2,3)(5,6), (1,5)(2,4)(3,6)]
            [(1,2,6)(3,4,5), (1,6,2)(3,5,4)]

            sage: W = ReflectionGroup((1,1,4))
            sage: sum(len(C) for C in W.conjugacy_classes()) == W.cardinality()
            True

            sage: W = ReflectionGroup((3,1,2))
            sage: sum(len(C) for C in W.conjugacy_classes()) == W.cardinality()
            True

            sage: W = ReflectionGroup(23)
            sage: sum(len(C) for C in W.conjugacy_classes()) == W.cardinality()
            True
       """
        return Family(self.conjugacy_classes_representatives(),
                      lambda w: w.conjugacy_class())

    def rank(self):
        r"""
        Return the rank of ``self``.

        This is the dimension of the underlying vector space.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.rank()
            2
            sage: W = ReflectionGroup((2,1,3))
            sage: W.rank()
            3
            sage: W = ReflectionGroup((4,1,3))
            sage: W.rank()
            3
            sage: W = ReflectionGroup((4,2,3))
            sage: W.rank()
            3
        """
        return self._rank

    @cached_method
    def degrees(self):
        r"""
        Return the degrees of ``self`` ordered within each irreducible
        component of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,4))
            sage: W.degrees()
            [2, 3, 4]

            sage: W = ReflectionGroup((2,1,4))
            sage: W.degrees()
            [2, 4, 6, 8]

            sage: W = ReflectionGroup((4,1,4))
            sage: W.degrees()
            [4, 8, 12, 16]

            sage: W = ReflectionGroup((4,2,4))
            sage: W.degrees()
            [4, 8, 12, 8]

            sage: W = ReflectionGroup((4,4,4))
            sage: W.degrees()
            [4, 8, 12, 4]

        Examples of reducible types::

            sage: W = ReflectionGroup((1,1,4), (3,1,2)); W
            Reducible complex reflection group of rank 5 and type A3 x G(3,1,2)
            sage: W.degrees()
            [2, 3, 4, 3, 6]

            sage: W = ReflectionGroup((1,1,4), (6,1,12), 23) # fails in GAP3
            sage: W.degrees()
            [2, 3, 4, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 2, 6, 10]
        """
        if self.is_irreducible():
            try:
                return self._gap_group.degrees.sage()
            except:
                return self._gap_group.ReflectionDegrees().sage()
        else:
            return flatten([comp.degrees() for comp in self.irreducible_components()])

    cardinality = ComplexReflectionGroups.Finite.ParentMethods.cardinality.__func__

    @cached_method
    def codegrees(self):
        r"""
        Return the codegrees of ``self`` ordered within each irreducible
        component of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,4))
            sage: W.codegrees()
            [2, 1, 0]

            sage: W = ReflectionGroup((2,1,4))
            sage: W.codegrees()
            [6, 4, 2, 0]

            sage: W = ReflectionGroup((4,1,4))
            sage: W.codegrees()
            [12, 8, 4, 0]

            sage: W = ReflectionGroup((4,2,4))
            sage: W.codegrees()
            [12, 8, 4, 0]

            sage: W = ReflectionGroup((4,4,4))
            sage: W.codegrees()
            [8, 4, 0, 8]

            sage: W = ReflectionGroup((1,1,4), (3,1,2))
            sage: W.codegrees()
            [2, 1, 0, 3, 0]

            sage: W = ReflectionGroup((1,1,4), (6,1,12), 23) # fails in GAP3
            sage: W.codegrees()
            [2, 1, 0, 66, 60, 54, 48, 42, 36, 30, 24, 18, 12, 6, 0, 8, 4, 0]
        """
        if self.is_irreducible():
            if self.is_well_generated():
                h = self.coxeter_number()
                return [h-d for d in self.degrees()]
            else:
                return sorted(self._gap_group.ReflectionCoDegrees().sage(), reverse=True)
        else:
            return flatten([comp.codegrees() for comp in self.irreducible_components()])

    @cached_method
    def reflection_eigenvalues_family(self):
        r"""
        Return the reflection eigenvalues of ``self`` as a finite family
        indexed by the class representatives of ``self``.

        OUTPUT:

        - list with entries `k/n` representing the eigenvalue `\zeta_n^k`.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.reflection_eigenvalues_family()
            Finite family {(): [0, 0], (1,4)(2,3)(5,6): [1/2, 0], (1,6,2)(3,5,4): [1/3, 2/3]}

            sage: W = ReflectionGroup((3,1,2))
            sage: reflection_eigenvalues = W.reflection_eigenvalues_family()
            sage: for elt in sorted(reflection_eigenvalues.keys()):
            ....:     print('%s %s'%(elt, reflection_eigenvalues[elt]))
            () [0, 0]
            (1,3,9)(2,4,10)(6,11,17)(8,12,18)(14,19,23)(15,16,20)(21,22,24) [1/3, 0]
            (1,3,9)(2,16,24)(4,20,21)(5,7,13)(6,12,23)(8,19,17)(10,15,22)(11,18,14) [1/3, 1/3]
            (1,5)(2,6)(3,7)(4,8)(9,13)(10,14)(11,15)(12,16)(17,21)(18,22)(19,20)(23,24) [1/2, 0]
            (1,7,3,13,9,5)(2,8,16,19,24,17)(4,14,20,11,21,18)(6,15,12,22,23,10) [1/6, 2/3]
            (1,9,3)(2,10,4)(6,17,11)(8,18,12)(14,23,19)(15,20,16)(21,24,22) [2/3, 0]
            (1,9,3)(2,20,22)(4,15,24)(5,7,13)(6,18,19)(8,23,11)(10,16,21)(12,14,17) [1/3, 2/3]
            (1,9,3)(2,24,16)(4,21,20)(5,13,7)(6,23,12)(8,17,19)(10,22,15)(11,14,18) [2/3, 2/3]
            (1,13,9,7,3,5)(2,14,24,18,16,11)(4,6,21,23,20,12)(8,22,17,15,19,10) [1/3, 5/6]

            sage: W = ReflectionGroup(23)
            sage: reflection_eigenvalues = W.reflection_eigenvalues_family()
            sage: for elt in sorted(reflection_eigenvalues.keys()):
            ....:     print('%s %s'%(elt, reflection_eigenvalues[elt]))
            () [0, 0, 0]
            (1,8,4)(2,21,3)(5,10,11)(6,18,17)(7,9,12)(13,14,15)(16,23,19)(20,25,26)(22,24,27)(28,29,30) [1/3, 2/3, 0]
            (1,16)(2,5)(4,7)(6,9)(8,10)(11,13)(12,14)(17,20)(19,22)(21,24)(23,25)(26,28)(27,29) [1/2, 0, 0]
            (1,16)(2,9)(3,18)(4,10)(5,6)(7,8)(11,14)(12,13)(17,24)(19,25)(20,21)(22,23)(26,29)(27,28) [1/2, 1/2, 0]
            (1,16)(2,17)(3,18)(4,19)(5,20)(6,21)(7,22)(8,23)(9,24)(10,25)(11,26)(12,27)(13,28)(14,29)(15,30) [1/2, 1/2, 1/2]
            (1,19,20,2,7)(3,6,11,13,9)(4,5,17,22,16)(8,12,15,14,10)(18,21,26,28,24)(23,27,30,29,25) [1/5, 4/5, 0]
            (1,20,7,19,2)(3,11,9,6,13)(4,17,16,5,22)(8,15,10,12,14)(18,26,24,21,28)(23,30,25,27,29) [2/5, 3/5, 0]
            (1,23,26,29,22,16,8,11,14,7)(2,10,4,9,18,17,25,19,24,3)(5,21,27,30,28,20,6,12,15,13) [1/10, 1/2, 9/10]
            (1,24,17,16,9,2)(3,12,13,18,27,28)(4,21,29,19,6,14)(5,25,26,20,10,11)(7,23,30,22,8,15) [1/6, 1/2, 5/6]
            (1,29,8,7,26,16,14,23,22,11)(2,9,25,3,4,17,24,10,18,19)(5,30,6,13,27,20,15,21,28,12) [3/10, 1/2, 7/10]
            """
        class_representatives = self.conjugacy_classes_representatives()
        Ev_list = self._gap_group.ReflectionEigenvalues().sage()
        return Family(class_representatives,
                      lambda w: Ev_list[class_representatives.index(w)])

    @cached_method
    def reflection_eigenvalues(self, w, test_class_repr=True):
        r"""
        Return the reflection eigenvalue of ``w`` in ``self``.

        INPUT:

        - ``test_class_repr`` -- boolean (default ``True``) whether to
          compute instead on the conjugacy class representative

        .. SEEALSO:: :meth:`reflection_eigenvalues_family`

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: for w in W:
            ....:     print('%s %s'%(w.reduced_word(), W.reflection_eigenvalues(w)))
            [] [0, 0]
            [2] [1/2, 0]
            [1] [1/2, 0]
            [1, 2] [1/3, 2/3]
            [2, 1] [1/3, 2/3]
            [1, 2, 1] [1/2, 0]
        """
        if test_class_repr:
            w_repr = w.conjugacy_class_representative()
        else:
            w_repr = w
        return self.reflection_eigenvalues_family()[w_repr]

    @cached_method
    def simple_roots(self):
        r"""
        Return the simple roots of ``self``.

        These are the roots corresponding to the simple reflections.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.simple_roots()
            [(1, 0), (0, 1)]

            sage: W = ReflectionGroup((1,1,4), (2,1,2))
            sage: W.simple_roots()
            [(1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 1, 0, 0),
             (0, 0, 0, 1, 0), (0, 0, 0, 0, 1)]

            sage: W = ReflectionGroup((3,1,2))
            sage: W.simple_roots()
            [(1, 0), (-1, 1)]

            sage: W = ReflectionGroup((1,1,4), (3,1,2))
            sage: W.simple_roots()
            [(1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 1, 0, 0),
             (0, 0, 0, 1, 0), (0, 0, 0, -1, 1)]
        """
        return self.roots()[:len(self.gens())]

    def simple_root(self, i):
        r"""
        Return the simple root with index ``i``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3])
            sage: W.simple_root(1)
            (1, 0, 0)
        """
        return self.simple_roots()[self._index_set_inverse[i]]

    @cached_method
    def simple_coroots(self):
        r"""
        Return the simple coroots of ``self``.

        These are the coroots corresponding to the simple reflections.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.simple_coroots()
            [[2, -1], [-1, 2]]

            sage: W = ReflectionGroup((1,1,4), (2,1,2))
            sage: W.simple_coroots()
            [[2, -1, 0, 0, 0],
             [-1, 2, -1, 0, 0],
             [0, -1, 2, 0, 0],
             [0, 0, 0, 2, -2],
             [0, 0, 0, -1, 2]]

            sage: W = ReflectionGroup((3,1,2))
            sage: W.simple_coroots()
            [[-2*E(3) - E(3)^2, 0], [-1, 1]]

            sage: W = ReflectionGroup((1,1,4), (3,1,2))
            sage: W.simple_coroots()
            [[2, -1, 0, 0, 0],
             [-1, 2, -1, 0, 0],
             [0, -1, 2, 0, 0],
             [0, 0, 0, -2*E(3) - E(3)^2, 0],
             [0, 0, 0, -1, 1]]
        """
        return self._gap_group.simpleCoroots.sage()

    @cached_method
    def independent_roots(self):
        r"""
        Return a collection of simple roots generating the underlying
        vector space of ``self``.

        For well-generated groups, these are all simple roots.
        Otherwise, a linearly independent subset of the simple roots is
        chosen.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.independent_roots()
            [(1, 0), (0, 1)]

            sage: W = ReflectionGroup((4,2,3))
            sage: W.simple_roots()
            [(1, 0, 0), (-E(4), 1, 0), (-1, 1, 0), (0, -1, 1)]
            sage: W.independent_roots()
            [(1, 0, 0), (-E(4), 1, 0), (0, -1, 1)]
        """
        Delta = self.simple_roots()
        if len(Delta) == self.rank():
            basis = Delta
        else:
            basis = []
            for alpha in Delta:
                if Matrix(basis+[alpha]).rank() == len(basis) + 1:
                    basis.append(alpha)
        return basis

    @cached_method
    def base_change_matrix(self):
        r"""
        Return the base change from the standard basis of the vector
        space of ``self`` to the basis given by the independent roots of
        ``self``.

        .. TODO::

            For non-well-generated groups there is a conflict with
            construction of the matrix for an element.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.base_change_matrix()
            [1 0]
            [0 1]

            sage: W = ReflectionGroup(23)
            sage: W.base_change_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: W = ReflectionGroup((3,1,2))
            sage: W.base_change_matrix()
            [1 0]
            [1 1]

            sage: W = ReflectionGroup((4,2,2))
            sage: W.base_change_matrix()
            [   1    0]
            [E(4)    1]
        """
        return Matrix( self.independent_roots() ).inverse()

    @cached_method
    def roots(self):
        r"""
        Return all roots corresponding to all reflections of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.roots()
            [(1, 0), (0, 1), (1, 1), (-1, 0), (0, -1), (-1, -1)]

            sage: W = ReflectionGroup((3,1,2))
            sage: W.roots()
            [(1, 0), (-1, 1), (E(3), 0), (-E(3), 1), (0, 1), (1, -1),
             (0, E(3)), (1, -E(3)), (E(3)^2, 0), (-E(3)^2, 1),
             (E(3), -1), (E(3), -E(3)), (0, E(3)^2), (1, -E(3)^2),
             (-1, E(3)), (-E(3), E(3)), (E(3)^2, -1), (E(3)^2, -E(3)),
             (E(3), -E(3)^2), (-E(3)^2, E(3)), (-1, E(3)^2),
             (-E(3), E(3)^2), (E(3)^2, -E(3)^2), (-E(3)^2, E(3)^2)]

            sage: W = ReflectionGroup((4,2,2))
            sage: W.roots()
            [(1, 0), (-E(4), 1), (-1, 1), (-1, 0), (E(4), 1), (1, 1),
             (0, -E(4)), (E(4), -1), (E(4), E(4)), (0, E(4)),
             (E(4), -E(4)), (0, 1), (1, -E(4)), (1, -1), (0, -1),
             (1, E(4)), (-E(4), 0), (-1, E(4)), (E(4), 0), (-E(4), E(4)),
             (-E(4), -1), (-E(4), -E(4)), (-1, -E(4)), (-1, -1)]

            sage: W = ReflectionGroup((1,1,4), (3,1,2))
            sage: W.roots()
            [(1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 1, 0, 0),
             (0, 0, 0, 1, 0), (0, 0, 0, -1, 1), (1, 1, 0, 0, 0),
             (0, 1, 1, 0, 0), (1, 1, 1, 0, 0), (-1, 0, 0, 0, 0),
             (0, -1, 0, 0, 0), (0, 0, -1, 0, 0), (-1, -1, 0, 0, 0),
             (0, -1, -1, 0, 0), (-1, -1, -1, 0, 0), (0, 0, 0, E(3), 0),
             (0, 0, 0, -E(3), 1), (0, 0, 0, 0, 1), (0, 0, 0, 1, -1),
             (0, 0, 0, 0, E(3)), (0, 0, 0, 1, -E(3)), (0, 0, 0, E(3)^2, 0),
             (0, 0, 0, -E(3)^2, 1), (0, 0, 0, E(3), -1), (0, 0, 0, E(3), -E(3)),
             (0, 0, 0, 0, E(3)^2), (0, 0, 0, 1, -E(3)^2), (0, 0, 0, -1, E(3)),
             (0, 0, 0, -E(3), E(3)), (0, 0, 0, E(3)^2, -1),
             (0, 0, 0, E(3)^2, -E(3)), (0, 0, 0, E(3), -E(3)^2),
             (0, 0, 0, -E(3)^2, E(3)), (0, 0, 0, -1, E(3)^2),
             (0, 0, 0, -E(3), E(3)^2), (0, 0, 0, E(3)^2, -E(3)^2),
             (0, 0, 0, -E(3)^2, E(3)^2)]
        """
        roots = [vector(sage_eval(str(root).replace("^","**")))
                 for root in self._gap_group.roots]
        for v in roots:
            v.set_immutable()
        return roots

    @cached_method
    def braid_relations(self):
        r"""
        Return the braid relations of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.braid_relations()
            [[[2, 1, 2], [1, 2, 1]]]

            sage: W = ReflectionGroup((2,1,3))
            sage: W.braid_relations()
            [[[2, 1, 2, 1], [1, 2, 1, 2]], [[3, 1], [1, 3]], [[3, 2, 3], [2, 3, 2]]]

            sage: W = ReflectionGroup((2,2,3))
            sage: W.braid_relations()
            [[[2, 1, 2], [1, 2, 1]], [[3, 1], [1, 3]], [[3, 2, 3], [2, 3, 2]]]
        """
        return self._gap_group.BraidRelations().sage()

    def fundamental_invariants(self):
        r"""
        Return the fundamental invariants of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: W.fundamental_invariants()
            [-2*x0^2 + 2*x0*x1 - 2*x1^2, 6*x0^2*x1 - 6*x0*x1^2]

            sage: W = ReflectionGroup((3,1,2))
            sage: W.fundamental_invariants()
            [x0^3 + x1^3, x0^3*x1^3]
        """
        import re
        from sage.rings.polynomial.all import PolynomialRing

        I = [ str(p) for p in gap3('List(Invariants(%s),x->ApplyFunc(x,List([0..%s],i->Mvp(SPrint("x",i)))))'%(self._gap_group._name,self.rank()-1)) ]
        P = PolynomialRing(QQ,['x%s'%i for i in range(0,self.rank())])
        x = P.gens()
        for i in range(len(I)):
            I[i] = I[i].replace('^','**')
            I[i] = re.compile('E(\d\d*)').sub(r'E(\1)', I[i])
            I[i] = re.compile('(\d)E\(').sub(r'\1*E(', I[i])
            for j in range(len(x)):
                I[i] = I[i].replace('x%s'%j,'*x[%s]'%j)
            I[i] = I[i].replace("+*","+").replace("-*","-").replace("ER(5)","*(E(5)-E(5)**2-E(5)**3+E(5)**4)").lstrip("*")
        # sage_eval is used since eval kills the rational entries!
        I = [ sage_eval(p, locals={'x':x}) for p in I ]
        return I

    def cartan_matrix(self):
        r"""
        Return the Cartan matrix associated with ``self``.

        If ``self`` is crystallographic, the returned Cartan matrix is
        an instance of :class:`CartanMatrix`, and a normal matrix
        otherwise.

        Let `s_1,\ldots,s_n` be a set of reflections which generate
        ``self`` with associated simple roots `s_1,\ldots,s_n` and
        simple coroots `s^\vee_i`. Then the Cartan matrix `C = (c_{ij})`
        is given by `s^\vee_i(s_j)`. The Cartan matrix completely
        determines the reflection representation if the `s_i` are
        linearly independent.

        EXAMPLES::

            sage: ReflectionGroup(['A',4]).cartan_matrix()
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [ 0  0 -1  2]

            sage: ReflectionGroup(['H',4]).cartan_matrix()
            [              2 E(5)^2 + E(5)^3               0               0]
            [E(5)^2 + E(5)^3               2              -1               0]
            [              0              -1               2              -1]
            [              0               0              -1               2]

            sage: ReflectionGroup(4).cartan_matrix()
            [-2*E(3) - E(3)^2           E(3)^2]
            [         -E(3)^2 -2*E(3) - E(3)^2]

            sage: ReflectionGroup((4,2,2)).cartan_matrix()
            [       2  -2*E(4)       -2]
            [    E(4)        2 1 - E(4)]
            [      -1 1 + E(4)        2]
        """
        if self.is_crystallographic():
            from sage.combinat.root_system.cartan_matrix import CartanMatrix as CartanMat
        else:
            from sage.matrix.all import Matrix as CartanMat
        return CartanMat(self._gap_group.CartanMat().sage())

    def invariant_form(self):
        r"""
        Return the form that is invariant under the action of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((3,1,2))
            sage: W.invariant_form()
            ?
        """
        C = self.cartan_matrix()
        n = self.rank()

        if self.is_crystallographic():
            ring = QQ
        else:
            from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField
            ring = UniversalCyclotomicField()

        from sage.matrix.constructor import zero_matrix
        form = zero_matrix(ring, n, n)

        # roots of unity of orders those of the simple reflections
        S = self.simple_reflections()
        exps = [1 - E(S[i].order()) for i in self.index_set()]

        for j in range(n):
            for i in range(j):
                if C[i,j] != 0:
                    form[j,j] = form[i,i].conjugate() * \
                                 ( C[i,j].conjugate() / C[j,i] ) * \
                                 ( exps[j] / exps[i].conjugate() )
            if form[j,j] == 0:
                form[j,j] = ring.one()
        for j in range(n):
            for i in range(j):
                form[i, j] = C[i, j] * form[i, i] / exps[j]
                form[j, i] = form[i, j].conjugate()

        form.set_immutable()
        return form

    def invariant_form_brute_force(self):
        r"""
        Return the form that is invariant under the action of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((3,1,2))
            sage: W.invariant_form_brute_force()
            ?
        """
        from sage.misc.cachefunc import cached_function

        Phi = self.roots()

        base_change = self.base_change_matrix()
        Delta = [self.simple_root(i) for i in self.index_set()]
        basis_is_Delta = base_change.is_one()
        if not basis_is_Delta:
            Delta = [beta * base_change for beta in Delta]

        S = self.simple_reflections()
        n = len(S)

        def act_on_root(w,beta):
            if basis_is_Delta:
                return w.act_on_root(beta)
            else:
                return beta * w.to_matrix()

        @cached_function
        def invariant_value(i,j):
            if i > j:
                return invariant_value(j,i).conjugate()
            val = sum((act_on_root(w,Delta[i])) * (act_on_root(w,Delta[j])).conjugate() for w in self)
            if val in QQ:
                val = QQ(val)
            return val

        coeffs = []
        for i in self.index_set():
            coeff = 1-E(S[i].order())
            if coeff in QQ:
                coeff = QQ(coeff)
            coeffs.append(coeff)

        return Matrix([ [ invariant_value(i,j)/self.cardinality() for j in range(n) ] for i in range(n) ])

    def set_reflection_representation(self,refl_repr=None):
        r"""
        Set the reflection representation of ``self``.

        INPUT:

        - ``refl_repr`` -- a dictionary representing the matrices of the
          generators of ``self`` with keys given by the index set, or
          ``None`` to reset to the default reflection representation

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: for w in W: w.to_matrix(); print("-----")
            [1 0]
            [0 1]
            -----
            [ 1  1]
            [ 0 -1]
            -----
            [-1  0]
            [ 1  1]
            -----
            [-1 -1]
            [ 1  0]
            -----
            [ 0  1]
            [-1 -1]
            -----
            [ 0 -1]
            [-1  0]
            -----

            sage: W.set_reflection_representation({1: matrix([[0,1,0],[1,0,0],[0,0,1]]), 2: matrix([[1,0,0],[0,0,1],[0,1,0]])})
            sage: for w in W: w.to_matrix(); print("-----")
            [1 0 0]
            [0 1 0]
            [0 0 1]
            -----
            [1 0 0]
            [0 0 1]
            [0 1 0]
            -----
            [0 1 0]
            [1 0 0]
            [0 0 1]
            -----
            [0 0 1]
            [1 0 0]
            [0 1 0]
            -----
            [0 1 0]
            [0 0 1]
            [1 0 0]
            -----
            [0 0 1]
            [0 1 0]
            [1 0 0]
            -----
            sage: W.set_reflection_representation()
        """
        self.one().to_matrix.clear_cache()
        if refl_repr is None or set(refl_repr.keys()) == set(self.index_set()):
            self._reflection_representation = refl_repr
        else:
            raise ValueError("the reflection representation must be defined for the complete index set")

    class Element(PermutationGroupElement):

        _reduced_word = None

        def apply_simple_reflection_right(self, i):
            r"""
            Return the product of ``self`` with the ``i``-th simple
            reflection.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: for w in W: w.apply_simple_reflection_right(1)
                (1,4)(2,3)(5,6)
                (1,2,6)(3,4,5)
                ()
                (1,5)(2,4)(3,6)
                (1,3)(2,5)(4,6)
                (1,6,2)(3,5,4)
            """
            gen = self.parent().gens()[self.parent()._index_set_inverse[i]]
            return self * gen

        def apply_simple_reflection_left(self, i):
            r"""
            Return the product of the ``i``-th simple reflection with
            ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: for w in W: w.apply_simple_reflection_left(1)
                (1,4)(2,3)(5,6)
                (1,6,2)(3,5,4)
                ()
                (1,3)(2,5)(4,6)
                (1,5)(2,4)(3,6)
                (1,2,6)(3,4,5)
            """
            gen = self.parent().gens()[self.parent()._index_set_inverse[i]]
            return gen * self

        @cached_in_parent_method
        def conjugacy_class_representative(self):
            r"""
            Return a representative of the conjugacy class of ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: for w in W:
                ....:     print('%s %s'%(w.reduced_word(), w.conjugacy_class_representative().reduced_word()))
                [] []
                [2] [1]
                [1] [1]
                [1, 2] [1, 2]
                [2, 1] [1, 2]
                [1, 2, 1] [1]
            """
            W = self.parent()
            for w in W._conjugacy_classes.keys():
                if self in W._conjugacy_classes[w]:
                    return w
            return W.conjugacy_classes_representatives()[ gap3("PositionClass(%s,%s)"%(W._gap_group._name,self)).sage()-1 ]

        def conjugacy_class(self):
            r"""
            Return the conjugacy class of ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: for w in W: sorted(w.conjugacy_class())
                [()]
                [(1,3)(2,5)(4,6), (1,4)(2,3)(5,6), (1,5)(2,4)(3,6)]
                [(1,3)(2,5)(4,6), (1,4)(2,3)(5,6), (1,5)(2,4)(3,6)]
                [(1,2,6)(3,4,5), (1,6,2)(3,5,4)]
                [(1,2,6)(3,4,5), (1,6,2)(3,5,4)]
                [(1,3)(2,5)(4,6), (1,4)(2,3)(5,6), (1,5)(2,4)(3,6)]
            """
            W = self.parent()
            if self not in W.conjugacy_classes_representatives():
                self = self.conjugacy_class_representative()
            if self in W._conjugacy_classes.keys():
                return W._conjugacy_classes[self]
            gens = W.simple_reflections()
            count = 0
            orbit = [self]
            orbit_set = set(orbit)
            while count < len(orbit):
                w = orbit[count]
                count += 1
                for s in gens:
                    w_new = s*w*s**-1
                    if w_new not in orbit_set:
                        orbit.append(w_new)
                        orbit_set.add(w_new)
            orbit_set = frozenset(orbit_set)
            W._conjugacy_classes[self] = orbit_set
            return orbit_set

        def reduced_word(self):
            r"""
            Return a word in the simple reflections to obtain ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup((5,1,1), index_set=['a'], hyperplane_index_set=['x'], reflection_index_set=['A','B','C','D'])
                sage: [w.reduced_word() for w in W]
                [[], ['a'], ['a', 'a'], ['a', 'a', 'a'], ['a', 'a', 'a', 'a']]

            .. SEEALSO:: :meth:`reduced_word_in_reflections`
            """
            I = self.parent()._index_set
            return [I[i] for i in self._reduced_word]

        @lazy_attribute
        def _reduced_word(self):
            r"""
            Computes a reduced word and stores it into ``self._reduced_word``.

            TESTS::

                sage: W = ReflectionGroup((5,1,1))
                sage: w = W.an_element()
                sage: w._reduced_word
                [0]
            """
            W = self.parent()
            gens = [W.simple_reflection(j) for j in W._index_set]
            return _gap_factorization(self, gens)

        @cached_in_parent_method
        def reduced_word_in_reflections(self):
            r"""
            Return a word in the reflections to obtain ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup((5,1,1), index_set=['a'], reflection_index_set=['A','B','C','D'])
                sage: [w.reduced_word_in_reflections() for w in W]
                [[], ['A'], ['B'], ['C'], ['D']]

                sage: W = ReflectionGroup(['A',2], index_set=['a','b'], reflection_index_set=['A','B','C'])
                sage: [(w.reduced_word(), w.reduced_word_in_reflections()) for w in W]
                [([], []),
                 (['b'], ['B']),
                 (['a'], ['A']),
                 (['a', 'b'], ['A', 'B']),
                 (['b', 'a'], ['A', 'C']),
                 (['a', 'b', 'a'], ['C'])]

            .. SEEALSO:: :meth:`reduced_word`
            """
            if self.is_one():
                return []

            # TODO: Move to the real implementation
            W = self.parent()
            if W.is_real():
                r = self.reflection_length()
                R = W.reflections()
                I = W.reflection_index_set()
                word = []
                while r > 0:
                    for i in I:
                        w = R[i]._mul_(self)
                        if w.reflection_length() < r:
                            word += [i]
                            r -= 1
                            self = w
                            break
                return word
            else:
                gens = [W.reflection(j) for j in W._reflection_index_set]
                word = _gap_factorization(self, gens)
                return [self.parent()._reflection_index_set[i] for i in word]

        def length(self):
            r"""
            Return the length of ``self`` in generating reflections.

            This is the minimal numbers of generating reflections needed
            to obtain ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup(4)
                sage: for w in W:
                ....:   print('%s %s'%(w.reduced_word(), w.length()))
                [] 0
                [1] 1
                [2] 1
                [1, 1] 2
                [1, 2] 2
                [2, 1] 2
                [2, 2] 2
                [1, 1, 2] 3
                [1, 2, 1] 3
                [1, 2, 2] 3
                [2, 1, 1] 3
                [2, 2, 1] 3
                [1, 1, 2, 1] 4
                [1, 1, 2, 2] 4
                [1, 2, 1, 1] 4
                [1, 2, 2, 1] 4
                [2, 1, 1, 2] 4
                [2, 2, 1, 1] 4
                [1, 1, 2, 1, 1] 5
                [1, 1, 2, 2, 1] 5
                [2, 1, 1, 2, 1] 5
                [1, 2, 2, 1, 1] 5
                [2, 1, 1, 2, 1, 1] 6
                [1, 1, 2, 2, 1, 1] 6
            """
            return len(self.reduced_word())

        @cached_in_parent_method
        def reflection_length(self, in_unitary_group=False):
            r"""
            Return the reflection length of ``self``.

            This is the minimal numbers of reflections needed to obtain
            ``self``.

            INPUT:

            - ``in_unitary_group`` -- (default:``False``) if ``True``,
              the reflection length is computed in the unitary group
              which is the dimension of the move space of ``self``

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: sorted([t.reflection_length() for t in W])
                [0, 1, 1, 1, 2, 2]

                sage: W = ReflectionGroup((2,1,2))
                sage: sorted([t.reflection_length() for t in W])
                [0, 1, 1, 1, 1, 2, 2, 2]

                sage: W = ReflectionGroup((2,2,2))
                sage: sorted([t.reflection_length() for t in W])
                [0, 1, 1, 2]

                sage: W = ReflectionGroup((3,1,2))
                sage: sorted([t.reflection_length() for t in W])
                [0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
            """
            W = self.parent()
            if self in W.conjugacy_classes_representatives():
                if in_unitary_group or W.is_real():
                    return W.rank() - self.reflection_eigenvalues(test_class_repr=False).count(0)
                else:
                    return len(self.reduced_word_in_reflections())
            else:
                w = self.conjugacy_class_representative()
                # the following assert a possible implementation bug and
                # is hopefully never needed
                assert w in self.parent().conjugacy_classes_representatives()
                return w.reflection_length(in_unitary_group=in_unitary_group)

        @cached_in_parent_method
        def to_matrix(self):
            r"""
            Return ``self`` as a matrix acting on the underlying vector
            space.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: for w in W:
                ....:     print(w.reduced_word())
                ....:     print(w.to_matrix())
                []
                [1 0]
                [0 1]
                [2]
                [ 1  1]
                [ 0 -1]
                [1]
                [-1  0]
                [ 1  1]
                [1, 2]
                [-1 -1]
                [ 1  0]
                [2, 1]
                [ 0  1]
                [-1 -1]
                [1, 2, 1]
                [ 0 -1]
                [-1  0]
            """
            W = self.parent()
            if W._reflection_representation is None:
                Delta = W.independent_roots()
                Phi = W.roots()
                M = Matrix([Phi[self(Phi.index(alpha)+1)-1] for alpha in Delta])
                return W.base_change_matrix() * M
            else:
                refl_repr = W._reflection_representation
                id_mat = identity_matrix(QQ,refl_repr[W.index_set()[0]].nrows())
                return prod([refl_repr[i] for i in self.reduced_word()], id_mat)

        @cached_in_parent_method
        def fix_space(self):
            r"""
            Return the fix space of ``self``.

            This is the sub vector space of the underlying vector space
            on which ``self`` acts trivially.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: for w in W:
                ....:     w.reduced_word()
                ....:     w.fix_space()
                []
                Vector space of degree 2 and dimension 2 over Rational Field
                Basis matrix:
                [1 0]
                [0 1]
                [2]
                Vector space of degree 2 and dimension 1 over Rational Field
                Basis matrix:
                [1 0]
                [1]
                Vector space of degree 2 and dimension 1 over Rational Field
                Basis matrix:
                [0 1]
                [1, 2]
                Vector space of degree 2 and dimension 0 over Rational Field
                Basis matrix:
                []
                [2, 1]
                Vector space of degree 2 and dimension 0 over Rational Field
                Basis matrix:
                []
                [1, 2, 1]
                Vector space of degree 2 and dimension 1 over Rational Field
                Basis matrix:
                [ 1 -1]

                sage: W = ReflectionGroup(23)
                sage: W.an_element().fix_space()
                Vector space of degree 3 and dimension 2 over Universal Cyclotomic Field
                Basis matrix:
                [0 1 0]
                [0 0 1]
            """
            I = identity_matrix(QQ, self.parent().rank())
            return (self.to_matrix() - I).right_kernel()

        @cached_in_parent_method
        def reflection_eigenvalues(self, test_class_repr=True):
            r"""
            Return the reflection eigenvalues of ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup(4)
                sage: for w in W: w.reflection_eigenvalues()
                [0, 0]
                [1/3, 0]
                [1/3, 0]
                [2/3, 0]
                [1/6, 1/2]
                [1/6, 1/2]
                [2/3, 0]
                [1/4, 3/4]
                [1/4, 3/4]
                [1/4, 3/4]
                [1/4, 3/4]
                [1/4, 3/4]
                [1/3, 0]
                [1/2, 5/6]
                [1/3, 0]
                [1/2, 5/6]
                [1/2, 5/6]
                [1/2, 5/6]
                [1/6, 1/2]
                [2/3, 0]
                [1/6, 1/2]
                [2/3, 0]
                [1/2, 1/2]
                [1/4, 3/4]
            """
            return self.parent().reflection_eigenvalues(self, test_class_repr=test_class_repr)

        @cached_in_parent_method
        def galois_conjugates(self):
            r"""
            Return all Galois conjugates of ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup(4)
                sage: for w in W: print(w.galois_conjugates())
                [[1 0]
                 [0 1]]
                [[   1    0]
                 [   0 E(3)], [     1      0]
                 [     0 E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
                 [ 4/3*E(3) + 2/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2], [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
                 [ 2/3*E(3) + 4/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]]
                [[     1      0]
                 [     0 E(3)^2], [   1    0]
                 [   0 E(3)]]
                [[ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2], [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
                 [ 2/3*E(3) - 2/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
                 [ 4/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2], [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) + 4/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
                 [ 2/3*E(3) + 4/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2], [ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
                 [ 4/3*E(3) + 2/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) - 4/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2], [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
                 [-4/3*E(3) - 2/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) + 2/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2], [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) - 2/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
                 [-4/3*E(3) - 2/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2], [ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) - 4/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
                 [ 4/3*E(3) + 2/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2], [-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) + 4/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) + 4/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2], [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
                 [ 4/3*E(3) + 2/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) - 4/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2], [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
                 [-4/3*E(3) - 2/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
                 [ 2/3*E(3) - 2/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2], [ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
                 [-2/3*E(3) + 2/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2], [-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) - 2/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
                 [-4/3*E(3) - 2/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2], [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
                 [-2/3*E(3) - 4/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]]
                [[   -1     0]
                 [    0 -E(3)], [     -1       0]
                 [      0 -E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) + 4/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2], [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
                 [ 4/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
                 [-2/3*E(3) - 4/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2], [-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
                 [-4/3*E(3) - 2/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) - 2/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2], [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
                 [-2/3*E(3) + 2/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]]
                [[     -1       0]
                 [      0 -E(3)^2], [   -1     0]
                 [    0 -E(3)]]
                [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
                 [-4/3*E(3) - 2/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2], [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) - 4/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]]
                [[-1  0]
                 [ 0 -1]]
                [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) - 2/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2], [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) + 2/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]]
            """
            rk = self.parent().rank()
            M = self.to_matrix().list()
            m = lcm([x.conductor() if hasattr(x,"conductor") else 1 for x in M])
            M_gals = [x.galois_conjugates(m) if hasattr(x,"galois_conjugates") else [x] for x in M]
            conjugates = []
            for i in range(len(M_gals[0])):
                conjugates.append(Matrix(rk, [X[i] for X in M_gals]))
            return conjugates

        def __cmp__(self, other):
            r"""
            Compare ``self`` with ``other``.

            Without this comparison method, the initialization of this
            permutation group fails ...

            TESTS::

                sage: W = ReflectionGroup((1,1,3))
                sage: a,b = W.gens()
                sage: a.__cmp__(b)
                1
                sage: b.__cmp__(a)
                -1
            """
            return super(ComplexReflectionGroup.Element, self).__cmp__(other)

class IrreducibleComplexReflectionGroup(ComplexReflectionGroup):

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3)); W
            Irreducible real reflection group of rank 2 and type A2
            sage: W = ReflectionGroup((3,1,4)); W
            Irreducible complex reflection group of rank 4 and type G(3,1,4)
        """
        type_str = self._irrcomp_repr_(self._type[0])
        return 'Irreducible complex reflection group of rank %s and type %s'%(self._rank,type_str)

    @cached_method
    def a_coxeter_element(self):
        r"""
        Return a Coxeter element of a well-generated, irreducible
        reflection group.

        This is an element having a regular eigenvector (a vector not
        contained in any reflecting hyperplane of ``self``).

        .. WARNING::

            ``self`` is assumed to be well-generated.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,4))
            sage: W.a_coxeter_element().reduced_word()
            [1, 2, 3]

            sage: W = ReflectionGroup((2,1,4))
            sage: W.a_coxeter_element().reduced_word()
            [1, 2, 3, 4]

            sage: W = ReflectionGroup((4,1,4))
            sage: W.a_coxeter_element().reduced_word()
            [1, 2, 3, 4]

            sage: W = ReflectionGroup((4,4,4))
            sage: W.a_coxeter_element().reduced_word()
            [1, 2, 3, 4]
        """
        if not self.is_irreducible() or not self.is_well_generated():
            raise ValueError("this method is available for irreducible, well-generated complex reflection groups")
        return self.from_word(self._index_set)

    @cached_method
    def coxeter_elements(self):
        r"""
        Return the (unique) conjugacy class in ``self`` containing all
        Coxeter elements.

        .. NOTE::

            - ``self`` is assumed to be well-generated.
            - This works even beyond real reflection groups, but the
              conjugacy class is not unique and we only obtain one such class.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))
            sage: sorted(c.reduced_word() for c in W.coxeter_elements())
            [[1, 2], [2, 1]]

            sage: W = ReflectionGroup((1,1,4))
            sage: sorted(c.reduced_word() for c in W.coxeter_elements())
            [[1, 2, 3], [1, 2, 3, 1, 2], [2, 3, 1],
             [2, 3, 1, 2, 1], [3, 1, 2], [3, 2, 1]]
        """
        return self.a_coxeter_element().conjugacy_class()

    @cached_method
    def standard_coxeter_elements(self):
        r"""
        Return all standard Coxeter elements in ``self``.

        This is the set of all elements in self obtained from any
        product of the simple reflections in ``self``.

        .. NOTE::

            - ``self`` is assumed to be well-generated.
            - This works even beyond real reflection groups, but the conjugacy
              class is not unique and we only obtain one such class.

        EXAMPLES::

            sage: W = ReflectionGroup(4)
            sage: sorted(W.standard_coxeter_elements())
            [(1,7,6,12,23,20)(2,8,17,24,9,5)(3,16,10,19,15,21)(4,14,11,22,18,13),
             (1,10,4,12,21,22)(2,11,19,24,13,3)(5,15,7,17,16,23)(6,18,8,20,14,9)]
        """
        if not self.is_irreducible() or not self.is_well_generated():
            raise ValueError("this method is available for irreducible, well-generated complex reflection groups")
        from sage.combinat.permutation import Permutations
        return set(self.from_word(w) for w in Permutations(self._index_set))

    def elements_below_coxeter_element(self, c=None):
        r"""
        Return all elements in ``self`` in the interval `[1,c]` in the
        absolute order of ``self``.

        This order is defined by

        .. MATH::

            \omega \leq_R \tau \Leftrightarrow \ell_R(\omega) +
            \ell_R(\omega^{-1} \tau) = \ell_R(\tau),

        where `\ell_R` denotes the reflection length.

        .. NOTE::

            ``self`` is assumed to be well-generated.

        INPUT:

        - ``c`` -- (default: ``None``) if an element ``c`` is given, it
          is used as the maximal element in the interval; if a list is
          given, the union of the various maximal elements is computed

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))

            sage: sorted( w.reduced_word() for w in W.elements_below_coxeter_element() )
            [[], [1], [1, 2], [1, 2, 1], [2]]

            sage: sorted( w.reduced_word() for w in W.elements_below_coxeter_element(W.from_word([2,1])) )
            [[], [1], [1, 2, 1], [2], [2, 1]]

            sage: sorted( w.reduced_word() for w in W.elements_below_coxeter_element(W.from_word([2])) )
            [[], [2]]
        """
        if c in self:
            cs = [c]
        elif c is None:
            cs = [self.a_coxeter_element()]
        else:
            cs = list(c)
        l = cs[0].reflection_length(in_unitary_group=True)
        f = lambda pi: any(pi.reflection_length(in_unitary_group=True)
                           + (c*pi**-1).reflection_length(in_unitary_group=True) == l
                           for c in cs)
        # first computing the conjugacy classes only needed if the interaction with gap3 is slow due to a bug
        #self.conjugacy_classes()
        return filter(f, self)

    @cached_method
    def noncrossing_partition_lattice(self, c=None, L=None, in_unitary_group=False):
        r"""
        Return the interval `[1,c]` in the absolute order of
        ``self`` as a finite lattice.

        .. SEEALSO:: :meth:`elements_below_coxeter_element`

        INPUT:

        - ``c`` -- (default: ``None``) if an element ``c`` in ``self`` is
          given, it is used as the maximal element in the interval

        - ``L`` -- (default: ``None``) if a subset ``L`` (must be hashable!)
          of ``self`` is given, it is used as the underlying set (only
          cover relations are checked)

        - ``in_unitary_group`` -- (default: ``False``) if ``False``, the
          relation is given by `\sigma \leq \tau` if
          `l_R(\sigma) + l_R(\sigma^{-1}\tau) = l_R(\tau)`;
          if ``True``, the relation is given by `\sigma \leq \tau` if
          `\dim(\mathrm{Fix}(\sigma)) + \dim(\mathrm{Fix}(\sigma^{-1}\tau))
          = \dim(\mathrm{Fix}(\tau))`

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))

            sage: sorted( w.reduced_word() for w in W.noncrossing_partition_lattice() )
            [[], [1], [1, 2], [1, 2, 1], [2]]

            sage: sorted( w.reduced_word() for w in W.noncrossing_partition_lattice(W.from_word([2,1])) )
            [[], [1], [1, 2, 1], [2], [2, 1]]

            sage: sorted( w.reduced_word() for w in W.noncrossing_partition_lattice(W.from_word([2])) )
            [[], [2]]
        """
        from sage.combinat.posets.all import Poset, LatticePoset
        if c is None:
            c = self.a_coxeter_element()

        smart_covers = not in_unitary_group

        if self.is_real():
            smart_covers = in_unitary_group = True

        R = self.reflections()
        if L is None:
            L = self.elements_below_coxeter_element(c=c)
            if c.is_coxeter_element():
                smart_covers = in_unitary_group = True
        rels = []
        ref_lens = {w: w.reflection_length(in_unitary_group=in_unitary_group)
                    for w in L}
        if smart_covers:
            for pi in L:
                for t in R:
                    tau = pi*t
                    if tau in L and ref_lens[pi] + 1 == ref_lens[tau]:
                        rels.append((pi,tau))
        else:
            rels = [(pi,tau) for pi in L for tau in L
                    if ref_lens[pi] + ref_lens[pi.inverse()*tau] == ref_lens[tau]]
        P = Poset((L,rels), cover_relations=smart_covers, facade=True)
        if P.is_lattice():
            return LatticePoset(P)
        else:
            return P

    def generalized_noncrossing_partitions(self, m, c=None, positive=False):
        r"""
        Return the set of all chains of length ``m`` in the noncrossing
        partition lattice of ``self``, see
        :meth:`noncrossing_partition_lattice`.

        .. NOTE::

            ``self`` is assumed to be well-generated

        INPUT:

        - ``c`` -- (default: ``None``) if an element ``c`` in ``self``
          is given, it is used as the maximal element in the interval

        - ``positive`` -- (default: ``False``) if ``True``, only those
          generalized noncrossing partitions of full support are returned

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))

            sage: sorted([w.reduced_word() for w in chain]
            ....:        for chain in W.generalized_noncrossing_partitions(2))
            [[[], [], [1, 2]],
             [[], [1], [2]],
             [[], [1, 2], []],
             [[], [1, 2, 1], [1]],
             [[], [2], [1, 2, 1]],
             [[1], [], [2]],
             [[1], [2], []],
             [[1, 2], [], []],
             [[1, 2, 1], [], [1]],
             [[1, 2, 1], [1], []],
             [[2], [], [1, 2, 1]],
             [[2], [1, 2, 1], []]]

            sage: sorted([w.reduced_word() for w in chain]
            ....:        for chain in W.generalized_noncrossing_partitions(2, positive=True))
            [[[], [1, 2], []],
             [[], [1, 2, 1], [1]],
             [[1], [2], []],
             [[1, 2], [], []],
             [[1, 2, 1], [], [1]],
             [[1, 2, 1], [1], []],
             [[2], [1, 2, 1], []]]
        """
        from sage.combinat.combination import Combinations
        NC = self.noncrossing_partition_lattice(c=c)
        one = self.one()
        if c is None:
            c = self.a_coxeter_element()
        chains = NC.chains()
        NCm = set()
        iter = chains.breadth_first_search_iterator()
        chain = next(iter)
        chain = next(iter)
        while len(chain) <= m:
            chain.append( c )
            for i in range(len(chain)-1, 0, -1):
                chain[i] = chain[i-1]**-1 * chain[i]
            k = m + 1 - len(chain)
            for positions in Combinations(range(m+1),k):
                ncm = []
                for l in range(m+1):
                    if l in positions:
                        ncm.append(one)
                    else:
                        l_prime = l - len([i for i in positions if i <= l])
                        ncm.append(chain[l_prime])
                if not positive or prod(ncm[:-1]).has_full_support():
                    NCm.add(tuple(ncm))
            try:
                chain = next(iter)
            except StopIteration:
                chain = range(m + 1)
        return NCm

    @cached_method
    def absolute_poset(self, in_unitary_group=False):
        r"""
        Return the poset induced by the absolute order of ``self`` as a
        finite lattice.

        INPUT:

        - ``in_unitary_group`` -- (default: ``False``) if ``False``, the
          relation is given by ``\sigma \leq \tau`` if ``l_R(\sigma) + l_R(\sigma^{-1}\tau) = l_R(\tau)``.
          If ``True``, the relation is given by ``\sigma \leq \tau`` if
          ``dim(Fix(\sigma)) + dim(Fix(\sigma^{-1}\tau)) = dim(Fix(\tau))``.

        .. SEEALSO:: :meth:`noncrossing_partition_lattice`

        EXAMPLES::

            sage: P = ReflectionGroup((1,1,3)).absolute_poset(); P
            Finite poset containing 6 elements

            sage: sorted(w.reduced_word() for w in P)
            [[], [1], [1, 2], [1, 2, 1], [2], [2, 1]]

            sage: W = ReflectionGroup(4); W
            Irreducible complex reflection group of rank 2 and type ST4
            sage: W.absolute_poset()
            Finite poset containing 24 elements
        """
        return self.noncrossing_partition_lattice(L=self)

    class Element(ComplexReflectionGroup.Element):

        @cached_in_parent_method
        def is_coxeter_element(self, which_primitive=1, test_class_repr=True):
            r"""
            Return ``True`` if ``self`` is a Coxeter element.

            .. SEEALSO:: :meth:`~IrreducibleComplexReflectionGroup.a_coxeter_element`

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: for w in W:
                ....:     print('%s %s'%(w.reduced_word(), w.is_coxeter_element()))
                [] False
                [2] False
                [1] False
                [1, 2] True
                [2, 1] True
                [1, 2, 1] False
            """
            if not self.parent().is_irreducible() or not self.parent().is_well_generated():
                raise ValueError("this method is available for elements in irreducible, well-generated complex reflection groups")
            h = self.parent().coxeter_number()
            # to check regularity for a Coxeter number h, we get that an eigenvector is regular for free
            return any(QQ(ev).denom() == h and QQ(ev).numer() == which_primitive
                       for ev in self.reflection_eigenvalues(test_class_repr=test_class_repr))

        @cached_in_parent_method
        def is_h_regular(self, test_class_repr=True):
            r"""
            Return ``True`` if ``self`` is regular.

            This is if ``self`` has an eigenvector with eigenvalue `h`
            and which does not lie in any reflecting hyperplane.
            Here, `h` denotes the Coxeter number.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: for w in W:
                ....:     print('%s %s'%(w.reduced_word(), w.is_h_regular()))
                [] False
                [2] False
                [1] False
                [1, 2] True
                [2, 1] True
                [1, 2, 1] False
            """
            if not self.parent().is_irreducible() or not self.parent().is_well_generated():
                raise ValueError("This method is available for elements in irreducible, well-generated complex reflection groups")
            h = self.parent().coxeter_number()
            # to check regularity for a Coxeter number h, we get that an eigenvector is regular for free
            return any(QQ(ev).denom() == h
                       for ev in self.reflection_eigenvalues(test_class_repr=test_class_repr))

        @cached_in_parent_method
        def is_regular(self, h, test_class_repr=True):
            r"""
            Return ``True`` if ``self`` is regular.

            This is, if ``self`` has an eigenvector with eigenvalue
            ``h`` and which does not lie in any reflecting hyperplane.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: for w in W:
                ....:     print('%s %s'%(w.reduced_word(), w.is_regular(W.coxeter_number())))
                [] False
                [2] False
                [1] False
                [1, 2] True
                [2, 1] True
                [1, 2, 1] False

                sage: W = ReflectionGroup(23)
                sage: for w in W:
                ....:     if w.is_regular(W.coxeter_number()):
                ....:         w.reduced_word()
                [1, 2, 3]
                [2, 3, 1]
                [3, 1, 2]
                [3, 2, 1]
                [2, 1, 2, 3, 2]
                [2, 3, 2, 1, 2]
                [1, 2, 1, 2, 3, 2, 1]
                [1, 2, 3, 2, 1, 2, 1]
                [1, 2, 3, 1, 2, 1, 2, 3, 2]
                [2, 1, 2, 3, 1, 2, 1, 2, 3]
                [2, 3, 1, 2, 1, 2, 3, 2, 1]
                [3, 1, 2, 1, 2, 3, 2, 1, 2]
                [3, 2, 1, 2, 3, 1, 2, 1, 2]
                [1, 2, 1, 2, 3, 1, 2, 1, 2]
                [3, 2, 1, 2, 3, 2, 1, 2, 1]
                [2, 1, 2, 3, 1, 2, 1, 2, 1]
                [3, 2, 1, 2, 3, 1, 2, 1, 2, 3, 1]
                [3, 2, 1, 2, 3, 1, 2, 1, 2, 3, 2]
                [2, 1, 2, 3, 1, 2, 1, 2, 3, 2, 1]
                [1, 2, 3, 1, 2, 1, 2, 3, 2, 1, 2]
                [3, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1]
                [2, 3, 1, 2, 1, 2, 3, 2, 1, 2, 1]
                [3, 2, 1, 2, 3, 1, 2, 1, 2, 3, 1, 2, 1]
                [3, 2, 1, 2, 3, 1, 2, 1, 2, 3, 2, 1, 2]
            """
            evs = self.reflection_eigenvalues(test_class_repr=test_class_repr)
            I = identity_matrix(self.parent().rank())
            for ev in evs:
                ev = QQ(ev)
                if h == ev.denom():
                    M = self.to_matrix() - E(ev.denom(),ev.numer()) * I
                    V = M.right_kernel()
                    if all(not V.is_subspace(H) for H in self.parent().reflecting_hyperplanes()):
                        return True
            return False

def _gap_factorization(w, gens):
    r"""
    Return a factorization of ``w`` using the generators ``gens``.

    .. WARNING::

        This is only available through GAP3 and chevie.

    EXAMPLES::

        sage: from sage.combinat.root_system.reflection_group_complex import _gap_factorization
        sage: W = ReflectionGroup((1,1,3))
        sage: gens = [ W.simple_reflection(i) for i in W.index_set() ]
        sage: for w in W: _gap_factorization(w,gens)
        []
        [1]
        [0]
        [0, 1]
        [1, 0]
        [0, 1, 0]
    """
    gap3.execute('W := GroupWithGenerators(%s)'%str(gens))
    gap3.execute(_gap_factorization_code)
    fac = gap3('MinimalWord(W,%s)'%str(w)).sage()
    return [i-1 for i in fac]

_gap_factorization_code = """
# MinimalWord(G,w)
# given a permutation group G find some expression of minimal length in the
# generators of G and their inverses of the element w (an inverse is
# representated by a negative index).
# To speed up  later calls to  the same function  the fields G.base, G.words,
# G.nbwordslength are kept.
MinimalWord:=function(G,w)
  local decode,i,p,g,h,n,bag,nbe,nbf,new,gens,inds;
  # to save space elements of G are represented as image of the base, and
  # words are represented as: index of previous elt, last generator applied;
  if not IsBound(G.base) then
    StabChain(G);g:=G; G.base:=[];
    while IsBound(g.orbit) do Add(G.base,g.orbit[1]); g:=g.stabilizer; od;
  fi;
  w:=OnTuples(G.base,w);
  if not IsBound(G.words) then
    G.words:=[G.base]; G.lastmult:=[[0,0]];
    G.nbwordslength:=[1];
  fi;
  gens:=ShallowCopy(G.generators);inds:=[1..Length(gens)];
  #  for g in G.generators do
  #    if g<>g^-1 then Add(gens,g^-1);Add(inds,-Position(gens,g));fi;
  #  od;
  bag:=Set(G.words);
  nbe:=0;nbf:=0;
  decode:=function(i)local w;w:=[];
    while i<>1 do Add(w,G.lastmult[i][2]); i:=G.lastmult[i][1];od;
    return Reversed(w);
  end;
  while true do
    if w in bag then return decode(Position(G.words,w));fi;
    new:=Length(G.words);
    for g in [1..Length(gens)] do
      for h in [1+Sum(G.nbwordslength{[1..Length(G.nbwordslength)-1]})..new] do
         n:=OnTuples(G.words[h],gens[g]);
         if n in bag then
           nbe:=nbe+1;# if nbe mod 500=1 then Print(".\c");fi;
         else
           nbf:=nbf+1;# if nbf mod 500=1 then Print("*\c");fi;
       Add(G.words,n);Add(G.lastmult,[h,inds[g]]);AddSet(bag,n);
         fi;
       od;
    od;
    Add(G.nbwordslength,Length(G.words)-new);
    Print("\n",G.nbwordslength[Length(G.nbwordslength)]," elements of length ",
      Length(G.nbwordslength)-1);
  od;
end;"""

def _gap_return(S, coerce_obj='self'):
    r"""
    Return the string ``S`` after a few modifications are done.

    This is a stupid internal function to take gap output as a string,
    replace a few things, to then turn it into a sage object.

    TESTS::

        sage: from sage.combinat.root_system.reflection_group_complex import _gap_return
        sage: _gap_return("[ (), (1,4)(2,3)(5,6), (1,6,2)(3,5,4) ]")
        "[self('()',check=False),self('(1,4)(2,3)(5,6)',check=False),self('(1,6,2)(3,5,4)',check=False)]"
    """
    S = S.replace(' ','').replace('\n','')
    S = S.replace(',(','\',check=False),%s(\'('%coerce_obj).replace('[','[%s(\''%coerce_obj).replace(']','\',check=False)]')
    return S
