r"""
Lie Conformal Algebras
AUTHORS

- Reimundo Heluani (08-09-2019): Initial implementation

.. include:: ../../../vertex_algebras/lie_conformal_algebra_desc.rst

EXAMPLES:

The base class for all Lie conformal algebras is :class:`LieConformalAlgebra`. 
All subclasses are called through its method ``__classcall_private__``. 
We provide some convenience functions to define named Lie conformal algebras
like :meth:`VirasoroLieConformalAlgebra` and :meth:`AffineLieConformalAlgebra`

- We construct the Virasoro Lie conformal algebra, its universal enveloping
  vertex algebra and lift some elements::

        sage: Vir = VirasoroLieConformalAlgebra(QQ)
        sage: Vir.inject_variables()
        Defining L, C
        sage: L.bracket(L)
        {0: TL, 1: 2*L, 3: 1/2*C}
        sage: cp = Family({C:1/2}) 
        sage: V = Vir.universal_enveloping_algebra(cp)
        sage: V
        The universal enveloping vertex algebra of Lie conformal algebra with generators (L, C) over Rational Field.
        sage: L.lift()
        L_-2|0>
        sage: L*L
        L_-2L_-2|0>
        sage: sorted(L.bracket(L*L).items())
        [(0, 2*L_-3L_-2|0>+L_-5|0>),
         (1, 4*L_-2L_-2|0>),
         (2, 3*L_-3|0>),
         (3, 17/2*L_-2|0>),
         (5, 3/2*|0>)]

- We construct the Current algebra for `\mathfrak{sl}_2`::

        sage: V = AffineLieConformalAlgebra(QQ, 'A1')
        sage: V.gens()
        (alpha[1], alphacheck[1], -alpha[1], K)
        sage: e = V.0; f = V.2; e.nproduct(f,1)
        K

- We construct the `\beta-\gamma` system by directly giving the
  `\lambda`-brackets of the generators::

        sage: betagamma_dict = {('b','a'):{0:{('K',0):1}}}
        sage: V = LieConformalAlgebra(QQ, betagamma_dict, names=('a','b'), weights=(1,0), central_elements=('K',))
        sage: V.category()
        Category of finitely generated H-graded Lie conformal algebras with basis over Rational Field
        sage: V.inject_variables()
        Defining a, b, K
        sage: a.bracket(b)
        {0: -K}

"""


#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.lie_conformal_algebras import LieConformalAlgebras
from sage.algebras.vertex_algebras.lie_conformal_algebra_element import \
        LCAStructureCoefficientsElement
from sage.sets.family import Family 
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.algebras.lie_algebras.lie_algebra import LieAlgebra
from sage.categories.cartesian_product import cartesian_product
from sage.combinat.free_module import CombinatorialFreeModule
from sage.functions.other import binomial
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.categories.commutative_rings import CommutativeRings
from sage.structure.indexed_generators import (IndexedGenerators, 
                                               standardize_names_index_set)


class LieConformalAlgebra(Parent, UniqueRepresentation):
    @staticmethod
    def __classcall_private__(cls, R=None, arg0=None, names=None, 
        central_elements=None, index_set=None, category=None, 
        weights=None, **kwds):
        
        if not R in CommutativeRings():
            raise ValueError("First argument must be a commutative ring" + 
                             " got {}".format(R))

        if isinstance(arg0, dict):
            graded=kwds.get("graded", False)
            if weights is not None or graded:
                #if not arg0:
                    #return AbelianLieConformalAlgebra(R, names, index_set,
                    #    weights, graded=True)
                return GradedLieConformalAlgebra(R, arg0,
                     names=names, index_set=index_set, 
                    central_elements=central_elements,
                    category=category,weights=weights, **kwds)
            else:
                #if not arg0:
                #    return AbelianLieConformalAlgebra(R, names, index_set)
                return LieConformalAlgebraWithStructureCoefficients(R, arg0,
                     names=names, index_set=index_set,
                     central_elements=central_elements,
                     category=category, **kwds)

        return NotImplementedError("Not implemented")

    def __init__(self, R, names=None, category=None):
        r"""
        Lie conformal algebras base class and factory

        INPUT:

        - ``R`` -- a commutative ring (Default: None); The base ring of
          this Lie conformal algebra. Behaviour is undefined if it is
          not a field of characteristic zero. 

        - ``arg0`` -- Dictionary (Default: None); 
          a dictionary containing the `\lambda` brackets of the
          generators of this Lie conformal algebra. The keys of this 
          dictionary are pairs of either names or indices of the 
          generators and the values are themselves dictionaries. For a 
          pair of generators `a` and `b`, the value of 
          ``arg0[('a','b')]`` is a dictionary whose keys are positive 
          integer numbers and the corresponding value for the
          key `j` is a dictionary itself representing the j-th product
          `a_{(j)}b`. Thus, for a positive integer number `j`, the 
          value of ``arg0[('a','b')][j]`` is a dictionary whose entries
          are pairs ``('c',n)`` where ``'c'`` is the name of a generator
          and `n` is a positive number. The value for this key is the 
          coefficient of `\frac{T^{n}}{n!} c` in `a_{(j)}b`. For 
          example the ``arg0`` for the *Virasoro* Lie conformal algebra
          is::

                {('L','L'):{0:{('L',1):1}, 1:{('L',0):2}, 3:{('C',0):1/2}}}


          Do not include central elements in this dictionary. Also, if 
          the key ``('a','b')`` is present, there is no need to include
          ``('b','a')`` as it is defined by skew-symmetry. Any missing 
          pair (besides the ones defined by skew-symmetry) is assumed 
          to have vanishing `\lambda`-bracket. 

        - ``names`` -- tuple of ``str`` (Default: None); The list of 
          names for generators of this Lie conformal algebra. Do not 
          include central elements in this list.

        - ``central_elements`` -- tuple of ``str`` (Default: None); 
          A list of names for central elements of this Lie conformal 
          algebra. 

        - ``index_set`` -- enumerated set (Default: None); an indexing 
          set for the generators of this Lie conformal algebra. Do not 
          include central elements in this list. 

        - ``weights`` -- tuple of non-negative integers (Default: None);
          a list of degrees for this Lie conformal algebra. The 
          returned Lie conformal algebra is H-Graded. This tuple needs
          to have the same cardinality as ``index_set`` or ``names``.
          Central elements are assumed to have weight `0`. 

        - ``category`` The category that this Lie conformal algebra 
          belongs to.

        In addition we accept the following keywords:

        - ``graded`` -- if present, the returned algebra is H-Graded.
          If ``weights`` is not specified, all non-central generators 
          are assigned degree `1`. This keyword is unnecessary if 
          ``weights`` is specified

        EXAMPLES:

        We construct the `\beta-\gamma` system or *Weyl* Lie conformal 
        algebra::

            sage: betagamma_dict = {('b','a'):{0:{('K',0):1}}}
            sage: V = LieConformalAlgebra(QQbar, betagamma_dict, names=('a','b'), weights=(1,0), central_elements=('K',))
            sage: V.category()
            Category of finitely generated H-graded Lie conformal algebras with basis over Algebraic Field
            sage: V.inject_variables()
            Defining a, b, K
            sage: a.bracket(b)
            {0: -K}

        We construct the current algebra for `\mathfrak{sl}_2`::

            sage: sl2dict = {('e','f'):{0:{('h',0):1}, 1:{('K',0):1}}, ('e','h'):{0:{('e',0):-2}}, ('f','h'):{0:{('f',0):2}}, ('h', 'h'):{1:{('K', 0):2}}}
            sage: V = LieConformalAlgebra(QQ, sl2dict, names=('e', 'h', 'f'), central_elements=('K',), graded=True)
            sage: V.inject_variables()
            Defining e, h, f, K
            sage: e.bracket(f)
            {0: h, 1: K}
            sage: h.bracket(e)
            {0: 2*e}
            sage: e.bracket(f.T())
            {0: Th, 1: h, 2: 2*K}
            sage: V.category()
            Category of finitely generated H-graded Lie conformal algebras with basis over Rational Field
            sage: e.degree()
            1

        """
        category = LieConformalAlgebras(R).or_subcategory(category)
        super(LieConformalAlgebra,self).__init__(base=R, names=names, 
                                                 category=category)

    def _element_constructor_(self,x):
        """
        Construct an element of this Lie conformal algebra.

        EXAMPLES::

            sage: R = VirasoroLieConformalAlgebra(QQbar)
            sage: R(0)
            0
            sage: R(1)
            Traceback (most recent call last):
            ...
            ValueError: can only convert the scalar 0 into a Lie conformal algebra element
        """
        if x in self.base_ring():
            if not self.base_ring()(x).is_zero():
                raise ValueError("can only convert the scalar 0 into a"\
                                 " Lie conformal algebra element")
            return self.zero()
        try: 
            x = self.module()(x)
        except (TypeError, ValueError) as e:
            raise ValueError("Don't know how to convert {0} into an element"\
                             " of {1}".format(x,self))
        return self.element_class(self,x)

class LieConformalAlgebraWithBasis(LieConformalAlgebra):
    def __init__(self,R, names=None, index_set=None, 
                 category=None, prefix='B', **kwds):
        """
        Base class for a Lie conformal algebra with a preferred basis.

        EXAMPLES::
            
            sage: R = VirasoroLieConformalAlgebra(QQbar);R
            Lie conformal algebra with generators (L, C) over Algebraic Field.
        """
        self._indices = index_set
        category = LieConformalAlgebras(R).WithBasis().or_subcategory(category)
        super(LieConformalAlgebraWithBasis,self).__init__(R, names, category)

    @cached_method
    def basis(self):
        r""" 
        The basis of this Lie conformal algebra.

        EXAMPLES::

            sage: Vir = VirasoroLieConformalAlgebra(QQ)
            sage: Vir.inject_variables()
            Defining L, C
            sage: B = Vir.basis(); B
            Lazy family (basis map(i))_{i in Disjoint union of Family (The Cartesian product of ({'L'}, Non negative integers), The Cartesian product of ({'C'}, {0}))}
            sage: B[('L',2)]
            T^(2)L
            sage: sorted(B[('L',3)].bracket(L).items())
            [(3, -TL), (4, -8*L), (6, -10*C)]
            sage: B[('L',0)] == L
            True

        """
        return Family(self._indices, self.monomial, name="basis map")

    def indices(self):
        """
        The index set that parametrizes the basis of this Lie conformal 
        algebra.

        EXAMPLES::

            sage: V = AffineLieConformalAlgebra(QQ, 'A2')
            sage: V.indices()
            Disjoint union of Family (The Cartesian product of ({alpha[2], alpha[1], alpha[1] + alpha[2], alphacheck[1], alphacheck[2], -alpha[2], -alpha[1], -alpha[1] - alpha[2]}, Non negative integers), The Cartesian product of ({'K'}, {0}))
            sage: Vir = VirasoroLieConformalAlgebra(QQ)
            sage: Vir.indices()
            Disjoint union of Family (The Cartesian product of ({'L'}, Non negative integers), The Cartesian product of ({'C'}, {0}))

        """
        return self._indices

    def module(self):
        r"""
        The underlying `R` module for this Lie conformal algebra.

        EXAMPLES::

            sage: Vir = VirasoroLieConformalAlgebra(QQ); Vir.module()
            Free module generated by Disjoint union of Family (The Cartesian product of ({'L'}, Non negative integers), The Cartesian product of ({'C'}, {0})) over Rational Field
            sage: V = AffineLieConformalAlgebra(QQ, 'A2'); V.module()
            Free module generated by Disjoint union of Family (The Cartesian product of ({alpha[2], alpha[1], alpha[1] + alpha[2], alphacheck[1], alphacheck[2], -alpha[2], -alpha[1], -alpha[1] - alpha[2]}, Non negative integers), The Cartesian product of ({'K'}, {0})) over Rational Field
    
        """
        return CombinatorialFreeModule(self.base_ring(), 
                                       basis_keys=self._indices,
                                       category=self.category())

    def monomial(self,i):
        """
        The monomial of this Lie conformal algebra parametrized by this index.

        EXAMPLES::

            sage: Vir = VirasoroLieConformalAlgebra(QQ)
            sage: Vir.monomial(('L',3))
            T^(3)L
            sage: L = Vir.0; Vir.monomial(('L',4)) == 1/24*L.T(4)
            True
        
        """
        B = self.module().basis()
        return self.element_class(self, self.module()(B[i]))
    
    def zero(self): 
        """
        The zero element of this Lie conformal algebra.

        EXAMPLES::

            sage: Vir = VirasoroLieConformalAlgebra(QQ)
            sage: Vir.zero()
            0
            sage: L = Vir.0; L.nproduct(L,2) == Vir.zero()
            True
        """
        return self.element_class(self, self.module().zero())

class LieConformalAlgebraWithGenerators(LieConformalAlgebraWithBasis):
    def __init__(self,R, names=None, index_set=None, central_elements=None,      
                 category=None, prefix='B', **kwds):
        """
        Base class for a Lie conformal algebra with distinguished
        generators.

        .. NOTE::

            We now only accept direct sums of free modules plus 
            finitely many central
            generators `C_i` such that `TC_i = 0`.
        """

        self._generators = tuple(index_set)

        E = cartesian_product([index_set, NonNegativeIntegers()])
        if central_elements is not None:
            self._generators = self._generators + tuple(central_elements)
            E = DisjointUnionEnumeratedSets((E, cartesian_product([
                tuple(central_elements), {Integer(0)}])))
    
        super(LieConformalAlgebraWithGenerators,self).__init__(
                R,names=names, index_set=E, category=category, 
                prefix=prefix, **kwds)

        self._central_elements = tuple(central_elements)

    def lie_conformal_algebra_generators(self):
        """
        The generators of this Lie conformal algebra.
        
        OUTPUT: a (possibly infinite) family of generators (as an 
        `R[T]`-module) of this Lie conformal algebra.

        EXAMPLES::

            sage: Vir = VirasoroLieConformalAlgebra(QQ)
            sage: Vir.lie_conformal_algebra_generators()
            Finite family {'L': L, 'C': C}
            sage: V = AffineLieConformalAlgebra(QQ,'A1')
            sage: V.lie_conformal_algebra_generators()
            Finite family {alpha[1]: alpha[1], alphacheck[1]: alphacheck[1], -alpha[1]: -alpha[1], 'K': K}
        """
        return Family(self._generators, 
                      lambda i: self.monomial((i,Integer(0))), 
                      name = "generator map")

    def central_elements(self):
        """
        The central generators of this Lie conformal algebra.

        EXAMPLES::

            sage: Vir = VirasoroLieConformalAlgebra(QQ)
            sage: Vir.central_elements()
            (C,)
            sage: V = AffineLieConformalAlgebra(QQ, 'A1')
            sage: V.central_elements()
            (K,)
        """
        G = self.lie_conformal_algebra_generators()
        return tuple(G[i] for i in self._central_elements)


class FinitelyGeneratedLieConformalAlgebra(LieConformalAlgebraWithGenerators):
    def __init__(self,R, names=None, index_set=None, central_elements=None,      
            category=None, prefix='B', **kwds):
        """
        Base class for finitely generated Lie conformal algebras.
        """
        super(FinitelyGeneratedLieConformalAlgebra,self).__init__(R,names=names,
            index_set=index_set,central_elements=central_elements, 
            category=category, prefix=prefix, **kwds)
        self.__ngens = len(self._generators)

    def _repr_(self):
        """
        The name of this Lie conformal algebra
        """
        if self.__ngens == 1:
            return "Lie conformal algebra generated by {0} over {1}".format(
                self.gen(0), self.base_ring())
        return "Lie conformal algebra with generators {0} over {1}.".format(
                                                 self.gens(), self.base_ring()) 
    def _an_element_(self):
        return self.sum(self.gens())

    def ngens(self):
        """
        The number of generators of this Lie conformal algebra

        EXAMPLES::

            sage: Vir = VirasoroLieConformalAlgebra(QQ); Vir.ngens()
            2
            sage: V = AffineLieConformalAlgebra(QQ, 'A1'); V.ngens()
            4
        """
        return self.__ngens

    @cached_method
    def gens(self):
        """
        The generators for this Lie conformal algebra.
        
        OUTPUT: This method returns a tuple with the (finite) generators
        of this Lie conformal algebra. 

        EXAMPLES::

            sage: Vir = VirasoroLieConformalAlgebra(QQ);
            sage: Vir.gens()
            (L, C)

        .. SEEALSO::

            :meth:`lie_conformal_algebra_generators<LieConformalAlgebraWithGenerators.lie_conformal_algebra_generators>`
        """
        G = self.lie_conformal_algebra_generators()
        try:
            return tuple(G[i] for i in self._generators)
        except (KeyError, ValueError):
            return tuple(G)
 
class LieConformalAlgebraWithStructureCoefficients(
            FinitelyGeneratedLieConformalAlgebra, IndexedGenerators):
    r"""
    A Lie conformal algebra with a set of specified structure
    coefficients.

    INPUT:
 
    - ``R`` -- a ring (Default: None); The base ring of this Lie 
      conformal algebra. Behaviour is undefined if it is not a field
      of characteristic zero. 

    - ``s_coeff`` -- Dictionary (Default: None); 
      a dictionary containing the `\lambda` brackets of the
      generators of this Lie conformal algebra. The keys of this 
      dictionary are pairs of either names or indices of the generators 
      and the values are themselves dictionaries. For a pair of 
      generators `a` and `b`, the value of ``s_coeff[('a','b')]`` is a 
      dictionary whose keys are positive integer numbers and the 
      corresponding value for the key `j` is a dictionary itself 
      representing the j-th product `a_{(j)}b`. 
      Thus, for a positive integer number `j`, the value of
      ``s_coeff[('a','b')][j]`` is a dictionary whose entries are pairs 
      ``('c',n)`` where ``'c'`` is the name of a generator and `n` is 
      a positive number. The value for this key is the coefficient of
      `\frac{T^{n}}{n!} c` in `a_{(j)}b`. For example the ``s_coeff`` 
      for the *Virasoro* Lie conformal algebra is::

            {('L','L'):{0:{('L',1):1}, 1:{('L',0):2}, 3:{('C',0):1/2}}}


      Do not include central elements in this dictionary. Also, if the 
      key ``('a','b')`` is present, there is no need to include
      ``('b','a')`` as it is defined by skew-symmetry.
      Any missing pair (besides the ones
      defined by skew-symmetry) is assumed to have vanishing
      `\lambda`-bracket. 

    - ``names`` -- tuple of ``str`` (Default: None); The list of names 
      for generators of this Lie conformal algebra. Do not include 
      central elements in this list.

    - ``central_elements`` -- tuple of ``str`` (Default: None); 
      A list of names for central
      elements of this Lie conformal algebra. 

    - ``index_set`` -- enumerated set (Default: None); 
      an indexing set for the generators of this Lie
      conformal algebra. Do not include central elements in this list. 
    """
    @staticmethod
    def __classcall__(cls,R,s_coeff, names=None,
            index_set=None, central_elements=None, **kwds):

        names, index_set = standardize_names_index_set(names,index_set)

        if names is not None and names != tuple(index_set):
            names2 = names + tuple(central_elements)
            index_set2 = DisjointUnionEnumeratedSets((index_set,
                Family(tuple(central_elements))))
            d = {x:index_set2[i] for i,x in enumerate(names2)}
            try:
                #If we are given a dictionary with names as keys, 
                #convert to index_set as keys
                s_coeff = {(d[k[0]],d[k[1]]):{a:{(d[x],y): 
                    s_coeff[k][a][(d[x],y)] for (x,y) in 
                    s_coeff[k][a]} for a in s_coeff[k]} for k in s_coeff}

            except KeyError:
                # We assume the dictionary was given with keys in the 
                # index_set
                pass


        s_coeff = LieConformalAlgebraWithStructureCoefficients\
                    ._standardize_s_coeff(s_coeff, index_set, central_elements)
   
 
        if names is not None and central_elements is not None:
            names += tuple(central_elements)

        return super(LieConformalAlgebraWithStructureCoefficients, cls)\
                .__classcall__(cls, R, s_coeff, names=names, 
                               index_set=index_set, 
                               central_elements=central_elements,
                               **kwds)

    @staticmethod
    def _standardize_s_coeff(s_coeff, index_set, ce):
        """
        Convert an input dictionary to structure constants of this 
        Lie conformal algebra.

        INPUT:

        - ``s_coeff`` -- a dictionary as in 
          :class:`LieConformalAlgebraWithStructureCoefficients`.
        - ``index_set` -- A finite enumerated set indexing the 
          generators not counting the central elements. 
        - ``ce`` -- a tuple of ``str``; a list of names for the central
          generators of this Lie conformal algebra

        OUTPUT: a Finite Family representing ``s_coeff`` in the input, 
        with the indexes ordered by applying skew-symmetry if necessary.
        """
        index_to_pos = {k:i for i,k in enumerate(index_set)}

        sc = {}
        #mypair has a pair of generators
        for mypair in s_coeff.keys():
            #e.g.  v = { 0: { (L,2):3, (G,3):1}, 1:{(L,1),2} }
            v = s_coeff[mypair]
            if index_to_pos[mypair[0]] > index_to_pos[mypair[1]]:
                key=(mypair[1],mypair[0])
                maxpole = max(v.keys())
                vals={} 
                for k in range(maxpole+1):
                    kth_product = {}
                    for j in range(maxpole+1-k):
                        if k+j in v.keys():
                            for i in v[k+j].keys():
                                if (i[0] not in ce) or (
                                    i[0] in ce and i[1] + j == 0):
                                    kth_product[(i[0],i[1]+j)] = \
                                            kth_product.get((i[0], i[1]+j), 0)
                                    kth_product[(i[0],i[1]+j)] += \
                                    v[k+j][i]*(-1)**(k+j+1)*binomial(i[1]+j,j)
                    kth_product = {k:v for k,v in kth_product.items() if v}
                    if kth_product:
                        vals[k]=kth_product
            else:
                if not index_to_pos[mypair[0]] < index_to_pos[mypair[1]]:
                    if mypair[0] == mypair[1]:
                        #Here I should check skew-symmetry!
                        pass
                key = tuple(mypair)
                vals={}
                for l in v.keys():
                    lth_product = {k:y for k,y in v[l].items() if y}
                    if lth_product:
                        vals[l]=lth_product

            myvals = tuple((i, tuple((x,v) for x,v in vals[i].items())) 
                            for i in vals.keys() if vals[i])
            
            if key in sc.keys():
                if sorted(sc[key]) != sorted(myvals): 
                        raise ValueError("two distinct values given for one "\
                                         "and the same bracket")
            if myvals:
                sc[key] = myvals
        return Family(sc)

    def __init__(self, R, s_coeff, names=None, index_set=None,
                 category=None, central_elements=None,  bracket=False,
                 latex_bracket=False, string_quotes=False, prefix='', **kwds):
        
        self._index_to_pos = {k: i for i,k in enumerate(index_set)}
        #Add central parameters to index_to_pos so that we can 
        #represent names
        for i,ce in enumerate(central_elements):
            self._index_to_pos[ce] = len(index_set)+i

        if "sorting_key" not in kwds:
            kwds["sorting_key"] = self._index_to_pos.__getitem__

        category = LieConformalAlgebras(R).WithBasis().FinitelyGenerated().\
                                                       or_subcategory(category)
        FinitelyGeneratedLieConformalAlgebra.__init__(
            self, R, names=names, index_set=index_set, category=category,
            central_elements=central_elements, prefix=prefix, **kwds)

        IndexedGenerators.__init__(self, self._indices, prefix=prefix,
                                   bracket=bracket, latex_bracket=latex_bracket,
                                   string_quotes=string_quotes)

        #at this stage we have access to self.module(). We convert the values
        #of s_coeff to be elements in self.
   
        self._s_coeff = Family({k: tuple((p[0], sum(c[1]*self.monomial(c[0]) 
                for c in p[1] )) for p in s_coeff[k]) for k in s_coeff.keys()})

    
    Element = LCAStructureCoefficientsElement
    _repr_term = IndexedGenerators._repr_generator
    _latex_term = IndexedGenerators._latex_generator

    def structure_coefficients(self):
        """
        The structure coefficients of this Lie conformal algebra.

        EXAMPLES::

            sage: R = AffineLieConformalAlgebra(QQbar, 'A1')
            sage: scoef = R.structure_coefficients()
            sage: sorted(scoef)
            [((0, -2*alpha[1]),),
             ((0, alphacheck[1]), (1, K)),
             ((0, -2*-alpha[1]),),
             ((1, 2*K),)]

            sage: Vir = VirasoroLieConformalAlgebra(AA)
            sage: scoef = Vir.structure_coefficients()
            sage: sorted(scoef)
            [((0, TL), (1, 2*L), (3, 1/2*C))]
        """
        return self._s_coeff

class GradedLieConformalAlgebra(LieConformalAlgebraWithStructureCoefficients):
    def __init__(self, R, s_coeff, names=None, index_set=None,
            category=None, central_elements=None,  bracket=False,
            latex_bracket=False, string_quotes=False, prefix='', 
            weights=None,**kwds):
        r"""
        An H-Graded Lie conformal algebra

        INPUT:

        - ``R`` -- a commutative ring (Default: None); The base ring of
          this Lie conformal algebra. Behaviour is undefined if it is 
          not a field of characteristic zero. 

        - ``s_coeff`` -- Dictionary (Default: None); As in the input of
          :class:`LieConformalAlgebraStructureCoefficients` 

        - ``names`` -- tuple of ``str`` (Default: None); as in the
          input of 
          :class:`LieConformalAlgebraStructureCoefficients` 

        - ``central_elements`` -- tuple of ``str`` (Default: None); as 
          in the input of 
          :class:`LieConformalAlgebraStructureCoefficients` 

        - ``index_set`` -- enumerated set (Default: None); As in the 
          input of 
          :class:`LieConformalAlgebraStructureCoefficients` 
          
        - ``weights`` -- tuple of non-negative integers 
          (Default: None); a list of degrees for this Lie conformal 
          algebra. The returned Lie conformal algebra is H-Graded. 
          This tuple needs to have the same cardinality as
          ``index_set`` or ``names``. Central elements are assumed 
          to have weight `0`. 

        - ``category`` The category that this Lie conformal algebra belongs 
          to.
        """
        category = LieConformalAlgebras(R).Graded().WithBasis()\
                   .FinitelyGenerated().or_subcategory(category)
        LieConformalAlgebraWithStructureCoefficients.__init__(
            self,R, s_coeff,
            names=names, category=category, central_elements=central_elements, 
            bracket=bracket, latex_bracket=latex_bracket, index_set=index_set,
            string_quotes=string_quotes, prefix=prefix, **kwds)
        if weights is None:
            weights = (1,)* (len(self._generators) - 
                             len(self.central_elements()))
        if len (weights) != (len(self._generators) - 
                                len(self.central_elements())):
            raise ValueError("weights and (non-central) generator lists "\
                             "must be of same length")
        self._weights = weights

    class Element(LCAStructureCoefficientsElement):
        def degree(self):
            """
            The degree of this element

            EXAMPLES::

                sage: V = VirasoroLieConformalAlgebra(QQ)
                sage: V.inject_variables()
                Defining L, C
                sage: C.degree()
                0
                sage: L.T(4).degree()
                6
            """
            if self.is_zero():
                return Infinity
            p = self.parent()
            ls = []
            for idx in self.value.monomial_coefficients().keys():
                ret = idx[1]
                if idx[0] not in p._central_elements:
                    ret+= p._weights[p._index_to_pos[idx[0]]] 
                ls.append(ret)
            if ls[1:] == ls[:-1]:
                return ls[0]
            raise ValueError("{} is not homogeneous!".format(self))
 
class VirasoroLieConformalAlgebra(GradedLieConformalAlgebra):
    #This is the problem of following the Lie algebra Factory design
    #To have a named algebra represented in a different way we need to 
    #go over all of this trouble instead of a simple def Viras...
    """
    The Virasoro Lie Conformal algebra over `R`

    INPUT:

    - ``R``: a commutative Ring. Behaviour is undefined if `R` is not 
      a Field of characteristic zero. 

    EXAMPLES::

            sage: Vir = VirasoroLieConformalAlgebra(QQ)
            sage: Vir.category()
            Category of finitely generated H-graded Lie conformal algebras with basis over Rational Field
            sage: Vir.gens()
            (L, C)
            sage: L = Vir.0
            sage: sorted(L.bracket(L).items())
            [(0, TL), (1, 2*L), (3, 1/2*C)]

    """
    @staticmethod
    def __classcall__(cls,R):
        virdict =  {('L','L'):{0:{('L',1):1}, 1:{('L',0): 2}, 
                    3:{('C', 0):R(2).inverse_of_unit()}}}
        return GradedLieConformalAlgebra.__classcall__(cls,R, virdict, 
            names = ('L',), central_elements = ('C'), weights = (2,))

    def __init__(self, R, s_coeff, names=None, index_set=None,
            category=None, central_elements=None, weights=None):
        GradedLieConformalAlgebra.__init__(self,R,s_coeff, 
            names=names, index_set=index_set,
            central_elements=central_elements,weights=weights)

    def _repr_(self):
        return "The Virasoro Lie conformal algebra over {}".format(
                                                            self.base_ring())



class AffineLieConformalAlgebra(GradedLieConformalAlgebra):
    r"""
    The current Lie conformal algebra

    INPUT:

    - ``R`` -- a Ring; The base ring for this Lie conformal algebra. 
    - ``ct`` -- a ``str`` or a ``CartanType``; The Cartan Type for the
      corresponding finite dimensional Lie algebra. It must correspond to a
      simple finite dimensional Lie algebra.

    EXAMPLES::

            sage: V = AffineLieConformalAlgebra(QQ, 'A1')
            sage: V
            Lie conformal algebra with generators (alpha[1], alphacheck[1], -alpha[1], K) over Rational Field.
            sage: AffineLieConformalAlgebra(QQ, CartanType('A1')) is V
            True

            sage: V = AffineLieConformalAlgebra(QQ, CartanType(["A",2,1]))
            Traceback (most recent call last):
            ...
            ValueError: Only affine algebras of simple finite dimensionalLie algebras are implemented

    """
    @staticmethod
    def __classcall__(cls,R,ct,**kwds):
        if type(ct) is str: 
            from sage.combinat.root_system.cartan_type import CartanType
            ct = CartanType(ct)
        if not ( ct.is_finite() and ct.is_irreducible ):
            raise ValueError("Only affine algebras of simple finite dimensional"
                "Lie algebras are implemented")
        hv = Integer(ct.dual_coxeter_number())
        g = LieAlgebra(R, cartan_type=ct)
        B = g.basis()
        S = B.keys()
        gdict = {}
        for k1 in S:
            for k2 in S:
                if S.rank(k2) < S.rank(k1):
                    myb = B[k1].bracket(B[k2]).monomial_coefficients()
                    myf = R(2).inverse_of_unit()*R(hv).inverse_of_unit()\
                          *g.killing_form(B[k1],B[k2])
                    if myb or myf:
                        gdict[(k1,k2)] = {}
                        if myb:
                            gdict[(k1,k2)][0] = {(nk,0):myb[nk] for nk in
                                                 myb.keys()} 
                        if myf:
                            gdict[(k1,k2)][1] = {('K',0):myf}

        weights = (1,)*B.cardinality()
        prefix = kwds.get('prefix','E')
        bracket = kwds.get('bracket','[')
        names = kwds.get('names', None)
        return GradedLieConformalAlgebra.__classcall__(cls,
                    R, gdict, index_set=B.keys(),  
                    central_elements=('K',), weights=weights, 
                    prefix=prefix, bracket=bracket, names=names, 
                    cartan_type=ct)

    def __init__(self, R, s_coeff, names=None, index_set=None,
            category=None, central_elements=None, weights=None,
            cartan_type=None, prefix=None, bracket=None):
        
        self._ct = cartan_type
        GradedLieConformalAlgebra.__init__(self,R,s_coeff, 
            names=names, index_set=index_set,
            central_elements=central_elements,weights=weights,
            bracket=bracket,prefix=prefix)

    def cartan_type(self):
        """
        The Cartan type of this Lie conformal algbera.

        EXAMPLES::
            sage: R = AffineLieConformalAlgebra(QQ, 'B3') 
            sage: R
            The affine Lie conformal algebra of type ['B', 3] over Rational Field
            sage: R.cartan_type()
            ['B', 3]     
        """
        return self._ct

    def _repr_(self):
        return "The affine Lie conformal algebra of "\
            "type {} over {}".format(self._ct,self.base_ring())

