"""
Lie Conformal Algebras With Structure Coefficients.

AUTHORS:

- Reimundo Heluani (08-09-2019): Initial implementation

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

from sage.functions.other import binomial
from sage.structure.indexed_generators import (IndexedGenerators, 
                                               standardize_names_index_set)
from sage.sets.family import Family
from .lie_conformal_algebra_element import LCAStructureCoefficientsElement
from sage.categories.lie_conformal_algebras import LieConformalAlgebras
from .finitely_freely_generated_lca import FinitelyFreelyGeneratedLCA
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.structure.indexed_generators import (IndexedGenerators,
                                               standardize_names_index_set)

class LieConformalAlgebraWithStructureCoefficients(
            FinitelyGeneratedLieConformalAlgebra, IndexedGenerators):
    @staticmethod
    def _standardize_s_coeff(s_coeff, index_set, ce, parity=None):
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
        - ``parity`` -- a tuple of `0` or `1` (Default: tuple of `0`); 
          this tuple specifies the parity of each non-central generator.

        OUTPUT: a Finite Family representing ``s_coeff`` in the input. 
        It contains superfluous information that can be obtained by 
        skew-symmetry 
        """
        if parity is None:
            parity = (0,)*index_set.cardinality()
        index_to_parity = {i:p for (i,p) in zip(index_set,parity)}
        sc = {}
        #mypair has a pair of generators
        for mypair in s_coeff.keys():
            #e.g.  v = { 0: { (L,2):3, (G,3):1}, 1:{(L,1),2} }
            v = s_coeff[mypair]
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
                        raise ValueError("two distinct values given for one "
                                         "and the same bracket, skew-symmetry" + 
                                         "is not satisfied?")
            if myvals:
                sc[key] = myvals
            
            #We now add the skew-symmetric part to optimize 
            #brackets computations later
            key=(mypair[1],mypair[0])
            if index_to_parity[mypair[0]]*index_to_parity[mypair[1]]:
                parsgn = -1
            else:
                parsgn = 1
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
                                kth_product[(i[0],i[1]+j)] += parsgn*\
                                v[k+j][i]*(-1)**(k+j+1)*binomial(i[1]+j,j)
                kth_product = {k:v for k,v in kth_product.items() if v}
                if kth_product:
                    vals[k]=kth_product
           
            myvals = tuple((i, tuple((x,v) for x,v in vals[i].items())) 
                            for i in vals.keys() if vals[i])

            if key in sc.keys():
                if sorted(sc[key]) != sorted(myvals): 
                        raise ValueError("two distinct values given for one "+
                                         "and the same bracket. "+ 
                                         "Skew-symmetry is not satisfied?")
            if myvals:
                sc[key] = myvals
        return Family(sc)

    def __init__(self, R, s_coeff, **kwds):
        r"""
        A Lie conformal algebra with a set of specified structure
        coefficients.

        INPUT:
     
        - ``R`` -- a ring (Default: None); The base ring of this Lie 
          conformal algebra. Behaviour is undefined if it is not a field
          of characteristic zero. 

        - ``s_coeff`` -- Dictionary (Default: None); 
          a dictionary containing the `\lambda` brackets of the
          generators of this Lie conformal algebra. The family encodes a 
          dictionary whose keys 
          are pairs of either names or indices of the generators 
          and the values are themselves dictionaries. For a pair of 
          generators `a` and `b`, the value of ``s_coeff[('a','b')]`` is
          a dictionary whose keys are positive integer numbers and the 
          corresponding value for the key `j` is a dictionary itself 
          representing the j-th product `a_{(j)}b`. 
          Thus, for a positive integer number `j`, the value of
          ``s_coeff[('a','b')][j]`` is a dictionary whose entries are 
          pairs ``('c',n)`` where ``'c'`` is the name of a generator 
          and `n` is a positive number. The value for this key is the 
          coefficient of `\frac{T^{n}}{n!} c` in `a_{(j)}b`. For example
          the ``s_coeff`` for the *Virasoro* Lie conformal algebra is::

                {('L','L'):{0:{('L',1):1}, 1:{('L',0):2}, 3:{('C',0):1/2}}}


          Do not include central elements in this dictionary. Also, if 
          the key ``('a','b')`` is present, there is no need to include
          ``('b','a')`` as it is defined by skew-symmetry.
          Any missing pair (besides the ones
          defined by skew-symmetry) is assumed to have vanishing
          `\lambda`-bracket. 

        - ``names`` -- tuple of ``str`` (Default: None); The list of 
          names for generators of this Lie conformal algebra. Do not 
          include central elements in this list.

        - ``central_elements`` -- tuple of ``str`` (Default: None); 
          A list of names for central
          elements of this Lie conformal algebra. 

        - ``index_set`` -- enumerated set (Default: None); 
          an indexing set for the generators of this Lie
          conformal algebra. Do not include central elements in this 
          list. 

        - ``parity`` -- tuple of `0` or `1` (Default: tuple of `0`);
           a tuple specifying the parity of each non-central generator.

        EXAMPLES:

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

        - We construct the centerless Virasoro Lie conformal algebra::

            sage: virdict =  {('L','L'):{0:{('L',1):1}, 1:{('L',0): 2}}}
            sage: R = LieConformalAlgebra(QQbar, virdict, names='L')
            sage: R.inject_variables()
            Defining L
            sage: L.bracket(L)
            {0: TL, 1: 2*L}

        - The construction checks that skew-symmetry is violated::

            sage: wrongdict =  {('L','L'):{0:{('L',1):2}, 1:{('L',0): 2}}}
            sage: LieConformalAlgebra(QQbar, wrongdict, names='L')
            Traceback (most recent call last):
            ...
            ValueError: two distinct values given for one and the same bracket. Skew-symmetry is not satisfied?
        """
        names = kwds.get('names', None)
        index_set = kwds.get('index_set', None)
        central_elements = kwds.get('central_elements', ())
        names, index_set = standardize_names_index_set(names,index_set)

        if names is not None and names != tuple(index_set):
            names2 = names + tuple(central_elements)
            index_set2 = DisjointUnionEnumeratedSets((index_set,
                Family(tuple(central_elements))))
            d = {x:index_set2[i] for i,x in enumerate(names2)}
            try:
                #If we are given a dictionary with names as keys, 
                #convert to index_set as keys
                s_coeff = {(d[k[0]],d[k[1]]):{a:{(d[x[1]],x[2]): 
                    s_coeff[k][a][x] for x in 
                    s_coeff[k][a]} for a in s_coeff[k]} for k in s_coeff.keys()}

            except KeyError:
                # We assume the dictionary was given with keys in the 
                # index_set
                pass

        parity = kwds.get('parity',None)
        if parity is None:
            issuper = kwds.get('super',False)
            parity = (0,)*index_set.cardinality()
        else:
            issuper = True
            
        try:
            assert len(parity) == index_set.cardinality()
        except AssertionError:
            raise ValueError("parity should have the same length as the "+
                             "number of generators, got {}".format(parity))

        s_coeff = LieConformalAlgebraWithStructureCoefficients\
                    ._standardize_s_coeff(s_coeff, index_set, central_elements,
                                          parity)
   
        if names is not None and central_elements is not None:
            names += tuple(central_elements)

        self._index_to_pos = {k: i for i,k in enumerate(index_set)}
        #Add central parameters to index_to_pos so that we can 
        #represent names
        for i,ce in enumerate(central_elements):
            self._index_to_pos[ce] = len(index_set)+i
        category = kwds.get('category', None)
        category = LieConformalAlgebras(R).WithBasis().FinitelyGenerated()\
                    .or_subcategory(category)
        if issuper:
            category = category.Super()

        FinitelyGeneratedLieConformalAlgebra.__init__(
            self, R, names=names, index_set=index_set, category=category,
            central_elements=central_elements)

        sorting_key = kwds.get("sorting_key", self._index_to_pos.__getitem__)
        prefix = kwds.get("prefix", '')
        bracket = kwds.get("bracket",'')
        string_quotes = kwds.get("string_quotes", None)
        latex_bracket = kwds.get("latex_bracket", None)
        IndexedGenerators.__init__(self, self._indices, prefix=prefix,
                                   bracket=bracket, latex_bracket=latex_bracket,
                                   string_quotes=string_quotes,
                                   sorting_key=sorting_key)

        #at this stage we have access to self.module(). We convert the values
        #of s_coeff to be elements in self.
        self._s_coeff = Family({k: tuple((p[0], sum(c[1]*self.monomial(c[0]) 
                for c in p[1] )) for p in s_coeff[k]) for k in s_coeff.keys()})
        self._parity = dict(zip(self.gens(),parity+(0,)*len(central_elements)))



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
            [((0, 2*alpha[1]),),
             ((0, -2*alpha[1]),),
             ((0, -alphacheck[1]), (1, K)),
             ((0, alphacheck[1]), (1, K)),
             ((0, 2*-alpha[1]),),
             ((0, -2*-alpha[1]),),
             ((1, 2*K),)]

            sage: Vir = VirasoroLieConformalAlgebra(AA)
            sage: scoef = Vir.structure_coefficients()
            sage: sorted(scoef)
            [((0, TL), (1, 2*L), (3, 1/2*C))]
        """
        return self._s_coeff

    def universal_enveloping_algebra(self, 
                                    central_parameters=None, 
                                    names=None):
        r"""
        The universal enveloping vertex algebra of this Lie 
        conformal algebra.

        INPUT:

        - ``central_parameters`` -- A family of constants in the 
          base ring of this Lie conformal algebra parametrized by 
          the central elements.

        - ``names`` -- The names of the generators of the universal
          enveloping vertex algebra.

        OUTPUT: 

        If `L` is a Lie conformal algebra over `R` with some 
        central elements `C_i \in L` indexed by a set `I`, 
        ``central_parameters`` is a family of elements `c_i \in R`
        indexed by the same set, then this method returns the 
        central quotient of the universal enveloping vertex algebra
        of `L` by the ideal generated by `C_i - c_i |0\rangle`. 
    
        EXAMPLES:

        We construct the universal enveloping vertex algebra of the
        Virasoro Lie conformal algebra at central charge `0` over 
        the complex numbers::

            sage: Vir = VirasoroLieConformalAlgebra(CC)
            sage: V = Vir.universal_enveloping_algebra(); V
            The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Complex Field with 53 bits of precision
            sage: V.0.bracket(V.0)
            {0: 1.00000000000000*1.00000000000000*L_-3|0>,
             1: 2.00000000000000*1.00000000000000*L_-2|0>}

        We construct the universal enveloping vertex algebra of the
        Virasoro Lie conformal algebra at central charge 2 over the
        rationals::

            sage: Vir = VirasoroLieConformalAlgebra(QQ); 
            sage: Vir.gens()
            (L, C)
            sage: cp = Family({Vir.1:2}); V = Vir.universal_enveloping_algebra(cp)
            sage: V.0.nproduct(V.0,3)
            |0>

        The Virasoro Algebra is not defined over `\ZZ`::

            sage: Vir = VirasoroLieConformalAlgebra(ZZ)
            Traceback (most recent call last):
            ...
            ArithmeticError: inverse does not exist

        The list of central parameters needs to be indexed by the
        central elements of the algebra::

            sage: Vir = VirasoroLieConformalAlgebra(QQ);
            sage: cp = Family({Vir.0 : 3})
            sage: V = Vir.universal_enveloping_algebra(cp)
            Traceback (most recent call last):
            ...
            ValueError: central_parameters must be parametrized by central elements
        """  
        if central_parameters == None:
            from sage.sets.family import Family
            central_parameters = Family({ce:0 for ce in
                                        self.central_elements()})
        if hasattr(self, 'lift'):
            cp = self.lift.codomain().central_parameters()
            if cp == central_parameters:
                return self.lift.codomain()
        V = self._construct_UEA(central_parameters, 
                                names=names)

        from sage.categories.homset import Hom
        self.lift = _LiftMorphism(Hom(self, V, category = 
                        LieConformalAlgebras(self.base_ring())))
        try: 
            self.lift.register_as_coercion()
        except AssertionError:
            #we already constructed this morphisms and its fine
            pass

        return V

    def _construct_UEA(self, central_parameters=None, names=None):
        """ 
        Returns the universal enveloping vertex algebra of ``self``. 

        see :meth:`universal_enveloping_algebra`
        """
        from sage.algebras.vertex_algebras.vertex_algebra\
                                                    import VertexAlgebra
        return VertexAlgebra(base_ring=self.base_ring(),
                            lie_conformal_algebra=self, 
                            central_parameters=central_parameters, 
                            names=names)

from sage.categories.morphism import Morphism
class _LiftMorphism(Morphism):
    """
    The lift morphism from this Lie conformal algebra to its universal
    enveloping vertex algebra

    EXAMPLES::
    
        sage: Vir = VirasoroLieConformalAlgebra(AA)
        sage: cp = Family({Vir.1:2})
        sage: V = Vir.universal_enveloping_algebra(cp)
        sage: Vir.lift
        Generic morphism:
          From: The Virasoro Lie conformal algebra over Algebraic Real Field
          To:   The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Algebraic Real Field
        sage: W = VirasoroVertexAlgebra(AA,1)
        sage: Vir.lift
        Generic morphism:
          From: The Virasoro Lie conformal algebra over Algebraic Real Field
          To:   The Virasoro vertex algebra of central charge 1 over Algebraic Real Field
    """
    def _call_(self,x):
        return x.lift()    

