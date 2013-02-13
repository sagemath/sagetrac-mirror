"""
Some basis for `AbstractPolynomialRing` and `FiniteAbstractPolynomialRing` :

`SchubertBasisOnVectors`, `DemazureBasisOnVectors`, `GrothendieckBasisOnVectors`

"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Viviane Pons < pons at univ-mlv.fr  >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.permutation import Permutation 
from sage.combinat.words.word import Word
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.combinat.sf.sf import SymmetricFunctions
from sage.structure.unique_representation import UniqueRepresentation

from basis import PolynomialRingWithBasisFromMorphism, FinitePolynomialRingWithBasisFromMorphism, PolynomialRingWithBasis

class LinearBasisOnVectors(PolynomialRingWithBasisFromMorphism):
    r"""
    Explain the use of this class...

    EXAMPLES::

        sage: # Fix a nice example

    TESTS::

        sage: # Fix a nice test
    """
    def __init__(self, abstract_polynomial_ring, ambient_space_basis, basis_name, basis_repr, on_basis_method =None, extra_parameters = (), variables_auto_coerce = False,  **keywords ):
        r"""
        TESTS::

            sage: # Fix a nice test
        """
        self._ambient_space_basis = ambient_space_basis
        get_morphism_on_basis = self._get_morphism_on_basis
        
        if(not keywords.has_key("triangular")): keywords["triangular"] = "upper"
        if(not keywords["triangular"] is None):
            if(not keywords.has_key("cmp")): keywords["cmp"] = self.cmp
        self._keywords = keywords
        
        self._extra_parameters = {}
        for (key,val) in extra_parameters: self._extra_parameters[key] = val
        self._on_basis_method = on_basis_method
        
        PolynomialRingWithBasisFromMorphism.__init__(
            self,
            abstract_polynomial_ring,
            1,
            ambient_space_basis,
            basis_name,
            basis_repr,
            self._get_basis_keys,
            get_morphism_on_basis,
            variables_auto_coerce = variables_auto_coerce,
            **keywords
        )
        
    def cmp(self, key1, key2):
        l = len(key1.parent()._basis_keys)
        d1 = sum( [key1[i] for i in xrange(l)]) 
        d2 = sum( [key2[i] for i in xrange(l)]) 
        if (d1 > d2): return -1
        if (d1 < d2): return 1
        for i in xrange(l-1,-1,-1):
            if (key1[i]>key2[i]):
                return 1
            if (key1[i]<key2[i]):
                return -1
        return 0   

    def _get_basis_keys(self, n):
        r"""
        TESTS::

            sage: # Fix a nice example
        """
        return self._ambient_space_basis.finite_basis(n)._basis_keys

    def ambient_space_basis(self):
        r"""
        Return the ambient space indexing the basis elements of ``self``.

        EXAMPLES::

            sage: # Fix a nice example
        """
        return self._ambient_space_basis

    def group_type(self):
        r"""
        Returns the group type of the group acting on ``self``.

        EXAMPLES::

            sage: # Fix a nice example
        """
        return self._ambient_space_basis.group_type() 

    def _finite_basis_instance(self, n):
        r"""
        TESTS::

            sage: # Fix a nice test
        """
        return self.abstract_algebra().finite_polynomial_ring(n).linear_basis_on_vectors(self, self._basis_name, self._basis_repr, **self._keywords) 

    def equivalent_basis(self, abstract_polynomial_ring):
        return abstract_polynomial_ring.linear_basis_on_vectors(self.group_type(), self._basis_name, self._basis_repr, self._on_basis_method, self._extra_parameters, **self._keywords)

        
    def morphism_method_wrapper(self, n):
        ambient_space_basis = self._ambient_space_basis.finite_basis(n)
        return self.MorphismMethodWrapper(ambient_space_basis, self._on_basis_method, **self._extra_parameters)
        
    def _get_morphism_on_basis(self, n):
        method_wrapper = self.morphism_method_wrapper(n)
        return method_wrapper.morphism_on_basis
        
    class MorphismMethodWrapper(UniqueRepresentation):
        
        def __init__(self, ambient_space_basis, on_basis_method, **keywords):
            self._parameters = keywords
            self._basis = ambient_space_basis
            self._call_back = self.morphism_on_list
            self._on_basis_method = on_basis_method
            
        def get_cached_method(self):
            return self.morphism_on_basis
            
        def clear_cache(self):
            self.get_cached_method().clear_cache()
            
        def get_cache(self):
            return self.get_cached_method().get_cache()
            
        def get_cached_keys(self):
             return [ i[0] for (i,j) in self.get_cache()]
             
        def morphism_on_list(self, x):
            return self.morphism_on_basis(self._basis.basis().keys()(x))
            
            
        @cached_method
        def morphism_on_basis(self, key):
            x = [key[i] for i in xrange(self._basis.nb_variables())]
            if(self._on_basis_method is None):
                raise NotImplementedError,"No compute method has been implemented"%()
            else:
                if (len(self._parameters)==0):
                    return self._on_basis_method(x, self._basis, self._call_back)
                else:
                    return self._on_basis_method(x, self._basis, self._call_back, **self._parameters)
            
            
        

class SchubertBasisOnVectors(LinearBasisOnVectors):
    r"""
    Explain this class

    EXAMPLES::

        sage: # Fix a nice example

    TESTS::

        sage: # Fix a nice test
    """
    
    def __init__(self, abstract_polynomial_ring, ambient_space_basis, basis_name, basis_repr):
        LinearBasisOnVectors.__init__(
            self,
            abstract_polynomial_ring,
            ambient_space_basis,
            basis_name,
            basis_repr,
            on_basis_method = self.on_basis_method,
            variables_auto_coerce = True
        )
        
        
    def equivalent_basis(self, abstract_polynomial_ring):
        return abstract_polynomial_ring.schubert_basis_on_vectors(self.group_type())

        
    def on_basis_method(self, x, basis, call_back):
        for i in xrange( len( x ) - 1 ):
            if( x[i]<x[i+1] ):
                x[i], x[i+1] = x[i+1]+1, x[i]
                return call_back(x).divided_difference(i+1)   
        return basis(x)
        
    
    def from_partition(self,part,nb_variables = None):
        r"""
        Returns the Schubert polynomial corresponding to the partition,
        i.e., the Schubert polynomial equal to the schur function indexed
        by the partition expanded in ``nb_variables`` variables. It is just the 
        Schubert element indexex by the reverse vector of the partition 
        (as a vector of size ``nb_variables``)
        
        INPUT: 
        - ``part``, a partition
        - ``nb_variables`` (optional) the number of variables to expand 
        the schur function (by default, the length of the partition)
        
        OUTPUT:
        - The Schubert polynomial corresponding to this partition
        
        EXAMPLES::
            sage: A = AbstractPolynomialRing(QQ)
            sage: Schub = A.schubert_basis_on_vectors()
            sage: l = Partition([3,2,1])
            sage: Schub.from_partition(l)
            Y(1, 2, 3)
            sage: Schub.from_partition(l,4)
            Y(0, 1, 2, 3)
            sage: Schub.from_partition(l,2)
            0
        
        """
        if nb_variables is None:
            nb_variables = len(part)
        if(nb_variables< len(part)):
            return self(0)
        key = [0 for i in xrange(nb_variables)]
        for i in xrange(len(part)): key[i] = part[i]
        key.reverse()
        return self(key)
    
        
    def from_schur(self, schur):
        r"""
        Returns the Schubert polynomial corresponding to the given 
        schur function
        
        INPUT:
        - a schur function
        
        OUTPUT:
        - the corresponding Schubert polynomial
        
        EXAMPLES::
            
            sage: A = AbstractPolynomialRing(QQ)
            sage: Schub = A.schubert_basis_on_vectors()
            sage: SF = SymmetricFunctions(QQ)
            sage: schur = SF.schur()
            sage: y = schur([3,2,1]) + schur([2,1])
            sage: Schub.from_schur(y)
            Y(1, 2, 3) + Y(0, 1, 2)
        
        """
        schur = list(schur)
        size = max([len(el[0]) for el in schur])
        return sum([el[1] * self.from_partition(el[0],size) for el in schur])
        
    class _divided_difference_wrapper(PolynomialRingWithBasis._divided_difference_wrapper):
        r"""
        This class is a wrapper for the on basis divided differences
        methods. 
        It contains a optimized version of the divided difference 
        for Schubert polynomials, it will be used instead of the
        default one on monomials.
                
        TESTS::
                sage: A = AbstractPolynomialRing(QQ)
                sage: Schub = A.schubert_basis_on_vectors()
                sage: Schub3 = Schub.finite_basis(3)
                sage: wrapp = Schub._divided_difference_wrapper(Schub3, 1)
        """
        def divided_difference_on_basis(self,key):
            r"""
            On basis action for the divided difference on the 
            Schubert basis.
            
            TESTS::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: Schub = A.schubert_basis_on_vectors()
                sage: Schub3 = Schub.finite_basis(3)
                sage: wrapp = Schub._divided_difference_wrapper(Schub3, 2)
                sage: pol = Schub[2,4,1]
                sage: key = list(pol)[0][0]
                sage: wrapp.divided_difference_on_basis(key)
                Y(2, 1, 3)
                sage: pol = Schub[2,1,4]
                sage: key = list(pol)[0][0]
                sage: wrapp.divided_difference_on_basis(key)
                0
                
            Test consistency..::
                sage: pol = Schub[2,4,1]
                sage: morph = Schub3._get_morphism_backup(1)
                sage: morph(pol) == pol.divided_difference(1)
                True
                sage: morph = Schub3._get_morphism_backup(2)
                sage: morph(pol) == pol.divided_difference(2)
                True
                sage: pol = Schub[2,2,1]
                sage: morph = Schub3._get_morphism_backup(1)
                sage: morph(pol) == pol.divided_difference(1)
                True
            """     
            i = self._i
            if(key[i-1] > key[i]):
                key2 = [key[j] for j in xrange(self._module.nb_variables())]
                key2[i-1], key2[i] = key2[i], key2[i-1]-1
                return self._module(key2)
            return self._module.zero()
                
                  


class GrothendieckPositiveBasisOnVectors(LinearBasisOnVectors):
    r"""
    Explain this class

    EXAMPLES::

        sage: # Fix a nice example

    TESTS::

        sage: # Fix a nice test
    """
    def __init__(self, abstract_polynomial_ring, ambient_space_basis, basis_name, basis_repr):
        LinearBasisOnVectors.__init__(
            self,
            abstract_polynomial_ring,
            ambient_space_basis,
            basis_name,
            basis_repr,
            on_basis_method = self.on_basis_method,
            variables_auto_coerce = True
        )
        
    def equivalent_basis(self, abstract_polynomial_ring):
        return abstract_polynomial_ring.grothendieck_positive_basis_on_vectors(self.group_type())

    def on_basis_method(self, x, basis, call_back):
        for i in xrange( len( x ) - 1 ):
            if( x[i]<x[i+1] ):
                x[i], x[i+1] = x[i+1]+1, x[i]
                res = call_back(x)
                fact = (basis.one() - basis.var(i+2))/basis.var(i+1)
                return (fact * res).isobaric_divided_difference(i+1)
        return basis(x)
 
            
class GrothendieckNegativeBasisOnVectors(LinearBasisOnVectors):
    r"""
    Explain this class

    EXAMPLES::

        sage: # Fix a nice example

    TESTS::

        sage: # Fix a nice test
    """
    def __init__(self, abstract_polynomial_ring, ambient_space_basis, basis_name, basis_repr):
        
        LinearBasisOnVectors.__init__(
            self,
            abstract_polynomial_ring,
            ambient_space_basis,
            basis_name,
            basis_repr,
            on_basis_method = self.on_basis_method,
            variables_auto_coerce = True,
            triangular = None
            
        )
            
    
            
    def equivalent_basis(self, abstract_polynomial_ring):
        return abstract_polynomial_ring.grothendieck_negative_basis_on_vectors(self.group_type())

    def on_basis_method(self, x, basis, call_back):
        for i in xrange( len( x ) - 1 ):
            if( x[i]<x[i+1] ):
                x[i], x[i+1] = x[i+1]+1, x[i]
                return call_back(x).isobaric_divided_difference(i+1)   
        prod = basis.one()
        for i in xrange(len(x)):
            inv_x_i = basis.var(i+1)**(-1)
            prod *= (basis.one() - inv_x_i)**x[i]
        return prod
        
    class _divided_difference_wrapper(PolynomialRingWithBasis._divided_difference_wrapper):
        r"""
        This class is a wrapper for the on basis isobaric divided differences
        methods. 
        It contains a optimized version of the divided difference 
        for Grothendieck polynomials.
                
        TESTS::
                sage: A = AbstractPolynomialRing(QQ)
                sage: Groth = A.grothendieck_negative_basis_on_vectors()
                sage: Groth3 = Groth.finite_basis(3)
                sage: wrapp = Groth._divided_difference_wrapper(Groth3, 1)

        """
        def isobaric_divided_difference_on_basis(self, key):
            r"""
            On basis action for the isobaric divided difference on the 
            Grothendieck basis.
            
            TESTS::
                sage: A = AbstractPolynomialRing(QQ)
                sage: Groth = A.grothendieck_negative_basis_on_vectors()
                sage: Groth3 = Groth.finite_basis(3)
                sage: wrapp = Groth._divided_difference_wrapper(Groth3, 1)
                sage: pol = Groth[3,2,3]
                sage: key = list(pol)[0][0]               
                sage: wrapp.isobaric_divided_difference_on_basis(key)
                G(2, 2, 3)
                sage: pol = Groth[2,3,3]
                sage: key = list(pol)[0][0]
                sage: wrapp.isobaric_divided_difference_on_basis(key)
                G(2, 3, 3)
                
            Test consistency..::
                
                sage: pol = Groth[3,2,3]
                sage: a = pol.isobaric_divided_difference(1).expand()
                sage: b = pol.expand().isobaric_divided_difference(1)
                sage: a == b
                True
                sage: a = pol.isobaric_divided_difference(2).expand()
                sage: b = pol.expand().isobaric_divided_difference(2)
                sage: a == b
                True
                sage: pol = Groth[2,3,3]
                sage: a = pol.isobaric_divided_difference(1).expand()
                sage: b = pol.expand().isobaric_divided_difference(1)
                sage: a == b
                True
                sage: a = pol.isobaric_divided_difference(2).expand()
                sage: b = pol.expand().isobaric_divided_difference(2)
                sage: a == b
                True
            """
            i = self._i
            if(key[i-1] > key[i]):
                key2 = [key[j] for j in xrange(self._module.nb_variables())]
                key2[i-1], key2[i] = key2[i], key2[i-1]-1
                return self._module(key2)
            else:
                return self._module(key)
                
        def hat_isobaric_divided_difference_on_basis(self, key):
            r"""
            On basis action for the hat isobaric divided difference on the 
            Grothendieck basis.
            
            TESTS::
                sage: A = AbstractPolynomialRing(QQ)
                sage: Groth = A.grothendieck_negative_basis_on_vectors()
                sage: Groth3 = Groth.finite_basis(3)
                sage: wrapp = Groth._divided_difference_wrapper(Groth3, 1)
                sage: pol = Groth[3,2,3]
                sage: key = list(pol)[0][0]               
                sage: wrapp.hat_isobaric_divided_difference_on_basis(key)
                G(2, 2, 3) - G(3, 2, 3)
                sage: pol = Groth[2,3,3]
                sage: key = list(pol)[0][0]
                sage: wrapp.hat_isobaric_divided_difference_on_basis(key)
                0
                
            Test consistency..::
                
                sage: pol = Groth[3,2,3]
                sage: a = pol.hat_isobaric_divided_difference(1).expand()
                sage: b = pol.expand().hat_isobaric_divided_difference(1)
                sage: a == b
                True
                sage: a = pol.hat_isobaric_divided_difference(2).expand()
                sage: b = pol.expand().hat_isobaric_divided_difference(2)
                sage: a == b
                True
                sage: pol = Groth[2,3,3]
                sage: a = pol.hat_isobaric_divided_difference(1).expand()
                sage: b = pol.expand().hat_isobaric_divided_difference(1)
                sage: a == b
                True
                sage: a = pol.hat_isobaric_divided_difference(2).expand()
                sage: b = pol.expand().hat_isobaric_divided_difference(2)
                sage: a == b
                True
            """
            return self.isobaric_divided_difference_on_basis(key) - self._module(key)
            
class MacdonaldBasisOnVectors(LinearBasisOnVectors):
    r"""
    Explain this class

    EXAMPLES::

        sage: # Fix a nice example

    TESTS::

        sage: # Fix a nice test
    """
    
    def __init__(self, abstract_polynomial_ring, ambient_space_basis, basis_name, basis_repr, t1, t2, q):
        r"""
        TESTS::

            sage: # Fix a nice test
        """
        extra_parameters = []
        extra_parameters.append(("t1",t1))
        extra_parameters.append(("t2",t2))
        extra_parameters.append(("q",q))
        extra_parameters = tuple(extra_parameters)
        LinearBasisOnVectors.__init__(
            self,
            abstract_polynomial_ring,
            ambient_space_basis,
            basis_name,
            basis_repr,
            on_basis_method = self.on_basis_method,
            extra_parameters = extra_parameters,
            cmp = self.cmp
        )

        
    def equivalent_basis(self, abstract_polynomial_ring):
        return abstract_polynomial_ring.macdonald_basis_on_vectors(self.group_type(), self._t1, self._t2, self._q)

    def cmp(self, key1, key2):
        l = len(key1.parent()._basis_keys)
            
        d1 = sum( [key1[i] for i in xrange(l) ] )
        d2 = sum( [key2[i] for i in xrange(l) ] )
        if(d1>d2): return 1
        if(d1<d2): return -1
        
        v1 = [key1[i] for i in xrange(l)]
        v2 = [key2[i] for i in xrange(l)]
        v1.sort()
        v2.sort()
        
        for i in xrange(l-1,-1,-1):
            if(v1[i] > v2[i]):
                return 1
            if(v1[i] < v2[i]):
                return -1
        
        for i in xrange(l):
            if (key1[i]>key2[i]):
                return 1
            if (key1[i]<key2[i]):
                return -1
        return 0

    
    def on_basis_method(self, u, basis, call_back, t1=1, t2=1, q=1):
        size = len(u)
        if(u[size-1] != 0):
            temp = u[size-1]
            for i in xrange(size-1,0,-1): u[i] = u[i-1]
            u[0]  = temp-1
            
            res = call_back(u)
                
            p = [i+2 for i in xrange(size)]
            p[size-1] = 1
            p = Permutation( p )
            res = res.perm_vars(p)
            replace = basis.var(size) * q**(-1)
            res = res.subs_var((size,replace))
            res*= basis.var(size) - (-t2)**(size-1)
            return res
            
        for i in xrange(size-1):
            if(u[i] > u[i+1]):
                u[i], u[i+1] = u[i+1], u[i]
                res = call_back(u)

                res1 = res.hecke_generator(i+1,t1, t2)
                w = Word( [u[j] for j in xrange(size-1,-1,-1)  ] )
                p = w.standard_permutation()
                u_spec = [ q**(u[j])*(-t1/t2)**(p(size-j) -1) for j in xrange(size) ]
                res2 = res * (t1  + t2)/(u_spec[i+1]/u_spec[i] -1) 
                return res1 + res2
                    
            
        return basis.one()
           

class DemazureBasisOnVectors(LinearBasisOnVectors):
    r"""
    Explain

    EXAMPLES::

        sage: A = AbstractPolynomialRing(QQ)
        sage: Dem = A.demazure_basis_on_vectors("B")
        sage: Dem( Dem.an_element().expand()) == Dem.an_element()
        True
        sage: Dem( Dem[2,-2,-1].expand() ) == Dem[2,-2,-1]
        True
        sage: Dem = A.demazure_basis_on_vectors("C")
        sage: Dem( Dem.an_element().expand()) == Dem.an_element()
        True
        sage: Dem( Dem[2,-2,-1].expand() ) == Dem[2,-2,-1]
        True
        sage: Dem = A.demazure_basis_on_vectors("D")
        sage: Dem( Dem.an_element().expand()) == Dem.an_element()
        True
        sage: Dem( Dem[2,-2,-1].expand() ) == Dem[2,-2,-1]
        True

    TESTS::

        sage: # Fix a nice test
    """    
    def __init__(self, abstract_polynomial_ring, ambient_space_basis, basis_name, basis_repr, method):
        r"""
        TESTS::

            sage: # Fix a nice test
        """
        extra_parameters = []
        extra_parameters.append(("method", method))
        extra_parameters = tuple(extra_parameters)
        self._method = method
        if(ambient_space_basis.group_type() != "A"):
            variables_auto_coerce = False
        else: 
            variables_auto_coerce = True
        LinearBasisOnVectors.__init__(
            self,
            abstract_polynomial_ring,
            ambient_space_basis,
            basis_name,
            basis_repr,
            on_basis_method = self.on_basis_method,
            variables_auto_coerce = variables_auto_coerce,
            extra_parameters = extra_parameters
        )
        if (self._method == "isobaric_divided_difference"):
            self._divided_difference_wrapper = self._divided_difference_wrapper_K
        elif(self._method == "hat_isobaric_divided_difference"):
            self._divided_difference_wrapper = self._divided_difference_wrapper_HatK
        

        
    def equivalent_basis(self, abstract_polynomial_ring):
        return abstract_polynomial_ring.demazure_basis_on_vectors(self.group_type(), method= self._method)

    def on_basis_method(self, x, basis, call_back, method = "isobaric_divided_difference"):
        otype = basis.group_type()
        if(otype=="D"):
            for i in xrange( len(x) - 1):
                if( x[i]<x[i+1] ):
                    x[i], x[i+1] = x[i+1], x[i]
                    return getattr(call_back(x), method)(i+1)
            i = len(x) - 2
            if(x[i] + x[i+1] <0):
                x[i], x[i+1] = -x[i+1], -x[i]
                return getattr(call_back(x),method)(i+2)
            return basis(x)    
        else:    
            for i in xrange( len( x ) - 1 ):
                if( x[i]<x[i+1] ):
                    x[i], x[i+1] = x[i+1], x[i]
                    return getattr(call_back(x), method)(i+1)   
            i = len(x) - 1
            if(x[i]<0):
                x[i] = -x[i]
                return getattr(call_back(x), method)(i+1)
            return basis(x)   
            

        
        
    class _divided_difference_wrapper_K(PolynomialRingWithBasis._divided_difference_wrapper):
        r"""
        This class is a wrapper for the on basis isobaric divided differences
        methods. 
        It contains a optimized version of the isobaric divided differences
        for K polynomials (Demazre basis), they will be used instead of the
        default ones on monomials.
        
        TESTS::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: K = A.demazure_basis_on_vectors()
                sage: K3 = K.finite_basis(3)
                sage: wrapp = K._divided_difference_wrapper_K(K3,1)
        
        """
        
        def isobaric_divided_difference_on_basis(self, key):
            r"""
            On basis action for the isobaric divided difference on the 
            demazure basis.
            
            TESTS::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: K = A.demazure_basis_on_vectors()
                sage: K3 = K.finite_basis(3)
                sage: wrapp = K._divided_difference_wrapper_K(K3,1)
                sage: pol = K[2,4,1]
                sage: key = list(pol)[0][0]
                sage: wrapp.isobaric_divided_difference_on_basis(key)
                K(2, 4, 1)
                sage: pol = K[4,2,1]
                sage: key = list(pol)[0][0]
                sage: wrapp.isobaric_divided_difference_on_basis(key)
                K(2, 4, 1)
            
            Test consistency..::
                sage: pol = K[2,4,1]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: pol = K[2,4,4]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: K = A.demazure_basis_on_vectors("B")
                sage: K3 = K.finite_basis(3)
                sage: pol = K[2,4,-1]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = K[2,-4,1]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = K[-2,-4,1]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = K[2,-4,-4]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: K = A.demazure_basis_on_vectors("C")
                sage: K3 = K.finite_basis(3)
                sage: pol = K[2,4,-1]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = K[2,-4,1]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = K[-2,-4,1]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = K[2,-4,-4]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: K = A.demazure_basis_on_vectors("D")
                sage: K3 = K.finite_basis(3)
                sage: pol = K[2,4,-1]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = K[2,-4,1]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = K[-2,-4,1]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = K[2,-4,-4]
                sage: morph = K3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
            """
            key2 = key.weyl_action([self._i])
            if(self._module.basis_tower().cmp(key,key2) < 0):
                return self._module(key2)
            else:
                return self._module(key)
                
        def hat_isobaric_divided_difference_on_basis(self, key):
            r"""
            On basis action for the hat isobaric divided difference on the 
            demazure basis.
            
            TESTS::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: K = A.demazure_basis_on_vectors()
                sage: K3 = K.finite_basis(3)
                sage: wrapp = K._divided_difference_wrapper_K(K3,1)
                sage: pol = K[2,4,1]
                sage: key = list(pol)[0][0]
                sage: wrapp.hat_isobaric_divided_difference_on_basis(key)
                0
                sage: pol = K[4,2,1]
                sage: key = list(pol)[0][0]
                sage: wrapp.hat_isobaric_divided_difference_on_basis(key)
                K(2, 4, 1) - K(4, 2, 1)
                
            Test consistency..::
                sage: pol = K[2,4,1]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: pol = K[2,4,4]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: K = A.demazure_basis_on_vectors("B")
                sage: K3 = K.finite_basis(3)
                sage: pol = K[2,4,-1]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = K[2,-4,1]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = K[-2,-4,1]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = K[2,-4,-4]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: K = A.demazure_basis_on_vectors("C")
                sage: K3 = K.finite_basis(3)
                sage: pol = K[2,4,-1]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = K[2,-4,1]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = K[-2,-4,1]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = K[2,-4,-4]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: K = A.demazure_basis_on_vectors("D")
                sage: K3 = K.finite_basis(3)
                sage: pol = K[2,4,-1]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = K[2,-4,1]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = K[-2,-4,1]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = K[2,-4,-4]
                sage: morph = K3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = K3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = K3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
            """
            return self.isobaric_divided_difference_on_basis(key) - self._module(key)

    class _divided_difference_wrapper_HatK(PolynomialRingWithBasis._divided_difference_wrapper):
        r"""
        This class is a wrapper for the on basis isobaric divided differences
        methods. 
        It contains a optimized version of the isobaric divided differences
        for hat K polynomials (hat Demazre basis), they will be used instead of the
        default ones on monomials.
        
        
        TESTS::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: HatK = A.demazure_hat_basis_on_vectors()
                sage: HatK3 = HatK.finite_basis(3)
                sage: wrapp = HatK._divided_difference_wrapper_HatK(HatK3,1)
        """
        def hat_isobaric_divided_difference_on_basis(self, key):
            r"""
            On basis action for the hat isobaric divided difference on the 
            hat demazure basis.
            
            TESTS::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: HatK = A.demazure_hat_basis_on_vectors()
                sage: HatK3 = HatK.finite_basis(3)
                sage: wrapp = HatK._divided_difference_wrapper_HatK(HatK3,1)
                sage: pol = HatK[2,4,1]
                sage: key = list(pol)[0][0]
                sage: wrapp.hat_isobaric_divided_difference_on_basis(key)
                -^K(2, 4, 1)
                sage: pol = HatK[4,2,1]
                sage: key = list(pol)[0][0]
                sage: wrapp.hat_isobaric_divided_difference_on_basis(key)
                ^K(2, 4, 1)
                
            Test consistency..::
                sage: pol = HatK[2,4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: pol = HatK[1,4,4]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: HatK = A.demazure_basis_on_vectors("B")
                sage: HatK3 = HatK.finite_basis(3)
                sage: pol = HatK[2,4,-1]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = HatK[2,-4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = HatK[-2,-4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = HatK[2,-4,-4]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: HatK = A.demazure_basis_on_vectors("C")
                sage: HatK3 = HatK.finite_basis(3)
                sage: pol = HatK[2,4,-1]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = HatK[2,-4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = HatK[-2,-4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = HatK[2,-4,-4]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: HatK = A.demazure_basis_on_vectors("D")
                sage: HatK3 = HatK.finite_basis(3)
                sage: pol = HatK[2,4,-1]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = HatK[2,-4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = HatK[-2,-4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
                sage: pol = HatK[2,-4,-4]
                sage: morph = HatK3._get_morphism_backup(1,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="hatpi")
                sage: morph(pol) == pol.hat_isobaric_divided_difference(3)
                True
            """
            key2 = key.weyl_action([self._i])
            c = self._module.basis_tower().cmp(key,key2)
            if(c < 0):
                return self._module(key2)
            elif(c ==0):
                return self._module.zero()
            else:
                return -self._module(key)
                
        def isobaric_divided_difference_on_basis(self, key):
            r"""
            On basis action for the hat isobaric divided difference on the 
            hat demazure basis.
            
            TESTS::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: HatK = A.demazure_hat_basis_on_vectors()
                sage: HatK3 = HatK.finite_basis(3)
                sage: wrapp = HatK._divided_difference_wrapper_HatK(HatK3,1)
                sage: pol = HatK[2,4,1]
                sage: key = list(pol)[0][0]
                sage: wrapp.isobaric_divided_difference_on_basis(key)
                0
                sage: pol = HatK[4,2,1]
                sage: key = list(pol)[0][0]
                sage: wrapp.isobaric_divided_difference_on_basis(key)
                ^K(2, 4, 1) + ^K(4, 2, 1)
                
            Test consistency..::
                sage: pol = HatK[2,4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: pol = HatK[1,4,4]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: HatK = A.demazure_basis_on_vectors("B")
                sage: HatK3 = HatK.finite_basis(3)
                sage: pol = HatK[2,4,-1]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = HatK[2,-4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = HatK[-2,-4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = HatK[2,-4,-4]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: HatK = A.demazure_basis_on_vectors("C")
                sage: HatK3 = HatK.finite_basis(3)
                sage: pol = HatK[2,4,-1]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = HatK[2,-4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = HatK[-2,-4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = HatK[2,-4,-4]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: HatK = A.demazure_basis_on_vectors("D")
                sage: HatK3 = HatK.finite_basis(3)
                sage: pol = HatK[2,4,-1]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = HatK[2,-4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = HatK[-2,-4,1]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True
                sage: pol = HatK[2,-4,-4]
                sage: morph = HatK3._get_morphism_backup(1,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(1)
                True
                sage: morph = HatK3._get_morphism_backup(2,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(2)
                True
                sage: morph = HatK3._get_morphism_backup(3,method="pi")
                sage: morph(pol) == pol.isobaric_divided_difference(3)
                True

            """
            
            return self.hat_isobaric_divided_difference_on_basis(key) + self._module(key)
            
            
                

class FiniteLinearBasisOnVectors(FinitePolynomialRingWithBasisFromMorphism):
    r"""
    Explains this class

    EXAMPLES::

        sage: # Fix a nice example

    TESTS::

        sage: # Fix a nice test
    """
    def __init__(self, abstract_polynomial_ring, polynomial_ring_tower, basis_name, basis_repr, **keywords):
        r"""
        TESTS::

            sage: # Fix a nice test
        """
                 
        
        FinitePolynomialRingWithBasisFromMorphism.__init__(
            self,
            abstract_polynomial_ring,
            polynomial_ring_tower,
            basis_name,
            basis_repr,
            **keywords
        )  
        self._ambient_space_basis = self._morphism_to_basis

    def ambient_space_basis(self):
        r"""
        Returns ...

        EXAMPLES::

            sage: # Fix a nice example
        """
        return self._ambient_space_basis

    def one_basis(self):
        r"""
        Returns the index of the one of ``self``.

        EXAMPLES::

            sage: # Fix a nice example
        """
        return self.basis().keys().zero() 

    def __getitem__(self, c, *rest):
        r"""
        TESTS::

            sage: # Fix a nice test
        """
        if len(rest) > 0 or type(c) is int or type(c) is Integer:
            c = tuple([c])+tuple(rest)
        return self.term( self._basis_keys( list(c) ) )    
        
    def __call__(self, obj):
        r"""
        TESTS::

            sage: # Fix a nice test
        """
        if( type(obj) is list ):
            return self.term( self._basis_keys( obj))
        else:
            return super(FiniteLinearBasisOnVectors, self).__call__(obj)
            
    def morphism_method_wrapper(self):
        return self.basis_tower().morphism_method_wrapper(self.nb_variables())
    
    def group_type(self):
        r"""
            Returns the basis group type
            
            EXAMPLES::
            sage: A = AbstractPolynomialRing(QQ)                                         
            sage: K = A.demazure_basis_on_vectors()
            sage: pol = K.an_element()
            sage: pol.parent().group_type()
            'A'
            sage: KB = A.demazure_basis_on_vectors("B")
            sage: pol = KB.an_element()
            sage: pol.parent().group_type()
            'B'

        """
        return self.basis_tower().group_type()
