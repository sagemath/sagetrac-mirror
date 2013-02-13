r"""
Multivariate Polynomials With Several Bases on two sets of variables
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Viviane Pons < pons at univ-mlv.fr  >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from multivariate_polynomials import AbstractPolynomialRing
from basis import Bases, Finite_bases, PolynomialRingWithBasis
from linear_basis_on_vectors import LinearBasisOnVectors, FiniteLinearBasisOnVectors

class DoubleAbstractPolynomialRing(AbstractPolynomialRing):
    """r
    This class implements the ring of abstract polynomials over a double set of variables
    This ring is actually an ``AbstractPolynomialRing`` over antoher ``AbstractPolynomialRing``
    but the class provides methods that are specific to double multivariate polynomials to make it easy 
    to use
    
    INPUT:
        - ``R``: the base ring of the algebra
        - ``repr_var1``, a string representing the main variable set, default is `x`
        - ``repr_var2``, a string representing the secondary variable set, default is `y`
        - ``inversed_ring``, the ring where the roles of the two sets of variables are inversed. 
        By default, nothing is sent and the ring is created, the field is then used to avoid infinite
        recursion
            
    OUTPUT:

        - The abstract ring of multivariate polynomials on ``repr_var1`` over te abstract ring
        of multivariate polynomials on ``repr_var2`` over ``R``
        
    EXAMPLES::
    
        sage: D = DoubleAbstractPolynomialRing(QQ);D
        The abstract ring of multivariate polynomials on x over The abstract ring of multivariate polynomials on y over Rational Field
        sage: D.an_element()
        y[]*x[1, 2, 3]

    All the bases and actions existing in ``AbstractPolynomialRing`` are available, they correpond to the 
    bases and actions on the `x` set of variables
    
        sage: M = D.monomial_basis(); M
        The ring of multivariate polynomials on x over The abstract ring of multivariate polynomials on y over Rational Field on the monomial basis
        sage: pol = M[1,2,3] + M[2,2,4]; pol
        y[]*x[1, 2, 3] + y[]*x[2, 2, 4]
        sage: pol.divided_difference(1)
        (-y[])*x[1, 1, 3]
        
    You can obtain the ring on the `y` and all the bases for the `y` polynomials
    
        sage: Y = D.coeffs_ring(); Y
        The abstract ring of multivariate polynomials on y over Rational Field
        sage: MY = Y.monomial_basis(); MY
        The ring of multivariate polynomials on y over Rational Field on the monomial basis
        sage: MY.an_element()
        y[1, 2, 3]
        sage: MY.an_element() * pol
        (y[1,2,3])*x[1, 2, 3] + (y[1,2,3])*x[2, 2, 4]
        
    You can change the bases of the `x` polynomial by the usual coercion system
    
        sage: Schub = D.schubert_basis_on_vectors(); Schub
        The ring of multivariate polynomials on x over The abstract ring of multivariate polynomials on y over Rational Field on the Schubert basis of type A (indexed by vectors)
        sage: Schub(pol)
        y[]*Y(1, 2, 3) + (-y[])*Y(1, 3, 2) + (-y[])*Y(2, 1, 3) + y[]*Y(2, 2, 4) + y[]*Y(2, 3, 1) + (-y[])*Y(2, 3, 3) + (-y[])*Y(2, 4, 2) + y[]*Y(3, 1, 2) + (-y[])*Y(3, 2, 1) + y[]*Y(3, 3, 2) + y[]*Y(4, 1, 1)

    Be carreful as the upper case `Y` here doesn't correspond to the set of variables but reprensents the
    Schubert polynomials. 
    You can also change the base for the `y`
    
        sage: YS = Y.schubert_basis_on_vectors(); YS
        The ring of multivariate polynomials on y over Rational Field on the Schubert basis of type A (indexed by vectors)
        sage: pol = MY.an_element() * pol; pol
        (y[1,2,3])*x[1, 2, 3] + (y[1,2,3])*x[2, 2, 4]
        sage: pol.change_coeffs_bases(YS)
        (Y(1,2,3)-Y(1,3,2)-Y(2,1,3)+Y(2,3,1)+Y(3,1,2)-Y(3,2,1)+Y(4,1,1))*x[1, 2, 3] + (Y(1,2,3)-Y(1,3,2)-Y(2,1,3)+Y(2,3,1)+Y(3,1,2)-Y(3,2,1)+Y(4,1,1))*x[2, 2, 4]

    Also, you can obtain the ring where the roles of `x` and `y` are exchanged
    
        sage: D2 = D.inversed_ring(); D2
        The abstract ring of multivariate polynomials on y over The abstract ring of multivariate polynomials on x over Rational Field
        sage: D2.an_element()
        x[0]*y[1, 2, 3]
        
    There is a coercion between `D` and `D2`
    
        sage: pol
        (y[1,2,3])*x[1, 2, 3] + (y[1,2,3])*x[2, 2, 4]
        sage: D2(pol)
        (x[1,2,3]+x[2,2,4])*y[1, 2, 3]
        
    But this coercion doesn't allow operations including polynomials from `D` and `D2` as the coercion 
    only exists between abstract polynomial rings but not between the concrete bases : sage coercion system
    doesn't look for a parent where the coercion could be made

    TESTS::

        sage: # WARNING: Fix _test_one and _test_zero when facade feature will be ready
        ssage: D = DoubleAbstractPolynomialRing(QQ);
        sage: TestSuite(D).run(skip=['_test_one','_test_zero'])
    


    """
    
    def __init__(self, R, repr_var1 = 'x', repr_var2 = 'y', inversed_ring = None):
        self._coeffs_ring = AbstractPolynomialRing(R, repr_var2, always_show_main_var = True)
        AbstractPolynomialRing.__init__(
            self,
            self._coeffs_ring,
            repr_var1,
            finite_bases_category_class = Finite_double_bases,
            always_show_main_var = True
        )
        self._repr_var1 = repr_var1
        self._repr_var2 = repr_var2
        self._coeffs_base_ring = R
        if(inversed_ring is None):
            self._inversed_ring = DoubleAbstractPolynomialRing(R, repr_var2, repr_var1, self)
        else:
            self._inversed_ring = inversed_ring
        
        m = SetMorphism( Hom(self, self.inversed_ring()), lambda x : x.swap_coeffs_elements())
        m.register_as_coercion()
       
   
    
    def coeffs_ring(self):
        """r
        returns the multivariate polynomial ring on the second set of variables
         used as coefficients of the main ring on the first set of variables
         
        OUPUT :
            - the ring of multivariate polynomials on the second set of variables 
         
        EXAMPLES::
        
            sage: D = DoubleAbstractPolynomialRing(QQ);D
            The abstract ring of multivariate polynomials on x over The abstract ring of multivariate polynomials on y over Rational Field
            sage: Y = D.coeffs_ring(); Y
            The abstract ring of multivariate polynomials on y over Rational Field
        
        """
        return self._coeffs_ring
        
    def inversed_ring(self):
        """r
        returns the ring of multivariate polynomials where the roles of the two sets of variables 
        are exchanged
        
        OUTPUT:
            - the ring of multivariate polynomials on the second set of variables 
            over the ring of multivariate polynomials on the first set of variables
        
        EXAMPLES::
            sage: D = DoubleAbstractPolynomialRing(QQ);D
            The abstract ring of multivariate polynomials on x over The abstract ring of multivariate polynomials on y over Rational Field
            sage: D2 = D.inversed_ring(); D2
            The abstract ring of multivariate polynomials on y over The abstract ring of multivariate polynomials on x over Rational Field
        
            There is a coercion between `D` and `D2`
        
            sage: M = D.monomial_basis()
            sage: pol = D.coeffs_ring().an_element() * (M[1,2,3] + M[2,2,4]);pol
            (y[1,2,3])*x[1, 2, 3] + (y[1,2,3])*x[2, 2, 4]
            sage: D2(pol)
            (x[1,2,3]+x[2,2,4])*y[1, 2, 3]
            
        But this coercion doesn't allow for operations including polynomials from both `D` and `D2` as the coercion 
        only exists between abstract polynomial rings but not between the concrete bases : sage coercion system
        doesn't look for a parent where the coercion could be made
        """
        return self._inversed_ring
        
    def double_schubert_basis_on_vectors(self, group_type = "A", basis_name = None, basis_repr = "YY"):
        if(basis_name is None):
            basis_name = "Double Schubert basis of type " +group_type+" (indexed by vectors)"  
        ambient_space_basis = self.ambient_space_basis(group_type)
        return DoubleSchubertBasisOnVectors(self, ambient_space_basis, basis_name, basis_repr)     
        
    def double_grothendieck_basis_on_vectors(self, group_type = "A", basis_name = None, basis_repr = "GG"):
        if(basis_name is None):
            basis_name = "Double Grothendieck basis of type " +group_type+" (indexed by vectors)"  
        ambient_space_basis = self.ambient_space_basis(group_type)
        return DoubleGrothendieckBasisOnVectors(self, ambient_space_basis, basis_name, basis_repr)  
        
class Finite_double_bases(Finite_bases):
    """r
    This category is a sub category of ``sage.combinat.multivariate_polynomials.basis.Finite_bases`` and
    is used for finite bases of ``DoubleAbstractPolynomialRing``, i.e. finite bases for polynomials on 
    two sets of variables. It contains element methds that are specific to these bases ant not to the 
    multivariate polynomials on one set of variables
    
    INPUT:
        - ``abstract_algebra``, the finite abstract algebra of which ``self`` is a base
        - ``algebra_tower``, the base of :class:``AbstractPolynomialRing`` that correspond to 
        this base on an unset number of variables 
        
    OUTPUT:
        the category of finite_double_bases of the abstract ring of multivariate polynomials  with `nb_variables` 
        variables
    
    EXAMPLES::
        sage: D = DoubleAbstractPolynomialRing(QQ)  
        sage: pol = D.an_element(); pol
        y[]*x[1, 2, 3]
        sage: pol.parent().category()
        The category of bases of The abstract ring of multivariate polynomials on x over The abstract ring of multivariate polynomials on y over Rational Field with 3 variables where algebra tower is The ring of multivariate polynomials on x over The abstract ring of multivariate polynomials on y over Rational Field on the monomial basis

    TESTS::
    
        sage: D = DoubleAbstractPolynomialRing(QQ)
        sage: F3 = D.finite_polynomial_ring(3)
        sage: M = D.monomial_basis()
        sage: from sage.combinat.multivariate_polynomials.double_multivariate_polynomials import Finite_double_bases
        sage: C = Finite_double_bases(F3,M)
        sage: TestSuite(C).run()
    """
    
    
    def super_categories(self):
        return [Finite_bases(self.base(), self.algebra_tower())]
    
    class ElementMethods:
        
        def change_coeffs_bases(self, new_base):
            """r
            This method changes the base of the coefficients of a given polynomial which are polynomials
            on the second set of variables
            
            INPUT:
                - ``new_base`` : a base of :class:`AbstractPolynomialRing` on the second set of variables
                
            OUTPUT: 
                - the polynomial where the base of the coefficients has been changed
                
            EXAMPLES::
            
                sage: D = DoubleAbstractPolynomialRing(QQ)
                sage: Y = D.coeffs_ring()
                sage: YM = Y.monomial_basis()
                sage: YS = Y.schubert_basis_on_vectors()
                sage: pol = YM.an_element() * D.an_element() + D.one(); pol
                y[0]*x[0, 0, 0] + (y[1,2,3])*x[1, 2, 3]
                sage: pol.change_coeffs_bases(YS)
                Y(0)*x[0, 0, 0] + (Y(1,2,3)-Y(1,3,2)-Y(2,1,3)+Y(2,3,1)+Y(3,1,2)-Y(3,2,1)+Y(4,1,1))*x[1, 2, 3]

            """
            return sum( [ self.parent().term(x,new_base(y)) for (x,y) in self ] )
            
        def inversed_ring(self):
            return self.parent().abstract_algebra().polynomial_ring_tower().inversed_ring()
            
        def swap_coeffs_elements(self):
            inversed_ring = self.inversed_ring()
            coeffs_ring = inversed_ring.coeffs_ring()
            base_x = self.parent().basis_tower().equivalent_basis(coeffs_ring)
            base_y = None
            default_base_y = None
            result = None
            for (x,y) in self:
                if(base_y is None):
                    default_base_y = y.parent().basis_tower()
                    base_y = default_base_y.equivalent_basis(inversed_ring)
                    result = base_y.zero()
                y = default_base_y( y )
                for (key, coeff) in y:
                    result += base_y.term( key, base_x.term(x,coeff))
            return result            
            


class DoubleSchubertBasisOnVectors(LinearBasisOnVectors):
    
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
            cmp = self.cmp
        )
   
    def cmp(self, key1, key2):
        l = len(key1.parent()._basis_keys)
        for i in xrange(l-1,-1,-1):
            if (key1[i]>key2[i]):
                return 1
            if (key1[i]<key2[i]):
                return -1
        return 0    

    def on_basis_method(self, x, basis, call_back):
        for i in xrange( len( x ) - 1 ):
            if( x[i]<x[i+1] ):
                x[i], x[i+1] = x[i+1]+1, x[i]
                return call_back(x).divided_difference(i+1)   
        basex = basis
        basey = basex.basis_tower().equivalent_basis(basex.base_ring())
        return prod( [basex.var(i+1) - basey.var(j+1) for i in xrange(len(x)) for j in xrange(x[i])], basex.one())
        
        

class DoubleGrothendieckBasisOnVectors(LinearBasisOnVectors):
    
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
        

    def on_basis_method(self, x, basis, call_back):
        for i in xrange( len( x ) - 1 ):
            if( x[i]<x[i+1] ):
                x[i], x[i+1] = x[i+1]+1, x[i]
                return call_back(x).isobaric_divided_difference(i+1)   
        prod = basis.one()
        basex = basis
        basey = basex.basis_tower().equivalent_basis(basex.base_ring())
        for i in xrange(len(x)):
            inv_x_i = basex.var(i+1)**(-1)
            for j in xrange(x[i]):
                prod *= (basis.one() - basey.var(j+1) * inv_x_i)
        return prod
        
    class _divided_difference_wrapper(PolynomialRingWithBasis._divided_difference_wrapper):
        r"""
        This class is a wrapper for the on basis isobaric divided differences
        methods. 
        It contains a optimized version of the divided difference 
        for Grothendieck polynomials.
                
        TESTS::
            sage: D = DoubleAbstractPolynomialRing(QQ)
            sage: DG = D.double_grothendieck_basis_on_vectors()
            sage: DG3 = DG.finite_basis(3)
            sage: wrapp = DG._divided_difference_wrapper(DG3,1)


        """
        def isobaric_divided_difference_on_basis(self, key):
            r"""
            On basis action for the isobaric divided difference on the 
            Grothendieck basis.
            
            TESTS::
                sage: D = DoubleAbstractPolynomialRing(QQ)
                sage: DG = D.double_grothendieck_basis_on_vectors()
                sage: DG3 = DG.finite_basis(3)
                sage: wrapp = DG._divided_difference_wrapper(DG3,1)
                sage: pol = DG[3,2,3]
                sage: key = list(pol)[0][0]               
                sage: wrapp.isobaric_divided_difference_on_basis(key)
                y[0]*GG(2, 2, 3)
                sage: pol = DG[2,3,3]
                sage: key = list(pol)[0][0]
                sage: wrapp.isobaric_divided_difference_on_basis(key)
                y[0]*GG(2, 3, 3)
                
            Test consistency..::
                
                sage: pol = DG[3,2,3]
                sage: a = pol.isobaric_divided_difference(1).expand()
                sage: b = pol.expand().isobaric_divided_difference(1)
                sage: a == b
                True
                sage: a = pol.isobaric_divided_difference(2).expand()
                sage: b = pol.expand().isobaric_divided_difference(2)
                sage: a == b
                True
                sage: pol = DG[2,3,3]
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
                sage: pol = DG[3,2,3]
                sage: key = list(pol)[0][0]               
                sage: wrapp.hat_isobaric_divided_difference_on_basis(key)
                y[0]*GG(2, 2, 3) + (-y[0])*GG(3, 2, 3)
                sage: pol = DG[2,3,3]
                sage: key = list(pol)[0][0]
                sage: wrapp.hat_isobaric_divided_difference_on_basis(key)
                0
                
            Test consistency..::
                
                sage: pol = DG[3,2,3]
                sage: a = pol.hat_isobaric_divided_difference(1).expand()
                sage: b = pol.expand().hat_isobaric_divided_difference(1)
                sage: a == b
                True
                sage: a = pol.hat_isobaric_divided_difference(2).expand()
                sage: b = pol.expand().hat_isobaric_divided_difference(2)
                sage: a == b
                True
                sage: pol = DG[2,3,3]
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
        
