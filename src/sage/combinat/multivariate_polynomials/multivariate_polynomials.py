r"""
Multivariate Polynomials With Several Bases
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Viviane Pons < pons at univ-mlv.fr  >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.categories.all import Category
from sage.categories.all import GradedModules, CommutativeAlgebras, Realizations, GradedAlgebras
from sage.categories.all import Rings
from sage.categories.category_types import Category_over_base_ring
from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.all import var

# Temporary!!!!
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


class Facade_polynomial_ring(GradedAlgebras):
    r"""
    This category is used only for (see :class:``AbstractPolynomialRing``)
    it is both a facade for its realizations ``PolynomialRingWithBasis`` but also
    for the ``FiniteAbstratPolynomialRing``. The ``WithRealizations`` Method adds the
    facades categories and methods. We only have to override the parent method
     ``facade_for`` for it to get both realizations and non realizations.

    EXAMPLES::
    
        sage: A = AbstractPolynomialRing(ZZ)
        sage: m = A.monomial_basis()
        sage: m3 = m.finite_basis(3)
        sage: A.realizations()
        [The ring of multivariate polynomials on x over Integer Ring on the monomial basis]
        sage: A.facade_for()
        [The ring of multivariate polynomials on x over Integer Ring on the monomial basis, The abstract ring of multivariate polynomials on x over Integer Ring with 3 variables, The abstract ring of multivariate polynomials on x over Integer Ring with 1 variable, The ring of multivariate polynomials on x over Integer Ring with 3 variables on the monomial basis]


    Actual elements from ``A`` are not directly elements 
    from its realizations but elements from realizations 
    of ``FinitePolynomialRing``, we have both these shemas:

    A -> facade for F3 -> realization FiniteMonomialBasis
    A -> realization MonomialBasis -> facade for FintieMonomialBasis

    TESTS::

        sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import Facade_polynomial_ring
        sage: C = Facade_polynomial_ring(QQ)
        sage: TestSuite(C).run()
    """
    def super_categories(self):
        r"""
        
        OUTPUT:
            - a list of the super categories of ``self``
            
        EXAMPLES::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import Facade_polynomial_ring
            sage: C = Facade_polynomial_ring(QQ)
            sage: C.super_categories()
            [Join of Category of graded algebras over Rational Field and Category of commutative additive monoids with realizations and Category of monoids with realizations, Category of commutative algebras over Rational Field]
        """
        R = self.base_ring()
        return [GradedAlgebras(R).WithRealizations(), CommutativeAlgebras(R)]
        
    class ParentMethods:
        
        def facade_for(self):
            r"""
            Returns all the parents this set is a facade for

            EXAMPLES::

                sage: A = AbstractPolynomialRing(CC)                                         
                sage: A.facade_for()
                []

            The facade_for parents are added dynamically.
            
            """
            l = [r for r in self.realizations()]
            l.extend(self._finite_rings)
            for r in self.realizations():
                l.extend(r.facade_for())
            return l



class AbstractPolynomialRing(UniqueRepresentation, Parent):
    r"""
    This class implements the ring of abstract polynomials.
    The number of variables is not set : this ring is the projective limit
    of all polynomial rings with a finite number of variables.

    INPUT:

    - ``R``: the base ring of the algebra
    - ``main_repr_var``, the letter corresponding to the set of variables, 
    it is used to represent several bases, default is ``x``
    - ``bases_category_class``, the class to use for the category of the 
    ring bases, default is basis.Bases
    - ``finite_bases_category_class``, the class to use for the category 
    of finite bases (bases of ``FiniteAbstractPolynomialRing``), default 
    is ``basis.Finite_bases``
    - ``always_show_main_var``, if True ``main_repr_var`` will be displayed 
    on elements of every basis, even the ones that don't use it directly 
    (Schubert basis, Demazure basis, ...), false by default, used on 
    ``DoubleAbstractPolynomialRing`` to differentiate the two sets of
    variables
    

    OUTPUT:

    - The abstract ring of multivariate polynomials on ``main_repr_var`` over ``R``

    EXAMPLES:: 

        sage: A = AbstractPolynomialRing(QQ); A
        The abstract ring of multivariate polynomials on x over Rational Field

    The abstract ring contains methods to obtain its basis::
    
        sage: m = A.monomial_basis(); m
        The ring of multivariate polynomials on x over Rational Field on the monomial basis
        sage: ma = A.ambient_space_basis("A"); ma
        The ring of multivariate polynomials on x over Rational Field on the Ambient space basis of type A

    From these basis, you can create elements with any number of variables and
    make simple operations::

        sage: pol1 = m[1,2,3]; pol1
        x[1, 2, 3]
        sage: pol2 = ma( [2,2,3,5] );pol2
        x(2, 2, 3, 5)
        sage: pol1 + pol2
        x(1, 2, 3, 0) + x(2, 2, 3, 5)
        sage: pol1*pol2
        x(3, 4, 6, 5)
        
    The coercion between elements with a different number of variables and
    between the monomial basis and ambient space basis are done automatically.

    TESTS::

        sage: A = AbstractPolynomialRing(QQ)
        sage: TestSuite(A).run()
    """    
    def __init__(self, R, main_repr_var = 'x', bases_category_class = None, finite_bases_category_class = None, always_show_main_var = False):
        r"""
        TESTS::

            sage: A = AbstractPolynomialRing(QQ)
        """
        assert(R in Rings())
        Parent.__init__(
            self,
            base = R,
            category = Facade_polynomial_ring(R)
        )
        self._finite_rings = set([])
        self._main_repr_var = main_repr_var
        self._show_main_var = always_show_main_var
        from basis import Bases, Finite_bases
        if(bases_category_class is None): self._bases_category_class = Bases
        else: self._bases_category_class = bases_category_class
        
        if(finite_bases_category_class is None): self._finite_bases_category_class = Finite_bases
        else: self._finite_bases_category_class = finite_bases_category_class
    
    def bases_category_class(self):
        r"""
        Returns the class to use for the category of the ring bases
        
        EXAMPLES::
        
        sage: A = AbstractPolynomialRing(QQ);
        sage: A.bases_category_class()
        <class 'sage.combinat.multivariate_polynomials.basis.Bases'>
        
        
        """
        return self._bases_category_class
        
    def finite_bases_category_class(self):
        r"""
        Returns the class to use for the category of finite bases (bases
        of ``FiniteAbstractPolynomialRing``)
        
        EXAMPLES::
        
        sage: A = AbstractPolynomialRing(QQ);
        sage: A.finite_bases_category_class()
        <class 'sage.combinat.multivariate_polynomials.basis.Finite_bases'>

        
        """
        return self._finite_bases_category_class

    def _repr_(self):
        r"""
        
        Print the name of ``self``
        
        EXAMPLES::

            sage: AbstractPolynomialRing(QQ)
            The abstract ring of multivariate polynomials on x over Rational Field
        """
        return "The abstract ring of multivariate polynomials on %s over %s"%(
            self._main_repr_var,
            self.base_ring()
        )

    def _repr_with_basis(self):
        r"""
        
        This methods is used by bases of the abstract ring to print its name. 
        As it is used as part of the name of concrete bases the word "abstract"
        has been removed from the usual ``_repr_`` method
        
        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: A._repr_with_basis()
            'The ring of multivariate polynomials on x over Rational Field'
        """
        return  "The ring of multivariate polynomials on %s over %s"%(
            self._main_repr_var,
            self.base_ring()
        )

    def a_realization(self):
        r"""
        Returns a default realization of ``self`` : the monomial basis
        
        EXAMPLES:
            sage: A = AbstractPolynomialRing(QQ)
            sage: A.a_realization()
            The ring of multivariate polynomials on x over Rational Field on the monomial basis

        """
        return self.monomial_basis()

    def an_element(self):
        r"""
        Returns an element of ``self``. By default, this element lies in the
        monomial basis.

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: A.an_element()
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
        """
        return self.monomial_basis().an_element()

    def _element_constructor_(self, element):
        r"""
        As ``self`` is an abstract algebra, this method will 
        just check if ``element`` belongs to ``self``
        
        INPUT:
        
        -``element`` the element to be contructed from
        
        OUTPUT:
        
        - The element itself if it belongs to ``self``, if not, 
        it raises a TypeError Exception
        
        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: p = A.an_element(); p 
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
            sage: A._element_constructor_(p)
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
            sage: A._element_constructor_(1)
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make '1' an element of 'The abstract ring of multivariate polynomials on x over Rational Field'

        """
        if any (parent.is_parent_of(element) for parent in self.facade_for()):
            return element
        raise TypeError, "do not know how to make '%s' an element of '%s'"%(element, self)

    def var(self, i, nb_variables = 0):
        r"""
        Returns the i_th variable as a monomial base element

        INPUT:

        - ``i``: the index of the variable to return
        - ``nb_variables``: the number of variables of the result, 
        default is ``i``, if ``nb_variables`` is lower than ``i`` it is
        ignored and changed to ``i``

        OUTPUT:

        - the ith variable as a monomial element

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ);
            sage: A.var(1)
            x[1]
            sage: A.var(3)
            x[0, 0, 1]
            sage: A.var(1,3)
            x[1, 0, 0]
            sage: A.var(4,3)
            x[0, 0, 0, 1]

        """
        if(nb_variables ==0 or nb_variables < i ): nb_variables =i
        vect = [0 for j in xrange(nb_variables)]
        vect[i-1] = 1
        return self.monomial_basis()( vect )

    def from_expr(self, expr, alphabet = None, second_alphabet = None):
        r"""
        Method delegated to the monomial basis.
        
        Constructs a polynomial from a symbolic expression.
        
        INPUT:
            - ``expr`` a symbolic expression, it must be a polynomial
            in all variables of ``variables``
            - ``alphabet`` (optional), a list of symbolic variables. 
            If not set, it takes ``expr.variables()``. The variables are matched
            to the vector key of the monomials by the order of the list.
            - ``second_alphabet`` (optional) a list of symbolic variables
            when working on a polynomial on 2 sets of variables
            
        EXAMPLES::
            
            sage: A = AbstractPolynomialRing(QQ)
            sage: var('x1,x2,x3')
            (x1, x2, x3)
            sage: expr = 3*x3 + x2^2 - x1*x3      
            sage: A.from_expr(expr)
            -x[1, 0, 1] + x[0, 2, 0] + 3*x[0, 0, 1] 
            sage: var('t1,t2')   
            (t1, t2)
            sage: K.<t1,t2> = QQ[]
            sage: K = K.fraction_field()
            sage: A = AbstractPolynomialRing(K) 
            sage: expr = t1*t2*x2 - x3^4*x1*4*t2^2
            sage: A.from_expr(expr,[x1,x2,x3])
            (-4*t2^2)*x[1, 0, 4] + t1*t2*x[0, 1, 0]
            
        Works with polynomials in two sets of variables::
        
            sage: D = DoubleAbstractPolynomialRing(QQ)
            sage: var('x1,x2,x3,y1,y2,y3')
            (x1, x2, x3, y1, y2, y3)
            sage: expr = x1*y1 +(y2*y3^2 - y1)*x3*x1^4
            sage: D.from_expr(expr,[x1,x2,x3],[y1,y2,y3]) 
            (y[1,0,0])*x[1, 0, 0] + (-y[1,0,0]+y[0,1,2])*x[4, 0, 1]

        """
        return self.monomial_basis().from_expr(expr, alphabet, second_alphabet)

    def finite_polynomial_ring(self, nb_variables, basis_repr = None):
        r"""
        This method returns the polynomial ring with a given finite number of variables
        
        INPUT:

        - ``nb_variables``: the number of variables
        - ``basis_repr``, the representation letter for the elements of the base, 
        by default, it is the main representation for the set of variable : 
        ``self._main_repr_var``

        OUTPUT:

        - The abstract ring of multivariate polynomials on x over ``R`` with 
          ``nb_variables`` variables where ``R`` is the algebra base ring

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ); A
            The abstract ring of multivariate polynomials on x over Rational Field
            sage: F3 = A.finite_polynomial_ring(3); F3
            The abstract ring of multivariate polynomials on x over Rational Field with 3 variables

        The finite number of variables ring contains method to obtain the ring
        basis on a finite number of variables::

            sage: m3 = F3.monomial_basis(); m3
            The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis
            sage: ma3 = F3.ambient_space_basis("A"); ma3
            The ring of multivariate polynomials on x over Rational Field with 3 variables on the Ambient space basis of type A

        Coercions between rings with a different number of variables are created
        dynamically::

            sage: m = A.monomial_basis()
            sage: pol1 = m[1,2,3]; pol1
            x[1, 2, 3]
            sage: pol1.parent()
            The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis
            sage: pol2 = m[1,1]; pol2
            x[1, 1]
            sage: pol2.parent()
            The ring of multivariate polynomials on x over Rational Field with 2 variables on the monomial basis
            sage: pol1 + pol2
            x[1, 1, 0] + x[1, 2, 3]
            sage: (pol1 + pol2).parent()
            The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis
        """
        if(basis_repr is None): basis_repr = self._main_repr_var
        F = FiniteAbstractPolynomialRing(self,nb_variables, basis_repr, bases_category_class = self.finite_bases_category_class())
        return F
        
    def _register_finite_ring(self,F):
        r"""
        Adds the finite ring F as one of ``self`` finite rings. It is called
        by the init function of ``F``
        
        EXAMPLES::
        sage: var('t')
        t
        sage: A = AbstractPolynomialRing(QQ[t])
        sage: A
        The abstract ring of multivariate polynomials on x over Univariate Polynomial Ring in t over Rational Field
        sage: A._finite_rings
        set([])
        sage: F3 = A.finite_polynomial_ring(3)
        sage: A._finite_rings
        set([The abstract ring of multivariate polynomials on x over Univariate Polynomial Ring in t over Rational Field with 3 variables])

        """
        if(not(F in self._finite_rings)):
            for ring in self._finite_rings:
                if( ring.nb_variables() != F.nb_variables() ): self._create_morphism(ring, F)
            self._finite_rings.add(F)

    def monomial_basis(self, basis_repr = None):
        r"""
        Returns the monomial basis of the polynomial ring.
        
        INPUT:
        - ``basis_repr``, the representation letter for the elements of the base, 
          by default, it is the main representation for the set of variable : 
          ``self._main_repr_var``

        OUTPUT:

        - The ring of multivariate polynomials on ``x`` over ``R`` on the 
          monomial basis where R is the algebra base ring defined in the abstract algebra

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis();m
            The ring of multivariate polynomials on x over Rational Field on the monomial basis
            sage: m[1,2,3]
            x[1, 2, 3]
            sage: m( [2,2] )
            x[2, 2]
        """
        if(basis_repr is None): basis_repr = self._main_repr_var
        from monomial import MonomialBasis
        return MonomialBasis(self, basis_repr)

    def ambient_space_basis(self, group_type, basis_repr = None):
        r"""
        Returns the basis where polynomials are indexed by a root system lattice

        INPUT:

        - ``group_type``: the letter that represents type of the weyl group
        - ``basis_repr``, the representation letter for the elements of the base, 
          by default, it is the main representation for the set of variable : 
          ``self._main_repr_var``

        OUTPUT:

        - The ring of multivariate polynomials on x over ``R`` on the Ambient space
          basis of type ``group_type`` where ``R`` is the algebra base ring
          defined in the abstract algebra

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: ma = A.ambient_space_basis("A"); ma
            The ring of multivariate polynomials on x over Rational Field on the Ambient space basis of type A
            sage: mb = A.ambient_space_basis("B"); mb
            The ring of multivariate polynomials on x over Rational Field on the Ambient space basis of type B

        Default coercions are created between the ambient space basis and
        monomial basis::

            sage: m = A.monomial_basis()
            sage: ma( m( [1,2,3] ) )
            x(1, 2, 3)
            sage: m( ma( [2,4] ) )
            x[2, 4]
            sage: mb( m( [1,2,3] ) )
            x(1, 2, 3)
            sage: m( mb( [2,4] ) )
            x[2, 4]
        """    
        if(basis_repr is None): basis_repr = self._main_repr_var
        from ambient_space_basis import PolynomialRingWithBasisFromAmbientSpace
        return PolynomialRingWithBasisFromAmbientSpace(self,group_type,"Ambient space basis of type " + group_type, basis_repr)    

    def from_morphism_basis(self, neutral_nb_variables, morphism_to_basis, get_basis_keys, get_morphism_on_basis, basis_name, basis_repr, variables_auto_coerce =False, **keywords):
        r"""
        Creates a basis defined by its morphism to another basis

        INPUT:

        - ``neutral_nb_variables``: the default number of variables to get the
            one element
        - ``morphism_to_basis``: the basis of the abstract polynomial ring on
            which the morphism will be defined
        - ``get_basis_keys``: a function with :
                input:
                    - ''nb_variables'' the number of variables
                output:
                    - the set of indexes that will be used to index elements 
                      of the basis on the given 
                      number of variables
        - ``get_morphism_on_basis``: a function with :
                input:
                    -''nb_variables'', the number of variables
                output:
                    - the function that will be used to create the module morphims on 
                      basis on the given number of variables
        - ``basis_name``: the name of the basis (used in repr)
        - ``basis_repr``: the basis representation for elements (exemple ``x``)
        - ``variables_auto_coerce``: if set to ``True``, a coercion will be 
        created between elements of the basis indexed by vectors of size
        n to basis on m>n variables by extending the vectors with zeros 
        (example: x[2,2,1] -> x[2,2,1,0,0]. Default is ``False``.
        - ``**keywords`` : the keywords sent to the ``CombinatorialFreeModule``
        morphism.

        OUTPUT:

        - the basis of which elements are indexed by the sets return
          by ``get_basis_keys`` and can be coerced on
          the ``morphims_to_basis`` basis

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: def get_basis_keys(n): code = "A" + str(n-1); return RootSystem(code).ambient_space(QQ)
            sage: def get_morphism_on_basis(n): return lambda key: m( [key[i] for i in xrange(n)])
            sage: MyBasis = A.from_morphism_basis(1,m,get_basis_keys,get_morphism_on_basis,"My Basis", "x"); MyBasis
            The ring of multivariate polynomials on x over Rational Field on the My Basis
            sage: MyBasis.an_element()
            x(2, 2, 3)
            sage: m( MyBasis.an_element() )
            x[2, 2, 3]

        We have recreated the basis on ambient space.
        
        """
        from basis import PolynomialRingWithBasisFromMorphism
        return PolynomialRingWithBasisFromMorphism(self, neutral_nb_variables, morphism_to_basis, basis_name, basis_repr, get_basis_keys, get_morphism_on_basis,variables_auto_coerce, **keywords)

    def linear_basis_on_vectors(self, group_type, basis_name, basis_repr, on_basis_method, extra_parameters = (), **keywords):
        r"""
        Creates a linear basis on objects inedexed by vectors based on an operation 
        to convert each object (through its vector) into a ambient space basis polynomial. 
         
        - ``group_type``: the letter that represents the type of the weyl group that
           will be used for the ambient space basis
        - ``basis_name``: the name of the basis (used in repr)
        - ``basis_repr``: the basis representation for elements 
        - ``on_basis_method``: a method that takes a vector (python list) and returns 
          the converted polynomial associated with it
          The ``on_basis_method`` should have this signature :
            Input :
                - ``v`` a python list representing the vector
                - ``basis`` the ambient space basis used to make the conversion 
                - ``call_back`` a call_back method to use the conversion recursively
                - ``**keywords`` extra parameters that could be used for convertion
                
            Output :
                - a polynomial expanded into the sent ``basis`` and corresponding to 
                  the objected indexed by the sent vector ``v``
        - ``extra_parameters``: (default : empty) a tuple containing the extra parameters 
        to be sent to the ``on_basis_method`` as tuples ``(key,val)``
        
        - ``**keyword`` : parameters used to create the morphism to the ambient space basis, 
          sent to ``CombinatorialFreeModule.module_morphism``. By default, ``triangular`` is
          set to ``upper`` : change it explicitly to ``None`` if you're basis is not. 
          A default ``cmp`` method is also used to order keys : it compares degrees first and 
          then orders vectors by inversed lexical order. This comparasion method is the classic 
          one, used by Schubert, Demazure and Grothendieck polynomials. 
         
        OUTPUT :
        
            - a basis named ``basis_name`` and defined by its conversion to the ambient space basis of type ``group_type``
            
        EXAMPLES::
        
            sage: A = AbstractPolynomialRing(QQ)
            sage: def schubert_on_basis(v, basis, call_back):
            ...     for i in xrange(len(v)-1):
            ...         if(v[i]<v[i+1]):
            ...             v[i], v[i+1] = v[i+1] + 1, v[i]
            ...             return call_back(v).divided_difference(i+1)
            ...     return basis(v)
            sage: myBasis = A.linear_basis_on_vectors("A","MySchub","Y",schubert_on_basis)
            sage: pol = myBasis[2,1,3];pol
            Y(2, 1, 3)
            sage: pol.expand()
            x(2, 1, 3) + x(2, 2, 2) + x(2, 3, 1) + x(3, 1, 2) + x(3, 2, 1) + x(4, 1, 1)
            sage: myBasis(A.an_element())
            Y(0, 0, 0) + 2*Y(1, 0, 0) + Y(1, 2, 3) - Y(1, 3, 2) + 3*Y(2, 0, 0) - Y(2, 1, 3) + Y(2, 3, 1) + Y(3, 1, 2) - Y(3, 2, 1) + Y(4, 1, 1)
        
        We have recreated the Schubert basis. Let's see an example with parameters::
        
            sage: def t_inverse(v, basis, call_back, t1=1, t2=1): return t1/(t2*basis(v))
            sage: tInverse = A.linear_basis_on_vectors("A","tInverse","T",t_inverse, (("t1",2), ("t2",4)), triangular = None)                    
            sage: pol = tInverse[1,1,2]; pol    
            T(1, 1, 2)
            sage: pol.expand()
            1/2*x(-1, -1, -2)

        """
        from linear_basis_on_vectors import LinearBasisOnVectors
        ambient_space_basis = self.ambient_space_basis(group_type)
        return LinearBasisOnVectors(self, ambient_space_basis, basis_name, basis_repr, on_basis_method, extra_parameters, **keywords)

    def schubert_basis_on_vectors(self, basis_name = None, basis_repr = "Y"):
        r"""
        Creates the simple Schubert basis where schubert polynomials are indexed
        by vectors.
        
        Here is the definition we use. For $v = (v_1, \cdots, v_n) \in \mathbb{N}^n$, we define
        
        $Y_v = x_1^{v_1}x_2^{v_2}\cdotsx_n^{v_n}$ if $v$ is a partition, i.e, if $v_1 \geq v_2 \geq \cdots \geq v_n$.
        
        Otherwise, we have for $ v_i > v_{i+1}$
        
        $Y_{\cdots v_{i+1} v_i-1 \cdots} = Y_v \partial_i$ where $\partial_i$ is the ith divided difference. 
        
        The vectors indexing the Schubert polynomials can as well been seen as 
        lehmer codes.
        
        INPUT:

        - ``basis_name``: (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``: (defaul: ``Y``) the basis representation for elements 

        OUTPUT:

        - The ring of multivariate polynomials on x over ``R`` on the Schubert basis
          of type ``group_type`` index by vectors where ``R`` is the algebra base
          ring defined in the abstract algebra.

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: Schub = A.schubert_basis_on_vectors(); Schub
            The ring of multivariate polynomials on x over Rational Field on the Schubert basis of type A (indexed by vectors)
            sage: Schub.an_element()
            Y(2, 2, 3)
            sage: Schub[1,2,3]
            Y(1, 2, 3)

        Let us see the coercions::

            sage: y = Schub.an_element(); y
            Y(2, 2, 3)
            sage: y.expand()
            x(2, 2, 3) + x(2, 3, 2) + x(3, 2, 2)
            sage: ma = A.ambient_space_basis("A")
            sage: ma( y )
            x(2, 2, 3) + x(2, 3, 2) + x(3, 2, 2)
            sage: m = A.monomial_basis()
            sage: m( y )
            x[2, 2, 3] + x[2, 3, 2] + x[3, 2, 2]
            sage: ma.an_element(); Schub( ma.an_element())
            x(2, 2, 3)
            Y(2, 2, 3) - Y(2, 3, 2)
            sage: m.an_element(); Schub(m.an_element())
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
            Y(0, 0, 0) + 2*Y(1, 0, 0) + Y(1, 2, 3) - Y(1, 3, 2) + 3*Y(2, 0, 0) - Y(2, 1, 3) + Y(2, 3, 1) + Y(3, 1, 2) - Y(3, 2, 1) + Y(4, 1, 1)


        Let us see some operations::

            sage: Schub[1,2] + Schub[3,0,0]
            Y(1, 2, 0) + Y(3, 0, 0)
            sage: Schub.an_element() * Schub.an_element()
            Y(4, 4, 6) + Y(4, 5, 5)
            sage: Schub[1,2] * Schub[3,0,0]
            Y(4, 2, 0) + Y(5, 1, 0)
        """
        from linear_basis_on_vectors import SchubertBasisOnVectors
        if(basis_name is None):
            basis_name = "Schubert basis of type A (indexed by vectors)"  
        if(self._show_main_var): basis_repr+= self._main_repr_var
        ambient_space_basis = self.ambient_space_basis("A")
        return SchubertBasisOnVectors(self, ambient_space_basis, basis_name, basis_repr)     

    def demazure_basis_on_vectors(self, group_type ="A", basis_name = None, basis_repr = "K"):
        r"""
        Creates the Demazure basis where demazure / key polynomials are indexed
        by vectors.
        
        Here is the definition we use for type A. For $v = (v_1, \cdots, v_n) \in \mathbb{N}^n$, 
        we define
        
        $K_v = x_1^{v_1}x_2^{v_2}\cdotsx_n^{v_n}$ if $v$ is a partition, i.e, if $v_1 \geq v_2 \geq \cdots \geq v_n$.
        
        Otherwise, we have for $ v_i > v_{i+1}$
        
        $K_{\cdots v_{i+1} v_i \cdots} = K_v \pi_i$ where $\pi_i$ is the ith isobar divided difference. 
        
        The vectors indexing the key polynomials can as well been seen as lehmer codes.

        INPUT:

        - ``group_type``: (default: ``A``) the letter that represents the type of the weyl group
        - ``basis_name``: (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``: (default: ``K``) the basis representation for elements 

        OUTPUT:

        - The ring of multivariate polynomials on x over ``R`` on the Demazure basis
          of type ``group_type`` indexed by vectors where ``R`` is the algebra
          base ring defined in the abstract algebra.

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: Dem = A.demazure_basis_on_vectors("A"); Dem
            The ring of multivariate polynomials on x over Rational Field on the Demazure basis of type A (indexed by vectors)
            sage: Dem.an_element()
            K(2, 2, 3)
            sage: Dem[1,2,3]
            K(1, 2, 3)

        Let us see some coercions::

            sage: k = Dem.an_element(); k
            K(2, 2, 3)
            sage: k.expand()
            x(2, 2, 3) + x(2, 3, 2) + x(3, 2, 2)
            sage: ma = A.ambient_space_basis("A")
            sage: ma( k )
            x(2, 2, 3) + x(2, 3, 2) + x(3, 2, 2)
            sage: m = A.monomial_basis()
            sage: m( k )
            x[2, 2, 3] + x[2, 3, 2] + x[3, 2, 2]
            sage: ma.an_element(); Dem( ma.an_element() )
            x(2, 2, 3)
            K(2, 2, 3) - K(2, 3, 2)
            sage: m.an_element(); Dem( m.an_element() )
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
            K(0, 0, 0) + 2*K(1, 0, 0) + K(1, 2, 3) - K(1, 3, 2) + 3*K(2, 0, 0) - K(2, 1, 3) + K(2, 3, 1) + K(3, 1, 2) - K(3, 2, 1)


        Let us see some operations::

            sage: Dem[1,2] + Dem[3,0,0]
            K(1, 2, 0) + K(3, 0, 0)
            sage: Dem.an_element() * Dem.an_element()
            K(4, 4, 6) + K(4, 5, 5)
            sage: Dem[1,2] * Dem[3,0,0]
            K(4, 2, 0) + K(5, 1, 0)
            
        We can also have type B, C or D key polynomials::
        
            sage: DemB = A.demazure_basis_on_vectors("B"); DemB
            The ring of multivariate polynomials on x over Rational Field on 
            the Demazure basis of type B (indexed by vectors)
            sage: pol = DemB[2,1,-2]; pol
            K(2, 1, -2)
            sage: pol.expand()
            x(2, 1, 0) + x(2, 1, -2) + x(2, 1, -1) + x(2, 1, 1) + x(2, 1, 2) + x(2, 2, 0) + x(2, 2, -1) + x(2, 2, 1)

        """
        from linear_basis_on_vectors import DemazureBasisOnVectors
        if(basis_name is None):
            basis_name = "Demazure basis of type " +group_type+" (indexed by vectors)"  
        if(self._show_main_var): basis_repr+= self._main_repr_var
        ambient_space_basis = self.ambient_space_basis(group_type)
        return DemazureBasisOnVectors(self, ambient_space_basis, basis_name, basis_repr, "isobaric_divided_difference") 
    
    def demazure_hat_basis_on_vectors(self, group_type ="A", basis_name = None, basis_repr = "^K"):
        r"""
        Creates the Demazure hat basis where demazure polynomials are indexed
        by vectors.
        
        Here is the definition we use for type A. For $v = (v_1, \cdots, v_n) \in \mathbb{N}^n$, we define
        
        $K_v = x_1^{v_1}x_2^{v_2}\cdotsx_n^{v_n}$ if $v$ is a partition, i.e, if $v_1 \geq v_2 \geq \cdots \geq v_n$.
        
        Otherwise, we have for $ v_i > v_{i+1}$
        
        $K_{\cdots v_{i+1} v_i \cdots} = K_v \pi_i$ where $\pi_i$ is the ith isobar hat divided difference. 
        
        The vectors indexing the key polynomials can as well been seen as lehmer codes.

        INPUT:

        - ``group_type``: (default: ``A``) the letter that represents the type of the weyl group
        - ``basis_name``: (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``: (default: ``^K``) the basis representation for elements 

        OUTPUT:

        - The ring of multivariate polynomials on x over ``R`` on the Demazure hat basis
          of type ``group_type`` indexed by vectors where ``R`` is the algebra
          base ring defined in the abstract algebra.    
          
        EXAMPLES::
            
            sage: A = AbstractPolynomialRing(QQ);
            sage: Demh = A.demazure_hat_basis_on_vectors("A"); Demh
            The ring of multivariate polynomials on x over Rational Field on the Demazure hat basis of type A (indexed by vectors)
            sage: Demh.an_element()
            ^K(2, 2, 3)
            
        Let us see some coercions::

            sage: kh = Demh[1,2,4]; kh
            ^K(1, 2, 4)
            sage: kh.expand()
            x(1, 2, 4) + x(1, 3, 3) + x(2, 2, 3)
            sage: ma = A.ambient_space_basis("A")
            sage: ma( kh )
            x(1, 2, 4) + x(1, 3, 3) + x(2, 2, 3)
            sage: m = A.monomial_basis()
            sage: m( kh )
            x[1, 2, 4] + x[1, 3, 3] + x[2, 2, 3]
            sage: Demh( ma[1,2,4] )
            ^K(1, 2, 4) - ^K(1, 3, 3) + ^K(2, 3, 2)
            sage: Demh( m[1,2,4] )
            ^K(1, 2, 4) - ^K(1, 3, 3) + ^K(2, 3, 2)

        Let us see some operations::

            sage: Demh[1,2] + Demh[3,0,0]
            ^K(1, 2, 0) + ^K(3, 0, 0)
            sage: Demh.an_element() * Demh.an_element()
            ^K(4, 4, 6) - ^K(4, 5, 5) - ^K(5, 4, 5)
            sage: Demh[1,2] * Demh[3,0,0]              
            ^K(4, 2, 0)
            
        We can also have type B, C or D hat key polynomials::
        
            sage: DemhB = A.demazure_hat_basis_on_vectors("B"); Demh
            The ring of multivariate polynomials on x over Rational Field on the Demazure hat basis of type A (indexed by vectors)
            sage: pol = DemhB[2,1,-2];pol
            ^K(2, 1, -2)
            sage: pol.expand()
            x(2, 1, 0) + x(2, 1, -2) + x(2, 1, -1) + x(2, 1, 1)

        """
        from linear_basis_on_vectors import DemazureBasisOnVectors
        if(basis_name is None):
            basis_name = "Demazure hat basis of type " +group_type+" (indexed by vectors)"  
        if(self._show_main_var): basis_repr+= self._main_repr_var
        ambient_space_basis = self.ambient_space_basis(group_type)
        return DemazureBasisOnVectors(self, ambient_space_basis, basis_name, basis_repr, "hat_isobaric_divided_difference") 
    

    def grothendieck_negative_basis_on_vectors(self, basis_name = None, basis_repr = "G"):
        r"""
        Creates the simple Grothendieck basis where Grothendieck polynomials are indexed
        by vectors. The Negative stands for that we use a definition of the basis where 
        the variable have negative exposants.
        
        Here is the definition we use. For $v = (v_1, \cdots, v_n) \in \mathbb{N}^n$, we define
        
        $G_v = \prod_{1 \cdots n} (1 - \frac{1}{x_i})^{v_i}$ if $v$ is a partition, i.e, if $v_1 \geq v_2 \geq \cdots \geq v_n$.
        
        Otherwise, we have for $ v_i > v_{i+1}$
        
        $G_{\cdots v_{i+1} v_i-1 \cdots} = G_v \pi_i$ where $\pi_i$ is the ith isobar divided difference. 
        
        The vectors indexing the Grothendieck polynomials can as well been seen as lehmer codes.
        
         

        INPUT:

        - ``basis_name``: (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``: (defaul: ``G``) the basis representation for elements 

        OUTPUT:

        - The ring of multivariate polynomials on x over ``R`` on the 
          Grothendieck basis of type ``group_type`` with negative exposants 
          indexd by vectors where ``R`` is the algebra base ring defined in 
          the abstract algebra.

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: Groth = A.grothendieck_negative_basis_on_vectors(); Groth
            The ring of multivariate polynomials on x over Rational Field on the Grothendieck basis of type A with negative exposants (indexed by vectors)
            sage: Groth.an_element()
            G(2, 2, 3)
            sage: Groth[1,2,3] 
            G(1, 2, 3)

        We can convert a Grothendieck polynomial into the Monomial or 
        Ambient Space basis but not the other way around as Grothendieck 
        polynomials are with negative exposants. Note that conversion 
        from monomials with negative exposants into Grothendieck polynomials 
        is NOT implemented ::

            sage: g = Groth[0,1];g
            G(0, 1)
            sage: g.expand()
            x(0, 0) - x(-1, -1)
            sage: ma = A.ambient_space_basis("A")
            sage: ma( g )
            x(0, 0) - x(-1, -1)
            sage: m = A.monomial_basis()
            sage: m( g )
            x[0, 0] - x[-1, -1]
            sage: pol = m[0,0] - m[-1,-1]  
            sage: Groth( pol)
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= x[0, 0] - x[-1, -1]) an element of self (=The ring of multivariate polynomials on x over Rational Field with 2 variables on the Grothendieck basis of type A with negative exposants (indexed by vectors))

        We can add Grothendieck polynomials but not multiply them as this 
        would use conversion from monomials into Grothendieck polynomials ::

            sage: Groth[1,2] + Groth[3,0,0]
            G(1, 2, 0) + G(3, 0, 0)
            sage: Groth[1,2] * Groth[3,0,0]
            Traceback (most recent call last):
            ...
            NotImplementedError: The product is not implemented for this basis
        """
        from linear_basis_on_vectors import GrothendieckNegativeBasisOnVectors
        if(basis_name is None):
            basis_name = "Grothendieck basis of type A with negative exposants (indexed by vectors)"  
        if(self._show_main_var): basis_repr+= self._main_repr_var
        ambient_space_basis = self.ambient_space_basis("A")
        return GrothendieckNegativeBasisOnVectors(self, ambient_space_basis, basis_name, basis_repr) 
        
    def grothendieck_positive_basis_on_vectors(self, basis_name = None, basis_repr = "G"):
        r"""
        Creates the simple Grothendieck basis where Grothendieck polynomials 
        are indexed by vectors. The positive stands for that we use a definition
        of the basis where variables have positive exposants. 
        
        It corresponds to the basis given by ``grothendieck_negative_basis_on_vectors`` 
        by a change of variables :
        
        $x_i = 1 - \fract{1}{x_i}$
        
        Here is the definition we use. For $v = (v_1, \cdots, v_n) \in \mathbb{N}^n$, we define
        
        $G_v = x_1^{v_1}x_2^{v_2}\cdotsx_n^{v_n}$ if $v$ is a partition, i.e, if $v_1 \geq v_2 \geq \cdots \geq v_n$.
        
        Otherwise, we have for $ v_i > v_{i+1}$
        
        $G_{\cdots v_{i+1} v_i-1 \cdots} = G_v \pi_i$ where $\pi_i$ is the ith isobar divided difference. 
        
        The vectors indexing the Grothendieck polynomials can as well been seen as lehmer codes.
       
        INPUT:

        - ``basis_name``: (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``: (default: ``G``) the basis representation for elements 

        OUTPUT:

        - The ring of multivariate polynomials on x over ``R`` on the Grothendieck basis
          of type ``group_type`` with positive exposants indexed by vectors where ``R``
          is the algebra base ring defined in the abstract algebra.

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: Groth = A.grothendieck_positive_basis_on_vectors(); Groth
            The ring of multivariate polynomials on x over Rational Field on the Grothendieck basis of type A, with positive exposants (indexed by vectors) 
            sage: Groth.an_element()
            G(2, 2, 3)
            sage: Groth[1,2,3]
            G(1, 2, 3)

        Let us see some coercions::

            sage: g = Groth.an_element(); g
            G(2, 2, 3)
            sage: g.expand()
            x(2, 2, 3) + x(2, 3, 2) - x(2, 3, 3) + x(3, 2, 2) - x(3, 2, 3) - x(3, 3, 2) + x(3, 3, 3)
            sage: ma = A.ambient_space_basis("A")
            sage: ma(g)
            x(2, 2, 3) + x(2, 3, 2) - x(2, 3, 3) + x(3, 2, 2) - x(3, 2, 3) - x(3, 3, 2) + x(3, 3, 3)
            sage: m = A.monomial_basis()
            sage: m(g)
            x[2, 2, 3] + x[2, 3, 2] - x[2, 3, 3] + x[3, 2, 2] - x[3, 2, 3] - x[3, 3, 2] + x[3, 3, 3]
            sage: ma.an_element(); Groth( ma.an_element())
            x(2, 2, 3)
            G(2, 2, 3) - G(2, 3, 2) + G(2, 3, 3) - G(3, 3, 2) + G(3, 3, 3)
            sage: m.an_element(); Groth( m.an_element())  
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
            G(0, 0, 0) + 2*G(1, 0, 0) + G(1, 2, 3) - G(1, 3, 2) + G(1, 3, 3) + 3*G(2, 0, 0) - G(2, 1, 3) + G(2, 2, 3) + G(2, 3, 1) - 2*G(2, 3, 2) + G(2, 3, 3) + G(3, 1, 2) - G(3, 1, 3) - G(3, 2, 1) + G(3, 2, 2) - G(3, 3, 2) + G(3, 3, 3) + G(4, 1, 1) - G(4, 1, 3)



        Let us see some operations::

            sage: Groth[1,2] + Groth[3,0,0]
            G(1, 2, 0) + G(3, 0, 0)
            sage: Groth.an_element() * Groth.an_element()
            G(4, 4, 6) + G(4, 5, 5) - G(4, 5, 6)
            sage: Groth[1,2] * Groth[3,0,0]
            G(4, 2, 0) + G(5, 1, 0) - G(5, 2, 0)
        """
        
        from linear_basis_on_vectors import GrothendieckPositiveBasisOnVectors
        if(basis_name is None):
            basis_name = "Grothendieck basis of type A, with positive exposants (indexed by vectors) "  
        if(self._show_main_var): basis_repr+= self._main_repr_var
        ambient_space_basis = self.ambient_space_basis("A")
        return GrothendieckPositiveBasisOnVectors(self, ambient_space_basis, basis_name, basis_repr) 
        
        
    def macdonald_basis_on_vectors(self, t1 =None , t2=None, q=None , basis_name = None, basis_repr = "M"):
        r"""
        Creates the the basis of non symmetric Macdonald polynomials indexed by vectors. 
        
        INPUT:

        - ``t1``: (default: symbolic variable t1) the first parameter for the Hecke algebra operator
        - ``t2``: (default: symbolic variable t2) the second parameter for the Hecke algebra operator
        - ``q``: (default: symbolic variable q) the specific q parmater of the polynomials
        - ``basis_name``: (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``: (default: ``M``) the basis representation for elements
        
        OUTPUT:

        - The ring of multivariate polynomials on x over ``R`` on the Macdolald basis
          of type A indexed by vectors where ``R ``is the algebra base ring defined 
          in the abstract algebra. 
          
        EXAMPLES::
        
            sage: var('t1 t2 q')
            (t1, t2, q)
            sage: K.<t1,t2,q> = QQ[]
            sage: K = K.fraction_field()
            sage: A = AbstractPolynomialRing(K)
            sage: Mac = A.macdonald_basis_on_vectors(); Mac
            The ring of multivariate polynomials on x over Fraction Field of Multivariate Polynomial Ring in t1, t2, q over Rational Field on the Macdonald basis of type A (indexed by vectors) 
            sage: Mac.an_element()
            M(2, 2, 3)
            sage: Mac[1,2]
            M(1, 2)
            
        Let us see some coercions::
        
            sage: pol = Mac[1,2];pol
            M(1, 2)
            sage: pol.expand()
            t2^3*x(0, 0) + t2^2*x(1, 0) + ((t2*q+t2)/q)*x(1, 1) + 1/q*x(1, 2) + ((t2^2*q+t2^2)/q)*x(0, 1) + t2/q*x(0, 2)
            sage: ma = A.ambient_space_basis("A")
            sage: ma( pol )
            t2^3*x(0, 0) + t2^2*x(1, 0) + ((t2*q+t2)/q)*x(1, 1) + 1/q*x(1, 2) + ((t2^2*q+t2^2)/q)*x(0, 1) + t2/q*x(0, 2)
            sage: m = A.monomial_basis()
            sage: m( pol )
            t2^3*x[0, 0] + t2^2*x[1, 0] + ((t2*q+t2)/q)*x[1, 1] + 1/q*x[1, 2] + ((t2^2*q+t2^2)/q)*x[0, 1] + t2/q*x[0, 2]
            sage: Mac( m[1,0] + m[0,1])
            (t1-t2)*M(0, 0) + (1/(-t2))*M(1, 0) + ((t1*q-t1)/(t1*q+t2))*M(0, 1)
            
        Let us see some operations::
            sage: Mac[1,2] + Mac[1,0]
            M(1, 0) + M(1, 2)
            sage: Mac[1,2] * Mac[1,0]
            ((t1^2*t2*q^2-t1^2*t2*q-t2^3*q+t2^3)/(-t1*q-t2))*M(1, 2) + ((t1*q^2+t2*q^2)/(t1*q+t2))*M(1, 3) + ((t1^2*t2*q^3-t1^2*t2*q^2-t2^3*q^2+t2^3*q)/(-t1^2*q^2-2*t1*t2*q-t2^2))*M(2, 2)
        """     
        if(t1 is None): t1 = self.base_ring()(var('t1'))
        if(t2 is None): t2 = self.base_ring()(var('t2'))
        if(q is None): q = self.base_ring()(var('q'))
        from linear_basis_on_vectors import MacdonaldBasisOnVectors
        if(basis_name is None):
            basis_name = "Macdonald basis of type A (indexed by vectors) "  
        if(self._show_main_var): basis_repr+= self._main_repr_var
        ambient_space_basis = self.ambient_space_basis("A")
        return MacdonaldBasisOnVectors(self, ambient_space_basis, basis_name, basis_repr, t1, t2,q)
        
    def _create_morphism(self,f1,f2):
        r"""
        Creates a morphism between two `FiniteAbstractPolynomialRing` on their `FiniteMonomialBasis`
        by adding extra variables to the elements of the one with least variables. 
        The morphism is then registered as a coercion.
        
        This method is called by `finite_polynomial_ring` each time a new `FinitePolynomialRing` is created
        
        INPUT:
            - ``f1`` a `FiniteAbstractPolynomialRing`
            - ``f2`` another `FiniteAbstractPolynomialRing` with a different number of variables
        
        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: F1 = A.finite_polynomial_ring(1); F1
            The abstract ring of multivariate polynomials on x over Rational Field with 1 variable
            sage: F2 = A.finite_polynomial_ring(2); F2
            The abstract ring of multivariate polynomials on x over Rational Field with 2 variables
            sage: M1 = F1.monomial_basis();M1
            The ring of multivariate polynomials on x over Rational Field with 1 variable on the monomial basis
            sage: M2 = F2.monomial_basis(); M2
            The ring of multivariate polynomials on x over Rational Field with 2 variables on the monomial basis
            sage: M2(M1.an_element())
            x[0, 0] + 3*x[1, 0] + 3*x[2, 0]
            
            by creating ``F1`` and ``F2`` through ``A`` a coercion between their monomial basis ``M1`` and ``M2``
            has been created

        """
        if(f1.nb_variables() > f2.nb_variables()):
            temp = f1
            f1 = f2
            f2 = temp
            
            

    def change_nb_variables(self, pol, nb_variables):
        r"""
        Forcing the addition of variables 

        INPUT:

        - ``pol``: a polynomial expressed in any basis
        - ``nb_variables``: the new number of variables

        OUTPUT:

        - the polynomial seen as a polynomial of ``nb_variables`` variables

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: M = A.monomial_basis()
            sage: pol = M.an_element(); pol
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
            sage: A.change_nb_variables(pol, 5)
            x[0, 0, 0, 0, 0] + 2*x[1, 0, 0, 0, 0] + x[1, 2, 3, 0, 0] + 3*x[2, 0, 0, 0, 0]
            sage: MA = A.ambient_space_basis("A")
            sage: pol = MA.an_element(); pol
            x(2, 2, 3)
            sage: A.change_nb_variables(pol, 5)
            x(2, 2, 3, 0, 0)
        """
        if(nb_variables < pol.nb_variables()): 
            raise NotImplementedError, "This method doesn't reduce the number of variables, use reduce_nb_variables"%()
        basis = pol.parent().basis_tower().finite_basis(nb_variables)
        return basis( pol )   

    def reduce_nb_variables(self, pol):
        """
        Creates a polynomial by removing all last variables with exposant 0 of the given polynomial

        INPUT:

        - ``pol``: the polynomial to be reduced

        OUTPUT:

        - a polynomial equal to ``pol`` and without all the last variables
          with exposant 0

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: M = A.monomial_basis()
            sage: pol = M[1,2,2,0,0] + M[2,5];pol
            x[1, 2, 2, 0, 0] + x[2, 5, 0, 0, 0]
            sage: pol.nb_variables()
            5
            sage: red = A.reduce_nb_variables(pol); red
            x[1, 2, 2] + x[2, 5, 0]
            sage: red.nb_variables()
            3
            sage: red == pol
            True
        """
        max = 1
        for ind, coeff in pol:
            for i in xrange(pol.nb_variables()-1,-1,-1):
                if(ind[i]!=0):
                    if(i+1>max): max = i+1
                    break
        if(max==pol.nb_variables()): return pol
        codomain = pol.parent().basis_tower().finite_basis(max)
        return sum( [coeff * codomain([ind[i] for i in xrange(0,max)]) for ind, coeff in pol] )

    def maxDiffDiv(self, pol):
        """
        Apply the maximum divided difference to the polynomial.
        As the result is a symmetrical function, it is writen as a symmetrical
        Schubert polynomial (same as Schur functions). Giving the result in this
        basis makes the algorithm faster and the result compact.

        INPUT:

        - ``pol``: the polynomial to apply the maximum divided difference on

        OUTPUT:

        - the result polynomial in Schubert basis after applying maximum
          divided difference

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: M = A.monomial_basis()
            sage: pol = M[1,2,3]
            sage: A.maxDiffDiv(pol)
            -Y(1, 1, 1)
        """
        return pol.maxDiffDiv()

    def maxPi(self, pol):
        """
        Apply the maximum isobaric divided difference to the polynomial.
        As the result is a symmetrical function, it is writen as a symmetrical
        Schubert polynomial (same as Schur functions). Giving the result in this
        basis makes the algorithm faster and the result compact.

        INPUT:

        - ``pol``: the polynomial to apply the maximum isobaric divided 
          difference on

        OUTPUT:

        - the result polynomial in Schubert basis after applying maximum
          isobaric divided difference

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: M = A.monomial_basis()
            sage: pol = M[2,1,3]
            sage: A.maxPi(pol)
            -Y(2, 2, 2)
            """
        return pol.maxPi()

class FiniteAbstractPolynomialRing(UniqueRepresentation, Parent):
    r"""
    This class implements the ring of abstract polynomial on a given number
    of variables. It is obtained by calling AbstractPolynomialRing.finite_polynmial_ring
    See the documentation of the above method for more information
    
    INPUT:
    
    - ``polynomial_ring_tower``: the class for Abstract polynomial ring on an unset number of variables 
    from which the ``FiniteAbstractPolynomialRing`` comes from. A ``FiniteAbstractPolynomialRing`` always 
    comes from a ``AbstractPolynomialRing`` which contains general informations linke the base ring
    - ``nb_variables`` : the number of variables for the polynomials
    - ``main_repr_var``, the letter corresponding to the set of variables, 
    it is used to represent several bases, default is ``x``
    - ``bases_category_class`` : the class to use for the category 
    of concrete bases, default is ``basis.Finite_bases``
    

    TESTS::

        sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteAbstractPolynomialRing
        sage: A = AbstractPolynomialRing(QQ)
        sage: B = FiniteAbstractPolynomialRing(A,4)
        sage: B
        The abstract ring of multivariate polynomials on x over Rational Field with 4 variables
        sage: TestSuite(B).run()
    """
    def __init__(self, polynomial_ring_tower, nb_variables, main_repr_var = 'x', bases_category_class = None):
        r"""
        TESTS::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteAbstractPolynomialRing
            sage: A = AbstractPolynomialRing(QQ)
            sage: B = FiniteAbstractPolynomialRing(A,4)
        """
        self._polynomial_ring_tower = polynomial_ring_tower
        self._nb_variables = nb_variables
        Parent.__init__(
            self,
            base = polynomial_ring_tower.base_ring(),
            category = GradedAlgebras(polynomial_ring_tower.base_ring()).WithRealizations()
        )
        self._main_repr_var = main_repr_var
        self._polynomial_ring_tower._register_finite_ring(self)
        m = SetMorphism( Hom(self, polynomial_ring_tower), lambda x: x)
        m.register_as_coercion()
        
        from basis import Finite_bases
        if(bases_category_class is None): self._bases_category_class = Finite_bases
        else: self._bases_category_class = bases_category_class

    def bases_category_class(self):
        r"""
        Returns the class to use for he category of concrete bases
        
        OUTPUT:
        
        - the class to use for he category of concrete bases
        
        EXAMPLES::
        
        sage: A = AbstractPolynomialRing(QQ)
        sage: F3 = A.finite_polynomial_ring(3)
        sage: F3.bases_category_class()
        <class 'sage.combinat.multivariate_polynomials.basis.Finite_bases'>   
        """
        
        return self._bases_category_class

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteAbstractPolynomialRing
            sage: A = AbstractPolynomialRing(QQ)
            sage: FiniteAbstractPolynomialRing(A,4)
            The abstract ring of multivariate polynomials on x over Rational Field with 4 variables
        """
        if(self.nb_variables()>1):
            variables_str = "variables"
        else:
            variables_str = "variable"
        return "%s with %s %s"%(self.polynomial_ring_tower(),self.nb_variables(),variables_str)

    def _repr_with_basis(self):
        r"""
        
        This methods is used by bases of the abstract ring to print its name. 
        As it is used as part of the name of concrete bases the word "abstract"
        has been removed from the usual ``_repr_`` method
        
        TESTS::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteAbstractPolynomialRing
            sage: A = AbstractPolynomialRing(QQ)
            sage: B = FiniteAbstractPolynomialRing(A,8)
            sage: B._repr_with_basis()
            'The ring of multivariate polynomials on x over Rational Field with 8 variables'
        """
        if(self.nb_variables()>1):
            variables_str = "variables"
        else:
            variables_str = "variable"
        return "%s with %s %s"%(self.polynomial_ring_tower()._repr_with_basis(),self.nb_variables(),variables_str)

    def nb_variables(self):
        r"""
        Returns the number of variables of ``self``.

        EXAMPLES::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteAbstractPolynomialRing
            sage: A = AbstractPolynomialRing(QQ)
            sage: B = FiniteAbstractPolynomialRing(A,8)
            sage: B.nb_variables()
            8
            sage: B = FiniteAbstractPolynomialRing(A,0)
            sage: B.nb_variables()
            0
        """
        return self._nb_variables

    def polynomial_ring_tower(self):
        r"""
        Returns the polynomial ring tower given to define ``self``.

        EXAMPLES::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteAbstractPolynomialRing
            sage: K = CyclotomicField(3)
            sage: A = AbstractPolynomialRing(K)
            sage: B = FiniteAbstractPolynomialRing(A,3)
            sage: B.polynomial_ring_tower()
            The abstract ring of multivariate polynomials on x over Cyclotomic Field of order 3 and degree 2
        """
        return self._polynomial_ring_tower


    def _element_constructor_(self, element):
        r"""
        
        As ``self`` is an abstract algebra, this method will 
        just check if ``element`` belongs to ``self``
        
        INPUT:
        
        -``element`` the element to be contructed from
        
        OUTPUT:
        
        - The element itself if it belongs to ``self``, if not, 
        it raises a TypeError Exception
        
        TESTS::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteAbstractPolynomialRing
            sage: A = AbstractPolynomialRing(QQ)
            sage: B = FiniteAbstractPolynomialRing(A,3)
            sage: p = B.an_element()
            sage: B._element_constructor_(p)
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
            sage: B._element_constructor_(1)
            Traceback (most recent call last):
            ...
            ValueError: '1' is not an element of 'The abstract ring of multivariate polynomials on x over Rational Field with 3 variables'
        """
        if self.is_parent_of(element):
            return element
        raise ValueError, "'%s' is not an element of '%s'"%(element, self)


    def an_element(self):
        r"""
        Returns an element of ``self``. By default, this element lies in
        the monomial basis.

        EXAMPLES::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteAbstractPolynomialRing
            sage: A = AbstractPolynomialRing(QQ)
            sage: B = FiniteAbstractPolynomialRing(A,7)
            sage: B.an_element()
            x[0, 0, 0, 0, 0, 0, 0] + 2*x[1, 0, 0, 0, 0, 0, 0] + x[1, 2, 3, 0, 0, 0, 0] + 3*x[2, 0, 0, 0, 0, 0, 0]

        """
        return self.a_realization().an_element()
        
    def a_realization(self):
        r"""
        Returns a default realization of ``self``, the monomial basis
        
        EXAMPLES::
            
            sage: A = AbstractPolynomialRing(QQ)                                         
            sage: F = A.finite_polynomial_ring(3)
            sage: F.a_realization()
            The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis
        """
        return self.monomial_basis()
        
        
    def monomial_basis(self, basis_repr = None):
        r"""
        Returns the algebra ``self`` view in the monomials basis.
        
        INPUT:
            - ``basis_repr``, the representation letter for the elements of the base, by default, it is the main representation for
        the set of variable : ``self._main_repr_var``

        EXAMPLES::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteAbstractPolynomialRing
            sage: A = AbstractPolynomialRing(QQ)
            sage: B = FiniteAbstractPolynomialRing(A,5)
            sage: B.monomial_basis()
            The ring of multivariate polynomials on x over Rational Field with 5 variables on the monomial basis
        """
        if(basis_repr is None): basis_repr = self._main_repr_var
        from monomial import FiniteMonomialBasis
        return FiniteMonomialBasis(self, basis_repr)

    def ambient_space_basis(self, letter, basis_repr = None):
        r"""
        Returns the algebra ``self`` view in the proper ambient space of the
        root system design by ``letter``.

        EXAMPLES::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteAbstractPolynomialRing
            sage: A = AbstractPolynomialRing(QQ)
            sage: B = FiniteAbstractPolynomialRing(A,2)
            sage: B.ambient_space_basis("B")
            The ring of multivariate polynomials on x over Rational Field with 2 variables on the Ambient space basis of type B
        """
        if(basis_repr is None): basis_repr = self._main_repr_var
        from ambient_space_basis import FinitePolynomialRingWithBasisFromAmbientSpace
        if(letter == "A"): number = self.nb_variables()-1
        else: number = self.nb_variables()
        code = str(letter) + str(number)
        basis = FinitePolynomialRingWithBasisFromAmbientSpace(self,code,letter,"Ambient space basis of type " + letter, basis_repr)
        return basis

    def from_morphism_basis(self, polynomial_ring_tower, basis_name, basis_repr):
        r"""
        Creates a basis defined by its morphism to another basis

        INPUT:

        - ``polynomial_ring_tower``: the basis of ``AbsractPolynomialRing`` which is a facade to this basis and represents 
        it on a undefined number of variables. It must have a `get_morphism_on_basis` method and a `get_basis_keys` method
        as well as a ``morphism_to_basis``
        - ``basis_name``: the name of the basis (used in repr)
        - ``basis_repr``: the basis representation for elements (exemple "x")

        OUTPUT:

        - the basis of which elements are indexed by the set return
          by ``polynomial_ring_tower.get_basis_keys(self.nb_variables())`` and can be coerced on
          the ``morphims_to_basis`` basis of the ``polynomial_ring_tower`` on the right number of variables (``self.nb_variables()``)

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: M = A.monomial_basis()
            sage: def get_basis_keys(n): code = "A" + str(n-1); return RootSystem(code).ambient_space(QQ)
            sage: def get_morphism_on_basis(n): return lambda key: M( [key[i] for i in xrange(n)])
            sage: MyBasis = A.from_morphism_basis(1,M,get_basis_keys,get_morphism_on_basis,"My Basis", "X"); MyBasis
            The ring of multivariate polynomials on x over Rational Field on the My Basis
            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteAbstractPolynomialRing
            sage: F2 = FiniteAbstractPolynomialRing(A,2)
            sage: MyFiniteBasis = F2.from_morphism_basis(MyBasis,"MyBasis", "X"); MyFiniteBasis
            The ring of multivariate polynomials on x over Rational Field with 2 variables on the MyBasis
            sage: MyFiniteBasis.an_element()
            X(2, 2)
            sage: M( MyFiniteBasis.an_element())
            x[2, 2]


        We have recreated the basis on ambient space.
        
        """

        from basis import FinitePolynomialRingWithBasisFromMorphism
        return FinitePolynomialRingWithBasisFromMorphism(self, polynomial_ring_tower, basis_name, basis_repr)


    def linear_basis_on_vectors(self, polynomial_ring_tower, basis_name, basis_repr, **keywords):
        r"""
        Creates a linear basis on objects inedexed by vectors based on an operation 
        to convert each object (through its vector) into a ambient space basis polynomial. 
        The type of the ambient space basis and the method of conversion are all contained into the
        ``polynomial_ring_tower``
         
        - ``polynomial_ring_tower``: the basis of ``AbsractPolynomialRing`` which is a facade to this basis and 
        represents it on a undefined number of variables. It should inherit from a basis.LinearBasisOnVectors
        - ``basis_name``: the name of the basis (used in repr)
        - ``basis_repr``: the basis representation for elements 
        - ``**keyword`` : parameters used to create the morphism to the ambient space basis, 
          sent to ``CombinatorialFreeModule.module_morphism``. 
         
        OUTPUT :
        
            - a basis named ``basis_name`` and defined by its conversion to an ambient space basis,
             the type of the ambient space basis and the method of conversion are all contained into the
             ``polynomial_ring_tower``
         

        EXAMPLES::
        
            sage: A = AbstractPolynomialRing(QQ)
            sage: def schubert_on_basis(v, basis, call_back):
            ...     for i in xrange(len(v)-1):
            ...         if(v[i]<v[i+1]):
            ...             v[i], v[i+1] = v[i+1] + 1, v[i]
            ...             return call_back(v).divided_difference(i+1)
            ...     return basis(v)
            sage: myBasis = A.linear_basis_on_vectors("A","MySchub","Y",schubert_on_basis)
            sage: F3 = A.finite_polynomial_ring(3)
            sage: myFiniteBasis = F3.linear_basis_on_vectors(myBasis,"MySchub","Y")
            sage: myFiniteBasis
            The ring of multivariate polynomials on x over Rational Field with 3 variables on the MySchub
            
        """
        from linear_basis_on_vectors import FiniteLinearBasisOnVectors
        return FiniteLinearBasisOnVectors(self, polynomial_ring_tower, basis_name, basis_repr, **keywords)    

   

