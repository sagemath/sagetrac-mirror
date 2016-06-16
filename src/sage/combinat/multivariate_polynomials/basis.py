"""
categories and technichal upperclasses for basis of
of the `AbstractPolynomialRing` and `FiniteAbstractPolynomialRing`
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Viviane Pons < pons at univ-mlv.fr  >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.all import var
from sage.categories.all import GradedAlgebrasWithBasis, CommutativeAlgebras, GradedAlgebras
from sage.categories.all import Realizations
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.facade_sets import FacadeSets
from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutation
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation



class Bases(Category_realization_of_parent):
    r"""
    The category of bases of (see :class:``AbstractPolynomialRing``). These
    classes are both realizations of ``AbstractPolynomialRing`` and facades
    for their finite equivalent.

    EXAMPLES::

        sage: A = AbstractPolynomialRing(QQ)
        sage: M = A.monomial_basis()
        sage: M.category()
        Category of bases of The abstract ring of multivariate polynomials on x over Rational Field
        sage: M.realization_of()
        The abstract ring of multivariate polynomials on x over Rational Field
        sage: pol = M[1,2,3]
        sage: M.facade_for()
        [The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis]

    Even so :meth:`__contains__<sage.categories.facade_sets.FacadeSets.ParentMethods.facade_for>`
    is defined in FacadeSets it does not override the default Parent method and
    so we have unexpectd behaviour::

        sage: M = A.monomial_basis()
        sage: Schub = A.schubert_basis_on_vectors("A")
        sage: pol = M[1,2,3]
        sage: M.is_parent_of(pol)
        True
        sage: Schub.is_parent_of(pol)
        False
        sage: pol in M
        True
        sage: pol in Schub
        True
        sage: M3 = M.finite_basis(3)
        sage: Schub3 = Schub.finite_basis(3)
        sage: M3.is_parent_of(pol)
        True
        sage: Schub3.is_parent_of(pol)
        False
        sage: pol in M3
        True
        sage: pol in Schub3
        False

    The default :meth:`__contains__<sage.categories.facade_sets.FacadeSets.ParentMethods.facade_for>`
    method says that a parent contains an element if it can construct it :
    both Schub and M can construct pol. But the combinatorialFreeModule
    :meth:`__contains__<sage.combinat.free_module.CombinatorialFreeModule.__contains__>`
    method returns true if the element's parents is actually the
    combinatorialFreeModule itself. That is why pol is considered to be part of
    Schub but not part of any of its finite equivalent. There is a contradiction
    here and a choice has to be make. It would seem logical that "pol in Schub"
    returns False, and that is what would happen if the category method was
    overritgthing the default Parent method.

    TESTS::

        sage: from sage.combinat.multivariate_polynomials.basis import Bases
        sage: A = AbstractPolynomialRing(QQ)
        sage: C = Bases(A)
        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        r"""
        OUTPUT:
            - a list of the super categories of ``self``
        
        EXAMPLES::
            
            sage: from sage.combinat.multivariate_polynomials.basis import Bases
            sage: A = AbstractPolynomialRing(QQ)
            sage: C = Bases(A)
            sage: C.super_categories()
            [Category of realizations of The abstract ring of multivariate polynomials on x over Rational Field, Category of graded algebras over Rational Field, Category of commutative algebras over Rational Field, Category of facade sets]

        """
        return [Realizations(self.base()),GradedAlgebras(self.base().base_ring()), CommutativeAlgebras(self.base().base_ring()), FacadeSets() ]

    class ParentMethods:
        def abstract_algebra(self):
            r"""
            Returns the abstract algebra over ``self``.

            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: m.abstract_algebra()
                The abstract ring of multivariate polynomials on x over Rational Field

            """
            return self.category().base()

        def facade_for(self):
            r"""
            Returns all the parents this set is a facade for. ``self`` is a bases where the number
            of variables is not defined, it is a facade for all equivalent bases with a fixed
            number of variables

            EXAMPLES::

                sage: A = AbstractPolynomialRing(ZZ.quotient_ring(2*ZZ))  
                sage: m = A.monomial_basis()
                sage: m.facade_for()
                []
                sage: m1 = m.finite_basis(1)
                sage: m2 = m.finite_basis(2)
                sage: m3 = m.finite_basis(3)
                sage: #random m.facade_for()
                [The ring of multivariate polynomials on x over Ring of integers modulo 2 with 3 variables on the monomial basis, The ring of multivariate polynomials on x over Ring of integers modulo 2 with 2 variables on the monomial basis, The ring of multivariate polynomials on x over Ring of integers modulo 2 with 1 variable on the monomial basis]

            TESTS::
                sage: A = AbstractPolynomialRing(ZZ.quotient_ring(2*ZZ))  
                sage: m = A.monomial_basis()
                sage: m1 = m.finite_basis(1)
                sage: m2 = m.finite_basis(2)
                sage: m3 = m.finite_basis(3)
                sage: set(m.facade_for()) == set([m1,m2,m3])
                True
            """
            return list(self._facade_for)
            
        def var(self, i, nb_variables = 0):
            r"""
            Returns the i_th variable as an element of ``self``
            
            INPUT:

            - ``i``: the index of the variable to return
            - ``nb_variables``: the number of variables of the result, 
            default is ``i``, if ``nb_variables`` is lower than ``i`` it is
            ignored and changed to ``i``
            
            OUTPUT:
            - the i_th variable as an element of ``self``
            
            EXAMPLES::
            
                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: m.var(3)
                x[0, 0, 1]
                sage: Schub = A.schubert_basis_on_vectors();     
                sage: Schub.var(1)
                Y(1)
                sage: Schub.var(2)
                -Y(1, 0) + Y(0, 1)
                
            """
            return self( self.abstract_algebra().var(i, nb_variables))      
            
        
        def one(self):
            r"""
            OUPUT:
                - the ``one`` of the basis.
                
            As ``self`` is a facade for the finie bases, the ``one`` returned comes
            from the method ``_default_finite_basis`` 
            
            EXAMPLES::
            
            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: m.one()
            x[0]

            """
            return self._default_finite_basis().one()
        
        def from_expr(self, expr, alphabet = None, second_alphabet = None):
            r"""
            Constructs a polynomial from a symbolic expression.
            
            INPUT:
                - ``expr`` a symbolic expression, it must be a polynomial
                in all variables of ``variables``
                - ``alphabet`` (optional), a list of symbolic variables. 
                If not set, it takes ``expr.variables()``. The variables are matched
                to the vector key of the monomials by the order of the list.
                - ``second_alphabet`` (optional) a list of symbolic variables
                when working on a polynomial on 2 sets of variables
                
            EXAMPLES:
                
                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: var('x1,x2,x3')
                (x1, x2, x3)
                sage: expr = 3*x3 + x2^2 - x1*x3
                sage: m.from_expr(expr)
                -x[1, 0, 1] + x[0, 2, 0] + 3*x[0, 0, 1]
                sage: K.<t1,t2> = QQ[]  
                sage: K = K.fraction_field()
                sage: A = AbstractPolynomialRing(K) 
                sage: m = A.monomial_basis()
                sage: expr = t1*t2*x2 - x3^4*x1*4*t2^2
                sage: m.from_expr(expr,[x1,x2,x3])
                (-4*t2^2)*x[1, 0, 4] + t1*t2*x[0, 1, 0]
            
            Works with polynomials in two sets of variables::
            
                sage: D = DoubleAbstractPolynomialRing(QQ)
                sage: dm = D.monomial_basis()
                sage: var('x1,x2,x3,y1,y2,y3')
                (x1, x2, x3, y1, y2, y3)
                sage: expr = x1*y1 +(y2*y3^2 - y1)*x3*x1^4
                sage: dm.from_expr(expr,[x1,x2,x3],[y1,y2,y3])
                (y[1,0,0])*x[1, 0, 0] + (-y[1,0,0]+y[0,1,2])*x[4, 0, 1]

            """
            if alphabet is None: alphabet = expr.variables()
            expr = expr.expand()
            
            #Test polynomial
            for v in alphabet:
                if(not expr.is_polynomial(v)):
                    raise ValueError, "The expression is not a polynomial of the variable %s"%(v)
            
            from sage.symbolic.expression import operator
            if(expr.operator() == operator.add):
                monomials = expr.operands()
            else:
                monomials = [expr]
            basis = self.abstract_algebra().monomial_basis().finite_basis(len(alphabet))
            result = basis(0)
            #Summing, monomial by monomial
            for monomial in monomials:
                key = []
                for v in alphabet: 
                    d = monomial.degree(v) 
                    if d != 0:
                        key.append(d)
                        monomial = monomial / v**d
                    else:
                        key.append(0)
                if hasattr(self.base_ring(),'from_expr'):
                    result = result + basis.base_ring().from_expr(monomial, second_alphabet) * basis(key)
                else:
                    result = result + basis.base_ring()(monomial) * basis(key)
            result = self(result)
            return result
            

    class ElementMethods:
        pass


class Finite_bases(Category_realization_of_parent):
    r"""
    This category is the default catergry for bases of :class:``FiniteAbstractPolynomialRing``,
    i.e. bases of :class:``AbstractPolynomialRing`` with a finite number of variables. It is really
    specific to multivariate polynomials and should not be used for other purpose. Nervertheless,
    it can be inherited by subcategories that would correspond to more specific multivariate polynomials,
    see :class:``Finite_double_bases``
    
    INPUT:
        - ``abstract_algebra``, the finite abstract algebra of which ``self`` is a base
        - ``algebra_tower``, the base of :class:``AbstractPolynomialRing`` that correspond to 
        this base on an unset number of variables 
        
    OUTPUT:
        the category of finite_bases of the abstract ring of multivariate polynomials  with `nb_variables` 
        variables
    
    EXAMPLES::
        sage: A = AbstractPolynomialRing(QQ)
        sage: pol = A.an_element()
        sage: pol.parent().category()
        The category of bases of The abstract ring of multivariate polynomials on x over Rational Field with 3 variables where algebra tower is The ring of multivariate polynomials on x over Rational Field on the monomial basis

    TESTS::
        
        sage: A = AbstractPolynomialRing(QQ);      
        sage: F3 = A.finite_polynomial_ring(3)
        sage: M = A.monomial_basis()
        sage: from sage.combinat.multivariate_polynomials.basis import Finite_bases
        sage: C = Finite_bases(F3, M)
        sage: TestSuite(C).run()
    """
    def __init__(self, abstract_algebra, algebra_tower):
        r"""
        TESTS::
            sage: A = AbstractPolynomialRing(QQ);      
            sage: F3 = A.finite_polynomial_ring(3)
            sage: M = A.monomial_basis()
            sage: from sage.combinat.multivariate_polynomials.basis import Finite_bases
            sage: C = Finite_bases(F3, M); C
            The category of bases of The abstract ring of multivariate polynomials on x over Rational Field with 3 variables where algebra tower is The ring of multivariate polynomials on x over Rational Field on the monomial basis

        """
        self._algebra_tower = algebra_tower
        Category_realization_of_parent.__init__(
            self,
            abstract_algebra
        )
        
        
        
    def _repr_(self):
        r"""
        
        Print the name of ``self``
        
        EXAMPLES::
        
            sage: A = AbstractPolynomialRing(QQ);      
            sage: F3 = A.finite_polynomial_ring(3)
            sage: M = A.monomial_basis()
            sage: from sage.combinat.multivariate_polynomials.basis import Finite_bases
            sage: C = Finite_bases(F3, M); C
            The category of bases of The abstract ring of multivariate polynomials on x over Rational Field with 3 variables where algebra tower is The ring of multivariate polynomials on x over Rational Field on the monomial basis

        """
        return "The category of bases of %s where algebra tower is %s"%(self.base(), self.algebra_tower())
    
    
    def algebra_tower(self):
        r"""
        OUTPUT:
        
        - the base of :class:``AbstractPolynomialRing`` that correspond to 
        this base on an unset number of variables 
        
        EXAMPLES::
        
            sage: A = AbstractPolynomialRing(QQ);      
            sage: F3 = A.finite_polynomial_ring(3)
            sage: M = A.monomial_basis()
            sage: from sage.combinat.multivariate_polynomials.basis import Finite_bases
            sage: C = Finite_bases(F3, M)
            sage: C.algebra_tower()
            The ring of multivariate polynomials on x over Rational Field on the monomial basis

        """
        return self._algebra_tower

    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ);      
            sage: F3 = A.finite_polynomial_ring(3)
            sage: M = A.monomial_basis()
            sage: from sage.combinat.multivariate_polynomials.basis import Finite_bases
            sage: C = Finite_bases(F3, M)
            sage: C.super_categories()
            [Category of realizations of The abstract ring of multivariate polynomials on x over Rational Field with 3 variables, Category of graded algebras with basis over Rational Field]
        """
        R = self.base().base_ring()
        return [Realizations(self.base()),GradedAlgebrasWithBasis(R)]

    class ParentMethods:
        def abstract_algebra(self):
            r"""
            Return the abstract algebra over ``self``, i.e. the abstract algebra with a fixed number of 
            variables

            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                    
                sage: m = A.monomial_basis()
                sage: m3 = m.finite_basis(3)
                sage: m3.abstract_algebra()
                The abstract ring of multivariate polynomials on x over Rational Field with 3 variables
                sage: F3 = A.finite_polynomial_ring(3)
                sage: F3
                The abstract ring of multivariate polynomials on x over Rational Field with 3 variables
                sage: F3 == m3.abstract_algebra()
                True

            """
            return self.category().base()
            
        def var(self, i):
            r"""
            Returns the i_th variable as an element of ``self`` if ``i`` is lower or equal 
            than the number of variables of ``self``, if not the result is given in the 
            corresponding basis.
            
            INPUT:

            - ``i``: the index of the variable to return
            
            OUTPUT:
            - the i_th variable as an element of ``self``
            
            EXAMPLES::
            
                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: m3 = m.finite_basis(3)
                sage: m3.var(1)
                x[1, 0, 0]
                sage: m3.var(2)
                x[0, 1, 0]
                sage: Schub = A.schubert_basis_on_vectors()
                sage: Schub3 = Schub.finite_basis(3)
                sage: Schub3.var(1)
                Y(1, 0, 0)
                sage: Schub3.var(2)
                -Y(1, 0, 0) + Y(0, 1, 0)
                sage: m3.var(4)
                x[0, 0, 0, 1]
                sage: Schub3.var(4)
                -Y(0, 0, 1, 0) + Y(0, 0, 0, 1)

            """
            return self.basis_tower().var(i, self.nb_variables())
            

    class ElementMethods:
        
        def maxDiffDiv(self):
            r"""
            Apply the maximum divided difference to the polynomial.
            As the result is a symemtrical function, it is writen as a
            symmetrical Schubert polynomial (same as Schur functions).
            Giving the result in this basis makes the algorithm faster and
            the result compact.

            OUTPUT:

            - the result polynomial in Schubert basis after applying maximum divided difference

            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: M = A.monomial_basis()
                sage: pol = M[1,2,3]
                sage: pol.maxDiffDiv()
                -Y(1, 1, 1)
            """
            M = self.parent().abstract_algebra().polynomial_ring_tower().monomial_basis()
            res = M( self )
            coeffs = list(res)
            return sum( [coeffs[i][1] * M._maxDiffDiv_on_basis(coeffs[i][0]) for i in xrange(len(coeffs)) ])

        def maxPi(self):
            r"""
            Apply the maximum isobaric divided difference to the polynomial.
            As the result is a symmetrical function, it is writen as a
            symmetrical Schubert polynomial (same as Schur functions). Giving
            the result in this basis makes the algorithm faster and the result
            compact.

            OUTPUT:

            - the result polynomial in Schubert basis after applying maximum isobaric divided difference

            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: M = A.monomial_basis()
                sage: pol = M[2,1,3]
                sage: pol.maxPi()
                -Y(2, 2, 2)
            """
            M = self.parent().abstract_algebra().polynomial_ring_tower().monomial_basis()
            max = [self.nb_variables() -i -1 for i in xrange(self.nb_variables())]
            return (M(max) * self).maxDiffDiv()

        def nb_variables(self):
            r"""
            Returns the number of variables of the polynomial ``self``.
            
            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: M = A.monomial_basis()
                sage: pol1 = M[1,2,3]; pol1
                x[1, 2, 3]
                sage: pol1.nb_variables()
                3
                sage: pol2 = M[2,2]; pol2
                x[2, 2]
                sage: pol2.nb_variables()
                2
                sage: sum = pol1 + pol2; sum
                x[1, 2, 3] + x[2, 2, 0]
                sage: sum.nb_variables()
                3
            """
            return self.parent().nb_variables()

        def change_nb_variables(self, nb_variables):
            r"""
            Creates a new polynomial with ``nb_variables`` variables by adding
            new variables to the given polynomial with exposant `0`.

            INPUT:

            - ``nb_variables``: the number of variables for the result polynomial

            OUTPUT:

            - A polynomomial with ``nb_variables`` variables and equal the the ``self``

            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: M = A.monomial_basis()
                sage: pol = M[2,2,3] + M[1,3,2]; pol
                x[1, 3, 2] + x[2, 2, 3]
                sage: pol.nb_variables()
                3
                sage: pol_5var = pol.change_nb_variables(5); pol_5var
                x[1, 3, 2, 0, 0] + x[2, 2, 3, 0, 0]
                sage: pol_5var == pol
                True

            .. warning::

                You can't use this method to remove variables.
            """
            return self.parent().abstract_algebra().polynomial_ring_tower().change_nb_variables(self,nb_variables)

        def reduce_nb_variables(self):
            r"""
            Creates a polynomial by removing all last variables with
            exposant `0` of the given polynomial.

            OUTPUT:

            - a polynomial equal to ``self`` and without all the last variables
              with exposant `0`

            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)   
                sage: M = A.monomial_basis()
                sage: pol = M[1,2,2,0,0] + M[2,5]; pol
                x[1, 2, 2, 0, 0] + x[2, 5, 0, 0, 0]
                sage: pol.nb_variables()
                5
                sage: red = pol.reduce_nb_variables(); red
                x[1, 2, 2] + x[2, 5, 0]
                sage: red.nb_variables()
                3
                sage: red == pol
                True
            """
            return self.parent().abstract_algebra().polynomial_ring_tower().reduce_nb_variables(self)
            
        def elements(self):
            r"""
            Returns a list of elements of ``self`` as tuple ``(key, coeff)``
            
            OUTPUT::
            
            - a list of tuples ``(key, coeff)`` where ``key`` is of type list and ``coeff`` is the coefficient 
            of the specific key
            
            Note: this is different of ``list(self)`` because the key are given as lists and don't depend 
            on the actual object we use in the Combinatorial free module 
            
            EXAMPLES::
            
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: m = A.monomial_basis()
                sage: pol = m[0,0,1] + m[2,1]; pol
                x[2, 1, 0] + x[0, 0, 1]
                sage: pol.elements()
                [([2, 1, 0], 1), ([0, 0, 1], 1)]
                sage: ma = A.ambient_space_basis("A")
                sage: pol = ma[2,2,1,0]           
                sage: pol = 3*ma[2,2,1,0] + 2*ma[1,1,3];pol
                2*x(1, 1, 3, 0) + 3*x(2, 2, 1, 0)
                sage: pol.elements()
                [([2, 2, 1, 0], 3), ([1, 1, 3, 0], 2)]
                sage: Schub = A.schubert_basis_on_vectors()
                sage: pol = Schub[1,2,1] + 3*Schub[3,3,1]; pol
                Y(1, 2, 1) + 3*Y(3, 3, 1)
                sage: pol.elements()
                [([3, 3, 1], 3), ([1, 2, 1], 1)]

            """
            kkeys = self.parent()._basis_keys._basis_keys # keys indexing the polynomials are themselves free modules
            return [([key[kk] for kk in kkeys],c) for (key,c) in self]
       
        

        def swap_vars(self, i, j):
            r"""
            Exchanges variables in the polynomial
            
            INPUT:

            - `i`,`j`: the index of the two variables to exchange

            OUTPUT:

            - the polynomial obtainend by exchanging the variables

            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: M = A.monomial_basis()
                sage: pol = M[1,2,3] + M[3,2,1]; pol
                x[1, 2, 3] + x[3, 2, 1]
                sage: res = pol.swap_vars(1,2); res
                x[2, 1, 3] + x[2, 3, 1]
                sage: Schub = A.schubert_basis_on_vectors("A")
                sage: y = Schub(pol); y
                Y(1, 2, 3) - Y(1, 3, 2) - Y(2, 1, 3) + Y(2, 3, 1) + Y(3, 1, 2) + Y(4, 1, 1)
                sage: y.swap_vars(1,2)
                Y(2, 1, 3) - Y(2, 2, 2) - Y(3, 1, 2) - Y(4, 1, 1)
            """
            m = self.parent().abstract_algebra().monomial_basis()
            pol = m(self)
            res = m.zero()
            for key, coeff in pol:
                l = [key[k] for k in xrange(self.nb_variables())]
                temp = l[i-1]
                l[i-1] = l[j-1]
                l[j-1] = temp
                res += coeff * m(l)
            return self.parent()(res)
            
        def perm_vars(self, p):
            r"""
            Action of the permutation ``p`` on the polynomial ``self``:
            
            `x^{(v1,v2,v3)} . p = x^{(p(v1),p(v2),p(v3))}` 
            
            INPUT:
            - A permutation ``p`` to act on the polynomial
            
            OUTPUT:
            - the polynomial given by the action of ``p`` on ``self``
            
            EXAMPLES::
            
                sage: A = AbstractPolynomialRing(QQ)
                sage: pol = A.an_element(); pol
                x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
                sage: p = Permutation([3,2,1])
                sage: pol.perm_vars(p)
                x[0, 0, 0] + x[3, 2, 1] + 2*x[0, 0, 1] + 3*x[0, 0, 2]
                sage: Schub = A.schubert_basis_on_vectors()
                sage: pol = Schub.an_element();pol
                Y(2, 2, 3)
                sage: pol.perm_vars(p)
                Y(2, 2, 3)
                sage: pol.expand()
                x(2, 2, 3) + x(2, 3, 2) + x(3, 2, 2)

            """
            size = len(p)
            pol = self
            parent = pol.parent()
            if(size> self.nb_variables()):
                pol = pol.change_nb_variables(len(p))
                parent = pol.parent()
            if(size<pol.nb_variables()):
                l = [i+1 for i in xrange(pol.nb_variables())]
                for i in xrange(size): l[i] = p[i]
                p = Permutation( l )
            m = pol.parent().abstract_algebra().monomial_basis()
            res = m.zero()
            pol = m(pol)
            for key, coeff in pol:
                l = [key[k] for k in xrange(pol.nb_variables())]
                l = p.action( l )
                res += coeff * m(l)
            return parent(res)  

        def subs_var(self, args):
            r"""
            Replacing the given variables by the given values

            INPUT:

            - A tuple or a list of tuple containing the number of the variable
              to replace and the value to replace it with

            OUTPUT:

            - the polynomial obtained by replacing the given variables by the
              given values

            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: M = A.monomial_basis()
                sage: pol = M[1,2,3] + M[2,2,2];pol
                x[1, 2, 3] + x[2, 2, 2]

            replacing `x_1` by `x_1^2`::

                sage: pol.subs_var((1,A.var(1)^2))
                x[2, 2, 3] + x[4, 2, 2]

            replacing `x_1` by `x_2`::

                sage: pol.subs_var((1,2))
                4*x[0, 2, 2] + 2*x[0, 2, 3]

            replacing `x_i` by `x_i^2` for all vars::

                sage: pol.subs_var([ (i,A.var(i)^2) for i in xrange(1,4)])
                x[2, 4, 6] + x[4, 4, 4]
            """
            if(type(args) is not list):
                args = [args]

            m = self.parent().abstract_algebra().monomial_basis()
            pol = m(self)

            for t in args:
                res = m.zero()
                for key, coeff in pol:
                    l = [key[k] for k in xrange(self.nb_variables())]
                    exp = l[t[0]-1]
                    l[t[0]-1] = 0
                    res += m(l) * t[1]**exp * coeff
                pol = res
            
            return self.parent()(res)
            
        def subs_basis(self, new_basis):
            r"""
            Substitutes the basis of ``self`` by ``new_basis`` - this is not
            a change of basis
            
            INPUT:
            - ``self`` a polynomial
            - ``new_basis`` a combinatorial free module that accepts lists
            as inputs or its keys
            
            EXAMPLES::
            
                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: Schub = A.schubert_basis_on_vectors()
                sage: Schub3 = Schub.finite_basis(3)
                sage: pol = m[0,0,1] + m[2,1]; pol
                x[2, 1, 0] + x[0, 0, 1]
                sage: pol.subs_basis(Schub3)
                Y(2, 1, 0) + Y(0, 0, 1)

            
            """
            keys = new_basis._basis_keys
            d = { keys(key):c for (key,c) in self.elements()}
            return new_basis._from_dict(d)
            
        def subs_on_keys(self, d):
            r"""
            Replace the given keys by their value in ``d``
            
            INPUT:
            - a dictonary of tuples
            
            EXAMPLES::
            
                sage: A = AbstractPolynomialRing(QQ)                                      
                sage: m = A.monomial_basis()
                sage: pol = m[0,0,1] + m[2,1]; pol
                x[2, 1, 0] + x[0, 0, 1]
                sage: d = {(2,1,0):0}     
                sage: pol.subs_on_keys(d)
                x[0, 0, 1]

            """
            res = []
            for (key,c) in self.elements():
                key = tuple(key)
                if(d.has_key(key)):
                    res.append(c*d[key])
                else:
                    res.append(c*self.parent()(key))
            return sum(res)
            
            
        def subs_on_coeffs(self, in_dict=None, **kwds):
            r"""
            Apply the ``subs`` method on all coefficients of ``self`` 
            
            INPUT:
                - ``in_dict`` - (optional) dict with variable:value pairs
                - ``**kw`` - names parameters

            OUPUT
                - the polynomial where ``subs`` has been applied on all coefficients
                
            EXAMPLES::
            
                sage: var('t1,t2,q')
                (t1, t2, q)
                sage: K.<t1,t2,q> = QQ[]    
                sage: K = K.fraction_field()
                sage: A = AbstractPolynomialRing(K)
                sage: pol = t1 + (t1 +t2 +q) * A.var(2);pol
                t1*x[0, 0] + (t1+t2+q)*x[0, 1]
                sage: pol.subs_on_coeffs(t1=t2)
                t2*x[0, 0] + (2*t2+q)*x[0, 1]

            """
            f = lambda x: x.subs(in_dict, **kwds)
            return self.map_coefficients(f)
            
        def to_expr(self, alphabet = None, alphabety = None):
            r"""
            Returns ``self`` as a symbolic expression on a given alphabet.
            
            If the polynomial basis has a method ``_to_expr_on_basis`` that
            converts each basis element into a symbolic expression, this method
            will be used.
            
            Otherwise, the polynomial will be converted into the monomial basis
            and use the ``MonomialBasis._to_expr_on_basis``method. The returned expression
            is then a sum of monomial factors.
            
            INPUT:
                - ``alphabet`` an optional alphabet for the variables,  
                if let to ``None``, it uses basis_repr + i (by default, 
                'x1' to 'xn')
                - ``alphabety`` an optional second alphabet for the 
                coeffients. it is used if the coeffincents have a 
                ``to_expr`` method.
            
            EXAMPLES::
                
                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: pol = m.an_element()
                sage: pol
                x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
                sage: pol.to_expr()
                x1*x2^2*x3^3 + 3*x1^2 + 2*x1 + 1
                sage: var('a,b,c')
                (a, b, c)
                sage: alphabet = [a,b,c]
                sage: pol.to_expr(alphabet = alphabet)
                a*b^2*c^3 + 3*a^2 + 2*a + 1
                sage: Schub = A.schubert_basis_on_vectors()
                sage: pol = Schub.an_element(); pol
                Y(2, 2, 3)
                sage: m(pol)
                x[2, 2, 3] + x[2, 3, 2] + x[3, 2, 2]
                sage: pol.to_expr()
                x1^3*x2^2*x3^2 + x1^2*x2^3*x3^2 + x1^2*x2^2*x3^3
                sage: pol.to_expr(alphabet = alphabet)
                a^3*b^2*c^2 + a^2*b^3*c^2 + a^2*b^2*c^3

            It also works with polynomials on two sets of variables::
            
                sage: D = DoubleAbstractPolynomialRing(QQ)  
                sage: dpol = D.an_element() * D.base_ring().an_element()
                sage: dpol
                (y[0,0,0]+2*y[1,0,0]+y[1,2,3]+3*y[2,0,0])*x[0, 0, 0] + (2*y[0,0,0]+4*y[1,0,0]+2*y[1,2,3]+6*y[2,0,0])*x[1, 0, 0] + (y[0,0,0]+2*y[1,0,0]+y[1,2,3]+3*y[2,0,0])*x[1, 2, 3] + (3*y[0,0,0]+6*y[1,0,0]+3*y[1,2,3]+9*y[2,0,0])*x[2, 0, 0]
                sage: dpol.to_expr()
                (y1*y2^2*y3^3 + 3*y1^2 + 2*y1 + 1)*x1*x2^2*x3^3 + y1*y2^2*y3^3 + 3*(y1*y2^2*y3^3 + 3*y1^2 + 2*y1 + 1)*x1^2 + 2*(y1*y2^2*y3^3 + 3*y1^2 + 2*y1 + 1)*x1 + 3*y1^2 + 2*y1 + 1
                sage: dpol.to_expr(alphabet = alphabet)
                (y1*y2^2*y3^3 + 3*y1^2 + 2*y1 + 1)*a*b^2*c^3 + y1*y2^2*y3^3 + 3*(y1*y2^2*y3^3 + 3*y1^2 + 2*y1 + 1)*a^2 + 2*(y1*y2^2*y3^3 + 3*y1^2 + 2*y1 + 1)*a + 3*y1^2 + 2*y1 + 1
                sage: var('d,e,f')
                (d, e, f)
                sage: alphabety = [d,e,f]
                sage: dpol.to_expr(alphabet = alphabet, alphabety = alphabety)
                (d*e^2*f^3 + 3*d^2 + 2*d + 1)*a*b^2*c^3 + d*e^2*f^3 + 3*(d*e^2*f^3 + 3*d^2 + 2*d + 1)*a^2 + 2*(d*e^2*f^3 + 3*d^2 + 2*d + 1)*a + 3*d^2 + 2*d + 1
            """
            from sage.symbolic.ring import SymbolicRing
            if(hasattr(self.parent().basis_tower(),'_to_expr_on_basis')):
                elt = list(self)
                if(hasattr(elt[0][1],'to_expr')):
                    return sum([coeff.to_expr(alphabet = alphabety)*self.parent().basis_tower()._to_expr_on_basis(key, alphabet = alphabet) for (key,coeff) in self])
                return sum([coeff*self.parent().basis_tower()._to_expr_on_basis(key, alphabet = alphabet) for (key,coeff) in self])

            else:
                pol = self.parent().abstract_algebra().monomial_basis()(self)
                return pol.to_expr(alphabet = alphabet, alphabety = alphabety)
        
        def divided_difference(self, i, otype=None):
            r"""
            Applies the ith divided difference on ``self``
            
            `\partial_i^A = (1 - s_i) \frac{1}{x_i - x_{i+1}}`
            
            `\partial_i^B = (1 - s_i^B) \frac{1}{x_i^{1/2} - x_i^{-1/2}}`
            
            `\partial_i^C = (1 - s_i^C) \frac{1}{x_i - x_i^{-1}}`
            
            `\partial_i^D = (1 - s_i^C) \frac{1}{x_{i-1}^{-1} - x_i}`
            

            INPUT: 

            - ``i``: the number of the divided difference
            - ``otype``: the type of divided difference for untyped basis,
            default is ``A``

            OUTPUT:

            - the polynomial obtained by performing the ith divided 
              difference
              

            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: pol = m[1,2,2]
                sage: pol.divided_difference(1)
                -x[1, 1, 2]
           
            The operation type is given either by the ``otype`` argument or 
            by the type of polynomial.
            
            The monomial basis is an untyped basis, so a type can be specified or
            is ``A`` by default. The applied operation will be `\partial_i^X`
            where ``X`` is the sent type regardless of the number of variables.
              
            ..::  
                sage: pol.divided_difference(1,"B")
                x[-1, 2, 2] + x[0, 2, 2]
             
            If the operation number is greater than the number of variables,
            the polynomial is changed so that operation can apply (only on 
            untyped basis).
            ..:
            
                sage: pol.divided_difference(3)
                x[1, 2, 1, 0] + x[1, 2, 0, 1]
                
            On a typed basis, the polynomial type is used and giving a 
            value to the ``otype`` will raise an exception. The operation 
            that is applied is the one that make sense for the Weyl Group 
            linked to the polynomial. As an example, in type ``B``, with 
            a polynomial in ``n`` variables, `\partial_i^A` will be used 
            for `1 \leq i < n` and `\partial_i^B` for `i=n`.
            
            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb.an_element(); pol
                x(2, 2, 3)
                sage: pol.divided_difference(2)
                -x(2, 2, 2)
                sage: pol.divided_difference(3)
                x(2, 2, 0) + x(2, 2, -3) + x(2, 2, -2) + x(2, 2, -1) + x(2, 2, 1) + x(2, 2, 2)
                sage: pol.divided_difference(3, "C")
                Traceback (most recent call last):
                ...
                TypeError: The type argument is not valid on this basis
                sage: pol.divided_difference(4)     
                Traceback (most recent call last):
                ...
                ValueError: 4 is not a valid operation number
                
            Some other examples in different bases::
            
                sage: K = A.demazure_basis_on_vectors()
                sage: pol = K[1,4,2]
                sage: pol.divided_difference(2)
                K(1, 2, 3)
                sage: KB = A.demazure_basis_on_vectors("B")
                sage: pol = KB[1,4,2]
                sage: pol.divided_difference(3)
                K(1, 4, -2) - K(1, 4, 2)
                sage: Schub = A.schubert_basis_on_vectors()
                sage: pol = Schub[1,4,2]
                sage: pol
                Y(1, 4, 2)
                sage: pol.divided_difference(2)
                Y(1, 2, 3)

            """
            return self.apply_morphism(i, otype, "divided_difference")
        
        def isobaric_divided_difference(self, i, otype=None):
            r"""
            Applies the ith isobaric divided difference on ``self``
            
            `\pi_i = x_i \partial_i`
            
            `\pi_i^B = x_i^{1/2}.\partial_i^B`
            
            `\pi_i^C = x_i.\partial_i^C`
            
            `\pi_i^D = (1 - s_i^D \frac{1}{x_{i-1}x_i})\frax{1}{1- \frac{1}{x_{i-1}x_i}}`
            
            INPUT:

            - ``i``: the number of the divided difference
            - ``otype``: the type of divided difference for untyped basis,
            default is ``A``

            OUTPUT:

            - the polynomial obtained by performing the ith isobaric divided
              difference.

            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: pol = m[2,1,1]
                sage: pol.isobaric_divided_difference(1)
                x[1, 2, 1] + x[2, 1, 1]
           
            The operation type is given either by the ``otype`` argument or 
            by the type of polynomial.
            
            The monomial basis is an untyped basis, so a type can be specified or
            is ``A`` by default. The applied operation will be `\pi_i^X`
            where ``X`` is the sent type regardless of the number of variables.
            
            ..::
                sage: pol.isobaric_divided_difference(1,"B")
                x[-2, 1, 1] + x[-1, 1, 1] + x[1, 1, 1] + x[2, 1, 1] + x[0, 1, 1]               
                
            If the operation number is greater than the number of variables,
            the polynomial is changed so that operation can apply (only on 
            untyped basis).
            
            ..::
                sage: pol.isobaric_divided_difference(3)
                x[2, 1, 1, 0] + x[2, 1, 0, 1]
                
            On a typed basis, the polynomial type is used and giving a 
            value to the ``otype`` will raise an exception. The operation 
            that is applied is the one that make sense for the Weyl Group 
            linked to the polynomial. As an example, in type ``B``, with 
            a polynomial in ``n`` variables, `\pi_i^A` will be used 
            for `1 \leq i < n` and `\pi_i^B` for ``i=n``.
            
            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb[2,1,1]; pol
                x(2, 1, 1)
                sage: pol.isobaric_divided_difference(1)
                x(1, 2, 1) + x(2, 1, 1)
                sage: pol.isobaric_divided_difference(3)
                x(2, 1, 0) + x(2, 1, -1) + x(2, 1, 1)
                sage: pol.isobaric_divided_difference(3, "C")
                Traceback (most recent call last):
                ...
                TypeError: The type argument is not valid on this basis
                sage: pol.isobaric_divided_difference(4)     
                Traceback (most recent call last):
                ...
                ValueError: 4 is not a valid operation number
                
            Some other examples in different bases::
            
                sage: K = A.demazure_basis_on_vectors()
                sage: pol = K[1,4,2]
                sage: pol.isobaric_divided_difference(2)
                K(1, 2, 4)
                sage: KB = A.demazure_basis_on_vectors("B")
                sage: pol = KB[1,4,2]
                sage: pol.isobaric_divided_difference(3)
                K(1, 4, -2)
                sage: Schub = A.schubert_basis_on_vectors()
                sage: pol = Schub[1,4,2]
                sage: pol
                Y(1, 4, 2)
                sage: pol.isobaric_divided_difference(2)
                Y(1, 2, 4)

            """
            return self.apply_morphism(i, otype, "isobaric_divided_difference")
        
        def hat_isobaric_divided_difference(self, i, otype=None):
            r"""
            Applies the ith hat isobaric divided difference on ``self``
            
            `\hat{\pi}_i = \partial_i x_{i+1}`
        
            `\hat{\pi}_i^B = \partial_i^B . x_i^{-1/2}.`
            
            `\hat{\pi}_i^C = \partial_i^C.x_i^{-1}`
            
            `\hat{\pi}_i^D = (1 - s_i^D)\frac{1}{x_{i-1}x_i} - 1}`
            
            INPUT:

            - ``i``: the number of the divided difference
            - ``otype``: the type of divided difference for untyped basis,
            default is ``A``

            OUTPUT:

            - the polynomial obtained by performing the ith hat isobaric divided
              difference.

            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: pol = m[2,1,1]
                sage: pol.hat_isobaric_divided_difference(1)
                x[1, 2, 1]
           
            The operation type is given either by the ``otype`` argument or 
            by the type of polynomial.
            
            The monomial basis is an untyped basis, so a type can be specified or
            is ``A`` by default. The applied operation will be `\hat{\pi}_i^X`
            where ``X`` is the sent type regardless of the number of variables.
            
            ..::
                sage: pol.hat_isobaric_divided_difference(1,"B")
                x[-2, 1, 1] + x[-1, 1, 1] + x[1, 1, 1] + x[0, 1, 1]
                
            If the operation number is greater than the number of variables,
            the polynomial is changed so that operation can apply (only on 
            untyped basis).
            
            ..::
                sage: pol.hat_isobaric_divided_difference(3)
                x[2, 1, 0, 1]
                
            On a typed basis, the polynomial type is used and giving a 
            value to the ``otype`` will raise an exception. The operation 
            that is applied is the one that make sense for the Weyl Group 
            linked to the polynomial. As an example, in type ``B``, with 
            a polynomial in ``n`` variables, `\hat{\pi}_i^A` will be used 
            for `1 \leq i < n` and `\hat{\pi}_i^B` for ``i=n``.
            
            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb[2,1,1]; pol
                x(2, 1, 1)
                sage: pol.hat_isobaric_divided_difference(1)
                x(1, 2, 1)
                sage: pol.hat_isobaric_divided_difference(3)
                x(2, 1, 0) + x(2, 1, -1)
                sage: pol.hat_isobaric_divided_difference(3, "C")
                Traceback (most recent call last):
                ...
                TypeError: The type argument is not valid on this basis
                sage: pol.hat_isobaric_divided_difference(4)     
                Traceback (most recent call last):
                ...
                ValueError: 4 is not a valid operation number
                
            Some other examples in different bases::
            
                sage: K = A.demazure_basis_on_vectors()
                sage: pol = K[1,4,2]
                sage: pol.hat_isobaric_divided_difference(2)
                K(1, 2, 4) - K(1, 4, 2)
                sage: KB = A.demazure_basis_on_vectors("B")
                sage: pol = KB[1,4,2]
                sage: pol.hat_isobaric_divided_difference(3)
                K(1, 4, -2) - K(1, 4, 2)
                sage: Schub = A.schubert_basis_on_vectors()
                sage: pol = Schub[1,4,2]
                sage: pol
                Y(1, 4, 2)
                sage: pol.hat_isobaric_divided_difference(2)
                Y(1, 2, 4) - Y(1, 4, 2)

            """
            return self.apply_morphism(i, otype, "hat_isobaric_divided_difference")
            
        def si(self, i, otype=None):
            r"""
            Applies the ``si`` action on ``self``.
            ``si`` is an action on the exponent vectors of a monomial.
        
            `v.s_i = [v_1, ..., v_{i+1}, v_{i}, ..., v_n]`
            
            `v.s_i^B = v.s_i^C =  [v_1, ..., -v_i, ..., v_n]`
            
            `v.s_i^D = [v_1, ..., -v_{i}, -v_{i-1}, ..., v_n]`
            
            
            INPUT:

            - ``i``: the number of the divided difference
            - ``otype``: the type of the action for untyped basis,
            default is ``A``

            OUTPUT:

            - the polynomial obtained by performing the ``si`` action.

            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: pol = m[2,1,1]
                sage: pol.si(1)
                x[1, 2, 1]
           
            The operation type is given either by the ``otype`` argument or 
            by the type of polynomial.
            
            The monomial basis is an untyped basis, so a type can be specified or
            is ``A`` by default. The applied operation will be `s_i^X`
            where ``X`` is the sent type regardless of the number of variables.
            
            ..::
                sage: pol.si(1,"B")
                x[-2, 1, 1]
                
            If the operation number is greater than the number of variables,
            the polynomial is changed so that operation can apply (only on 
            untyped basis).
            
            ..::
                sage: pol.si(3)
                x[2, 1, 0, 1]
                
            On a typed basis, the polynomial type is used and giving a 
            value to the ``otype`` will raise an exception. The operation 
            that is applied is the Weyl Group action linked to the polynomial. 
            As an example, in type ``B``, with a polynomial in ``n`` variables,
             `si_i^A` will be used for `1 \leq i < n` and `s_i^B` for ``i=n``.
            
            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb[2,1,1]; pol
                x(2, 1, 1)
                sage: pol.si(1)
                x(1, 2, 1)
                sage: pol.si(3)
                x(2, 1, -1)
                sage: pol.si(3, "C")
                Traceback (most recent call last):
                ...
                TypeError: The type argument is not valid on this basis
                sage: pol.hat_isobaric_divided_difference(4)     
                Traceback (most recent call last):
                ...
                ValueError: 4 is not a valid operation number
            
            Some other examples in different bases::
            
                sage: K = A.demazure_basis_on_vectors()
                sage: pol = K[1,4,2]
                sage: pol = K[1,4,2] 
                sage: pol.si(2)
                K(1, 2, 4) - K(1, 3, 3) - K(1, 4, 2) + K(3, 1, 3) + K(4, 1, 2)
                sage: KB = A.demazure_basis_on_vectors("B")
                sage: pol = KB[1,4,2]
                sage: pol.si(3)
                K(1, 4, -2) - K(1, 4, -1) - K(1, 4, 2) - K(2, 4, 0)
                sage: Schub = A.schubert_basis_on_vectors()
                sage: pol = Schub[1,4,2]
                sage: pol
                Y(1, 4, 2)
                sage: pol.si(2)
                Y(1, 2, 4) - Y(1, 3, 3) - Y(1, 4, 2) + Y(3, 1, 3)

            """
            return self.apply_morphism(i, otype, "si")
            
        def hecke_generator(self, i, t1=None, t2=None):
            r"""
            Applies the ith Hecke algebra generator.
            
            `T_i = \pi_i (t_1 + t_2) - s_i t_2`
            
            This is only implemented in type ``A``.
            
            INPUT:

            - ``i``: the number of the operation
            - ``t_1`` the first Hecke algebra parameter
            - ``t_2`` the second Hecke algebra parameter

            OUTPUT:

            - the polynomial obtained by performing the ith ``T_i`` operator

            EXAMPLES::

                sage: K.<t1,t2> = QQ[]
                sage: A = AbstractPolynomialRing(K)                                         
                sage: pol = A.an_element(); pol
                x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
                sage: pol.hecke_generator(2)
                t1*x[0, 0, 0] + 2*t1*x[1, 0, 0] + (-t2)*x[1, 3, 2] + 3*t1*x[2, 0, 0]
                
            If the operation number is greater than the number of variables,
            the polynomial is changed so that operation can apply (only on 
            untyped basis).
            
            ..::
                sage: pol.hecke_generator(3)
                t1*x[0, 0, 0, 0] + 2*t1*x[1, 0, 0, 0] + (t1+t2)*x[1, 2, 1, 2] + (t1+t2)*x[1, 2, 2, 1] + (t1+t2)*x[1, 2, 3, 0] + t1*x[1, 2, 0, 3] + 3*t1*x[2, 0, 0, 0]
                
            The operation is only defined in type ``A``, if used on typed
            basis different from ``A``, an exception is raised.
            
            ..:
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb.an_element()
                sage: pol.hecke_generator(1)
                Traceback (most recent call last):
                ...
                NotImplementedError: The hecke algebra operator is only implemented in type A
                
            Some other examples in different bases:
            
                sage: Dem = A.demazure_basis_on_vectors()
                sage: pol = Dem[1,4,2] 
                sage: pol.hecke_generator(2)
                t1*K(1, 2, 4) + t2*K(1, 3, 3) + t2*K(1, 4, 2) + (-t2)*K(3, 1, 3) + (-t2)*K(4, 1, 2)
                sage: Schub = A.schubert_basis_on_vectors()
                sage: pol = Schub[1,4,2]
                sage: pol.hecke_generator(2)
                t1*Y(1, 2, 4) + t2*Y(1, 3, 3) + t2*Y(1, 4, 2) + (-t2)*Y(3, 1, 3)
            """
            if(t1 is None): t1 = self.base_ring()(var('t1'))
            if(t2 is None): t2 = self.base_ring()(var('t2'))
            return self.apply_morphism(i, method="hecke_generator",t1=t1,t2=t2)
                
        def apply_morphism(self, i, otype=None, method="divided_difference", **keywords):
            r"""
            Uses the method ``get_morphism`` on the polynomial parent to apply
            the divided difference sent in ``method`` on ``self``
            
            INPUT: 

            - ``i``: the number of the operation, 
            on the monomial basis, it should only be greater than 0
            on typed basis like the ambient space basis, it should be greater
            than 0 and smaller than the number of variables.
            - ``otype``: the type of divided difference if needed by the 
            parent. On the monomial basis, the default value ``A`` will be
           used.
            - ``method``, the method name to apply, here are the names that are
            accepted:
                - ``divided_difference`` or ``d`` for the usual divided difference
                - ``isobaric_divided_difference`` or ``pi`` for the isobaric 
                divided difference
                - ``hat_isobaric_divided_difference`` or ``hatpi`` for the
                hat isobaric divided difference
                - ``si``  or ``s`` for the action of the Weyl group generator
                - ``product_variable`` or ``x`` to apply a product by the ``ith``
                variable
                - ``hecke_generator`` for the hecke algebra generator
                
            OUPUT:
            The result of applying the method on ``self``
            
            If the method is not implemented by the basis of ``self``, then
            a ``NotImplementedError`` exception is raised.
            
            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: pol = A.an_element() 
                sage: pol.apply_morphism(1)           
                2*x[0, 0, 0] + 3*x[1, 0, 0] - x[1, 1, 3] + 3*x[0, 1, 0]
                sage: pol.apply_morphism(1, method="pi")
                x[0, 0, 0] + 2*x[1, 0, 0] + 3*x[1, 1, 0] + 3*x[2, 0, 0] + 2*x[0, 1, 0] + 3*x[0, 2, 0]
                sage: pol.apply_morphism(1, method="hatpi")
                3*x[1, 1, 0] - x[1, 2, 3] + 2*x[0, 1, 0] + 3*x[0, 2, 0]

            The ``otype`` argument can be used on untyped polynomials (like
            elements of the monomial basis).
            
            If used on a typed element, a ``TypeError`` Exception is raised.
            
            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: pol = A.an_element()
                sage: pol.apply_morphism(1, otype="B")
                5*x[0, 0, 0] + 3*x[-2, 0, 0] + 5*x[-1, 0, 0] + x[-1, 2, 3] + 3*x[1, 0, 0] + x[0, 2, 3]
                sage: ma = A.ambient_space_basis("A")
                sage: pol = ma(pol); pol
                x(0, 0, 0) + 2*x(1, 0, 0) + x(1, 2, 3) + 3*x(2, 0, 0)
                sage: pol.apply_morphism(1, otype="B")
                Traceback (most recent call last):
                ...
                TypeError: The type argument is not valid on this basis
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb(pol)
                sage: pol.apply_morphism(1)          
                2*x(0, 0, 0) + 3*x(1, 0, 0) - x(1, 1, 3) + 3*x(0, 1, 0)
                sage: pol.apply_morphism(3)
                x(1, 2, 0) + x(1, 2, -3) + x(1, 2, -2) + x(1, 2, -1) + x(1, 2, 1) + x(1, 2, 2)
            """
            pol = self._right_number_of_variables_pol(i,otype)
            morph = pol.parent().get_morphism(i,otype,method, **keywords)
            return morph(pol)

        def apply_composed_morphism(self, operation_list):
            r"""
            Applies a composed morhism on ``self`` from an operation list.
            
            The operations are given through a list of tuples (or list) of the 
            form ``(operation,number,otype)`` 
            ``operation`` coressponds to a keyword string:
            - ``"x"`` : product by a variable
            - ``"d"`` : divided difference
            - ``"s"`` : action of ``si``
            - ``"pi"``: isobaric divided difference
            - ``"hatpi"``: hat isobaric divided difference
            
            ``number`` is the argument of the method, ``otype`` is the method
            type (optional and only for non typed basis).
            
            A endomorphism is created by composing all the operations (acting
            on the left).
            
            INPUT:
            - ``operation_list`` a list of tuples representing the operations
            
            OUTPUT:
            the polynomial optained by applying the composed operation on
            ``self``. Each operation is acting  on the left, and so the 
            operations are applied from the left to the right.
            
            EXAMPLES::
                
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: m = A.monomial_basis()
                sage: pol = m[1,2,3]
                sage: pol.apply_composed_morphism([("d",1),("d",2)])
                x[1, 1, 2] + x[1, 2, 1]
                sage: pol = m[2,1,1]
                sage: pol.apply_composed_morphism([("x",1),("d",1)])
                x[1, 2, 1] + x[2, 1, 1]
                sage: pol.apply_composed_morphism([("pi",1)])       
                x[1, 2, 1] + x[2, 1, 1]
                sage: pol.apply_composed_morphism([("pi",2,"B"),("s",1,"C"),("hatpi",2,"D")])
                -x[-2, -1, 1] - x[-2, 1, 1] - x[-2, 0, 1] - x[-1, 1, 1] - x[-1, 0, 1] - x[0, 1, 1]
                sage: pol.apply_composed_morphism([("pi",3),("d",2,"C")])
                x[2, 0, 1, 0] + x[2, 0, 0, 1]
                
            On a typed basis, no type parameter should be sent.
            ..::
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb.an_element()
                sage: pol                                
                x(2, 2, 3)
                sage: pol.apply_composed_morphism([("d",2),("pi",1),("hatpi",3)])
                -x(2, 2, 0) - x(2, 2, -2) - x(2, 2, -1) - x(2, 2, 1)
                sage: pol.apply_composed_morphism([("d",2,"C"),("pi",1),("hatpi",3)])
                Traceback (most recent call last):
                ...
                TypeError: The type argument is not valid on this basis
             """
            pol = self
            for op in operation_list:
                if(len(op)>2):
                    pol = pol._right_number_of_variables_pol(op[1],otype = op[2])
                else:
                    pol = pol._right_number_of_variables_pol(op[1])

            morph = pol.parent().get_composed_morphism(operation_list)
            return morph(pol)
            
        def apply_reduced_word(self, word, method="divided_difference"):
            r"""
            Apply an operation indexed by a reduced word. 
            
            INPUT:
            - ``word``, a reduced word 
            - ``method`` the method to apply : 
                - ``divided_difference`` or ``d`` for the usual divided difference
                - ``isobaric_divided_difference`` or ``pi`` for the isobaric 
                divided difference
                - ``hat_isobaric_divided_difference`` or ``hatpi`` for the
                hat isobaric divided difference
                - ``si``  or ``s`` for the action of the Weyl group generator
                - ``product_variable`` or ``x`` to apply a product by the ``ith``
                variable
                - ``hecke_generator`` for the hecke algebra generator
                
            On a untyped basis (like the monomial basis), the reduced word
            should be of type ``A``. On typed basis, it a reduced word of the
            corresponding Weyl group.
            
            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: pol = A.an_element();pol
                x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
                sage: pol.apply_reduced_word(Permutation([2,3,1]).reduced_word())      
                3*x[0, 0, 0] + x[1, 1, 2] + x[1, 2, 1]
                sage: pol.apply_reduced_word(Permutation([2,3,1]).reduced_word(), method="pi")
                x[0, 0, 0] + 2*x[1, 0, 0] + 3*x[1, 1, 0] + 3*x[1, 0, 1] + 3*x[2, 0, 0] + 2*x[0, 1, 0] + 3*x[0, 1, 1] + 3*x[0, 2, 0] + 2*x[0, 0, 1] + 3*x[0, 0, 2]
                sage: pol.apply_reduced_word(Permutation([2,3,1]).reduced_word(), method="si")
                x[0, 0, 0] + x[2, 3, 1] + 2*x[0, 0, 1] + 3*x[0, 0, 2]
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb(pol)
                sage: W = WeylGroup("B3")
                sage: word = W.an_element().reduced_word()
                sage: word
                [1, 2, 3]
                sage: pol.apply_reduced_word(word)                                            
                x(1, 1, 0) + x(1, 1, -2) + x(1, 1, -1) + x(1, 1, 1) + x(1, 2, 0) + x(1, 2, -1)
                sage: pol.apply_reduced_word(word, method="pi")
                6*x(0, 0, 0) + 5*x(1, 0, 0) + 3*x(1, 1, 0) + 3*x(1, 0, -1) + 3*x(1, 0, 1) + 3*x(2, 0, 0) + 5*x(0, 1, 0) + 3*x(0, 1, -1) + 3*x(0, 1, 1) + 3*x(0, 2, 0) + 3*x(0, 0, -2) + 5*x(0, 0, -1) + 5*x(0, 0, 1) + 3*x(0, 0, 2)
                sage: pol.apply_reduced_word(word, method="si")
                x(0, 0, 0) + x(2, 3, -1) + 3*x(0, 0, -2) + 2*x(0, 0, -1)
            
            """
            list = [(method,i) for i in word]
            return self.apply_composed_morphism(list)
       
        def _right_number_of_variables_pol(self, i, otype=None):
            r"""
            This methods verifies that the number of variables of a polynomial
            is appropriate for an operation.
            
            The generique implementation only raises an Exception if this number
            is not appropriate, but it can be overriten by subclasses to actually
            return a different polynomial by changing the number of variables.
            
            INPUT:
            - ``i``, the operation number
            - ``otype``, the operation type. By default this is the polynomial 
            group type or ``A`` for untyped polynomials
            
            OUTPUT:
            the polynomial itself, or a equivalent polynomial in a different
            number of variables
            
            EXAMPLES::
            
                sage: A = AbstractPolynomialRing(QQ)  
                sage: ma = A.ambient_space_basis("A")
                sage: pol = ma.an_element()
                sage: pol._right_number_of_variables_pol(2)
                x(2, 2, 3)
                sage: pol._right_number_of_variables_pol(0)
                Traceback (most recent call last):
                ...
                ValueError: 0 is not a valid operation number
                sage: pol._right_number_of_variables_pol(3)
                Traceback (most recent call last):
                ...
                ValueError: 3 is not a valid operation number
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb(pol)
                sage: pol._right_number_of_variables_pol(2)
                x(2, 2, 3)
                sage: pol._right_number_of_variables_pol(0)
                Traceback (most recent call last):
                ...
                ValueError: 0 is not a valid operation number
                sage: pol._right_number_of_variables_pol(3)
                x(2, 2, 3)
                
                TESTS::
                sage: mc = A.ambient_space_basis("C")
                sage: pol = mc(pol)
                sage: pol._right_number_of_variables_pol(2)
                x(2, 2, 3)
                sage: pol._right_number_of_variables_pol(0)
                Traceback (most recent call last):
                ...
                ValueError: 0 is not a valid operation number
                sage: pol._right_number_of_variables_pol(3)
                x(2, 2, 3)
                sage: md = A.ambient_space_basis("D")
                sage: pol = md(pol)
                sage: pol._right_number_of_variables_pol(2)
                x(2, 2, 3)
                sage: pol._right_number_of_variables_pol(0)
                Traceback (most recent call last):
                ...
                ValueError: 0 is not a valid operation number
                sage: pol._right_number_of_variables_pol(3)
                x(2, 2, 3)
                sage: pol._right_number_of_variables_pol(1)
                x(2, 2, 3)
            """
            min = self.parent().basis_tower()._right_number_of_variables(i, otype)
            if(self.nb_variables() < min):
                raise ValueError, "%s is not a valid operation number"%(i)
            return self
        
        def group_type(self):
            r"""
            Returns the polynomial group type, or an exception in case
            of untyped polynomials.
            
            EXAMPLES::
            
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: m = A.monomial_basis()
                sage: pol = m.an_element()
                sage: pol.group_type()
                Traceback (most recent call last):
                ...
                AttributeError: This polynomial has no group type
                sage: ma = A.ambient_space_basis("A")
                sage: pol = ma.an_element()
                sage: pol.group_type()
                'A'                                     
                sage: K = A.demazure_basis_on_vectors()
                sage: pol = K.an_element()
                sage: pol.group_type()
                'A'

            """
            try:
                return self.parent().group_type()
            except:
                raise AttributeError, "This polynomial has no group type"%()


class PolynomialRingWithBasis(UniqueRepresentation, Parent):
    """
    This class is the upperclass of all AbstractPolynomialRing bases. 
    It is not supposed to be directly called, its main subclasses are :
     - ``sage.combinat.multivariate_polynomials.basis.PolynomialRingWithBasisFromMorphism``
     - ``sage.combinat.multivariate_polynomials.monomial.MonomialBasis``
     - ``sage.combinat.multivariate_polynomials.ambient_space_basis.PolynomialRingWithBasisFromAmbientSpace``
     - ``sage.combinat.multivariate_polynomials.linear_basis_on_vectors.LinearBasisOnVectors``
     
    The number of variables is not set : to create an actual polynomial, the
    class is using the ``finite_basis`` method to create its finite basis on a
    given number of variables which inherit from ``FinitePolynomialRingWithBasis``

    INPUT:

    - ``abstract_polynomial_ring``: the abstract algebra of type ``AbstractPolynomialRing`` which is
    a facade for all polynomial bases
    - ``basis_name``: The name of the basis (example : "Monomial basis")
    - ``neutral_nb_variables``: the number of variables used to create the one element,
    - ``basis_repr`` : the string used to represent element of the basis (example : "x" for the Monomial Basis)
    - ``category`` : category of ``self``, by default, ``abstract_polynomial_ring.bases_category_class()`` 
    is used
    - ``variable_auto_coerce``, if True, a trivial embeding will be created between polynomial of this basis 
    having different number of variables (by adding zeros to the index vector)

    EXAMPLES::

        sage: A = AbstractPolynomialRing(QQ)
        sage: m = A.monomial_basis()
        sage: m
        The ring of multivariate polynomials on x over Rational Field on the monomial basis
        sage: Schub = A.schubert_basis_on_vectors()
        sage: Schub
        The ring of multivariate polynomials on x over Rational Field on the Schubert basis of type A (indexed by vectors)
    
    ``Schub`` and ``m`` are instances of subclasses of ``PolynomialRingWithBasis`` with abstract 
    polynomial ring ``A``

    """   
    def __init__(self, abstract_polynomial_ring, basis_name, neutral_nb_variables, basis_repr = None, category = None, variables_auto_coerce = True):
        r"""
        TESTS::
            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: m
            The ring of multivariate polynomials on x over Rational Field on the monomial basis
        """
        self._basis_name = basis_name
        self._neutral_nb_variables = neutral_nb_variables
        self._variables_auto_coerce = variables_auto_coerce
        if(category is None): category = abstract_polynomial_ring.bases_category_class()(abstract_polynomial_ring)
        Parent.__init__(
            self,
            base = abstract_polynomial_ring.base(),
            category = category
        )
        self._facade_for = set([])
        self._basis_repr = basis_repr
        m = SetMorphism( Hom(self, abstract_polynomial_ring), lambda x: x)
        m.register_as_coercion()
        

    def _repr_(self):
        r"""
        
        Print the name of ``self``
        
        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: m
            The ring of multivariate polynomials on x over Rational Field on the monomial basis
        """
        return "%s on the %s"%(self.abstract_algebra()._repr_with_basis(), self._basis_name)
    
    def basis_repr(self): 
        r"""
        OUTPUT:
            - the string representing the elements of the basis
        
        EXAMPLES::
            
            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: m.basis_repr()
            'x'
            sage: m.an_element()
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
            sage: Schub = A.schubert_basis_on_vectors()
            sage: Schub.basis_repr()
            'Y'
            sage: Schub.an_element()
            Y(2, 2, 3)
            
        """
        return self._basis_repr

    def __call__(self, obj):
        r"""
        The call function is used to create elements of the basis from 
        an object.
        If a polynomial from another basis is given, it is sent to the finite
        basis corresponding to the number of variables in the polynomial, this way
        the coercion system that exist between finite bases can work. (There is no
        automatic coercion between bases with an unset number of variables)
        
        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: Schub = A.schubert_basis_on_vectors()
            sage: m(Schub.an_element())
            x[2, 2, 3] + x[2, 3, 2] + x[3, 2, 2]

        The coercion exists only between the monomial basis on 3 variables and the
        Schubert basis on 3 variables. So when the monomial basis on all variables
        receive the polynomial, it gives it to its finite basis.
        
        If the given element type is the key type of the ``CombinatorialFreeModule`` 
        used in the finite basis, the element is created. As the non finite basis is 
        not a ``CombinatorialFreeModule`` itself, this would not without this test.
        
        EXAMPLES::
        
            sage: pol = m.an_element(); pol
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
            sage: key = list(pol)[0][0]
            sage: key
            [1, 2, 3]
            sage: m(key)
            x[1, 2, 3]
            
        As ``m`` is not itself a ``CombinatorialFreeModule``, it needs to identify ``key``
        as a key of its finite basis.
        
        If a list is sent, then the method uses the __getitem__method.
        
        EXAMPLES::
        
            sage: A = AbstractPolynomialRing(QQ)                                       
            sage: m = A.monomial_basis()  
            sage: m([2,2,3])
            x[2, 2, 3]

        """
        if( hasattr(obj, 'nb_variables')):
            return self.finite_basis(obj.nb_variables())(obj)     
        if( hasattr(obj, 'parent') and type( obj.parent()) == type(self._default_finite_basis()._basis_keys) ):
            return self.term(obj)
        if(type(obj) is list or type(obj) is tuple):
            return self.__getitem__(obj) 
        return super(PolynomialRingWithBasis, self).__call__(obj)
        
    def __getitem__(self, vect):
        r"""
        Allows the creation of elements with [ ]. The method computes the
        length of the sent vector and gives it to a finite basis of the right
        size.

        TESTS::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis();
            sage: m[2,2,1]
            x[2, 2, 1]

        """
        if( hasattr(vect, '__iter__')):
            size = len(vect)
        else:
            size=1
            vect = [vect]
        return self.finite_basis(size)( list(vect) )

    def finite_basis(self, nb_variables):
        r"""
        Creates the basis corresponding to ``self`` with ``nb_variables`` variables
        
        INPUT:
            - ``nb_variables``, the number of variables
            
        OUTPUT:
            - the basis corresponding to ``self`` with ``nb_variables`` variables
        
        TESTS::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: m3 = m.finite_basis(3)
            sage: m
            The ring of multivariate polynomials on x over Rational Field on the monomial basis
            sage: m3
            The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis

        """
        fb = self._finite_basis_instance(nb_variables)
        
        return fb

    def _default_finite_basis(self):
        r"""
        
        OUTPUT:
        - a default finite basis used to create the ``one`` element 
        
        TESTS::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: m._default_finite_basis()
            The ring of multivariate polynomials on x over Rational Field with 1 variable on the monomial basis

        """
        return self.finite_basis(self._neutral_nb_variables) 
        
    def _right_number_of_variables(self, i, otype=None):
        r"""
            This methods returns the minimal number of variables needed to apply
            an operation. If the operation number is not correct for any 
            number of variables, then an exception is raised.
            
            INPUT:
            - ``i``, the operation number
            - ``otype``, the operation type, this should be send only on 
            untyped basis, otherwise it will be ignored. 
            If set to ``None`` on a untyped basis, the default value is 
            ``A``
            
            OUTPUT:
            
            the minimal number of variables needed to apply the ``ith`` 
            divided difference opperator (or other predefined operations)
            
            Here is an example on a untyped basis:
            
            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                                             
                sage: m = A.monomial_basis()
                sage: m._right_number_of_variables(2)
                3
                sage: m._right_number_of_variables(2,"B")
                2
                sage: m._right_number_of_variables(0)
                Traceback (most recent call last):
                ...
                ValueError: 0 is not a valid operation number
                sage: m._right_number_of_variables(1,"D")
                Traceback (most recent call last):
                ...
                ValueError: 1 is not a valid operation number
                
            And now, on typed basis:
            ..::
                sage: ma = A.ambient_space_basis("A")
                sage: ma._right_number_of_variables(2)
                3
                sage: ma._right_number_of_variables(0)
                Traceback (most recent call last):
                ...
                ValueError: 0 is not a valid operation number
                sage: mb = A.ambient_space_basis("B")
                sage: mb._right_number_of_variables(2)
                2
                sage: md = A.ambient_space_basis("D")
                sage: md._right_number_of_variables(2)
                2
                sage: md._right_number_of_variables(1)
                2
                
            TESTS::
                sage: ma = A.ambient_space_basis("A")
                sage: ma._right_number_of_variables(0)
                Traceback (most recent call last):
                ...
                ValueError: 0 is not a valid operation number
                sage: ma._right_number_of_variables(3)
                4
                sage: mb = A.ambient_space_basis("B")
                sage: mb._right_number_of_variables(0)
                Traceback (most recent call last):
                ...
                ValueError: 0 is not a valid operation number
                sage: mb._right_number_of_variables(3)
                3
                sage: mc = A.ambient_space_basis("C")
                sage: mc._right_number_of_variables(0)
                Traceback (most recent call last):
                ...
                ValueError: 0 is not a valid operation number
                sage: mc._right_number_of_variables(3)
                3
                sage: md = A.ambient_space_basis("D")  
                sage: md._right_number_of_variables(1)
                2
                sage: md._right_number_of_variables(0)
                Traceback (most recent call last):
                ...
                ValueError: 0 is not a valid operation number
                sage: md._right_number_of_variables(3)
                3
        """
        if(hasattr(self, "group_type")):
            otype = self.group_type()
            if(i<=0):
                raise ValueError, "%s is not a valid operation number"%(i)   
            if(otype=="D" and i==1): return 2
            if(otype=="A"): return i+1
            return i
        else:
            if(otype is None): otype ="A"
            if(i<=0 or (otype=="D" and i==1)):
                raise ValueError, "%s is not a valid operation number"%(i)            
            if(otype=="A"): return i+1
            return i
   
        
    def _register_finite_basis(self, finite_basis):
        r"""
        This method registers a finite basis as a facade of the main basis.
        It also creates the coercion between polynomials on different number 
        of variables. It is called by the initialisation method of 
        ``FinitePolynomialRingWithBasis``.
        
        The method uses the parameter ``self._variables_auto_coerce``. It is
        set to ``True`` when indexes of the basis elements can be extended with 
        zeros to give indexes for a basis on a greater number of variables.
        
        As an example, on the monomial basis, ``x[2,2,1]`` is equal to ``x[2,2,1,0]``.
        When such an embeding exists, trivial coercions are created between 
        polynomials on different number of variables.
        
        On some basis (like the Macdonald basis), this embeding is not true, then
        ``_variables_auto_coerce`` is set to False, and the coercion between 
        polynomials on different number of variables is done throught the Monomial
        Basis.
        
        INPUT:
        
        - the final basis to register, this must be a finite basis for ``self``
        
        EXAMPLES::
        
            sage: K.<t1,t2,q> = QQ[]
            sage: K = K.fraction_field()
            sage: A = AbstractPolynomialRing(K) 
            sage: m = A.monomial_basis()
            sage: pol = m[2,3]; pol
            x[2, 3]
            sage: m3 = m.finite_basis(3)
            sage: m3(pol)
            x[2, 3, 0]
            sage: Mac = A.macdonald_basis_on_vectors()
            sage: pol = Mac[1]; pol
            M(1)
            sage: Mac2 = Mac.finite_basis(2)
            sage: Mac2(pol)
            (t1-1)*M(0, 0) + (1/(-t2))*M(1, 0) + ((-t1-t2)/(t1*q+t2))*M(0, 1)
        
        """
        if(self._variables_auto_coerce):
            for f2 in self._facade_for:
                self._create_morphism(finite_basis, f2)
        
                
        self._facade_for.add(finite_basis)
                
    def _create_morphism(self,f1,f2):
        r"""
        Creates a trivial module morphism between two realizations of a basis
        on different number of variables. The morphism only extend the index 
        vector with zeros. This method is called by ``_register_finite_basis``
        when needed.
        
        INPUT:
        - ``f1``,``f2`` two finite bases of ``self`` with a different number of
        variables. The number of variables should only be different, the order 
        between ``f1`` and ``f2`` doesn't matter.
        
        EXAMPLES::
        
            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: m2 = m.finite_basis(2)
            sage: m3 = m.finite_basis(3)
            sage: m2.an_element()
            x[0, 0] + 2*x[1, 0] + x[1, 2] + 3*x[2, 0]
            sage: m3(m2.an_element())
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 0] + 3*x[2, 0, 0]
        
        """
        if(f1.nb_variables() > f2.nb_variables()):
            f1, f2 = f2, f1
        
        morph = f1._module_morphism(
            lambda key : f2([key[i] for i in xrange(f2.nb_variables())]),
            codomain = f2
        )
        
        morph.register_as_coercion()

    def an_element(self):
        r"""
        Returns an element of ``self``. By default, it takes an element
        of the finite basis on 3 variables.

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: m.an_element()
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]

        """
        return self.finite_basis(3).an_element()    

    def term(self, key, coeff = None):
        r"""
        Return the term of ``self`` indexed by ``key``. ``self`` is not
        a ``CominatorialFreeModule``, it has to compute the size of ``key``
        to know which finite basis corresponds to ``key``.

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: pol = m.an_element()
            sage: pol
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
            sage: list(pol)
            [([1, 2, 3], 1), ([1, 0, 0], 2), ([2, 0, 0], 3), ([0, 0, 0], 1)]
            sage: key = list(pol)[0][0]; key
            [1, 2, 3]
            sage: m.term(key)
            x[1, 2, 3]

        """
        nb_variables = len(key.parent()._basis_keys)
        return self.finite_basis(nb_variables).term(key, coeff)
    
    def add_operator(self, name, method):
        r"""
        Adds dynamically a operator method to the basis. Then the operator
        can be used as a morphism by the polynomial through the methods 
        ``apply_morphism`` or ``apply_composed_morphism``.
        
        INPUT:
            - ``name`` the name of the operator
            - ``method``, an on basis method to create the morphism
            the method should have the following signature:
            ``def method(self, key)``
            the ``self`` corresponds to the inner wrapper of the basis,
            it contains at least two attributes :
                - ``self._i`` a integer between 1 and ``n`` to specify which
                divided difference is applying
                - ``self._module`` the basis it is applying to.
            Depending on the basis, more parameters can be available, like
            type, or parameters of the Hecke algebra operator. See the 
            documentation of the specific basis for more details.
            ``self`` also contains on basis methods for other operators
            that are defined on the basis.
            
            The parameter ``key`` is an index of an element of the finite 
            basis. Its parent and type depends on the basis.
            
        EXAMPLES::
            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()                                    
            sage: def affine(self,key): return self.divided_difference_on_basis(key) - self._module(key)
            sage: m.add_operator("affDiff",affine)
            sage: pol = m.an_element(); pol
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
            sage: pol.isobaric_divided_difference(2)
            x[0, 0, 0] + 2*x[1, 0, 0] + 3*x[2, 0, 0]
            sage: pol.divided_difference(1) - pol
            x[0, 0, 0] + x[1, 0, 0] - x[1, 1, 3] - x[1, 2, 3] - 3*x[2, 0, 0] + 3*x[0, 1, 0]
            sage: pol.divided_difference(1,otype="C")
            2*x[0, 0, 0] + 3*x[-1, 0, 0] + 3*x[1, 0, 0] + x[0, 2, 3]
            sage: pol2 = pol.divided_difference(1,otype="C") - pol
            sage: pol2
            x[0, 0, 0] + 3*x[-1, 0, 0] + x[1, 0, 0] - x[1, 2, 3] - 3*x[2, 0, 0] + x[0, 2, 3]
            sage: pol.apply_morphism(1, method="affDiff")
            x[0, 0, 0] + x[1, 0, 0] - x[1, 1, 3] - x[1, 2, 3] - 3*x[2, 0, 0] + 3*x[0, 1, 0]
            sage: pol2 = pol.apply_morphism(1,method="affDiff", otype="C")
            sage: pol2
            x[0, 0, 0] + 3*x[-1, 0, 0] + x[1, 0, 0] - x[1, 2, 3] - 3*x[2, 0, 0] + x[0, 2, 3]
            sage: ma = A.ambient_space_basis("A")
            sage: pol2 = ma(pol2)
            sage: pol2
            x(0, 0, 0) + 3*x(-1, 0, 0) + x(1, 0, 0) - x(1, 2, 3) - 3*x(2, 0, 0) + x(0, 2, 3)
            sage: pol2 = pol2.isobaric_divided_difference(2)
            sage: pol2
            x(0, 0, 0) + 3*x(-1, 0, 0) + x(1, 0, 0) - 3*x(2, 0, 0)
            sage: pol2 = m(pol2)
            sage: pol2
            x[0, 0, 0] + 3*x[-1, 0, 0] + x[1, 0, 0] - 3*x[2, 0, 0]
            sage: pol.apply_composed_morphism([("affDiff",1,"C"),("pi",2)])            
            x[0, 0, 0] + 3*x[-1, 0, 0] + x[1, 0, 0] - 3*x[2, 0, 0]
            
        In a typed basis, like the ambient space basis, the key is an 
        element of the ambient space..::
        
            sage: A = AbstractPolynomialRing(QQ)                                         
            sage: ma = A.ambient_space_basis("A")
            sage: def affine(self,key): return self._module(key.weyl_action(reduced_word=[self._i])) - self._module(key)
            sage: ma.add_operator("affSi",affine)
            sage: pol = ma.an_element();pol
            x(2, 2, 3)
            sage: pol.apply_morphism(2, method="affSi")
            -x(2, 2, 3) + x(2, 3, 2)

        """
        setattr(self._divided_difference_wrapper, name +"_on_basis", method)
        
    class _divided_difference_wrapper(UniqueRepresentation):
        r"""
        This inner class is used to wrap the divided difference on basis
        methods. The morphism module on basis method of 
        ``CombinatorialFreeModule`` only takes the element key as a 
        parameter. But we also need to know the number of the divided difference
        as well the basis we are applying it to.
        
        This wrapper does not contain any methods by default, but they can 
        be added by the subclasses or dynamically through the method ``self.add_operator``
        
        """

        def __init__(self, module, i):
            r"""
            INPUT:
            - ``module``, a ``CombinatorialFreeModule``, the
            finite basis on which the morphism apply
            - ``i`` the number of the wanted divided difference
            
            TESTS::

                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: m3 = m.finite_basis(3)
                sage: wrapp = m._divided_difference_wrapper(m3,1,"A")
            """
            self._module = module
            self._i = i
            

class FinitePolynomialRingWithBasis(CombinatorialFreeModule):   
    r"""
    This class is the upperclass of all finite basis. It corresponds to 
    finite realizations of basis of ``AbstractPolynomialRing`` and inherit
    from ``CombinatorialFreeModule``
    
    The finie basis ``basisName`` in n variables is a realization of both 
    ``multivariate_polynomials.FiniteAbstractPolynomialRing`` in n variables
    anf the basis ``basisName`` on a unset number of variables.
    
    It is not supposed to be directly called, its main subclasses are :
     - ``sage.combinat.multivariate_polynomials.basis.FinitePolynomialRingWithBasisFromMorphism``
     - ``sage.combinat.multivariate_polynomials.monomial.FiniteMonomialBasis``
     - ``sage.combinat.multivariate_polynomials.ambient_space_basis.FinitePolynomialRingWithBasisFromAmbientSpace``
     - ``sage.combinat.multivariate_polynomials.linear_basis_on_vectors.FiniteLinearBasisOnVectors``
     

    INPUT:

    - ``abstract_polynomial_ring``, the abstract polynomial ring in n variables
    - ``basis_tower`` the basis of ``AbstractPolynomialRing`` corresponding 
    to ``self`` in an unset number of variables
    - ``basis_name``, the name of the basis
    - ``basis_repr`` the letter representing the basis elements
    - ``category`` (optional), the category of ``self``, if not set,
    the category class ``abstract_polynomial_ring.basis_category_class()``
    is taken

    EXAMPLES::

        sage: A = AbstractPolynomialRing(QQ)
        sage: F3 = A.finite_polynomial_ring(3)
        sage: m3 = F3.monomial_basis(); m3
        The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis
        sage: Schub = A.schubert_basis_on_vectors()
        sage: Schub3 = Schub.finite_basis(3); Schub3
        The ring of multivariate polynomials on x over Rational Field with 3 variables on the Schubert basis of type A (indexed by vectors)

    Schub3 and m3 are subclasses of FinitePolynomialRingWithBasis

    """
    def __init__(self, abstract_polynomial_ring, basis_tower, basis_keys, basis_name, basis_repr = "x", category = None, base_ring = None):
        r"""
        TESTS::
            
            sage: A = AbstractPolynomialRing(QQ)
            sage: F3 = A.finite_polynomial_ring(3)
            sage: m3 = F3.monomial_basis(); m3
            The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis

        """
        self._basis_tower = basis_tower
        self._basis_name = basis_name
        self._nb_variables = abstract_polynomial_ring.nb_variables()
        self._basis_repr = basis_repr
        if(base_ring is None): base_ring = abstract_polynomial_ring.base_ring()
        if(category is None): category = abstract_polynomial_ring.bases_category_class()(abstract_polynomial_ring, basis_tower)
        CombinatorialFreeModule.__init__(
            self,
            base_ring,
            basis_keys,
            category = category
        )
        m1 = SetMorphism( Hom(self, abstract_polynomial_ring), lambda x: x)
        m2 = SetMorphism( Hom(self, basis_tower), lambda x: x)
        m3 = SetMorphism( Hom(self, abstract_polynomial_ring.polynomial_ring_tower()), lambda x: x)
        m1.register_as_coercion()
        m2.register_as_coercion()
        m3.register_as_coercion()
        basis_tower._register_finite_basis(self)

    def _repr_(self):
        r"""
        
        Print the name of ``self``
        
        TESTS::

            sage: A = AbstractPolynomialRing(QQ)
            sage: F3 = A.finite_polynomial_ring(3)
            sage: m3 = F3.monomial_basis(); m3
            The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis

        """
        return "%s on the %s"%(self.abstract_algebra()._repr_with_basis(), self._basis_name)

    def _repr_term(self, c):
        r"""
        Returns the string representation of a term of the 
        ``CombinatorialFreeModule``.
        
        TESTS::
            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: m3 = m.finite_basis(3)  
            sage: key = list(m3.an_element())[0][0]
            sage: m3._repr_term(key)
            'x[1, 2, 3]'
        """
        return self._basis_repr+str(c)

    def basis_tower(self):
        r"""
        Returns the polynomial ring tower given to define ``self``. If the
        finite basis is the Monomial basis on 3 variables, then the ring 
        tower is the Monomial basis where the number of variables is not set.

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: m3 = m.finite_basis(3)
            sage: m3.basis_tower()
            The ring of multivariate polynomials on x over Rational Field on the monomial basis
            sage: m3.basis_tower() == m
            True

        """
        return self._basis_tower    

    def nb_variables(self):
        r"""
        Returns the number of variables in ``self``.

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: m3 = m.finite_basis(3)
                sage: m3.nb_variables()
                3

        """
        return self._nb_variables 
    
    def divided_difference_morphism(self, i, otype=None):
        r"""
            Returns a morphism for the ith divided difference 
            
            `\partial_i^A = (1 - s_i) \frac{1}{x_i - x_{i+1}}`
            
            `\partial_i^B = (1 - s_i^B) \frac{1}{x_i^{1/2} - x_i^{-1/2}}`
            
            `\partial_i^C = (1 - s_i^C) \frac{1}{x_i - x_i^{-1}}`
            
            `\partial_i^D = (1 - s_i^C) \frac{1}{x_{i-1}^{-1} - x_i}`
            

            INPUT: 

            - ``i``: the number of the divided difference
            - ``otype``: the type of divided difference for untyped basis,
            default is ``A``

            OUTPUT:

            - a morphism applying the ith divided 
              difference
              

            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: pol = m[1,2,2]
                sage: m3 = pol.parent()
                sage: morph = m3.divided_difference_morphism(1)
                sage: morph
                Generic endomorphism of The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis
                sage: morph(pol)
                -x[1, 1, 2]
           
            The operation type is given either by the ``otype`` argument or 
            by the type of polynomial.
            
            The monomial basis is an untyped basis, so a type can be specified or
            is ``A`` by default. The applied operation will be `\partial_i^X`
            where ``X`` is the sent type regardless of the number of variables.
              
            ..::
                sage: morph = m3.divided_difference_morphism(1,"B")
                sage: morph(pol)
                x[-1, 2, 2] + x[0, 2, 2]
            
            On a typed basis, the polynomial type is used and giving a 
            value to the ``otype`` will raise an exception. The operation 
            that is applied is the one that make sense for the Weyl Group 
            linked to the polynomial. As an example, in type ``B``, with 
            a polynomial in ``n`` variables, `\partial_i^A` will be used 
            for `1 \leq i < n` and `\partial_i^B` for `i=n`.
            
            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb.an_element(); pol
                x(2, 2, 3)
                sage: mb3 = pol.parent()
                sage: morph = mb3.divided_difference_morphism(2)
                sage: morph(pol)
                -x(2, 2, 2)
                sage: morph = mb3.divided_difference_morphism(3)
                sage: morph(pol)
                x(2, 2, 0) + x(2, 2, -3) + x(2, 2, -2) + x(2, 2, -1) + x(2, 2, 1) + x(2, 2, 2)
                sage: morph = mb3.divided_difference_morphism(3,"C")
                Traceback (most recent call last):
                ...
                TypeError: The type argument is not valid on this basis
                sage: morph = mb3.divided_difference_morphism(4)
                Traceback (most recent call last):
                ...
                ValueError: 4 is not a valid operation number
            """
        return self.get_morphism(i, otype, "divided_difference")
        
    def isobaric_divided_difference_morphism(self, i, otype=None):
        r"""
            Returns a morphism for the ith isobaric divided difference 
            
            `\pi_i = x_i \partial_i`
            
            `\pi_i^B = x_i^{1/2}.\partial_i^B`
            
            `\pi_i^C = x_i.\partial_i^C`
            
            `\pi_i^D = (1 - s_i^D \frac{1}{x_{i-1}x_i})\frax{1}{1- \frac{1}{x_{i-1}x_i}}`
            
            INPUT:

            - ``i``: the number of the divided difference
            - ``otype``: the type of divided difference for untyped basis,
            default is ``A``

            OUTPUT:

            - a morphism applying the ith isobaric divided 
              difference
              
            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: pol = m[2,1,1]
                sage: m3 = pol.parent()
                sage: morph = m3.isobaric_divided_difference_morphism(1)
                sage: morph
                Generic endomorphism of The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis
                sage: morph(pol)
                x[1, 2, 1] + x[2, 1, 1]
           
            The operation type is given either by the ``otype`` argument or 
            by the type of polynomial.
            
            The monomial basis is an untyped basis, so a type can be specified or
            is ``A`` by default. The applied operation will be `\pi_i^X`
            where ``X`` is the sent type regardless of the number of variables.
            
            ..::
                sage: morph = m3.isobaric_divided_difference_morphism(1,"B")
                sage: morph(pol)
                x[-2, 1, 1] + x[-1, 1, 1] + x[1, 1, 1] + x[2, 1, 1] + x[0, 1, 1]
                
                
            On a typed basis, the polynomial type is used and giving a 
            value to the ``otype`` argument will raise an exception. The operation 
            that is applied is the one that make sense for the Weyl Group 
            linked to the polynomial. As an example, in type ``B``, with 
            a polynomial in ``n`` variables, `\pi_i^A` will be used 
            for `1 \leq i < n` and `\pi_i^B` for `i=n`.
            
            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb[2,1,1]; pol
                x(2, 1, 1)
                sage: mb3 = pol.parent()
                sage: morph = mb3.isobaric_divided_difference_morphism(1)
                sage: morph(pol)
                x(1, 2, 1) + x(2, 1, 1)
                sage: pol.isobaric_divided_difference(3)
                x(2, 1, 0) + x(2, 1, -1) + x(2, 1, 1)
                sage: morph = mb3.isobaric_divided_difference_morphism(3, "C")
                Traceback (most recent call last):
                ...
                TypeError: The type argument is not valid on this basis
                sage: morph = mb3.isobaric_divided_difference_morphism(4)     
                Traceback (most recent call last):
                ...
                ValueError: 4 is not a valid operation number

            """
        return self.get_morphism(i, otype, "isobaric_divided_difference")
        
    def hat_isobaric_divided_difference_morphism(self, i, otype=None):
        r"""
            Returns a morphism for the ith hat isobaric divided difference 
            
            `\hat{\pi}_i = \partial_i x_{i+1}`
        
            `\hat{\pi}_i^B = \partial_i^B . x_i^{-1/2}.`
            
            `\hat{\pi}_i^C = \partial_i^C.x_i^{-1}`
            
            `\hat{\pi}_i^D = (1 - s_i^D)\frac{1}{x_{i-1}x_i} - 1}`
            
            INPUT:

            - ``i``: the number of the divided difference
            - ``otype``: the type of divided difference for untyped basis,
            default is ``A``

            OUTPUT:

            - a morphism applying the ith hat isobaric divided 
              difference

            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: pol = m[2,1,1]
                sage: m3 = pol.parent()
                sage: morph = m3.hat_isobaric_divided_difference_morphism(1)
                sage: morph
                Generic endomorphism of The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis
                sage: morph(pol)
                x[1, 2, 1]
           
            The operation type is given either by the ``otype`` argument or 
            by the type of polynomial.
            
            The monomial basis is an untyped basis, so a type can be specified or
            is ``A`` by default. The applied operation will be `\hat{\pi}_i^X`
            where ``X`` is the sent type regardless of the number of variables.
            
            ..::
                sage: morph = m3.hat_isobaric_divided_difference_morphism(1,"B")
                sage: morph(pol)
                x[-2, 1, 1] + x[-1, 1, 1] + x[1, 1, 1] + x[0, 1, 1]
                
            On a typed basis, the polynomial type is used. Giving a 
            value to the ``otype`` argument will raise an exception. The operation 
            that is applied is the one that make sense for the Weyl Group 
            linked to the polynomial. As an example, in type ``B``, with 
            a polynomial in ``n`` variables, `\hat{\pi}_i^A` will be used 
            for `1 \leq i < n` and `\hat{\pi}_i^B` for ``i=n``.
            
            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb[2,1,1]; pol
                x(2, 1, 1)
                sage: mb3 = pol.parent()
                sage: morph = mb3.hat_isobaric_divided_difference_morphism(1)
                sage: morph(pol)
                x(1, 2, 1)
                sage: morph = mb3.hat_isobaric_divided_difference_morphism(3)
                sage: morph(pol)
                x(2, 1, 0) + x(2, 1, -1)
                sage: morph = mb3.hat_isobaric_divided_difference_morphism(3, "C")
                Traceback (most recent call last):
                ...
                TypeError: The type argument is not valid on this basis
                sage: morph = mb3.hat_isobaric_divided_difference_morphism(4)     
                Traceback (most recent call last):
                ...
                ValueError: 4 is not a valid operation number

            """

        return self.get_morphism(i, otype, "hat_isobaric_divided_difference")
        
    def hecke_generator_morphism(self, i, t1=None, t2=None):
        r"""
            Returns a morphism for the ith Hecke algebra generator 
            
            `T_i = \pi_i (t_1 + t_2) - s_i t_2`
            
            This is only implemented in type ``A``.
            
            INPUT:

            - ``i``: the number of the operation
            - ``t_1`` the first Hecke algebra parameter
            - ``t_2`` the second Hecke algebra parameter

            OUTPUT:

            - a morphism applying the ith Hecke algebra generator, `T_i`

            EXAMPLES::

                sage: K.<t1,t2> = QQ[]
                sage: A = AbstractPolynomialRing(K)                                         
                sage: pol = A.an_element(); pol
                x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
                sage: m3 = pol.parent()
                sage: morph = m3.hecke_generator_morphism(2)
                sage: morph
                Generic endomorphism of The ring of multivariate polynomials on x over Multivariate Polynomial Ring in t1, t2 over Rational Field with 3 variables on the monomial basis
                sage: morph(pol)
                t1*x[0, 0, 0] + 2*t1*x[1, 0, 0] + (-t2)*x[1, 3, 2] + 3*t1*x[2, 0, 0]
                
            """
        if(t1 is None): t1 = self.base_ring()(var('t1'))
        if(t2 is None): t2 = self.base_ring()(var('t2'))
        return self.get_morphism(i, method="hecke_generator",t1=t1,t2=t2)
    
    def get_composed_morphism(self, operation_list):
        r"""
            Returns a composed morphism from operation list.
            
            The operations are given through a list of tuples (or list) of the 
            form ``(operation,number,otype)`` 
            ``operation`` coressponds to a keyword string:
            - ``"x"`` : product by a variable
            - ``"d"`` : divided difference
            - ``"s"`` : action of ``si``
            - ``"pi"``: isobaric divided difference
            - ``"hatpi"``: hat isobaric divided difference
            
            ``number`` is the argument of the method, ``otype`` is the method
            type (optional and only for non typed basis).
            
            A endomorphism is created by composing all the operations (acting
            on the left).
            
            INPUT:
            - ``operation_list`` a list of tuples representing the operations
            
            OUTPUT:
            A composed morphism from the operation list.
            Each operation is acting  on the left, and so the 
            operations are applied from the left to the right.
            
            EXAMPLES::
                
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: m = A.monomial_basis()
                sage: pol = m[1,2,3]
                sage: m3 = pol.parent()
                sage: morph = m3.get_composed_morphism([("d",1),("d",2)])
                sage: morph(pol)
                x[1, 1, 2] + x[1, 2, 1]
                sage: pol = m[2,1,1]
                sage: morph = m3.get_composed_morphism([("x",1),("d",1)])
                sage: morph(pol)
                x[1, 2, 1] + x[2, 1, 1]
                sage: morph = m3.get_composed_morphism([("pi",1)])   
                sage: morph(pol)    
                x[1, 2, 1] + x[2, 1, 1]
                sage: morph = m3.get_composed_morphism([("pi",2,"B"),("s",1,"C"),("hatpi",2,"D")])
                sage: morph(pol)
                -x[-2, -1, 1] - x[-2, 1, 1] - x[-2, 0, 1] - x[-1, 1, 1] - x[-1, 0, 1] - x[0, 1, 1]
                
            On a typed basis, no type parameter should be sent.
            ..::
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb.an_element()
                sage: pol                                
                x(2, 2, 3)
                sage: mb3 = pol.parent()
                sage: morph = mb3.get_composed_morphism([("d",2),("pi",1),("hatpi",3)])
                sage: morph(pol)
                -x(2, 2, 0) - x(2, 2, -2) - x(2, 2, -1) - x(2, 2, 1)
                sage: morph = mb3.get_composed_morphism([("d",2,"C"),("pi",1),("hatpi",3)])
                Traceback (most recent call last):
                ...
                TypeError: The type argument is not valid on this basis
             """
        
        composed_morph = None
        for op in operation_list:
            
            if(len(op) >2):
                morph = self.get_morphism(op[1],otype = op[2], method= op[0])
            else:
                morph = self.get_morphism(op[1], method= op[0])
            if(composed_morph is None): composed_morph = morph 
            else: composed_morph = morph * composed_morph
        
        return composed_morph
        
    def get_morphism(self, i, otype=None, method="divided_difference", **keywords):
        r"""
            Returns a morphism that applies a divided difference method on a
            polynomial.
            
            INPUT: 

            - ``i``: the number of the operation, it should be greater
            than 0 and smaller than the number of variables.
            - ``otype``: the type of divided difference if needed by the 
            parent. On the monomial basis, the default value ``A`` will be
           used.
            - ``method``, the method name to apply, here are the names that are
            accepted:
                - ``divided_difference`` or ``d`` for the usual divided difference
                - ``isobaric_divided_difference`` or ``pi`` for the isobaric 
                divided difference
                - ``hat_isobaric_divided_difference`` or ``hatpi`` for the
                hat isobaric divided difference
                - ``si``  or ``s`` for the action of the Weyl group generator
                - ``product_variable`` or ``x`` to apply a product by the ``ith``
                variable
                - ``hecke_generator`` for the hecke algebra generator
            - ``**keywords``, extra parameters needed for the method
                
            OUTPUT:
            A morphim on polynomials appling the sent method.
                        
            If the method is not implemented by the basis of ``self``, then
            a ``NotImplementedError`` exception is raised.
            
            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: pol = A.an_element();pol
                x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
                sage: m3 = pol.parent()
                sage: morph = m3.get_morphism(1); morph
                Generic endomorphism of The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis
                sage: morph(pol)
                2*x[0, 0, 0] + 3*x[1, 0, 0] - x[1, 1, 3] + 3*x[0, 1, 0]
                sage: morph = m3.get_morphism(1, method="pi"); morph
                Generic endomorphism of The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis
                sage: morph(pol)
                x[0, 0, 0] + 2*x[1, 0, 0] + 3*x[1, 1, 0] + 3*x[2, 0, 0] + 2*x[0, 1, 0] + 3*x[0, 2, 0]

            The ``otype`` argument can be used on untyped polynomials (like
            elements of the monomial basis).
            
            If used on a typed element, a ``TypeError`` Exception is raised.
            
            EXAMPLES::
                sage: morph = m3.get_morphism(1, otype="B"); morph   
                Generic endomorphism of The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis
                sage: morph(pol)
                5*x[0, 0, 0] + 3*x[-2, 0, 0] + 5*x[-1, 0, 0] + x[-1, 2, 3] + 3*x[1, 0, 0] + x[0, 2, 3]
                sage: mb = A.ambient_space_basis("B")
                sage: pol = mb(pol)
                sage: mb3 = pol.parent()
                sage: morph = mb3.get_morphism(1,otype="C")
                Traceback (most recent call last):
                ...
                TypeError: The type argument is not valid on this basis
                sage: morph = mb3.get_morphism(1)         
                sage: morph(pol)
                2*x(0, 0, 0) + 3*x(1, 0, 0) - x(1, 1, 3) + 3*x(0, 1, 0)
                sage: morph = mb3.get_morphism(3)
                sage: morph(pol)
                x(1, 2, 0) + x(1, 2, -3) + x(1, 2, -2) + x(1, 2, -1) + x(1, 2, 1) + x(1, 2, 2)
                
        The method actually looks into the class for a divided difference 
        wrapper class in which it looks for an on basis method corresponding
        to the divided difference. If it fails, then it uses the 
        ``_get_morphism_backup`` method.
        
        TESTS::
            sage: A = AbstractPolynomialRing(QQ)                                         
            sage: pol = A.an_element()    
            sage: m3 = pol.parent()
            sage: m3.get_morphism(1,method="bla")
            Traceback (most recent call last):
            ...
            NotImplementedError: the method bla is not implemented of this basis
        """
        #method dict
        methods = {}
        methods["x"] = "product_variable"
        methods["d"] = "divided_difference"
        methods["s"] = "si"
        methods["pi"] = "isobaric_divided_difference"
        methods["hatpi"] = "hat_isobaric_divided_difference"
        
        if(methods.has_key(method)):
            method = methods[method]
        if(self.basis_tower()._right_number_of_variables(i,otype) > self.nb_variables()):
            raise ValueError, "%i is not a valid operation number"%(i)
        try:
            if(otype is None):
                wrapper = self.basis_tower()._divided_difference_wrapper(self,i, **keywords)
            else:
                wrapper = self.basis_tower()._divided_difference_wrapper(self,i, otype = otype, **keywords)
            
            on_basis = getattr(wrapper, method + "_on_basis")
        except AttributeError:
            return self._get_morphism_backup(i , otype, method, **keywords)
        except TypeError:
            if(otype is not None and hasattr(self,"group_type")):    
                raise TypeError, "The type argument is not valid on this basis"%()
            return self._get_morphism_backup(i , otype, method, **keywords)
        return self._module_morphism(
            on_basis,
            codomain = self
        )

                
    def _get_morphism_backup(self, i, otype = None, method = "divided_difference", **keywords):
        r"""
        This method is used as a "backup" when the ``get_morphism`` method fails.
        By default, it only raises a ``NotImplementedError`` but it can be overwritten
        bu a subclass to actually do something.
        
        INPUT: 

            - ``i``: the number of the operation, it should be greater
            than 0 and smaller than the number of variables.
            - ``otype``: the type of divided difference if needed by the 
            parent. On the monomial basis, the default value ``A`` will be
           used.
            - ``method``, the method name to apply.
            - ``**keywords``, extra parameters needed for the method
                
        OUTPUT:
            A morphim on polynomials applying the sent method.
            
        TESTS::
            sage: A = AbstractPolynomialRing(QQ)                                         
            sage: m3 = A.an_element().parent()
            sage: m3._get_morphism_backup(1,method="bla")
            Traceback (most recent call last):
            ...
            NotImplementedError: the method bla is not implemented of this basis

        """
        raise NotImplementedError, "the method %s is not implemented of this basis"%(method)
    
            
    class Element(CombinatorialFreeModule.Element):
        r"""
        The upperclass for all polynomials elements : they are subclasses
        ``CombinatorialFreeModule`` elements
        """
        def __div__(self, element):
            """
                Division by a coefficient or by an invertible element
                
                INPUT:
                    - ``element`` the element to divide ``self`` with
                    
                OUTPUT:
                    - the result of the division
                    
                EXAMPLES::
                    sage: A = AbstractPolynomialRing(QQ);
                    sage: M = A.monomial_basis()
                    sage: pol = M[1,2,3];pol
                    x[1, 2, 3]
                    sage: 1/pol
                    x[-1, -2, -3]
                    sage: pol/2
                    1/2*x[1, 2, 3]
                    sage: pol2 = M[1,1,1];pol2
                    x[1, 1, 1]
                    sage: pol/pol2
                    x[0, 1, 2]
                    sage: (pol+pol2)/pol2
                    x[0, 0, 0] + x[0, 1, 2]
                    sage: pol/(pol+pol2)
                    Traceback (most recent call last):
                    ...
                    ValueError: can not divide x[1, 2, 3] by x[1, 1, 1] + x[1, 2, 3]
                    sage: Schub = A.schubert_basis_on_vectors("A");
                    sage: y = Schub.an_element(); y
                    Y(2, 2, 3)
                    sage: y/pol
                    Y(1, 0, 0) + Y(1, 1, -1) + Y(2, 0, -1)
                    
                
            """
            if(element in self.base_ring()):
                return super(FinitePolynomialRingWithBasis.Element, self).__div__(element)
            try:
                element = element**-1
            except ValueError:
                raise ValueError,"can not divide %s by %s"%(self, element)
            return self*element
            
        def __invert__(self):
            r"""
            A default inversion function that convert into the monomial 
            basis to invert the element. The monomial basis can invert only 
            element that are a single monomial.
            
            The result is given into the monomial basis as most of the basis
            only make sense with positive exponents.
            
            EXAMPLES::
            
                sage: A = AbstractPolynomialRing(QQ)
                sage: m = A.monomial_basis()
                sage: Schub = A.schubert_basis_on_vectors()
                sage: pol = Schub(m[1,2,3])
                sage: pol
                Y(1, 2, 3) - Y(1, 3, 2) - Y(2, 1, 3) + Y(2, 3, 1) + Y(3, 1, 2) - Y(3, 2, 1) + Y(4, 1, 1)
                sage: pol^-1
                x[-1, -2, -3]
            
            """
            m = self.parent().abstract_algebra().monomial_basis()
            return m(self)**-1
        
        def __getitem__(self,c):
            r"""
            Returns the coefficient of the input in ``self`` the input can be a list
            or an element of the basis keys.
            
            EXAMPLES::
                sage: A = AbstractPolynomialRing(QQ)                                         
                sage: m = A.monomial_basis()
                sage: pol = m[0,0,1] + m[2,1]; pol
                x[2, 1, 0] + x[0, 0, 1]
                sage: pol[2,1,0]
                1
                sage: ma = A.ambient_space_basis("A")
                sage: pol = 3*ma[2,2,1,0] + 2*ma[1,1,3];pol
                2*x(1, 1, 3, 0) + 3*x(2, 2, 1, 0)
                sage: pol[2,2,1,0]
                3
                sage: key = list(pol)[0][0]; key     
                (2, 2, 1, 0)
                sage: pol[key]
                3


            """
            if not c in self.parent()._basis_keys:
                c = self.parent()._basis_keys(c)
            
            return self._coefficient_fast(c)
       
       
        def subs(self, **args):
            r"""
            The inherited ``subs`` method is not implemented for these kinds
            of polynomials. If you want to use a subs method, you can either:
            
            - use the ``subs_on_coeffs`` if your substitution is to be applied
            on coefficient
            
            - use the ``subs_var`` method if you want to replace any `x_i`
            by a value
            
            - use the ``subs_basis`` method if you want to substitue the 
            basis and keeping the indices
            
            - use the ``subs_on_keys`` is you want to replace given keys
            by a new value 
            
            - use ``swap_vars`` or ``perm_vars`` to swap or permute vars
            
            - use ``to_expr`` to obtain an actual symbolic expression
            
            - use the ``elements`` method to obtain a list of coefficients / keys to work on
            
            """
            raise NotImplementedError, "The subs method is not implemented for elements of %s"%(self.parent())
            
    
        
        
class PolynomialRingWithBasisFromMorphism(PolynomialRingWithBasis):    
    r"""
        Implements a basis defined by its morphism to another basis. It is
        mainlu used as an upperclass of ``LinearBasisOnVectors``.

        INPUT:

        - ``abstract_polynomial_ring``, the abstract polynomial ring
        of which ``self`` is a basis.
        - ``neutral_nb_variables``: the default number of variables to get the
            one element
        - ``morphism_to_basis``: the basis of the abstract polynomial ring on
            which the morphism will be defined
        - ``basis_name``: the name of the basis (used in repr)
        - ``basis_repr``: the basis representation for elements (exemple ``x``)
        - ``get_basis_keys``: a function with :
                input:
                    - ``nb_variables`` the number of variables
                output:
                    - the set of indexes that will be used to index elements 
                      of the basis on the given 
                      number of variables
        - ``get_morphism_on_basis``: a function with :
                input:
                    -``nb_variables``, the number of variables
                output:
                    - the function that will be used to create the module morphims on 
                      basis on the given number of variables
        - ``variables_auto_coerce``: if set to ``True``, a coercion will be 
        created between elements of the basis indexed by vectors of size
        n to basis on m>n variables by extending the vectors with zeros 
        (example: x[2,2,1] -> x[2,2,1,0,0]. Default is ``False``.
        - ``**keywords`` : the keywords sent to the ``CombinatorialFreeModule``
        morphism.

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
    def __init__(self, abstract_polynomial_ring, neutral_nb_variables, morphism_to_basis, basis_name, basis_repr, get_basis_keys, get_morphism_on_basis, variables_auto_coerce = False, **keywords):
        r"""
        TESTS::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: def get_basis_keys(n): code = "A" + str(n-1); return RootSystem(code).ambient_space(QQ)
            sage: def get_morphism_on_basis(n): return lambda key: m( [key[i] for i in xrange(n)])
            sage: MyBasis = A.from_morphism_basis(1,m,get_basis_keys,get_morphism_on_basis,"My Basis", "x"); MyBasis
            The ring of multivariate polynomials on x over Rational Field on the My Basis
        """
        self._morphism_to_basis = morphism_to_basis
        self._get_morphism_on_basis = get_morphism_on_basis
        self._get_basis_keys = get_basis_keys
        self._keywords = keywords
        PolynomialRingWithBasis.__init__(
            self,
            abstract_polynomial_ring,
            basis_name,
            neutral_nb_variables,
            basis_repr,
            variables_auto_coerce = variables_auto_coerce
        )

    def _finite_basis_instance(self, n):
        r"""
        INPUT: 
        - ``n`` : the number of variables
        
        OUTPUT:
            - the Polynomial ring with basis from morphism on ``n`` variables
        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: def get_basis_keys(n): code = "A" + str(n-1); return RootSystem(code).ambient_space()
            sage: def get_morphism_on_basis(n): return lambda key: m( [key[i] for i in xrange(n)])
            sage: MyBasis = A.from_morphism_basis(1,m,get_basis_keys,get_morphism_on_basis,"My Basis", "x")
            sage: MyBasis._finite_basis_instance(3)
            The ring of multivariate polynomials on x over Rational Field with 3 variables on the My Basis
        """
        return self.abstract_algebra().finite_polynomial_ring(n).from_morphism_basis(self, self._basis_name, self.basis_repr(), **self._keywords)

    def morphism_to_basis(self):
        r"""
        The basis of the polynomial ring where the morphism of ``self`` is 
        defined to.
        
        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: def get_basis_keys(n): code = "A" + str(n-1); return RootSystem(code).ambient_space()
            sage: def get_morphism_on_basis(n): return lambda key: m( [key[i] for i in xrange(n)])
            sage: MyBasis = A.from_morphism_basis(1,m,get_basis_keys,get_morphism_on_basis,"My Basis", "x")
            sage: MyBasis.morphism_to_basis()
            The ring of multivariate polynomials on x over Rational Field on the monomial basis

        """
        return self._morphism_to_basis    

class FinitePolynomialRingWithBasisFromMorphism(FinitePolynomialRingWithBasis):    
    r"""
    This class implements the basis obtained by a morphism to another basis
    on a given number of variables it is obtained automatically by 
    PolynomialRingWithBasisFromMorphism when a polynomial is created
    see AbastractPolynomialRing.from_morphism_basis for more information.
    
    INPUT:

        - ``abstract_polynomial_ring`` the abstract ring of polynomials
        in n variables, where n is wanted number of variables for the 
        basis
        - ``polynomial_ring_tower``: the basis of ``AbsractPolynomialRing`` which
         is a facade to this basis and represents it on a undefined number 
         of variables. It must have a ``get_morphism_on_basis`` method and 
         a ``get_basis_keys`` method as well as a ``morphism_to_basis``
        - ``basis_name``: the name of the basis (used in repr)
        - ``basis_repr``: the basis representation for elements (exemple "x")
        - ``**keywords`` : the keywords sent to the ``CombinatorialFreeModule``
        morphism.

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
    """    
    def __init__(self, abstract_polynomial_ring, polynomial_ring_tower, basis_name, basis_repr,**keywords):
        r"""
        TESTS::

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
        """
        self._morphism_on_basis = polynomial_ring_tower._get_morphism_on_basis(abstract_polynomial_ring.nb_variables())
        self._morphism_to_basis = polynomial_ring_tower.morphism_to_basis().finite_basis(abstract_polynomial_ring.nb_variables())
        if(not keywords.has_key("triangular") or keywords["triangular"] is None):
            self._invertible = False
        else:
            self._invertible = True
        
        FinitePolynomialRingWithBasis.__init__(
            self,
            abstract_polynomial_ring,
            polynomial_ring_tower,
            polynomial_ring_tower._get_basis_keys(abstract_polynomial_ring.nb_variables()),
            basis_name,
            basis_repr
        )         
        self._morphism = self._module_morphism(
            self._morphism_on_basis,
            codomain = self._morphism_to_basis,
            **keywords
        )
        self._morphism.register_as_coercion()
        if(self._invertible):
            #self._morphism.__invert__().register_as_coercion()
            #For some reason, the default inver uses a preimage on basis
            # which is very slow. We changed it to use the preimage directly
            
            inv = SetMorphism(Hom(self._morphism.codomain(),self._morphism.domain()),self._morphism.preimage)
            self._inverse = inv
            inv.register_as_coercion()
        else:
            self._inverse = None
            

    def morphism_to_basis(self):
        r"""
        The finite basis of the polynomial ring where the morphism of 
        ``self`` is defined to.
        
        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: m = A.monomial_basis()
            sage: def get_basis_keys(n): code = "A" + str(n-1); return RootSystem(code).ambient_space()
            sage: def get_morphism_on_basis(n): return lambda key: m( [key[i] for i in xrange(n)])
            sage: MyBasis = A.from_morphism_basis(1,m,get_basis_keys,get_morphism_on_basis,"My Basis", "x")
            sage: MyFiniteBasis = MyBasis.finite_basis(3)
            sage: MyFiniteBasis.morphism_to_basis()
            The ring of multivariate polynomials on x over Rational Field with 3 variables on the monomial basis

        """
        return self._morphism_to_basis    

    def product(self, elt1, elt2):
        r"""
        Returns the product of the elements ``elt1`` and ``elt2``. The 
        product can be made only if the morphism of ``self`` is invertible,
        i.e. if elements of ``morphism_to_basis`` can be converted into 
        ``self``. In this case, elements are converted into ``morphism_to_basis``
        then multiply, then converted back to ``self``. Otherwise, an 
        exception is raised.

        EXAMPLES::

            sage: A = AbstractPolynomialRing(QQ)
            sage: Schub = A.schubert_basis_on_vectors()
            sage: Schub.an_element() * Schub.an_element()
            Y(4, 4, 6) + Y(4, 5, 5)
            sage: Groth = A.grothendieck_negative_basis_on_vectors()
            sage: Groth.an_element() * Groth.an_element()
            Traceback (most recent call last):
            ...
            NotImplementedError: The product is not implemented for this basis

        """
        if(self._invertible):
            return self( self.morphism_to_basis()(elt1) * self.morphism_to_basis()(elt2) )
        else:
            raise NotImplementedError, "The product is not implemented for this basis"%()
            
    def _get_morphism_backup(self, i, otype=None, method="divided_difference", **keywords):
        r"""
        This is a redefinition of ``FinitePolynomialRingWithBasis._get_morphism_backup``.
        It is used as a "backup" when the ``get_morphism`` method fails.
        If the morphism that defines the basis is invertible and if the morphism codomain
        returns a morphism for the given method, then this method returns a composed morphism
        that converts into the morphism codomain, apply the operation and then converts 
        back into the basis. 
        
        As this method is only called as a backup of the ``get_morphism``
        method, it will only be used if no direct morphism has been implemented
        on the basis for the given method.
        
        INPUT: 

            - ``i``: the number of the operation, it should be greater
            than 0 and smaller than the number of variables.
            - ``otype``: the type of divided difference if needed by the 
            parent. On the monomial basis, the default value ``A`` will be
           used.
            - ``method``, the method name to apply.
            - ``**keywords``, extra parameters needed for the method
                
        OUTPUT:
            A morphim on polynomials applying the sent method.
            
        EXAMPLES::
            sage: A = AbstractPolynomialRing(QQ)                                         
            sage: K = A.demazure_basis_on_vectors()
            sage: pol = K[3,2,1]
            sage: K3 = pol.parent()
            sage: morph = K3._get_morphism_backup(1)
            sage: morph(pol)
            K(2, 2, 1)
        """
        if(self._inverse is None):
            return super(FinitePolynomialRingWithBasisFromMorphism, self)._get_morphism_backup(i, otype, method, **keywords)
        morph = self.morphism_to_basis().get_morphism(i, otype, method, **keywords)
        return self._inverse * morph * self._morphism
        

    class Element(FinitePolynomialRingWithBasis.Element):
        def expand(self):
            """
            Returns the polynomial seen in the ``morphism_to_basis`` of the
            parent basis

            EXAMPLES::

                sage: A = AbstractPolynomialRing(QQ)
                sage: Schub = A.schubert_basis_on_vectors()        
                sage: pol = Schub.an_element(); pol
                Y(2, 2, 3)
                sage: pol.expand()
                x(2, 2, 3) + x(2, 3, 2) + x(3, 2, 2)

            """
            return self.parent()._morphism(self)
            
        


