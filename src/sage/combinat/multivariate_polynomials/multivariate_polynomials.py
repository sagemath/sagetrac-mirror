r"""
Multivariate Polynomial Algebra

This modules implements the polynomial ring seen as an algebra with multibases.
Especially, it implements bases such as the Schubert, Grothendieck, and
Key polynomials, and any basis based on a divided difference type operation.

In the monomial basis, a multivariate polynomial is seen as a linear combination
of vectors. Where each vector represents the exponents of the given monomial.

The number of variables is not set :the algebra can be understood as the projective limit
of all polynomial rings with a finite number of variables.

EXAMPLES::

    sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
    sage: A
    The Multivariate polynomial algebra on x over Rational Field
    sage: A.an_element()
    x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
    sage: x
    The Multivariate polynomial algebra on x over Rational Field on the monomial basis
    sage: x[1,1,2] + x[3,2,4]
    x[1, 1, 2] + x[3, 2, 4]
    sage: x[1,1,2] * x[3,2,4]
    x[4, 3, 6]

"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Viviane Pons < pons at univ-mlv.fr  >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.graded_algebras import GradedAlgebras
from sage.categories.category import Category
from sage.categories.graded_modules import GradedModules
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.realizations import Realizations
from sage.categories.rings import Rings
from sage.categories.category_types import Category_over_base_ring
from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.all import var

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

def MultivariatePolynomialAlgebra(base_ring, names = None):
    r"""
    Return the multivariate polynomial algebra.

    This is a polynomial ring in one or two inifinite sets of variables,
    interpreted as a multibases algebra.

    OUTPUT:

    The multivariate polynomial algebra.

    EXAMPLES::

        sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
        sage: A
        The Multivariate polynomial algebra on x over Rational Field

    The generator ``x`` is not a single variable, it represents a infinite
    set of variables. More precisely, it is the monomial basis of the algebra.

    ::

        sage: x
        The Multivariate polynomial algebra on x over Rational Field on the monomial basis
        sage: x[1,1,2] + x[2,3,4]
        x[1, 1, 2] + x[2, 3, 4]
        sage: x == A.monomial_basis()
        True

    Here is how to access a single variable::

        sage: x1 = A.var(1)
        sage: x1
        x[1]
        sage: x2 = A.var(2)
        sage: x2
        x[0, 1]
        sage: x1 * x2
        x[1, 1]

    """
    assert(base_ring in Rings())
    if names is None:
        main_repr_var = 'x'
    elif len(names) > 1:
        raise ValueError("The multivariate polynomial algebra only allows one set of variables")
    else:
        main_repr_var = names[0]
    return MultivariatePolynomialAlgebra_generic(base_ring, main_repr_var)


class MultivariatePolynomialAlgebra_generic(UniqueRepresentation, Parent):
    r"""
    A class implementing the multivariate polynomial ring as multibases
    algebra.

    INPUT:

    - ``R``: the base ring of the algebra
    - ``main_repr_var``, the letter corresponding to the set of variables,
      it is used to represent several bases, default is ``x``
    - ``always_show_main_var``, if True ``main_repr_var`` will be displayed
      on elements of every basis, even the ones that don't use it directly
      (Schubert basis, Demazure basis, ...), false by default, used on
      ``DoubleMultivariatePolynomialAlgebra`` to differentiate the two sets of
      variables


    OUTPUT:

    - The abstract ring of multivariate polynomials on ``main_repr_var`` over ``R``

    EXAMPLES::

        sage: A.<x> = MultivariatePolynomialAlgebra(QQ); A
        The Multivariate polynomial algebra on x over Rational Field

    The monomial basis is given as the algebra generator::

        sage: x
        The Multivariate polynomial algebra on x over Rational Field on the monomial basis
        sage: A.monomial_basis() == x
        True

    You can use it to create polynomials.

    ::

        sage: pol = x[1,1,2] + x[3,4]; pol
        x[1, 1, 2] + x[3, 4, 0]
        sage: pol * x[2,3]
        x[3, 4, 2] + x[5, 7, 0]

    You can also access a single variable.

    ::

        sage: A.var(1)
        x[1]
        sage: A.var(2)
        x[0, 1]
        sage: A.var(3)
        x[0, 0, 1]

    The coercion between elements with a different number of variables is done
    automatically.

    ::

        sage: pol1 = x[1,1,2]
        sage: pol1.nb_variables()
        3
        sage: pol2 = x[2]
        sage: pol2.nb_variables()
        1
        sage: (pol1 * pol2).nb_variables()
        3

    TESTS::

        sage: A = MultivariatePolynomialAlgebra(QQ)
        sage: TestSuite(A).run()
    """
    def __init__(self, R, main_repr_var, always_show_main_var = False):
        r"""
        TESTS::

            sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
            sage: A
            The Multivariate polynomial algebra on x over Rational Field
            sage: x
            The Multivariate polynomial algebra on x over Rational Field on the monomial basis
        """
        Parent.__init__(
            self,
            base = R,
            category = [GradedAlgebras(R).WithRealizations(), CommutativeAlgebras(R)]
        )
        self._finite_rings = set([])
        self._main_repr_var = main_repr_var
        self._show_main_var = always_show_main_var

    def gens(self):
        r"""
        Return a tuple whose entries are the generators for this
        object.

        In the case of the multivariate polynomial algebra, the number
        of actual generators is potentatially infinite. So this method
        actually return a tuple containing the monomial basis which can
        be seen as a multivariate gerator.

        EXAMPLES::

            sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
            sage: x
            The Multivariate polynomial algebra on x over Rational Field on the monomial basis
            sage: x == A.gens()[0]
            True
            sage: x[1,2,3]
            x[1, 2, 3]
        """
        return [self.monomial_basis()]

    def facade_for(self):
        r"""
        Return all the parents ``self`` is a facade for

        EXAMPLES::

            sage: A = MultivariatePolynomialAlgebra(CC)
            sage: A.facade_for()
            []

        The facade_for parents are added dynamically.

        """
        l = [r for r in self.realizations()]
        l.extend(self._finite_rings)
        for r in self.realizations():
            l.extend(r.facade_for())
        return l

    def _repr_(self):
        r"""
        Return the string representation of ``self``

        EXAMPLES::

            sage: MultivariatePolynomialAlgebra(QQ)
            The Multivariate polynomial algebra on x over Rational Field
        """
        return "The Multivariate polynomial algebra on %s over %s"%(
            self._main_repr_var,
            self.base_ring()
        )


    def a_realization(self):
        r"""
        Return a default realization of ``self`` : the monomial basis

        EXAMPLES:
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: A.a_realization()
            The ring of multivariate polynomials on x over Rational Field on the monomial basis

        """
        return self.monomial_basis()

    def an_element(self):
        r"""
        Return an element of ``self``. By default, this element lies in the
        monomial basis.

        EXAMPLES::

            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: A.an_element()
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
        """
        return self.monomial_basis().an_element()

    def _element_constructor_(self, element):
        r"""
        Construct an element of ``self``.

        As ``self`` is an abstract algebra, this method will
        just check if ``element`` belongs to ``self``

        INPUT:

        -``element`` the element to be contructed from

        OUTPUT:

        The element itself if it belongs to ``self``, if not,
        it raises a TypeError Exception

        EXAMPLES::

            sage: A = MultivariatePolynomialAlgebra(QQ)
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
        Return the i_th variable as a monomial base element

        INPUT:

        - ``i``: the index of the variable to return
        - ``nb_variables``: the number of variables of the result,
          default is ``i``, if ``nb_variables`` is lower than ``i`` it is
          ignored and changed to ``i``

        OUTPUT:

        - the ith variable as a monomial element

        EXAMPLES::

            sage: A = MultivariatePolynomialAlgebra(QQ);
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
        Construct a polynomial from a symbolic expression.

        INPUT:
            - ``expr`` a symbolic expression, it must be a polynomial
              in all variables of ``variables``
            - ``alphabet`` (optional), a list of symbolic variables.
              If not set, it takes ``expr.variables()``. The variables are matched
              to the vector key of the monomials by the order of the list.
            - ``second_alphabet`` (optional) a list of symbolic variables
              when working on a polynomial on 2 sets of variables

        EXAMPLES::

            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: var('x1,x2,x3')
            (x1, x2, x3)
            sage: expr = 3*x3 + x2^2 - x1*x3
            sage: A.from_expr(expr)
            -x[1, 0, 1] + x[0, 2, 0] + 3*x[0, 0, 1]
            sage: var('t1,t2')
            (t1, t2)
            sage: K.<t1,t2> = QQ[]
            sage: K = K.fraction_field()
            sage: A = MultivariatePolynomialAlgebra(K)
            sage: expr = t1*t2*x2 - x3^4*x1*4*t2^2
            sage: A.from_expr(expr,[x1,x2,x3])
            (-4*t2^2)*x[1, 0, 4] + t1*t2*x[0, 1, 0]

        Works with polynomials in two sets of variables::

            sage: D = DoubleMultivariatePolynomialAlgebra(QQ)
            sage: var('x1,x2,x3,y1,y2,y3')
            (x1, x2, x3, y1, y2, y3)
            sage: expr = x1*y1 +(y2*y3^2 - y1)*x3*x1^4
            sage: D.from_expr(expr,[x1,x2,x3],[y1,y2,y3])
            (y[1,0,0])*x[1, 0, 0] + (-y[1,0,0]+y[0,1,2])*x[4, 0, 1]

        """
        return self.monomial_basis().from_expr(expr, alphabet, second_alphabet)

    def algebra_finite_nb_variables(self, nb_variables, basis_repr = None):
        r"""
        Return the realization of ``self`` in a given number of variables.

        INPUT:

        - ``nb_variables``: the number of variables
        - ``basis_repr``, the representation letter for the elements of the base,
          by default, it is the main representation for the set of variable :
          ``self._main_repr_var``

        OUTPUT:

        The multivariate polynomial algebra in ``nb_variables`` variables

        EXAMPLES::

            sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
            sage: A3 = A.algebra_finite_nb_variables(3); A3
            The Multivariate polynomial algebra on x over Rational Field with 3 variables.

        The finite number of variables ring contains method to obtain the algebra
        bases on a finite number of variables::

            sage: m3 = F3.monomial_basis(); m3
            The Multivariate polynomial algebra on x over Rational Field with 3 variables on the monomial basis
            sage: ma3 = F3.monomial_basis_with_type("A"); ma3
            The ring of multivariate polynomials on x over Rational Field with 3 variables on the Ambient space basis of type A

        Coercions between rings with a different number of variables are created
        dynamically::

            sage: x = A.monomial_basis()
            sage: pol1 = x[1,2,3]; pol1
            x[1, 2, 3]
            sage: pol1.parent()
            The Multivariate polynomial algebra on x over Rational Field with 3 variables on the monomial basis
            sage: pol2 = x[1,1]; pol2
            x[1, 1]
            sage: pol2.parent()
            The Multivariate polynomial algebra on x over Rational Field with 2 variables on the monomial basis
            sage: pol1 + pol2
            x[1, 1, 0] + x[1, 2, 3]
            sage: (pol1 + pol2).parent()
            The Multivariate polynomial algebra on x over Rational Field with 3 variables on the monomial basis
        """
        if(basis_repr is None): basis_repr = self._main_repr_var
        F = FiniteMultivariatePolynomialAlgebra(self,nb_variables, basis_repr)
        return F

    def _register_finite_ring(self,F):
        r"""
        Add ``F`` as one of ``self`` finite variables realizations. It is called
        by the init function of ``F``

        INPUT:

        - ``F`` -- a realization of ``self`` in a given number of variables

        EXAMPLES::

        sage: var('t')
        t
        sage: A = MultivariatePolynomialAlgebra(QQ[t])
        sage: A._finite_rings
        set([])
        sage: F3 = A.algebra_finite_nb_variables(3)
        sage: A._finite_rings
        {The Multivariate polynomial algebra on x over Univariate Polynomial Ring in t over Rational Field with 3 variables}
        """
        if(not(F in self._finite_rings)):
            for ring in self._finite_rings:
                if( ring.nb_variables() != F.nb_variables() ): self._create_morphism(ring, F)
            self._finite_rings.add(F)

    def monomial_basis(self, basis_repr = None):
        r"""
        Return the monomial basis of ``self``.

        INPUT:

        - ``basis_repr``, the representation letter for the elements of the base,
          by default, it is the main representation for the set of variable :
          ``self._main_repr_var``

        OUTPUT:

        The monomial basis of the multivariate polynomial algebra.

        EXAMPLES::

            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: x = A.monomial_basis();x
            The Multivariate polynomial algebra on x over Rational Field on the monomial basis.
            sage: x[1,2,3]
            x[1, 2, 3]
            sage: x( [2,2] )
            x[2, 2]
        """
        if(basis_repr is None): basis_repr = self._main_repr_var
        from monomial import MonomialBasis
        return MonomialBasis(self, basis_repr)

    def monomial_basis_with_type(self, group_type, basis_repr = None):
        r"""
        Return a typed monomial basis.

        Monomials are indexed by a root system lattice. They embed a group
        type and all divided difference are done within this group type.

        INPUT:

        - ``group_type``: -- the letter that represents type of the weyl group,
          can be either ``"A"``, ``"B"``, ``"C"``, or ``"D"``
        - ``basis_repr`` -- (optional) the representation letter for the elements of the base,
          by default, it is using both ``self._main_repr_var`` and ``group_type``

        OUTPUT:

        The monomial basis with type ``group_type``.

        EXAMPLES::

            sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
            sage: xA = A.monomial_basis_with_type("A"); xA
            The Multivariate polynomial algebra on x over Rational Field on the Ambient space basis of type A
            sage: xB
            The Multivariate polynomial algebra on x over Rational Field on the Ambient space basis of type B
            sage: xA == x
            False
            sage: xA == xB
            False
            sage: xA.group_type()
            'A'


        Default coercions are created between the typed and untyped basis::

            sage: xA( x[1,2,3])
            xA(1, 2, 3)
            sage: x( xA[2,4] )
            x[2, 4]
            sage: xB( x[1,2,3])
            xB(1, 2, 3)
            sage: x( xB[2,4])
            x[2, 4]

        """
        if(basis_repr is None): basis_repr = self._main_repr_var + group_type
        from ambient_space_basis import PolynomialRingWithBasisFromAmbientSpace
        return PolynomialRingWithBasisFromAmbientSpace(self,group_type,"Ambient space basis of type " + group_type, basis_repr)

    def from_morphism_basis(self, neutral_nb_variables, morphism_to_basis, get_basis_keys, get_morphism_on_basis, basis_name, basis_repr, variables_auto_coerce =False, **keywords):
        r"""
        Create a basis defined by its morphism to another basis

        INPUT:

        - ``neutral_nb_variables`` -- the default number of variables to get the
            one element
        - ``morphism_to_basis`` -- the basis of the polynomial algebra on
            which the morphism will be defined
        - ``get_basis_keys`` -- a function with :
                input:
                    - ``nb_variables`` -- the number of variables
                output:
                    - the set of indexes that will be used to index elements
                      of the basis on the given
                      number of variables
        - ``get_morphism_on_basis`` -- a function with :
                input:
                    -``nb_variables``, the number of variables
                output:
                    - the function that will be used to create the module morphims on
                      basis on the given number of variables
        - ``basis_name`` -- the name of the basis (used in repr)
        - ``basis_repr``-- the basis representation for elements (exemple ``x``)
        - ``variables_auto_coerce`` -- if set to ``True``, a coercion will be
          created between elements of the basis indexed by vectors of size
          n to basis on m>n variables by extending the vectors with zeros
          (example: x[2,2,1] -> x[2,2,1,0,0]. Default is ``False``.
        - ``**keywords`` -- the keywords sent to the ``CombinatorialFreeModule``
          morphism.

        OUTPUT:

        - the basis of which elements are indexed by the sets return
          by ``get_basis_keys`` and can be coerced on
          the ``morphims_to_basis`` basis

        EXAMPLES::

            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: m = A.monomial_basis()
            sage: def get_basis_keys(n): code = "A" + str(n-1); return RootSystem(code).ambient_space(QQ)
            sage: def get_morphism_on_basis(n): return lambda key: m( [key[i] for i in xrange(n)])
            sage: MyBasis = A.from_morphism_basis(1,m,get_basis_keys,get_morphism_on_basis,"My Basis", "x"); MyBasis
            The ring of multivariate polynomials on x over Rational Field on the My Basis
            sage: MyBasis.an_element()
            x(2, 2, 3)
            sage: m( MyBasis.an_element() )
            x[2, 2, 3]

        We have recreated the typed basis.
        """
        from basis import PolynomialRingWithBasisFromMorphism
        return PolynomialRingWithBasisFromMorphism(self, neutral_nb_variables, morphism_to_basis, basis_name, basis_repr, get_basis_keys, get_morphism_on_basis,variables_auto_coerce, **keywords)

    def linear_basis_on_vectors(self, group_type, basis_name, basis_repr, on_basis_method, extra_parameters = (), **keywords):
        r"""
        Create a linear basis on objects inedexed by vectors based on an operation
        to convert each object (through its vector) into a typed polynomial.

        INPUT:

        - ``group_type`` -- the letter that represents the type of the weyl group that
           will be used for the ambient space basis
        - ``basis_name`` -- the name of the basis (used in repr)
        - ``basis_repr``-- the basis representation for elements
        - ``on_basis_method`` -- a method that takes a vector (python list) and returns
          the converted polynomial associated with it
          The ``on_basis_method`` should have the following signature :
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

            sage: A = MultivariatePolynomialAlgebra(QQ)
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
        monomial_basis_with_type = self.monomial_basis_with_type(group_type)
        return LinearBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr, on_basis_method, extra_parameters, **keywords)

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

            sage: A = MultivariatePolynomialAlgebra(QQ)
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
            sage: ma = A.monomial_basis_with_type("A")
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
        monomial_basis_with_type = self.monomial_basis_with_type("A")
        return SchubertBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr)

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

            sage: A = MultivariatePolynomialAlgebra(QQ)
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
            sage: ma = A.monomial_basis_with_type("A")
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
        monomial_basis_with_type = self.monomial_basis_with_type(group_type)
        return DemazureBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr, "isobaric_divided_difference")

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

            sage: A = MultivariatePolynomialAlgebra(QQ);
            sage: Demh = A.demazure_hat_basis_on_vectors("A"); Demh
            The ring of multivariate polynomials on x over Rational Field on the Demazure hat basis of type A (indexed by vectors)
            sage: Demh.an_element()
            ^K(2, 2, 3)

        Let us see some coercions::

            sage: kh = Demh[1,2,4]; kh
            ^K(1, 2, 4)
            sage: kh.expand()
            x(1, 2, 4) + x(1, 3, 3) + x(2, 2, 3)
            sage: ma = A.monomial_basis_with_type("A")
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
        monomial_basis_with_type = self.monomial_basis_with_type(group_type)
        return DemazureBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr, "hat_isobaric_divided_difference")


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

            sage: A = MultivariatePolynomialAlgebra(QQ)
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
            sage: ma = A.monomial_basis_with_type("A")
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
        monomial_basis_with_type = self.monomial_basis_with_type("A")
        return GrothendieckNegativeBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr)

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

            sage: A = MultivariatePolynomialAlgebra(QQ)
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
            sage: ma = A.monomial_basis_with_type("A")
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
        monomial_basis_with_type = self.monomial_basis_with_type("A")
        return GrothendieckPositiveBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr)


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
            sage: A = MultivariatePolynomialAlgebra(K)
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
            sage: ma = A.monomial_basis_with_type("A")
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
        monomial_basis_with_type = self.monomial_basis_with_type("A")
        return MacdonaldBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr, t1, t2,q)

    def _create_morphism(self,f1,f2):
        r"""
        Creates a morphism between two `FiniteMultivariatePolynomialAlgebra` on their `FiniteMonomialBasis`
        by adding extra variables to the elements of the one with least variables.
        The morphism is then registered as a coercion.

        This method is called by `finite_polynomial_ring` each time a new `FinitePolynomialRing` is created

        INPUT:
            - ``f1`` a `FiniteMultivariatePolynomialAlgebra`
            - ``f2`` another `FiniteMultivariatePolynomialAlgebra` with a different number of variables

        EXAMPLES::

            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: F1 = A.algebra_finite_nb_variables(1); F1
            The abstract ring of multivariate polynomials on x over Rational Field with 1 variable
            sage: F2 = A.algebra_finite_nb_variables(2); F2
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

            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: M = A.monomial_basis()
            sage: pol = M.an_element(); pol
            x[0, 0, 0] + 2*x[1, 0, 0] + x[1, 2, 3] + 3*x[2, 0, 0]
            sage: A.change_nb_variables(pol, 5)
            x[0, 0, 0, 0, 0] + 2*x[1, 0, 0, 0, 0] + x[1, 2, 3, 0, 0] + 3*x[2, 0, 0, 0, 0]
            sage: MA = A.monomial_basis_with_type("A")
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

            sage: A = MultivariatePolynomialAlgebra(QQ)
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

            sage: A = MultivariatePolynomialAlgebra(QQ)
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

            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: M = A.monomial_basis()
            sage: pol = M[2,1,3]
            sage: A.maxPi(pol)
            -Y(2, 2, 2)
            """
        return pol.maxPi()

class FiniteMultivariatePolynomialAlgebra(UniqueRepresentation, Parent):
    r"""
    This class implements the polynomial algebra in a given number of variables.

    INPUT:

    - ``polynomial_ring_tower`` -- the class of the polynomial algebra in an
      unset number of variables. from which the ``FiniteMultivariatePolynomialAlgebra``
      comes from. A ``FiniteMultivariatePolynomialAlgebra`` always comes from a
      ``MultivariatePolynomialAlgebra`` which contains general informations like
       the base ring.
    - ``nb_variables`` -- the number of variables
    - ``main_repr_var`` -- the letter corresponding to the set of variables,
      it is used to represent several bases, default is ``x``


    TESTS::

        sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteMultivariatePolynomialAlgebra
        sage: A = MultivariatePolynomialAlgebra(QQ[t])
        sage: B = FiniteMultivariatePolynomialAlgebra(A,4)
        sage: B
        The Multivariate polynomial algebra on x over Univariate Polynomial Ring in t over Rational Field with 4 variables
        sage: TestSuite(B).run()
    """
    def __init__(self, polynomial_ring_tower, nb_variables, main_repr_var = 'x'):
        r"""
        TESTS::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteMultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: B = FiniteMultivariatePolynomialAlgebra(A,4)
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

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteMultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: FiniteMultivariatePolynomialAlgebra(A,4)
            The Multivariate polynomial algebra on x over Univariate Polynomial Ring in t over Rational Field with 4 variables
        """
        if(self.nb_variables()>1):
            variables_str = "variables"
        else:
            variables_str = "variable"
        return "%s with %s %s"%(self.polynomial_ring_tower(),self.nb_variables(),variables_str)


    def nb_variables(self):
        r"""
        Return the number of variables of ``self``.

        EXAMPLES::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteMultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: B = FiniteMultivariatePolynomialAlgebra(A,8)
            sage: B.nb_variables()
            8
            sage: B = FiniteMultivariatePolynomialAlgebra(A,0)
            sage: B.nb_variables()
            0
        """
        return self._nb_variables

    def polynomial_ring_tower(self):
        r"""
        Return the polynomial ring tower given to define ``self``.

        EXAMPLES::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteMultivariatePolynomialAlgebra
            sage: K = CyclotomicField(3)
            sage: A = MultivariatePolynomialAlgebra(K)
            sage: B = FiniteMultivariatePolynomialAlgebra(A,3)
            sage: B.polynomial_ring_tower()
            The Multivariate polynomial algebra on x over Cyclotomic Field of order 3 and degree 2
        """
        return self._polynomial_ring_tower


    def _element_constructor_(self, element):
        r"""
        Construct an element of ``self``.

        As ``self`` is an abstract algebra, this method will
        just check if ``element`` belongs to ``self``

        INPUT:

        -``element`` the element to be contructed from

        OUTPUT:

        The element itself if it belongs to ``self``, if not,
        it raises a TypeError Exception

        TESTS::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteMultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: B = FiniteMultivariatePolynomialAlgebra(A,3)
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
        Return an element of ``self``. By default, this element lies in
        the monomial basis.

        EXAMPLES::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteMultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: B = FiniteMultivariatePolynomialAlgebra(A,7)
            sage: B.an_element()
            x[0, 0, 0, 0, 0, 0, 0] + 2*x[1, 0, 0, 0, 0, 0, 0] + x[1, 2, 3, 0, 0, 0, 0] + 3*x[2, 0, 0, 0, 0, 0, 0]

        """
        return self.a_realization().an_element()

    def a_realization(self):
        r"""
        Returns a default realization of ``self``, the monomial basis

        EXAMPLES::

            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: F = A.algebra_finite_nb_variables(3)
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

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteMultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: B = FiniteMultivariatePolynomialAlgebra(A,5)
            sage: B.monomial_basis()
            The ring of multivariate polynomials on x over Rational Field with 5 variables on the monomial basis
        """
        if(basis_repr is None): basis_repr = self._main_repr_var
        from monomial import FiniteMonomialBasis
        return FiniteMonomialBasis(self, basis_repr)

    def monomial_basis_with_type(self, letter, basis_repr = None):
        r"""
        Return the algebra ``self`` view in the proper ambient space of the
        root system design by ``letter``.

        EXAMPLES::

            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteMultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: B = FiniteMultivariatePolynomialAlgebra(A,2)
            sage: B.monomial_basis_with_type("B")
            The ring of multivariate polynomials on x over Rational Field with 2 variables on the Ambient space basis of type B
        """
        if(basis_repr is None): basis_repr = self._main_repr_var + letter
        from ambient_space_basis import FinitePolynomialRingWithBasisFromAmbientSpace
        if(letter == "A"): number = self.nb_variables()-1
        else: number = self.nb_variables()
        code = str(letter) + str(number)
        basis = FinitePolynomialRingWithBasisFromAmbientSpace(self,code,letter,"Monomial basis of type " + letter, basis_repr)
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

            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: M = A.monomial_basis()
            sage: def get_basis_keys(n): code = "A" + str(n-1); return RootSystem(code).ambient_space(QQ)
            sage: def get_morphism_on_basis(n): return lambda key: M( [key[i] for i in xrange(n)])
            sage: MyBasis = A.from_morphism_basis(1,M,get_basis_keys,get_morphism_on_basis,"My Basis", "X"); MyBasis
            The ring of multivariate polynomials on x over Rational Field on the My Basis
            sage: from sage.combinat.multivariate_polynomials.multivariate_polynomials import FiniteMultivariatePolynomialAlgebra
            sage: F2 = FiniteMultivariatePolynomialAlgebra(A,2)
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

            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: def schubert_on_basis(v, basis, call_back):
            ...     for i in xrange(len(v)-1):
            ...         if(v[i]<v[i+1]):
            ...             v[i], v[i+1] = v[i+1] + 1, v[i]
            ...             return call_back(v).divided_difference(i+1)
            ...     return basis(v)
            sage: myBasis = A.linear_basis_on_vectors("A","MySchub","Y",schubert_on_basis)
            sage: F3 = A.algebra_finite_nb_variables(3)
            sage: myFiniteBasis = F3.linear_basis_on_vectors(myBasis,"MySchub","Y")
            sage: myFiniteBasis
            The ring of multivariate polynomials on x over Rational Field with 3 variables on the MySchub

        """
        from linear_basis_on_vectors import FiniteLinearBasisOnVectors
        return FiniteLinearBasisOnVectors(self, polynomial_ring_tower, basis_name, basis_repr, **keywords)



