r"""
Shortcuts for most used polynomial bases
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Viviane Pons < pons at univ-mlv.fr  >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from multivariate_polynomials import AbstractPolynomialRing

def SchubertPolynomialsOnVectors(R, basis_name = None, basis_repr= "Y", **keywords):
    r"""
    The Schubert polynomial ring indexed by vectors as a basis of 
    multivariate polynomials on different bases, see :class:`AbstractPolynomialRing`
    
    For double Schubert polynomials, see :class:`DoubleAbstractPolynomialRing`
    
    Here is the definition we use. For $v = (v_1, \cdots, v_n) \in 
    \mathbb{N}^n$, we define
        
    $Y_v = x_1^{v_1}x_2^{v_2}\cdotsx_n^{v_n}$ if $v$ is a partition, 
    i.e, if $v_1 \geq v_2 \geq \cdots \geq v_n$.
    
    Otherwise, we have for $ v_i > v_{i+1}$
    
    $Y_{\cdots v_{i+1} v_i-1 \cdots} = Y_v \partial_i$ where $\partial_i$ 
    is the ith divided difference. 
    
    The vectors indexing the Schubert polynomials can as well been seen as 
    lehmer codes.
    
    INPUT:
    
    - ``R`` -- a ring
    - ``basis_name`` -- (default: canonical name) the name of the basis
     (used in repr)
    - ``basis_repr`` -- (default: ``Y``) the basis representation for elements 
    - ``**keywords`` --  other keywords to send t the abstract polynomial ring

    OUTPUT:

    - The ring of multivariate polynomials on x over ``R`` on the Schubert basis
    
    EXAMPLES::
    
        sage: Schub = SchubertPolynomialsOnVectors(QQ)
        sage: Schub
        The ring of multivariate polynomials on x over Rational Field on the Schubert basis of type A (indexed by vectors)
        sage: Schub.an_element()
        Y(2, 2, 3)
        sage: Schub[1,2,3] + Schub[2,2]
        Y(1, 2, 3) + Y(2, 2, 0)
        
     some operations::
     
        sage: (Schub[2,2,3])^2
        Y(4, 4, 6) + Y(4, 5, 5)
        sage: Schub[2,3] * Schub[1,2]
        Y(3, 5) + Y(4, 4)
        sage: Schub[2,3].expand()
        x(2, 3) + x(3, 2)
        sage: Schub[2,3].to_expr()
        x1^3*x2^2 + x1^2*x2^3
        sage: pol = Schub[3,2,3] + Schub[3,1,1]
        sage: pol.divided_difference(1)
        Y(1, 2, 1) + Y(2, 2, 3)
        sage: pol.divided_difference(2)
        0
        sage: pol.isobaric_divided_difference(1) 
        Y(1, 3, 1) + Y(2, 3, 3) + Y(2, 4, 2)
        sage: pol.isobaric_divided_difference(2)
        Y(3, 1, 1) + Y(3, 2, 3)

        
    some coercions::
    
        sage: A = Schub.abstract_algebra();A
        The abstract ring of multivariate polynomials on x over Rational Field    
        sage: m = A.monomial_basis()
        sage: pol = m[2,2,3] + m[3,1,2];pol
        x[2, 2, 3] + x[3, 1, 2]
        sage: Schub(pol)
        Y(2, 2, 3) - Y(2, 3, 2) + Y(3, 1, 2) - Y(3, 2, 1)
        sage: var('x1,x2,x3')
        (x1, x2, x3)
        sage: expr = x1^2 * x2^2 * x3^3 + x1^3 * x2 *x3^2
        sage: Schub.from_expr(expr)
        Y(2, 2, 3) - Y(2, 3, 2) + Y(3, 1, 2) - Y(3, 2, 1)
        sage: Dem = DemazurePolynomials(QQ);
        sage: Schub(Dem[1,0,2])
        Y(1, 0, 2) - Y(3, 0, 0)

    """
    A = AbstractPolynomialRing(R, **keywords)
    return A.schubert_basis_on_vectors(basis_name=basis_name, basis_repr=basis_repr)


def DemazurePolynomials(R, group_type ="A", basis_name = None, basis_repr = "K", **keywords):
    r"""
        Creates the Demazure polynomials where demazure / key polynomials are indexed
        by vectors, as basis of multivariate polynomials on different bases,
         see :class:`AbstractPolynomialRing`
        
        Here is the definition we use for type A. For $v = (v_1, \cdots, v_n) \in \mathbb{N}^n$, 
        we define
        
        $K_v = x_1^{v_1}x_2^{v_2}\cdotsx_n^{v_n}$ if $v$ is a partition, i.e, if $v_1 \geq v_2 \geq \cdots \geq v_n$.
        
        Otherwise, we have for $ v_i > v_{i+1}$
        
        $K_{\cdots v_{i+1} v_i \cdots} = K_v \pi_i$ where $\pi_i$ is the 
        ith isobar divided difference. 
        
        The vectors indexing the key polynomials can as well been seen 
        as lehmer codes.

        INPUT:

        - ``R`` -- a ring
        - ``group_type`` -- (default: ``A``) the letter that represents the type of the weyl group
        - ``basis_name`` -- (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``-- (default: ``K``) the basis representation for elements 
        - ``**keywords`` --  other keywords to send t the abstract polynomial ring

        OUTPUT:

        - The ring of multivariate polynomials on x over ``R`` on the Demazure basis
          of type ``group_type`` 
        
        EXAMPLES::
        
            sage: Dem = DemazurePolynomials(QQ);Dem
            The ring of multivariate polynomials on x over Rational Field on the Demazure basis of type A (indexed by vectors)
            sage: Dem.an_element()
            K(2, 2, 3)
            sage: Dem[1,2,2] + Dem[2,3,1]
            K(1, 2, 2) + K(2, 3, 1)
        
        some operations::
        
            sage: Dem[2,2,3]^2
            K(4, 4, 6) + K(4, 5, 5)
            sage: Dem[2,3] * Dem[1,2]
            K(3, 5) + K(4, 4)
            sage: Dem[2,3].expand()
            x(2, 3) + x(3, 2)
            sage: Dem[2,3].to_expr()
            x1^3*x2^2 + x1^2*x2^3
            sage: pol = Dem[3,2,3] + Dem[3,1,1]
            sage: pol.divided_difference(1)
            K(1, 2, 1) + K(2, 2, 3) - K(2, 3, 2)
            sage: pol.divided_difference(2)
            0
            sage: pol.isobaric_divided_difference(1)
            K(1, 3, 1) + K(2, 3, 3)
            sage: pol.isobaric_divided_difference(2)
            K(3, 1, 1) + K(3, 2, 3)
            
        some coercions::
        
            sage: A = Dem.abstract_algebra();A  
            The abstract ring of multivariate polynomials on x over Rational Field
            sage: m = A.monomial_basis()
            sage: pol = m[2,2,3] + m[3,1,2];pol
            x[2, 2, 3] + x[3, 1, 2]
            sage: Dem(pol)
            K(2, 2, 3) - K(2, 3, 2) + K(3, 1, 2) - K(3, 2, 1)
            sage: var('x1,x2,x3')
            (x1, x2, x3)
            sage: expr = x1^2 * x2^2 * x3^3 + x1^3 * x2 *x3^2     
            sage: Dem.from_expr(expr)
            K(2, 2, 3) - K(2, 3, 2) + K(3, 1, 2) - K(3, 2, 1)
            sage: Schub = SchubertPolynomialsOnVectors(QQ)
            sage: Dem(Schub[1,0,2])
            K(1, 0, 2) + K(3, 0, 0)

    """
    A = AbstractPolynomialRing(R, **keywords)
    return A.demazure_basis_on_vectors(group_type=group_type, basis_name = basis_name, basis_repr=basis_repr)
    
def DemazureHatPolynomials(R, group_type ="A", basis_name = None, basis_repr = "^K", **keywords):
    r"""
        Creates the Demazure hat polynomials where demazure / key polynomials are indexed
        by vectors, as basis of multivariate polynomials on different bases,
         see :class:`AbstractPolynomialRing`
        
        Here is the definition we use for type A. For $v = (v_1, \cdots, v_n) 
        \in \mathbb{N}^n$, we define
        
        $\hat{K}_v = x_1^{v_1}x_2^{v_2}\cdotsx_n^{v_n}$ if $v$ is a partition, i.e,
         if $v_1 \geq v_2 \geq \cdots \geq v_n$.
        
        Otherwise, we have for $ v_i > v_{i+1}$
        
        $\hat{K}_{\cdots v_{i+1} v_i \cdots} = \hat{K}_v \hat{\pi}_i$ where
         $\hat{\pi}_i$ is the ith isobar hat divided difference. 
        
        The vectors indexing the key polynomials can as well been seen 
        as lehmer codes.

        INPUT:

        - ``R`` -- a ring
        - ``group_type`` -- (default: ``A``) the letter that represents the type of the weyl group
        - ``basis_name`` -- (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``-- (default: ``^K``) the basis representation for elements 
        - ``**keywords`` --  other keywords to send t the abstract polynomial ring

        OUTPUT:

        - The ring of multivariate polynomials on x over ``R`` on the Demazure basis
          of type ``group_type`` 
        
        EXAMPLES::
        
            sage: HatDem = DemazureHatPolynomials(QQ); HatDem
            The ring of multivariate polynomials on x over Rational Field on the Demazure hat basis of type A (indexed by vectors)
            sage: HatDem.an_element()
            ^K(2, 2, 3)
            sage: HatDem[1,2,2] + HatDem[2,3,1]
            ^K(1, 2, 2) + ^K(2, 3, 1)
        
        some operations::

            sage: HatDem[2,2,3]^2
            ^K(4, 4, 6) - ^K(4, 5, 5) - ^K(5, 4, 5)
            sage: HatDem[2,3] * HatDem[1,2]
            ^K(3, 5) - ^K(4, 4)
            sage: HatDem[2,3].expand()
            x(2, 3)
            sage: HatDem[2,3].to_expr()
            x1^2*x2^3
            sage: pol = HatDem[3,2,3] + HatDem[3,1,1]
            sage: pol.divided_difference(1)
            ^K(1, 2, 1) + ^K(2, 1, 1) + ^K(2, 2, 3)
            sage: pol.divided_difference(2)
            -^K(3, 2, 2)
            sage: pol.isobaric_divided_difference(1)
            ^K(1, 3, 1) + ^K(2, 3, 3) + ^K(3, 1, 1) + ^K(3, 2, 3)
            sage: pol.isobaric_divided_difference(2)
            ^K(3, 1, 1)
            sage: pol.hat_isobaric_divided_difference(1)
            ^K(1, 3, 1) + ^K(2, 3, 3)
            sage: pol.hat_isobaric_divided_difference(2)
            -^K(3, 2, 3)

            
        some coercions::
        
            sage: A = HatDem.abstract_algebra(); A
            The abstract ring of multivariate polynomials on x over Rational Field
            sage: m = A.monomial_basis()
            sage: pol = m[2,2,3] + m[3,1,2];pol
            x[2, 2, 3] + x[3, 1, 2]
            sage: HatDem(pol)
            ^K(2, 2, 3) + ^K(3, 1, 2)
            sage: var('x1,x2,x3')
            (x1, x2, x3)
            sage: expr = x1^2 * x2^2 * x3^3 + x1^3 * x2 *x3^2       
            sage: HatDem.from_expr(expr)
            ^K(2, 2, 3) + ^K(3, 1, 2)
            sage: Schub = SchubertPolynomialsOnVectors(QQ)
            sage: HatDem(Schub[1,0,2])
            ^K(1, 2, 0) + ^K(1, 0, 2) + ^K(2, 1, 0) + ^K(2, 0, 1) + ^K(3, 0, 0)

    """
    A = AbstractPolynomialRing(R, **keywords)
    return A.demazure_hat_basis_on_vectors(group_type=group_type, basis_name = basis_name, basis_repr=basis_repr)

def GrothendieckPolynomials(R, basis_name = None, basis_repr= "G", **keywords):
    r"""
    The Grothendieck polynomial ring indexed by vectors as a basis of 
    multivariate polynomials on different bases, see :class:`AbstractPolynomialRing`
    
    For double Grothendieck polynomials, see :class:`DoubleAbstractPolynomialRing`.
    
   Here is the definition we use. For $v = (v_1, \cdots, v_n) \in \mathbb{N}^n$, we define
        
    $G_v = x_1^{v_1}x_2^{v_2}\cdotsx_n^{v_n}$ if $v$ is a partition, i.e, if $v_1 \geq v_2 \geq \cdots \geq v_n$.
    
    Otherwise, we have for $ v_i > v_{i+1}$
    
    $G_{\cdots v_{i+1} v_i-1 \cdots} = G_v \pi_i$ where $\pi_i$ is the ith isobar divided difference. 
    
    The vectors indexing the Grothendieck polynomials can as well been seen as lehmer codes.
       
    
    The vectors indexing the Grothendieck polynomials can as well been seen as 
    lehmer codes.
    
    INPUT:
    
    - ``R`` -- a ring
    - ``basis_name`` -- (default: canonical name) the name of the basis
     (used in repr)
    - ``basis_repr`` -- (default: ``G``) the basis representation for elements 
    - ``**keywords`` --  other keywords to send to the abstract polynomial ring

    OUTPUT:

    - The ring of multivariate polynomials on x over ``R`` on the Grothendieck basis
    
    EXAMPLES::
    
        sage: Groth = GrothendieckPolynomials(QQ)
        sage: Groth
        The ring of multivariate polynomials on x over Rational Field on the Grothendieck basis of type A, with positive exposants (indexed by vectors) 
        sage: Groth.an_element()
        G(2, 2, 3)
        sage: Groth[1,2,3] + Groth[2,2]
        G(1, 2, 3) + G(2, 2, 0)

     some operations::
     
        sage: (Groth[2,2,3])^2
        G(4, 4, 6) + G(4, 5, 5) - G(4, 5, 6)
        sage: Groth[2,3] * Groth[1,2]
        G(3, 5) + G(4, 4) - G(4, 5)
        sage: Groth[2,3].expand()
        x(2, 3) + x(3, 2) - x(3, 3)
        sage: Groth[2,3].to_expr()
        -x1^3*x2^3 + x1^3*x2^2 + x1^2*x2^3
        sage: pol = Groth[3,2,3] + Groth[3,1,1]    
        sage: pol.divided_difference(1)
        G(1, 2, 1) + G(2, 2, 1) + G(2, 2, 3)
        sage: pol.divided_difference(2)
        0
        sage: pol.isobaric_divided_difference(1)
        G(1, 3, 1) + G(2, 3, 1) + G(2, 3, 3) + G(2, 4, 2) - G(2, 4, 3) + G(3, 3, 1) + G(3, 3, 3) + G(3, 4, 2) - G(3, 4, 3)
        sage: pol.isobaric_divided_difference(2)
        G(3, 1, 1) + G(3, 2, 3)
        
    some coercions::
    
        sage: A = Groth.abstract_algebra();A
        The abstract ring of multivariate polynomials on x over Rational Field
        sage: m = A.monomial_basis()
        sage: pol = m[2,2,3] + m[3,1,2];pol    
        x[2, 2, 3] + x[3, 1, 2]
        sage: Groth(pol)
        G(2, 2, 3) - G(2, 3, 2) + G(2, 3, 3) + G(3, 1, 2) - G(3, 2, 1) + G(3, 2, 2) - G(3, 3, 2) + G(3, 3, 3)
        sage: var('x1,x2,x3')
        (x1, x2, x3)
        sage: expr = x1^2 * x2^2 * x3^3 + x1^3 * x2 *x3^2   
        sage: Groth.from_expr(expr)
        G(2, 2, 3) - G(2, 3, 2) + G(2, 3, 3) + G(3, 1, 2) - G(3, 2, 1) + G(3, 2, 2) - G(3, 3, 2) + G(3, 3, 3)
        sage: Schub = SchubertPolynomialsOnVectors(QQ)
        sage: Groth(Schub[1,0,2])
        G(1, 1, 2) + G(1, 2, 2) + G(1, 0, 2) + G(2, 2, 2) + G(2, 0, 2) + G(3, 0, 2)

    """
    A = AbstractPolynomialRing(R, **keywords)
    return A.grothendieck_positive_basis_on_vectors(basis_name=basis_name, basis_repr=basis_repr)
    
