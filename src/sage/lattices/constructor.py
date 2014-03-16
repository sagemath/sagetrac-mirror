"""
Lattice constructor
      
Lattices are created using the `Lattice` factory function, given
either a basis, an inner product matrix, a quadratic form,
or an order in an algebraic number field.

AUTHORS:

- Jan Poeschko (2012-08-10): initial version
"""

#*****************************************************************************
#       Copyright (C) 2012 Jan Poeschko <jan@poeschko.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix, random_matrix
from lattice import Lattice_with_basis, Lattice_ZZ
        
def Lattice(basis=None, inner_product_matrix=None, quadratic_form=None, algebraic_order=None,
        base_ring=ZZ, **kwargs):
    """
    The :func:`Lattice` function creates lattices using a given base,
    an inner product matrix, an underlying quadratic form,
    or an order in an algebraic number field.
    
    INPUT:
    
    - ``basis``
    - ``inner_product_matrix``
    - ``quadratic_form``
    - ``algebraic_order``
    - ``base_ring``
    
    OUTPUT:
    
    A lattice.
    
    EXAMPLES::
    
        sage: Lattice([[2, 0, 0], [0, 1, 0]])
        ZZ-lattice of degree 3 and rank 2
        Inner product matrix:
        [1 0]
        [0 4]
        Basis matrix:
        [ 0  1  0]
        [-2  0  0]
        
    A lattice can be specified by a quadratic form::
        
        sage: Lattice(quadratic_form=QuadraticForm(ZZ, 3, [1,2,3,4,5,6]))
        Lattice of degree 3 and rank 3 over Integer Ring
        Inner product matrix:
        [  1   1 3/2]
        [  1   4 5/2]
        [3/2 5/2   6]
        Basis matrix:
        [                  1                   0                   0]
        [                  1  1.732050807568878?                   0]
        [                3/2 0.5773502691896258?  1.848422751068237?]
        
    It is also possible to specify an order in an algebraic number field::
    
        sage: K.<a> = NumberField(x ^ 2 + 1)
        sage: O = K.order(a)
        sage: Lattice(algebraic_order=O)
        ZZ-lattice of degree 2 and rank 2
        Inner product matrix:
        [1 0]
        [0 1]
        Basis matrix:
        [1 0]
        [0 1]
        
    It is an error to specify an empty basis::
    
        sage: Lattice([])
        Traceback (most recent call last):
        ...
        ValueError: basis must not be empty
    """
    if quadratic_form is not None:
        inner_product_matrix = quadratic_form.Gram_matrix_rational()
    if algebraic_order is not None:
        basis = algebraic_order.free_module().basis()
    if basis is not None:
        if not basis:
            raise ValueError("basis must not be empty")
        basis = matrix(basis)
        rank = basis.dimensions()[0]
        if inner_product_matrix is None:
            inner_product_matrix = basis * basis.transpose()
    elif inner_product_matrix is not None:
        inner_product_matrix = matrix(inner_product_matrix)
        rank = inner_product_matrix.dimensions()[0]
        basis = inner_product_matrix.cholesky()
    else:
        raise TypeError("basis or inner_product_matrix must be given")
    if base_ring == ZZ and inner_product_matrix.base_ring() == ZZ:
        return Lattice_ZZ(rank, inner_product_matrix, basis, **kwargs)
    else:
        return Lattice_with_basis(base_ring, rank, inner_product_matrix, basis)
    
def random_lattice(dimension):
    """
    Construct a random ZZ-lattice of a given dimension.
    
    INPUT:
    
    - ``dimension`` -- dimension of the constructed lattice.
    
    OUTPUT:
    
    A lattice with integer basis vectors.
    
    EXAMPLES::
    
        sage: random_lattice(3)
        ZZ-lattice of degree 3 and rank 3
        Inner product matrix:
        [   2    0    0]
        [   0   66  -22]
        [   0  -22 4246]
        Basis matrix:
        [ 0  1 -1]
        [-8  1  1]
        [14 45 45]
        sage: random_lattice(100).dimension()
        100
    """
    basis = random_matrix(ZZ, dimension, dimension)
    return Lattice(basis)
