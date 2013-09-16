"""
Examples for CHomP

These are mostly used for doctesting.
"""



def matrix_complex_circle():
    """
    Construct the chain complex for a circle

    EXAMPLES::

        sage: from sage.libs.chomp.examples import matrix_complex_circle
        sage: matrix_complex_circle()
        bd(index, dim) = boundary chain
        bd(0, 0) = 0
        bd(1, 0) = 0
        bd(0, 1) = 1[0] + 2[1]
        bd(1, 1) = 1[0] + 2[1]
    """
    from sage.libs.chomp.matrix_complex import MatrixComplex
    circle = [[0], [0, 0, 1], [1, 0, -1], [0, 1, 1], [1, 1, -1]]
    return MatrixComplex(circle)
    
