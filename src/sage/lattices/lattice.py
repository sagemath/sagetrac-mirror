"""
Lattices

This module provides the base class :class:`Lattice_with_basis` from which
all lattices in Sage derive, as well as a selection of more
specific base classes.

The class inheritance hierarchy is:

- :class:`FreeQuadraticModule_ambient_pid`

  - :class:`Lattice_with_basis`
    
    - :class:`Lattice_ZZ`
      
Lattices are created using the :func:`Lattice` factory function.

AUTHORS:

- Jan Poeschko (2012-05-26): initial version
"""

#*****************************************************************************
#       Copyright (C) 2012 Jan Poeschko <jan@poeschko.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.matrix.matrix_space

from sage.misc.latex import latex

from sage.modules.free_module import FreeModule_submodule_with_basis_pid, element_class
from sage.modules.free_module_element import vector
from sage.modules.free_quadratic_module import FreeQuadraticModule_ambient_pid
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.libs.pari.pari_instance import pari
from sage.matrix.matrix_space import MatrixSpace
from sage.symbolic.ring import SymbolicRing
        
from diamond_cutting import calculate_voronoi_cell

class Lattice_with_basis(FreeQuadraticModule_ambient_pid):
    """
    Construct a general lattice over a PID with an inner product matrix
    and a given basis of an embedding.
    
    INPUT:
    
    - ``base_ring`` -- the base ring of the lattice;
    
    - ``rank`` -- the rank of the lattice;
    
    - ``inner_product_matrix`` -- the inner product matrix of the lattice;

    - ``embedded_basis`` -- basis elements of the embedding of the lattice.
    
    OUTPUT:
    
    A lattice.
      
    EXAMPLES::
    
        sage: Lattice([[2, 0], [0, 1]], base_ring=Zp(5))
        Lattice of degree 2 and rank 2 over 5-adic Ring with capped relative precision 20
        Inner product matrix:
        [4 0]
        [0 1]
        Basis matrix:
        [2 0]
        [0 1]
        
    Real-valued (induced) inner product matrix::
    
        sage: L = Lattice([[1.0, 0, 0], [0, 1, 0]]); L
        Lattice of degree 3 and rank 2 over Integer Ring
        Inner product matrix:
        [ 1.00000000000000 0.000000000000000]
        [0.000000000000000  1.00000000000000]
        Basis matrix:
        [ 1.00000000000000 0.000000000000000 0.000000000000000]
        [0.000000000000000  1.00000000000000 0.000000000000000]
    """
    
    def __init__(self, base_ring, rank, inner_product_matrix, embedded_basis):
        """
        See :class:`Lattice_with_basis` for documentation.
        
        TESTS::
        
            sage: L = Lattice([[1, 0], [0, 1]], base_ring=Zp(5))
            sage: L.random_element() # random
            (-1, 2)
            sage: TestSuite(L).run()
        """
        super(Lattice_with_basis, self).__init__(base_ring, rank, inner_product_matrix)
        self._embedded_basis = matrix(embedded_basis).rows()
        if not self._embedded_basis:
            raise ValueError("basis must not be empty")
        
    def degree(self):
        """
        The degree of this lattice, i.e., the degree of the embedding vector space.
        
        EXAMPLES::
            
            sage: Lattice([[1, 0, 0], [0, 1, 0]]).degree()
            3
        """
        return len(self._embedded_basis[0])
        
    def embedded_basis(self):
        """
        The basis of the embedding of this lattice.
        
        EXAMPLES::
            
            sage: Lattice([[1, 0, 0], [0, 1, 0]]).embedded_basis()
            [(1, 0, 0), (0, 1, 0)]
        """
        return self._embedded_basis
        
    def embedded_basis_matrix(self):
        """
        The basis of the embedding of this lattice in matrix form
        
        EXAMPLES::
            
            sage: Lattice([[1, 0, 0], [0, 1, 0]]).embedded_basis_matrix()
            [1 0 0]
            [0 1 0]
        """
        return matrix(self._embedded_basis)
    
    def to_embedded(self, vector):
        """
        Convert a coordinate representation of a vector in this lattice
        to its embedded form.
        
        EXAMPLES::
            
            sage: L = Lattice([[2, 1], [0, 3]], reduce=False)
            sage: L.to_embedded((2, 3))
            (4, 11)
        """
        return sum(v * b for v, b in zip(vector, self.embedded_basis()))
    
    def determinant(self):
        """
        Calculate the determinant of this lattice, i.e.,
        the determinant of its Gram (inner product) matrix.
        
        EXAMPLES::
        
            sage: L = Lattice(inner_product_matrix=[[1, 0, 0], [0, 2, 0], [0, 0, 3]])
            sage: L.determinant()
            6
        """
        return self.inner_product_matrix().determinant()
    
    def discriminant(self):
        """
        Calculate the discriminant of this lattice, i.e.,
        the absolute value of its determinant.
        
        EXAMPLES::
        
            sage: L = Lattice(inner_product_matrix=[[1, 0, 0], [0, 2, 0], [0, 0, 3]])
            sage: L.discriminant()
            6
        """
        return abs(self.determinant())
    
    def is_unimodular(self):
        """
        Determines whether this lattice is unimodular, i.e.,
        whether its determinant is 1.
        
        EXAMPLES::
            
            sage: L = Lattice([[1, 0], [0, 1]])
            sage: L.is_unimodular()
            True
            sage: Lattice([[2, 0], [0, 3]]).is_unimodular()
            False
        """
        return self.determinant() == 1
        
    def _repr_(self):
        """
        Text representation of this lattice.
        
        TESTS::
        
            sage: Lattice([[2, 0], [0, 1]], base_ring=Zp(5))
            Lattice of degree 2 and rank 2 over 5-adic Ring with capped relative precision 20
            Inner product matrix:
            [4 0]
            [0 1]
            Basis matrix:
            [2 0]
            [0 1]        
            sage: Lattice([[1.0, 0, 0], [0, 1, 0]])
            Lattice of degree 3 and rank 2 over Integer Ring
            Inner product matrix:
            [ 1.00000000000000 0.000000000000000]
            [0.000000000000000  1.00000000000000]
            Basis matrix:
            [ 1.00000000000000 0.000000000000000 0.000000000000000]
            [0.000000000000000  1.00000000000000 0.000000000000000]
        """
        return "Lattice of degree %s and rank %s over %s\nInner product matrix:\n%s\nBasis matrix:\n%s"%(
            self.degree(), self.rank(), self.base_ring(), self.inner_product_matrix(), self.embedded_basis_matrix())
        
class Lattice_ZZ(Lattice_with_basis):
    """
    Construct a ZZ-lattice.
    
    INPUT:
    
    - ``rank``;
    - ``inner_product_matrix``;
    - ``embedded_basis``;
    - ``reduce`` -- whether to LLL-reduce the given lattice
      (both its inner product matrix and its embedded basis will be reduced if ``True``,
      which is the default).
      
    OUTPUT:
    
    A ZZ-lattice.
    
    EXAMPLES::
    
        sage: L = Lattice([[1, 0, 0], [0, 1, 0]]); L
        ZZ-lattice of degree 3 and rank 2
        Inner product matrix:
        [1 0]
        [0 1]
        Basis matrix:
        [1 0 0]
        [0 1 0]
        
    By default, the basis is reduced using the LLL algorithm::
    
        sage: L = Lattice([[6, 1], [9, 0]]); L
        ZZ-lattice of degree 2 and rank 2
        Inner product matrix:
        [ 9 -3]
        [-3 10]
        Basis matrix:
        [ 0  3]
        [ 3 -1]
        sage: L.discriminant()
        81
        
    However, you can prevent this::
    
        sage: Lattice([[6, 1], [9, 0]], reduce=False)
        ZZ-lattice of degree 2 and rank 2
        Inner product matrix:
        [37 54]
        [54 81]
        Basis matrix:
        [6 1]
        [9 0]
        sage: L.discriminant()
        81
    """
    def __init__(self, rank, inner_product_matrix, embedded_basis, reduce=True):
        """
        See :class:`Lattice_ZZ` for documentation.
        
        TESTS::
        
            sage: L = Lattice([[2, 1, 0], [0, 1, 0]])
            sage: TestSuite(L).run()
        """
        self.__reduced_inner_product_matrix = None
        self.__reduced_embedded_basis_matrix = None
        if reduce:
            inner_product_matrix, embedded_basis_matrix = self._reduce(inner_product_matrix,
                matrix(embedded_basis))
            embedded_basis = embedded_basis_matrix.rows()
        
        base_ring = ZZ
        super(Lattice_ZZ, self).__init__(base_ring, rank, inner_product_matrix, embedded_basis)
        
        self.__voronoi_cell = None      # cached result of voronoi_cell
        self.__shortest_vectors = None  # cached result of _shortest_vectors
        
    def _reduce(self, inner_product_matrix, embedded_basis_matrix):
        """
        Perform the LLL-reduction of inner product matrix and embedded basis.
        
        TESTS::
        
            sage: L = Lattice([[6, 1], [9, 0]], reduce=False)
            sage: L.reduced_embedded_basis_matrix()
            [ 0  3]
            [ 3 -1]
        """
        U = inner_product_matrix.LLL_gram()
        reduced_ipm = U.transpose() * inner_product_matrix * U
        reduced_ebm = U.transpose() * embedded_basis_matrix
        self.__reduced_inner_product_matrix = reduced_ipm
        self.__reduced_embedded_basis_matrix = reduced_ebm
        return reduced_ipm, reduced_ebm
    
    def reduced_inner_product_matrix(self):
        """
        Return an LLL-reduced inner product matrix for this lattice.
                
        EXAMPLES::
        
            sage: L = Lattice([[6, 1], [9, 0]], reduce=False); L
            ZZ-lattice of degree 2 and rank 2
            Inner product matrix:
            [37 54]
            [54 81]
            Basis matrix:
            [6 1]
            [9 0]
            sage: L.reduced_inner_product_matrix()
            [ 9 -3]
            [-3 10]
        """
        if self.__reduced_inner_product_matrix is not None:
            return self.__reduced_inner_product_matrix
        ipm, ebm = self._reduce(self.inner_product_matrix(), self.embedded_basis_matrix())
        return ipm
        
    def reduced_embedded_basis_matrix(self):
        """
        Return an LLL-reduced embedded basis for this lattice in matrix form.
        
        EXAMPLES::
        
            sage: L = Lattice([[6, 1], [9, 0]], reduce=False); L
            ZZ-lattice of degree 2 and rank 2
            Inner product matrix:
            [37 54]
            [54 81]
            Basis matrix:
            [6 1]
            [9 0]
            sage: L.reduced_embedded_basis_matrix()
            [ 0  3]
            [ 3 -1]
        """
        if self.__reduced_embedded_basis_matrix is not None:
            return self.__reduced_embedded_basis_matrix
        ipm, ebm = self._reduce(self.inner_product_matrix(), self.embedded_basis_matrix())
        return ebm
        
    def reduced_embedded_basis(self):
        """
        Return an LLL-reduced embedded basis for this lattice.
        
        EXAMPLES::
        
            sage: L = Lattice([[6, 1], [9, 0]], reduce=False)
            sage: L.reduced_embedded_basis()
            [(0, 3), (3, -1)]
        """
        return self.reduced_embedded_basis_matrix().rows()
        
    def _repr_(self):
        """
        Text representation of this lattice.
        
        TESTS::
        
            sage: Lattice([[1, 0, 0], [0, 1, 0]])
            ZZ-lattice of degree 3 and rank 2
            Inner product matrix:
            [1 0]
            [0 1]
            Basis matrix:
            [1 0 0]
            [0 1 0]
        """
        return "ZZ-lattice of degree %s and rank %s\nInner product matrix:\n%s\nBasis matrix:\n%s"%(
            self.degree(), self.rank(), self.inner_product_matrix(), self.embedded_basis_matrix())
        
    def _shortest_vectors(self, max_count=None, max_length=0):
        """
        Compute shortest vectors and their number and length (potentially).
        
        INPUT:
        
        - ``max_count`` -- maximum number of vectors to store in the result
          (default: all);
        - ``max_length`` -- maximal length of vectors to consider
          (default: 0, i.e., consider shortest vectors).
          
        OUTPUT:
        
        A triple consisting of the number of corresponding (shortest) vectors,
        their length, and a list of at most ``max_count`` such vectors. 
        
        TESTS::
        
            sage: L = Lattice([[1, 0], [0, 1]])
            sage: L.shortest_vectors_count()
            4
        """
        default = max_length == 0 and max_count == 0
        if default and self.__shortest_vectors is not None:
            return self.__shortest_vectors
        qf = self.reduced_inner_product_matrix()
        if max_count is None:
            # determine the number of shortest vectors
            max_count = self._shortest_vectors(max_length=max_length, max_count=0)[0] #self.degree()
        count, length, vectors = pari(qf).qfminim(0, max_count)
        result = count, length, vectors.python().columns()
        if default:
            self.__shortest_vectors = result
        return result
    
    def shortest_vectors_count(self, max_length=0):
        """
        Find the number of shortest vectors in the lattice.
        
        INPUT:
        
        - ``max_length`` -- the maximal length of vectors to consider.
        
        EXAMPLES::
        
            sage: L = Lattice([[1, 0], [0, 1]])
            sage: L.shortest_vectors_count()
            4
        """
        count, length, vectors = self._shortest_vectors(max_count=0, max_length=max_length)
        return count
    
    def shortest_vectors_length(self):
        """
        Find the length of the shortest vectors in the lattice.
        
        EXAMPLES::
        
            sage: L = Lattice([[1, 0], [0, 1]])
            sage: L.shortest_vectors_length()
            1
        """
        count, length, vectors = self._shortest_vectors(max_count=0)
        return length
    
    def shortest_vectors(self, max_count=None, max_length=0):
        """
        Find shortest vectors using Pari's Fincke-Pohst algorithm.
        
        INPUT:
          
        - ``max_count`` - limit on the number of vectors returned
          (default: no limit);
        
        - ``max_length`` - maximum length to consider in search
          (default of 0 means search for shortest vectors).
          
        OUTPUT:
        
        A list of at most ``max_count`` shortest vectors.
        Vectors are given in their integer representations with respect to the
        embedded basis.
        
        EXAMPLES::
        
            sage: L = Lattice([[2, 0], [0, 3]])
            sage: L.shortest_vectors()
            [(1, 0)]
            sage: map(L.to_embedded, L.shortest_vectors())
            [(2, 0)]
            
        Note that the given basis might be reduced, leading to unexpected results::
        
            sage: L = Lattice([[2, 1], [1, 1]])
            sage: L.shortest_vectors()
            [(0, 1), (1, 0)]
            sage: map(L.to_embedded, L.shortest_vectors())
            [(0, -1), (-1, 0)]
        """
        count, length, vectors = self._shortest_vectors(max_count=max_count,
            max_length=max_length)
        return vectors
    
    def voronoi_cell(self, radius=None):
        """
        Compute the Voronoi cell of a lattice, returning a Polyhedron.
        
        INPUT:
        
        - ``radius`` -- radius of ball containing considered vertices
          (default: automatic determination).
          
        OUTPUT:
        
        The Voronoi cell as a Polyhedron instance.
        
        The result is cached so that subsequent calls to this function
        return instantly.
        
        EXAMPLES::
        
            sage: L = Lattice([[1, 0], [0, 1]])
            sage: V = L.voronoi_cell()
            sage: V.Vrepresentation()
            (A vertex at (1/2, -1/2), A vertex at (1/2, 1/2), A vertex at (-1/2, 1/2), A vertex at (-1/2, -1/2))
            
        The volume of the Voronoi cell is the square root of the discriminant of the lattice::
        
            sage: L = random_lattice(4); L
                ZZ-lattice of degree 4 and rank 4
                Inner product matrix:
                [  2   1   0  -1]
                [  1   7   3   1]
                [  0   3  54   3]
                [ -1   1   3 673]
                Basis matrix:
                [  0   0   1  -1]
                [  1  -1   2   1]
                [ -6   0   3   3]
                [ -6 -24  -6  -5]
            sage: V = L.voronoi_cell()
            sage: V.volume()
            678.0
            sage: sqrt(L.discriminant())
            678
            
        Lattices not having full dimension are handled as well::
        
            sage: L = Lattice([[2, 0, 0], [0, 2, 0]])
            sage: V = L.voronoi_cell()
            sage: V.Hrepresentation()
            (An inequality (-1, 0, 0) x + 1 >= 0, An inequality (0, -1, 0) x + 1 >= 0, An inequality (1, 0, 0) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0)
        
        "Over-dimensional" lattices are reduced first::
        
            sage: L = Lattice([[1, 0], [2, 0], [0, 2]])
            sage: L.voronoi_cell().Vrepresentation()
            (A line in the direction (0, 0, 1), A vertex at (1, -1, 0), A vertex at (1, 1, 0), A vertex at (-1, 1, 0), A vertex at (-1, -1, 0))
            
        ALGORITHM:
        
        Uses parts of the algorithm from [Vit1996].
        
        REFERENCES:
        
        .. [Vit1996] E. Viterbo, E. Biglieri. Computing the Voronoi Cell
          of a Lattice: The Diamond-Cutting Algorithm.
          IEEE Transactions on Information Theory, 1996.
        """
        if self.__voronoi_cell is None:
            basis_matrix = self.reduced_embedded_basis_matrix()
            self.__voronoi_cell = calculate_voronoi_cell(basis_matrix, radius=radius)
        return self.__voronoi_cell
    
    def voronoi_relevant_vectors(self):
        """
        Compute the embedded vectors inducing the Voronoi cell.
        
        OUTPUT:
        
        The list of Voronoi relevant vectors.
        
        EXAMPLES::
        
            sage: L = Lattice([[3, 0], [4, 0]])
            sage: L.voronoi_relevant_vectors()
            [(-1, 0), (1, 0)]
        """
        V = self.voronoi_cell()
        
        def defining_point(ieq):
            """
            Compute the point defining an inequality.
            
            INPUT:
            
            - ``ieq`` - an inequality in the form [c, a1, a2, ...]
              meaning a1 * x1 + a2 * x2 + ... <= c
              
            OUTPUT:
            
            The point orthogonal to the hyperplane defined by ``ieq``
            in twice the distance from the origin.
            """
            c = ieq[0]
            a = ieq[1:]
            n = sum(y ** 2 for y in a)
            return vector([2 * y * c / n for y in a])
            
        return [defining_point(ieq) for ieq in V.inequality_generator()]
    
    def closest_vector(self, t):
        """
        Compute the closest vector in the embedded lattice to a given vector.
        
        INPUT:
        
        - ``t`` -- the target vector to compute the closest vector to.
        
        OUTPUT:
        
        The vector in the lattice closest to ``t``.
        
        EXAMPLES::
        
            sage: L = Lattice([[1, 0], [0, 1]])
            sage: L.closest_vector((-6, 5/3))
            (-6, 2)
            
        ALGORITHM:
        
        Uses the algorithm from [Mic2010].
        
        REFERENCES:
        
        .. [Mic2010] D. Micciancio, P. Voulgaris. A Deterministic Single
          Exponential Time Algorithm for Most Lattice Problems based on
          Voronoi Cell Computations.
          Proceedings of the 42nd ACM Symposium Theory of Computation, 2010.
        """
        voronoi_cell = self.voronoi_cell()
        
        def projection(M, v):
            Mt = M.transpose()
            P = Mt * (M * Mt) ** (-1) * M
            return P * v
        
        t = projection(matrix(self.basis()), vector(t))
        
        def CVPP_2V(t, V, voronoi_cell):
            t_new = t
            while not voronoi_cell.contains(t_new.list()):
                v = max(V, key=lambda v: t_new * v / v.norm() ** 2)
                t_new = t_new - v
            return t - t_new
            
        V = self.voronoi_relevant_vectors()
        t = vector(t)
        p = 0
        while not (ZZ(2 ** p) * voronoi_cell).contains(t):
            p += 1
        t_new = t
        i = p
        while i >= 1:
            V_scaled = [v * (2 ** (i - 1)) for v in V]
            t_new = t_new - CVPP_2V(t_new, V_scaled, ZZ(2 ** (i - 1)) * voronoi_cell)
            i -= 1
        return t - t_new
