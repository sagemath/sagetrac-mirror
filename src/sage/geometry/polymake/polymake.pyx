# distutils: language = c++
# distutils: libraries = polymake gmp xml2 perl
### distutils: include_dirs = /opt/food/include

###############################################################################
#       Copyright (C) 2011-2012 Burcin Erocal <burcin@erocal.org>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

from libc.stdlib cimport malloc

from sage.matrix.matrix_rational_dense cimport Matrix_rational_dense
from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.structure.sage_object cimport SageObject
from sage.rings.integer cimport Integer as SageInteger
from sage.rings.rational cimport Rational as SageRational
from sage.libs.gmp.mpz cimport mpz_set
from sage.libs.gmp.mpq cimport mpq_set
from sage.rings.all import QQ, ZZ
from sage.matrix.matrix_space import MatrixSpace

from defs cimport Main, PerlObject, MatrixRational, Rational, Integer, \
        VectorInteger
from defs cimport CallPolymakeFunction, CallPolymakeFunction1, \
        CallPolymakeFunction2, CallPolymakeFunction3, \
        new_PerlObject_from_PerlObject
from defs cimport pm_get, pm_get_MatrixRational, pm_get_PerlObject, \
        pm_get_VectorInteger, \
        pm_assign, get_element

# FIXME: pass user-settings parameter
cdef Main pm

cdef Vector_integer_dense pm_VectorInteger_to_sage(VectorInteger pm_vec):
    cdef Py_ssize_t size = pm_vec.size()
    # this seems to be a rather unreliable way to getting a new element
    # even though it is documented to work in the comments of zero_vector()
    cdef Vector_integer_dense svec = (ZZ**size).zero_vector()
    cdef Py_ssize_t i
    for i in range(size):
        mpz_set(svec._entries[i], pm_vec.get(i).get_rep())
    return svec

cdef Matrix_rational_dense pm_mat_to_sage(MatrixRational pm_mat):
    """
    Convert a polymake Matrix<Rational> to a Sage Matrix_rational_dense.
    """
    cdef Py_ssize_t nr = pm_mat.rows(), nc = pm_mat.cols()
    cdef Matrix_rational_dense smat = MatrixSpace(QQ, nr, nc)() #zero matrix
    cdef Py_ssize_t i, j
    for i in range(nr):
        for j in range(nc):
            mpq_set(smat._matrix[i][j],
                    get_element(pm_mat, i, j).get_rep())
    return smat

cdef MatrixRational* sage_mat_to_pm(Matrix_rational_dense mat):
    """
    Convert a Sage Matrix_rational_dense to a polymake Matrix<Rational>.
    """
    cdef Py_ssize_t nr = mat.nrows()
    cdef Py_ssize_t nc = mat.ncols()
    # create polymake matrix with dimensions of mat
    cdef MatrixRational* pm_mat = new MatrixRational(nr, nc)
    cdef Py_ssize_t i, j
    cdef Rational *tmp_rat
    # loop through the elements and assign values
    for i in range(nr):
        for j in range(nc):
            get_element(pm_mat[0], i, j).set(mat._matrix[i][j])
    return pm_mat

cdef class Polytope(SageObject):
    cdef PerlObject* pm_obj
    def __init__(self, prop_name, data):
        if prop_name not in ['VERTICES', 'POINTS', 'FACETS']:
            raise ValueError("property must be VERTICES, POINTS or FACETS")
        pm.set_application("polytope")
        self.pm_obj = new PerlObject("Polytope<Rational>")
        cdef MatrixRational* pm_mat = sage_mat_to_pm(data)
        pm_assign(self.pm_obj.take(prop_name), pm_mat[0])
        del pm_mat

    def __dealloc__(self):
        del self.pm_obj

    def _repr_(self):
        pass

    def _save(self, filename):
        """
        Saves this polytope to a file using polymake's representation.
        """
        self.pm_obj.save(filename)

    def _get_bool_property(self, prop):
        cdef Integer pm_res
        pm_get(self.pm_obj.give(prop),pm_res)
        return bool(pm_res.compare(0))

    def _get_integer_property(self, prop):
        cdef Integer pm_res
        pm_get(self.pm_obj.give(prop), pm_res)
        cdef SageInteger res = SageInteger.__new__(SageInteger)
        mpz_set(res.value, pm_res.get_rep())
        return res

    def _get_matrix_property(self, prop):
        cdef MatrixRational pm_mat
        pm_get_MatrixRational(self.pm_obj.give(prop), pm_mat)
        return pm_mat_to_sage(pm_mat)

    def _get_vector_property(self, prop):
        cdef VectorInteger pm_vec
        pm_get_VectorInteger(self.pm_obj.give(prop), pm_vec)
        return pm_VectorInteger_to_sage(pm_vec)

    def __add__(left, right):
        if not (isinstance(left, Polytope) and isinstance(right, Polytope)):
            raise TypeError("both arguments must be instances of Polytope")

    def f_vector(self):
        return self._get_vector_property("F_VECTOR")

    def h_star_vector(self):
        return self._get_vector_property("H_STAR_VECTOR")

    def num_facets(self):
        """
        EXAMPLES::

            sage: import sage.geometry.polymake as polymake  # optional - polymake
            sage: m = matrix(QQ,[[1,0,0,0], [1,0,0,1], [1,0,1,0], [1,0,1,1],  [1,1,0,0], [1,1,0,1], [1,1,1,0], [1,1,1,1]])
            sage: p = polymake.Polytope('POINTS',m)
            sage: p.num_facets()
            6


            sage: m = matrix(QQ,10,4,[1,0,0,0, 1,1/16,1/4,1/16, 1,3/8,1/4,1/32, 1,1/4,3/8,1/32, 1,1/16,1/16,1/4, 1,1/32,3/8,1/4, 1,1/4,1/16,1/16, 1,1/32,1/4,3/8, 1,3/8,1/32,1/4, 1,1/4,1/32,3/8])
            sage: p = polymake.Polytope('POINTS',m)
            sage: p.num_facets()
            13
        """
        return self._get_integer_property("N_FACETS")

    def num_points(self):
        """
        EXAMPLES::

            sage: import sage.geometry.polymake as polymake  # optional - polymake
            sage: m = matrix(QQ,[[1,0,0,0], [1,0,0,1], [1,0,1,0], [1,0,1,1],  [1,1,0,0], [1,1,0,1], [1,1,1,0], [1,1,1,1]])
            sage: p = polymake.Polytope('POINTS',m)
            sage: p.num_points()
            8

            sage: m = matrix(QQ,10,4,[1,0,0,0, 1,1/16,1/4,1/16, 1,3/8,1/4,1/32, 1,1/4,3/8,1/32, 1,1/16,1/16,1/4, 1,1/32,3/8,1/4, 1,1/4,1/16,1/16, 1,1/32,1/4,3/8, 1,3/8,1/32,1/4, 1,1/4,1/32,3/8])
            sage: p = polymake.Polytope('POINTS',m)
            sage: p.num_points()
            10
        """
        return self._get_integer_property("N_POINTS")

    def num_vertices(self):
        """
        EXAMPLES::

            sage:
        """
        return self._get_integer_property("N_VERTICES")

    def is_simple(self):
        """
        EXAMPLES::

            sage: import sage.geometry.polymake as polymake  # optional - polymake
            sage: m = matrix(QQ,[[1, 3, 0, 0], [1, 0, 3, 0], [1, 1, 1, 1], [1, 0, 0, 3]])
            sage: p = polymake.Polytope('POINTS',m)
            sage: p.is_simple()
            True

            sage: import sage.geometry.polymake as polymake  # optional - polymake
            sage: m = matrix(QQ,[[1,0,0,0], [1,0,0,1], [1,0,1,0], [1,0,1,1],  [1,1,0,0], [1,1,0,1], [1,1,1,0], [1,1,1,1]])
            sage: p = polymake.Polytope('POINTS',m)
            sage: p.is_simple()
            True

            sage: m = matrix(QQ,10,4,[1,0,0,0, 1,1/16,1/4,1/16, 1,3/8,1/4,1/32, 1,1/4,3/8,1/32, 1,1/16,1/16,1/4, 1,1/32,3/8,1/4, 1,1/4,1/16,1/16, 1,1/32,1/4,3/8, 1,3/8,1/32,1/4, 1,1/4,1/32,3/8])
            sage: p = polymake.Polytope('POINTS',m)
            sage: p.is_simple()
            False
        """
        return self._get_bool_property("SIMPLE")

    def is_simplicial(self):
        """
        EXAMPLES::

            sage: import sage.geometry.polymake as polymake  # optional - polymake
            sage: m = matrix(QQ,[[1, 3, 0, 0], [1, 0, 3, 0], [1, 1, 1, 1], [1, 0, 0, 3]])
            sage: p = polymake.Polytope('POINTS',m)
            sage: p.is_simplicial()
            True
        """
        return self._get_bool_property("SIMPLICIAL")

    def graph(self):
        cdef MatrixRational pm_mat
        cdef PerlObject *graph = new PerlObject("Graph<Undirected>")
        pm_get_PerlObject(self.pm_obj.give("GRAPH"), graph[0])
        pm_get_MatrixRational(graph[0].give("ADJACENCY"), pm_mat)
        # FIXME: this is broken
        # FIXME: how do we read the adjacency matrix?
        return pm_mat_to_sage(pm_mat)

    def visual(self):
        pm.set_preference("jreality")
        self.pm_obj.VoidCallPolymakeMethod("VISUAL")

    def vertices(self):
        """
        EXAMPLES::

            sage: import sage.geometry.polymake as polymake  # optional - polymake
            sage: cube = polymake.cube(3,0)
            sage: cube.vertices()
            [1 0 0 0]
            [1 1 0 0]
            [1 0 1 0]
            [1 1 1 0]
            [1 0 0 1]
            [1 1 0 1]
            [1 0 1 1]
            [1 1 1 1]
        """
        return self._get_matrix_property("VERTICES")

    def facets(self):
        """
        EXAMPLES::

            sage: import sage.geometry.polymake as polymake  # optional - polymake
            sage: cube = polymake.cube(3,0)
            sage: cube.facets()
            [ 0  1  0  0]
            [ 1 -1  0  0]
            [ 0  0  1  0]
            [ 1  0 -1  0]
            [ 0  0  0  1]
            [ 1  0  0 -1]
        """
        return self._get_matrix_property("FACETS")


def new_Polytope_from_function(name, *args):
    pm.set_application("polytope")
    cdef PerlObject pm_obj
    if len(args) == 0:
        pm_obj = CallPolymakeFunction(name)
    elif len(args) == 1:
        pm_obj = CallPolymakeFunction1(name, args[0])
    elif len(args) == 2:
        pm_obj = CallPolymakeFunction2(name, args[0], args[1])
    elif len(args) == 3:
        pm_obj = CallPolymakeFunction3(name, args[0], args[1], args[2])
    else:
        raise NotImplementedError("can only handle 1-3 arguments")
    cdef Polytope res = Polytope.__new__(Polytope)
    res.pm_obj = new_PerlObject_from_PerlObject(pm_obj)
    return res

def cube(dimension, scale=1):
    """
    Return a cube of given dimension.  +/-1-coordinates by default.

    EXAMPLES::

        sage: import sage.geometry.polymake as polymake # optional - polymake
        sage: cube = polymake.cube(3)
    """
    return new_Polytope_from_function("cube", dimension, scale)

def cell24():
    """
    EXAMPLES::

        sage: import sage.geometry.polymake as polymake # optional - polymake
        sage: c24 = polymake.cell24()
    """
    return new_Polytope_from_function("create_24_cell")

def birkhoff(n, even=False):
    """
    Birkhoff polytope

    EXAMPLES::

    """
    return new_Polytope_from_function("birkhoff", n, even)

def associahedron(dim):
    """
    EXAMPLES::

        sage: import sage.geometry.polymake as polymake # optional - polymake
        sage: a3 = polymake.associahedron(3)
    """
    return new_Polytope_from_function("associahedron", dim)

def rand_sphere(dim, npoints):
    """
    Return a random spherical polytope of given dimension and number of points.

    EXAMPLES::

        sage: import sage.geometry.polymake as polymake # optional - polymake
        sage: s3 = polymake.rand_sphere(3,20)
    """
    return new_Polytope_from_function("rand_sphere", dim, npoints)
