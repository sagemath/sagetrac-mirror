"""
This module provides access to polymake, which 'has been developed
since 1997 in the Discrete Geometry group at the Institute of
Mathematics of Technische Universitat Berlin. Since 2004 the
development is shared with Fachbereich Mathematik, Technische
Universitat Darmstadt. The system offers access to a wide variety
of algorithms and packages within a common framework. polymake is
flexible and continuously expanding. The software supplies C++ and
Perl interfaces which make it highly adaptable to individual
needs.'

AUTHORS:

- Ewgenij Gawrilow, Michael Joswig: main authors of polymake

- Burcin Erocal (2011-2012): total rewrite to 
  take advantage of library interface, called pypolymake
  
- Timo Kluck (2012-2013): integrate pypolymake fully into sage,
  making sure to keep backwards compatibility with older version of
  polymake interface

- William Stein: first Sage interface

EXAMPLES:

Use :func:`polytope` to construct cones. The easiest way it to define it by a 
set of affine inequalities::

    sage: import sage.libs.polymake.polymake as pm    # optional - polymake
    sage: p1 = pm.polytope([x < 3, x > 0], coordinates=(x,))  # optional - polymake
    sage: p2 = pm.polytope([x < 2, x > -3], coordinates=(x,))  # optional - polymake
    sage: p3 = p1.intersection(p2)  # optional - polymake
    sage: p3.get_defining_equations()  # optional - polymake
    [x > 0, -x + 2 > 0]
    sage: p3.volume()  # optional - polymake
    2

One can also define lower dimensional polytopes by specifying additional equations::

    sage: x,y,z = var('x,y,z')
    sage: p1 = pm.polytope([x+y+z==1, x> 1/10, y > 1/10, z > 1/10], coordinates = (x,y,z))  # optional - polymake
    sage: p1.get_defining_equations()  # optional - polymake
    [x + y + z - 1 == 0, x - 1/10 > 0, y - 1/10 > 0, z - 1/10 > 0]
    sage: p1.is_feasible()  # optional - polymake
    True
    sage: p1.volume()  # optional - polymake
    0

Doctests from the previous version of this interface::

    sage: P = pm.convex_hull([[1,0,0,0], [1,0,0,1], [1,0,1,0], [1,0,1,1],  [1,1,0,0], [1,1,0,1], [1,1,1,0], [1,1,1,1]])   # optional - polymake
    sage: P.facets()                            # optional - polymake
    [(0, 0, 0, 1), (0, 1, 0, 0), (0, 0, 1, 0), (1, 0, 0, -1), (1, 0, -1, 0), (1, -1, 0, 0)]
    sage: P.vertices()                            # optional - polymake
    [(1, 0, 0, 0), (1, 0, 0, 1), (1, 0, 1, 0), (1, 0, 1, 1), (1, 1, 0, 0), (1, 1, 0, 1), (1, 1, 1, 0), (1, 1, 1, 1)]
    sage: P.is_simple()                              # optional - polymake
    True
    sage: P.n_facets()                            # optional - polymake
    6
    sage: pm.cell24()            # optional - polymake
    The 24-cell
    sage: R.<x,y,z> = PolynomialRing(QQ,3)
    sage: f = x^3 + y^3 + z^3 + x*y*z
    sage: e = f.exponents()
    sage: a = [[1] + list(v) for v in e]
    sage: a
    [[1, 3, 0, 0], [1, 0, 3, 0], [1, 1, 1, 1], [1, 0, 0, 3]]
    sage: n = pm.convex_hull(a)       # optional - polymake
    sage: n                                 # optional - polymake
    Convex hull of points [[1, 0, 0, 3], [1, 0, 3, 0], [1, 1, 1, 1], [1, 3, 0, 0]]
    sage: n.facets()                        # optional - polymake
    [(0, 1, 0, 0), (3, -1, -1, 0), (0, 0, 1, 0)]
    sage: n.is_simple()                     # optional - polymake
    True
    sage: n.graph()                         # optional - polymake
    'GRAPH\n{1 2}\n{0 2}\n{0 1}\n\n'
    
CITE:

Ewgenij Gawrilow and Michael Joswig. polymake: a framework for analyzing
convex polytopes. Polytopes—combinatorics and computation (Oberwolfach,
1997), 43–73, DMV Sem., 29, Birkhäuser, Basel, 2000. MR1785292 (2001f:52033).

"""

###############################################################################
#       Copyright (C) 2011-2012 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2012-2013 Timo Kluck <tkluck@infty.nl>
#       Copyright (C) 2013 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

# This does some magic to ensure that __module__ attributes are right, thus
# ensuring that pickling is possible
include "../../ext/stdsage.pxi"

from libc.stdlib cimport malloc
from libcpp.string cimport string
import operator

from sage.matrix.matrix_rational_dense cimport Matrix_rational_dense
from sage.matrix.constructor import matrix
from sage.modules.all import vector as SageVector
from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.structure.sage_object cimport SageObject
from sage.rings.integer cimport Integer as SageInteger
from sage.rings.rational cimport Rational as SageRational
from sage.libs.gmp.mpz cimport mpz_set
from sage.libs.gmp.mpq cimport mpq_set
from sage.rings.all import QQ, ZZ
from sage.matrix.matrix_space import MatrixSpace
from sage.symbolic.expression import Expression
from sage.symbolic.ring import SR

include "defs.pxd"

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
    cdef PerlObject* _polymake_obj
    cdef object _coordinates
    cdef object _name
    def __init__(self, data=None, coordinates=None, name=None):
        """
        return a new polytope object

        """
        cdef MatrixRational* pm_mat
        cdef PerlObject perlobj_data
        self._coordinates = coordinates
        self._name = "A Polytope" if name is None else name
        if data in ['VERTICES', 'POINTS', 'FACETS']: 
            # compatibility with pypolymake: assume the
            # second argument is a matrix containing the data
            self._polymake_obj = new PerlObject("Polytope<Rational>")
            pm_mat = sage_mat_to_pm(coordinates)
            pm_assign(self._polymake_obj.take(data), pm_mat[0])
            del pm_mat
        elif isinstance(data, str):
            # backward compatibity, assume that we got a
            # (filename, description) pair
            self.load(data)
            if isinstance(coordinates, str):
                self._name = coordinates
                self._coordinates = None
        else:
            self._polymake_obj = new PerlObject("Polytope<Rational>")
            if data:
                if isinstance(data, list):
                    if isinstance(data[0], (Expression, bool)):
                        # Polytope([x < 2, x > 1], (1,x))
                        self.set_defining_equations(data, coordinates)
                    else:
                        # assume that it is numeric and that they represent inequalities
                        self._set_matrix_property("INEQUALITIES", data)
                    
    def set_defining_equations(self, equations, coordinates):
        self._coordinates = coordinates
            
        inequalities_coefficients = []
        equalities_coefficients = []
        for expr in equations:
            if expr is True:
                continue
            elif expr is False:
                inequalities_coefficients.append([-1] + [0 for _ in coordinates])
                continue
            elif expr.operator() in (operator.le, operator.lt):
                zero_expr = expr.rhs() - expr.lhs()
            elif expr.operator() in (operator.ge, operator.gt):
                zero_expr = expr.lhs() - expr.rhs()
            elif expr.operator() == operator.eq:
                zero_expr = expr.lhs() - expr.rhs()
            else:
                raise ValueError("equations should have operator >, <, >=, <= or ==")

            constant_part = QQ(zero_expr.subs({c:0 for c in coordinates}))
            coeffs = [constant_part] + [QQ(zero_expr.coefficient(var)) for var in coordinates]
            homogeneous_coordinates = [1] + list(coordinates)
            if sum(coeff*var for coeff, var in zip(coeffs, homogeneous_coordinates)) != zero_expr:
                raise ValueError("equation is not affine")

            if expr.operator() == operator.eq:
                equalities_coefficients.append(coeffs)
            else:
                inequalities_coefficients.append(coeffs)
        self._set_matrix_property("INEQUALITIES", matrix(QQ, inequalities_coefficients))
        self._set_matrix_property("EQUATIONS", matrix(QQ, equalities_coefficients))

    def get_defining_equations(self):
        """
        Return the defining equations.

        .. TODO::

            avoid using the symbolic ring

        EXAMPLES::
        """
        ineq_coeffs = self.facets()
        if self._coordinates is None:
            self._coordinates = tuple(SR.symbol() for _ in ineq_coeffs[0])
        homogeneous_coordinates = [1] + list(self._coordinates)
        ineqs = [sum(c * var
                     for c, var in zip(coeff, homogeneous_coordinates)) > 0
                 for coeff in ineq_coeffs]
        eq_coeffs = self.equations()
        eqs = [sum(c * var
                   for c, var in zip(coeff, homogeneous_coordinates)) == 0
               for coeff in eq_coeffs]
        return eqs + ineqs

    def __dealloc__(self):
        """
        De-allocate the object.
        """
        del self._polymake_obj

    def _repr_(self):
        """
        Return the string representation.

        EXAMPLES::

            sage: import sage.libs.polymake.polymake as pm    # optional - polymake
            sage: p1 = pm.polytope([x < 3, x > 0], coordinates=(x,))  # optional - polymake
            sage: p1._repr_()
            ?
        """
        return self._name if self._name else 'A Polytope'

    def load(self, filename):
        cdef char* fn = <char*>filename
        self._polymake_obj = new_PerlObject_from_PerlObject(load(fn))
        self._name = "A Polytope"

    def save(self, filename):
        """
        Save this polytope to a file using polymake's representation.

        INPUT:

        - filename

        EXAMPLES::
        """
        self._polymake_obj.save(filename)

    def _get_bool_property(self, prop):
        cdef Integer pm_res
        pm_get(self._polymake_obj.give(prop),pm_res)
        return bool(pm_res.compare(0))

    def _get_integer_property(self, prop):
        cdef Integer pm_res
        pm_get(self._polymake_obj.give(prop), pm_res)
        cdef SageInteger res = SageInteger.__new__(SageInteger)
        mpz_set(res.value, pm_res.get_rep())
        return res

    def _get_rational_property(self, prop):
        cdef Rational pm_res
        pm_get_Rational(self._polymake_obj.give(prop), pm_res)
        cdef SageRational res = SageRational.__new__(SageRational)
        mpq_set(res.value, pm_res.get_rep())
        return res

    def _get_matrix_property(self, prop):
        cdef MatrixRational pm_mat
        pm_get_MatrixRational(self._polymake_obj.give(prop), pm_mat)
        return pm_mat_to_sage(pm_mat)

    def _get_vector_list_property(self, prop):
        return [SageVector(row) for row in self._get_matrix_property(prop)]

    def _get_string_property(self, prop):
        cdef string pm_string
        pm_get_String(self._polymake_obj.give(prop), pm_string)
        return pm_string

    def _set_matrix_property(self, prop, value):
        cdef MatrixRational* pm_mat = sage_mat_to_pm(value)
        pm_assign(self._polymake_obj.take(prop), pm_mat[0])
        del pm_mat

    def _get_vector_property(self, prop):
        cdef VectorInteger pm_vec
        pm_get_VectorInteger(self._polymake_obj.give(prop), pm_vec)
        return pm_VectorInteger_to_sage(pm_vec)

    def __add__(left, right):
        if not (isinstance(left, Polytope) and isinstance(right, Polytope)):
            raise TypeError("both arguments must be instances of Polytope")

    def f_vector(self):
        """
        Return the f-vector.

        EXAMPLES::
        """
        return self._get_vector_property("F_VECTOR")

    def h_star_vector(self):
        """
        Return the h*-vector.

        EXAMPLES::
        """
        return self._get_vector_property("H_STAR_VECTOR")

    def num_facets(self):
        """
        Return the number of facets.

        EXAMPLES::

            sage: import sage.libs.polymake.polymake as pm  # optional - polymake
            sage: m = matrix(QQ,[[1,0,0,0], [1,0,0,1], [1,0,1,0], [1,0,1,1],  [1,1,0,0], [1,1,0,1], [1,1,1,0], [1,1,1,1]])
            sage: p = pm.Polytope('POINTS',m)  # optional - polymake
            sage: p.num_facets()  # optional - polymake
            6


            sage: m = matrix(QQ,10,4,[1,0,0,0, 1,1/16,1/4,1/16, 1,3/8,1/4,1/32, 1,1/4,3/8,1/32, 1,1/16,1/16,1/4, 1,1/32,3/8,1/4, 1,1/4,1/16,1/16, 1,1/32,1/4,3/8, 1,3/8,1/32,1/4, 1,1/4,1/32,3/8])
            sage: p = pm.Polytope('POINTS',m)  # optional - polymake
            sage: p.num_facets()  # optional - polymake
            13
        """
        return self._get_integer_property("N_FACETS")

    def num_points(self):
        """
        Return the number of points.

        EXAMPLES::

            sage: import sage.libs.polymake.polymake as pm  # optional - polymake
            sage: m = matrix(QQ,[[1,0,0,0], [1,0,0,1], [1,0,1,0], [1,0,1,1],  [1,1,0,0], [1,1,0,1], [1,1,1,0], [1,1,1,1]])
            sage: p = pm.Polytope('POINTS',m)  # optional - polymake
            sage: p.num_points()  # optional - polymake
            8

            sage: m = matrix(QQ,10,4,[1,0,0,0, 1,1/16,1/4,1/16, 1,3/8,1/4,1/32, 1,1/4,3/8,1/32, 1,1/16,1/16,1/4, 1,1/32,3/8,1/4, 1,1/4,1/16,1/16, 1,1/32,1/4,3/8, 1,3/8,1/32,1/4, 1,1/4,1/32,3/8])
            sage: p = pm.Polytope('POINTS',m)  # optional - polymake
            sage: p.num_points()  # optional - polymake
            10
        """
        return self._get_integer_property("N_POINTS")

    def num_vertices(self):
        """
        Return the number of vertices.

        EXAMPLES::

            sage: import sage.libs.polymake.polymake as pm  # optional - polymake
            sage: m = matrix(QQ,[[1,0,0,0], [1,0,0,1], [1,0,1,0], [1,0,1,1],  [1,1,0,0], [1,1,0,1], [1,1,1,0], [1,1,1,1]])
            sage: p = pm.Polytope('POINTS',m)  # optional - polymake
            sage: p.num_vertices()  # optional - polymake
            ?
        """
        return self._get_integer_property("N_VERTICES")

    # some useful aliases

    n_points = num_points
    n_vertices = num_vertices
    n_facets = num_facets

    def is_simple(self):
        """
        Return whether the polytope is simple.

        EXAMPLES::

            sage: import sage.libs.polymake.polymake as pm  # optional - polymake
            sage: m = matrix(QQ,[[1, 3, 0, 0], [1, 0, 3, 0], [1, 1, 1, 1], [1, 0, 0, 3]])
            sage: p = pm.Polytope('POINTS',m)  # optional - polymake
            sage: p.is_simple()  # optional - polymake
            True

            sage: m = matrix(QQ,[[1,0,0,0], [1,0,0,1], [1,0,1,0], [1,0,1,1],  [1,1,0,0], [1,1,0,1], [1,1,1,0], [1,1,1,1]])
            sage: p = pm.Polytope('POINTS',m)  # optional - polymake
            sage: p.is_simple()  # optional - polymake
            True

            sage: m = matrix(QQ,10,4,[1,0,0,0, 1,1/16,1/4,1/16, 1,3/8,1/4,1/32, 1,1/4,3/8,1/32, 1,1/16,1/16,1/4, 1,1/32,3/8,1/4, 1,1/4,1/16,1/16, 1,1/32,1/4,3/8, 1,3/8,1/32,1/4, 1,1/4,1/32,3/8])
            sage: p = pm.Polytope('POINTS',m)  # optional - polymake
            sage: p.is_simple()  # optional - polymake
            False
        """
        return self._get_bool_property("SIMPLE")

    def is_simplicial(self):
        """
        Return whether the polytope is simplicial.

        EXAMPLES::

            sage: import sage.libs.polymake.polymake as pm  # optional - polymake
            sage: m = matrix(QQ,[[1, 3, 0, 0], [1, 0, 3, 0], [1, 1, 1, 1], [1, 0, 0, 3]])
            sage: p = pm.Polytope('POINTS',m)  # optional - polymake
            sage: p.is_simplicial()  # optional - polymake
            True
        """
        return self._get_bool_property("SIMPLICIAL")

    def is_bounded(self):
        """
        Return whether ``self`` is bounded.

        EXAMPLES::
        """
        return self._get_bool_property("BOUNDED")

    def is_centered(self):
        """
        Return whether ``self`` is centered.

        EXAMPLES::
        """
        return self._get_bool_property("CENTERED")

    def is_feasible(self):
        """
        Return whether ``self`` is feasible.

        EXAMPLES::
        """
        return self._get_bool_property("FEASIBLE")

    def affine_hull(self):
        """
        Return the affine hull.

        EXAMPLES::
        """
        return self._get_vector_list_property("AFFINE_HULL")

    def gale_transform(self):
        """
        Return the Gale transform.

        EXAMPLES::
        """
        return self._get_matrix_property("GALE_TRANSFORM")

    def _points(self):
        return self._get_vector_list_property("POINTS")

    def graph(self):
        return self._get_string_property("GRAPH")
        #cdef MatrixRational pm_mat
        #cdef PerlObject *graph = new PerlObject("Graph<Undirected>")
        #pm_get_PerlObject(self._polymake_obj.give("GRAPH"), graph[0])
        #pm_get_MatrixRational(graph[0].give("ADJACENCY"), pm_mat)
        ## FIXME: this is broken
        ## FIXME: how do we read the adjacency matrix?
        #return pm_mat_to_sage(pm_mat)

    def volume(self):
        """
        Return the volume of ``self``.

        EXAMPLES::

            sage: import sage.libs.polymake.polymake as pm    # optional - polymake
            sage: p1 = pm.polytope([x < 3, x > 0], coordinates=(x,))  # optional - polymake
            sage: p1.volume()  # optional - polymake
            3?
        """
        return self._get_rational_property("VOLUME")

    def visual(self):
        main.set_preference("jreality")
        self._polymake_obj.VoidCallPolymakeMethod("VISUAL")

    def vertices(self):
        return self._get_vector_list_property("VERTICES")

    def vertices(self):
        """
        Return the vertices of ``self``.

        EXAMPLES::

            sage: import sage.libs.polymake.polymake as pm  # optional - polymake
            sage: cube = pm.cube(3,0)
            sage: cube.vertices()
            [(1, 0, 0, 0), (1, 1, 0, 0), (1, 0, 1, 0), (1, 1, 1, 0), (1, 0, 0, 1), (1, 1, 0, 1), (1, 0, 1, 1), (1, 1, 1, 1)]
        """
        return self._get_matrix_property("VERTICES")

    def equations(self):
        try:
            return self._get_vector_list_property("LINEAR_SPAN")
        except ValueError:
            return []

    def facets(self):
        return self._get_vector_list_property("FACETS")

    def facets(self):
        """
        Return the facets of ``self``.

        EXAMPLES::

            sage: import sage.libs.polymake.polymake as pm  # optional - polymake
            sage: cube = pm.cube(3,0)
            sage: cube.facets()
            [(0, 1, 0, 0), (1, -1, 0, 0), (0, 0, 1, 0), (1, 0, -1, 0), (0, 0, 0, 1), (1, 0, 0, -1)]
        """
        return self._get_matrix_property("FACETS")

    def intersection(self, other):
        """
        Return the intersection of ``self`` with ``other``.

        EXAMPLES::
        """
        cdef Polytope s = <Polytope?>self
        cdef Polytope o = <Polytope?>other
        cdef PerlObject polymake_obj = CallPolymakeFunction_PerlObject2("intersection", s._polymake_obj[0], o._polymake_obj[0])
        return new_Polytope_from_PerlObject(polymake_obj, coordinates=s._coordinates)

    def join(self, other):
        """
        Return the join of ``self`` with ``other``.

        EXAMPLES::
        """
        cdef Polytope s = <Polytope?>self
        cdef Polytope o = <Polytope?>other
        cdef PerlObject polymake_obj = CallPolymakeFunction_PerlObject2("join_polytopes", s._polymake_obj[0], o._polymake_obj[0])
        return new_Polytope_from_PerlObject(polymake_obj, coordinates=s._coordinates)

    def __richcmp__(self, other, operation):
        if operation != 2:
            raise NotImplementedError
        # per the Cython spec, it is possible for some symmetric
        # operations to be called on the second operand. This means
        # that we cannot even be sure that `self` is a Polytope, and
        # we need to cast it. For consistency, we do this cast even in
        # non-special methods.
        cdef Polytope s = <Polytope?>self
        cdef Polytope o = <Polytope?>other
        return BoolCallPolymakeFunction_PerlObject2("equal_polyhedra", s._polymake_obj[0], o._polymake_obj[0])

    def congruent(self, other):
        """
        Return whether ``self`` is congruent to ``other``.

        EXAMPLES::
        """
        cdef Polytope s = <Polytope?>self
        cdef Polytope o = <Polytope?>other
        # FIXME: we'd like to return the Rational that gives the
        # congruence scale factor, but I haven't found a way to
        # Cythonize a call to numerator(Rational) etc.
        return BoolCallPolymakeFunction_PerlObject2("congruent", s._polymake_obj[0], o._polymake_obj[0])

    def __reduce__(self):
        """
        EXAMPLES::

            sage: import sage.libs.polymake.polymake as pm
            sage: cube = pm.cube(3)
            sage: from cPickle import loads, dumps
            sage: cube2 = loads(dumps(cube))
            sage: cube2 == cube
            True
        """
        import tempfile
        import os
        # we need to add the extension '.poly', or polymake will do it for us
        with tempfile.NamedTemporaryFile(suffix='.poly', delete=False) as file:
            file.close()
            ## argghh race condition! Why doesn't polymake allow us to save to a stream?
            self.save(file.name)
            data = open(file.name,'r').read()
            os.unlink(file.name)
        return (new_Polytope_from_data, (data, self._coordinates, self._name))

properties_to_wrap = """MINIMAL_VERTEX_ANGLE: common::Float
 
ONE_VERTEX: common::Vector<Scalar>
 
STEINER_POINTS: common::Matrix<Scalar, NonSymmetric>
 
VALID_POINT: common::Vector<Scalar>
 
VERTEX_BARYCENTER: common::Vector<Scalar>
 
VERTEX_LABELS: common::Array<String>
 
VERTEX_NORMALS: common::Matrix<Scalar, NonSymmetric>
 
ZONOTOPE_INPUT_VECTORS: common::Matrix<Scalar, NonSymmetric>
    def _
    """

# Sage's convention is to (also) have a lowercase function for
# constructing objects
polytope = Polytope

cdef new_Polytope_from_PerlObject(PerlObject polymake_obj, coordinates=None,
                                  name=None):
    cdef Polytope res = Polytope.__new__(Polytope)
    res._polymake_obj = new_PerlObject_from_PerlObject(polymake_obj)
    res._coordinates = coordinates
    res._name = name
    return res

def new_Polytope_from_data(data, coordinates, name):
    cdef Polytope res = Polytope.__new__(Polytope)
    import tempfile
    import os
    with tempfile.NamedTemporaryFile(delete=False) as file:
        file.write(data)
        file.close()
        res.load(file.name)
        os.unlink(file.name)
    res._coordinates = coordinates
    res._name = name
    return res

def new_Polytope_from_function(name, function_name, *args):
    cdef PerlObject polymake_obj
    if len(args) == 0:
        polymake_obj = CallPolymakeFunction(function_name)
    elif len(args) == 1:
        polymake_obj = CallPolymakeFunction1(function_name, args[0])
    elif len(args) == 2:
        polymake_obj = CallPolymakeFunction2(function_name, args[0], args[1])
    elif len(args) == 3:
        polymake_obj = CallPolymakeFunction3(function_name, args[0], args[1], args[2])
    else:
        raise NotImplementedError("can only handle 1-3 arguments")
    return new_Polytope_from_PerlObject(polymake_obj, coordinates=None, name=name)

def cube(dimension, scale=1):
    """
    Return a cube of given dimension.

    The cube has +/-1-coordinates by default.

    INPUT:

    - ``scale`` -- default 1

    EXAMPLES::

        sage: import sage.libs.polymake.polymake as pm  # optional - polymake
        sage: cube = pm.cube(3)  # optional - polymake
    """
    return new_Polytope_from_function("Cube of dimension %s (scale %s)" % (dimension, scale), "cube", dimension, scale)

def cell24():
    """
    Return the 24-cell.

    EXAMPLES::

        sage: import sage.libs.polymake.polymake as pm  # optional - polymake
        sage: c24 = pm.cell24()  # optional - polymake
    """
    return new_Polytope_from_function("The 24-cell", "create_24_cell")

def birkhoff(n, even=False):
    """
    Return the Birkhoff polytope.

    EXAMPLES::

        sage: import sage.libs.polymake.polymake as pm  # optional - polymake
        sage: c24 = pm.birkhoff()  # optional - polymake
    """
    return new_Polytope_from_function("%s Birkhoff %s" % ("even" if even else "odd", n), "birkhoff", n, even)

def associahedron(dim):
    """
    Return the associahedron.

    EXAMPLES::

        sage: import sage.libs.polymake.polymake as pm  # optional - polymake
        sage: a3 = pm.associahedron(3)  # optional - polymake
    """
    return new_Polytope_from_function("%s-dimensional associahedron" % dim, "associahedron", dim)

def rand_sphere(dim, npoints):
    """
    Return a random spherical polytope of given dimension and number of points.

    EXAMPLES::

        sage: import sage.libs.polymake.polymake as pm  # optional - polymake
        sage: s3 = pm.rand_sphere(3,20)  # optional - polymake
    """
    return new_Polytope_from_function("Random spherical polytope", "rand_sphere", dim, npoints)

def convex_hull(points=[]):
    """
    Return the convex hull of the given points.

    EXAMPLES::
    """
    points.sort()
    return Polytope("POINTS", matrix(QQ, points), name="Convex hull of points %s" % points)

def from_data(data):
    raise NotImplementedError

def rand01(d, n, seed=None):
    raise NotImplementedError

# "none" signifies that we don't want polymake to save user configuration data
cdef Main* main
main = new Main(<char*>"none")
main.set_application("polytope")
main.set_custom("$Polymake::User::Verbose::external", 0)
main.set_custom("$Polymake::User::Verbose::credits", 0)

# TODO: delete main upon unloading this module
