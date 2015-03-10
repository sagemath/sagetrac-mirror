"""
This package provides a method of representing Jacobi forms and
computing Jacobi theta series. (See page 1 of [EZ1985]_ for a
definition of Jacobi theta series.)

EXAMPLES:

We can create the Jacobi-theta series for the `E_8` lattice with the
vector (0,1,0,0,0,0,0,0)::

    sage: ME8 = [[2, 0, 0, 0, 0, 0, 0, 1],\
                 [0, 2, 1, 1, 1, 1, 1, 1],\
                 [0, 1, 2, 1, 1, 1, 1, 1],\
                 [0, 1, 1, 2, 1, 1, 1, 1],\
                 [0, 1, 1, 1, 2, 1, 1, 1],\
                 [0, 1, 1, 1, 1, 2, 1, 1],\
                 [0, 1, 1, 1, 1, 1, 2, 0],\
                 [1, 1, 1, 1, 1, 1, 0, 2]]
    sage: QE8 = QuadraticForm(matrix(QQ,8,8,ME8))
    sage: Jacobi_Form(QE8,[0,1,0,0,0,0,0,0],4)
    Jacobi Form of weight 4 and index 1 with Fourier expansion:
    1 + (Z^-2 + 56*Z^-1 + 126 + 56*Z + Z^2)*q + (126*Z^-2 + 576*Z^-1 +
    756 + 576*Z + 126*Z^2)*q^2 + (56*Z^-3 + 756*Z^-2 + 1512*Z^-1 +
    2072 + 1512*Z + 756*Z^2 + 56*Z^3)*q^3 + (Z^-4 + 576*Z^-3 +
    2072*Z^-2 + 4032*Z^-1 + 4158 + 4032*Z + 2072*Z^2 + 576*Z^3 +
    Z^4)*q^4 + O(q^5)

Note that the Jacobi theta series for the zero vector recovers the
classical theta series::

    sage: D = DiagonalQuadraticForm(ZZ,[1,2,3,4])
    sage: Jacobi_Form(D,[0,0,0,0],9)
    Jacobi Form of weight 2 and index 0 with Fourier expansion:
    1 + 2*q + 2*q^2 + 6*q^3 + 8*q^4 + 8*q^5 + 16*q^6 + 16*q^7 + 14*q^8
    + 22*q^9 + O(q^10)
    sage: D.theta_series()
    1 + 2*q + 2*q^2 + 6*q^3 + 8*q^4 + 8*q^5 + 16*q^6 + 16*q^7 + 14*q^8
    + 22*q^9 + O(q^10)

ALGORITHMS:

Jacobi forms are created using the function
``compute_theta_series``. This algorithm is described in Section 5.2
of [deQ2010]_. It uses the method of Lagrange multipliers to find
the vectors x in the lattice for which the quadratic form
q(x) <= the precision which we want for the Fourier expansion.

..TODO::

    - Throw exceptions.
    - Change matrix in ``compute_theta_series`` to an integral matrix.
    - Consider other functions.
    - Add documentation for all the functions.
    - Copywrite statement.

REFERENCES:

.. [deQ2010] Victoria de Quehen. Jacobi Forms. Master's thesis, Department of
    Mathematics, McGill University, 2010.

.. [EZ1985] Martin Eichler and Don Zagier. The Theory of Jacobi
    Forms, Progress in Mathematics, 55, Boston, MA: Birkhäuser Boston,
    1985.

AUTHORS:

- Victoria de Quehen
- Andrew Fiori
- David Roe
"""
cdef extern from "math.h":
    double sqrt(double)
    double floor(double)
    double ceil(double)
    double round(double)
include "../../ext/stdsage.pxi"

# TODO -- Copywrite statement

from warnings import warn
from copy import copy

from sage.matrix.matrix import Matrix
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import IntegerRing, ZZ
from sage.rings.rational_field import RationalField, QQ
from sage.rings.ring import Ring
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.integer import Integer
from sage.rings.rational import Rational

#cdef compute_theta_series_precomp_nomaloc(self,mtx, mtx_inv, vec, index, LRing, PRing, double* mtx_final_col, double* adjust_radius, double * mtxs_inv_tip, double** center):
#        """
#           Secretly an internal function of jacobi forms for the purpose of using precomputed values
#
#           I would prefer to have a forward declaration and then put this code below the one inside to ease commenting
#           the setup of this function, and how it interacts with compute_all_theta is very hackish.
#
#           Functions the same as Jacobi_Form::compute_theta_series except that it takes as parameters values that
#               are purely precomputed
#        """
        # #### Initialize Various Values ####
        # cdef int z_prec = Integer(4*index*self._prec).isqrt() + 1
        # q = self._q
        # cdef int dim = q.dim()
        # cdef int dim_minus_one = dim - 1
        # cdef double error=.01
        # cdef double radius=<double>(self._prec)+.01
        # cdef double orig_radius = radius

        # #### Precomputations ####
        # _bil_mtx = 2*mtx*vec
        # #### Convert Python Variables to C Variables ####
        # cdef int i, k
        # cdef int* bil_mtx = <int*>sage_malloc(dim*sizeof(int))
        # for i from 0 <= i < dim:
        #     bil_mtx[i] = _bil_mtx[i]

        # cdef double* boundary = <double*>sage_malloc(dim*sizeof(double))
        # for i from 0 <= i < dim:
        #     boundary[i] = 0
        # boundary[0]=sqrt(<double>radius*(mtx_inv[0,0]))
        # cdef int* x = <int*>sage_malloc(dim*sizeof(int))
        # for i from 0 <= i < dim:
        #     x[i] = 0
        # x[0] = <int>ceil(-boundary[0])
        # cdef int** coeff = <int**>sage_malloc((self._prec+1)*sizeof(int*))
        # for i from 0 <= i <= self._prec:
        #     coeff[i] = <int*>sage_malloc((z_prec*2)*sizeof(int))
        #     for k from 0 <= k < z_prec*2:
        #         coeff[i][k] = 0
        # cdef double* cen = <double*>sage_malloc(dim*sizeof(double))
        # for i from 0 <= i < dim:
        #     cen[i] = 0

        # #### Initialize ####
        # cdef int num_hyperpl=0
        # cdef double eval_bil = bil_mtx[0]*x[0] # eval_bil = <x, vec>, calculated iteratively
        # cdef int eval_q # eval_q = q(x), calculated iteratively
        # cdef int adjust_eval_q = <int>(2*mtx_final_col[0]*x[0])
        # cdef int top
        # cdef double diff

        # #### Main Iteration ####
        # while 1:
        #     top = <int>floor(boundary[num_hyperpl] + cen[num_hyperpl])
        #     while num_hyperpl < dim - 1 and x[num_hyperpl]<=top:
        #         num_hyperpl=num_hyperpl+1

        #         ## Compute new boundary while we advance in hyperplanes
        #         radius=radius+adjust_radius[num_hyperpl]*(x[num_hyperpl-1]-cen[num_hyperpl-1])**2
        #         if -error<radius and radius<0:
        #             # correct rounding errors
        #             radius=0
        #         boundary[num_hyperpl]=sqrt(radius*mtxs_inv_tip[num_hyperpl])
        #         cen[num_hyperpl]=0
        #         ## Compute the new adjust_center
        #         for i from 1 <= i <= num_hyperpl:
        #             cen[num_hyperpl]=cen[num_hyperpl]+adjust_center[i][num_hyperpl]*(x[i-1]-cen[i-1])
        #         ## Initialize next iteration
        #         x[num_hyperpl] = <int>ceil(-boundary[num_hyperpl]+cen[num_hyperpl])
        #         eval_bil += bil_mtx[num_hyperpl]*x[num_hyperpl]
        #         adjust_eval_q += 2*x[num_hyperpl]*mtx_final_col[num_hyperpl]
        #         top=<int>floor( boundary[num_hyperpl]+ cen[num_hyperpl])

        #     ## Check the next vector
        #     if x[num_hyperpl]<=top:
        #         ## In order to iteratively compute quadratic form we
        #         ## precompute its value at first integral point
        #         diff = boundary[num_hyperpl]-cen[num_hyperpl]+x[num_hyperpl]
        #         eval_q= <int>(orig_radius + adjust_eval_q*diff - diff**2 * mtx_final_col[dim-1] + .5)
        #         # .5 so that it rounds to the nearest integer.  We shouldn't need nearly this much.
        #         ## Check the first one before iterating
        #         coeff[eval_q][<int>(eval_bil+z_prec)] = coeff[eval_q][<int>(eval_bil+z_prec)]+ 1

        #     while x[num_hyperpl] < top:
        #         ## Tight iteration over final dimension
        #         x[num_hyperpl]=x[num_hyperpl]+1
        #         eval_bil +=  bil_mtx[num_hyperpl]
        #         ## q(x+[0,..,1]) = q(x)+B(x,[0,...,0,1])+q([0,..,0,1])
        #         eval_q = <int>(eval_q + adjust_eval_q + mtx_final_col[dim-1])
        #              ## B(x+1,1) = B(x,1) + 2q(1)
        #         adjust_eval_q += 2*mtx_final_col[num_hyperpl]
        #         coeff[eval_q][<int>(eval_bil+z_prec)] = coeff[eval_q][<int>(eval_bil+z_prec)]+ 1

        #     while x[num_hyperpl] >=floor(boundary[num_hyperpl]+cen[num_hyperpl]) - 1 + error:
        #         ## Step back through to higher dimensional hyperplanes we have
        #         ## completed
        #         radius = radius - adjust_radius[num_hyperpl]*(x[num_hyperpl-1]-cen[num_hyperpl-1])**2
        #         eval_bil -= bil_mtx[num_hyperpl]*x[num_hyperpl]
        #         adjust_eval_q -= 2*mtx_final_col[num_hyperpl]*x[num_hyperpl]
        #         num_hyperpl=num_hyperpl-1
        #         if num_hyperpl < 0:
        #             ## We are done now!
        #             coeff_list = []
        #             for i from 0 <= i <= self._prec:
        #                 coeff_list.append([])
        #                 for k from 0 <= k < 2*z_prec:
        #                     coeff_list[i].append(coeff[i][k])
        #             tmp = [ LRing(coeff_list[i],-z_prec) for i in range(self._prec+1) ]
        #             returnresult = PRing(tmp, self._prec+1)
        #             #### Clean up C variables ####
        #             sage_free(bil_mtx)
        #             sage_free(boundary)
        #             sage_free(x)
        #             for i from 0 <= i <= self._prec:
        #                 sage_free(coeff[i])
        #             sage_free(coeff)
        #             sage_free(cen)
        #             return Jacobi_Form(q, vec, self._prec,self._weight, index, self._level, returnresult)

        #     #### Advance the next iteration ####
        #     x[num_hyperpl]=x[num_hyperpl]+1
        #     eval_bil += bil_mtx[num_hyperpl]
        #     adjust_eval_q += 2*mtx_final_col[num_hyperpl]


class Jacobi_Form():
    """
    A generalization of modular forms which can be represented by a
    Fourier expansion. Currently we can only create Jacobi theta
    series.

    .. WARNING::

        If the precision is too small the Jacobi form could be represented
        by a Fourier expansion that does not determine the Jacobi form
        uniquely.

    .. TODO::

        - fill in all the functions.
        - add examples for all the functions.

    AUTHORS:

    - Victoria de Quehen
    - Andrew Fiori
    - David Roe
    """
    def _reduce_basis(self, q, vec):
       """
       Changes basis to try to minimize radius

       TODO - enable once basis of short vectors works
       (ie fix basis of short vectors)
       """
#       M = Matrix(q.basis_of_short_vectors());
#       self._q = quadraticForm(M*q.matrix()*M.transpose())
#       self._vec = M.transpose().inverse()*vec
       self._q = q
       self._vec = vec
       return

    def _compute_radius(self):
        """
        Compute the radius of a circle that encloses the ellipse our
        quadratic form defines.
        """
        from sage.functions.all import sqrt
        MinEigen = min( self._q.matrix().eigenvalues() )
        radius = floor(sqrt(self._prec / MinEigen)) + 1
        return radius

    def _compute_rectangle(self):
        """
        Computes a bounding rectangle for our ellipse.
        """
        from sage.functions.all import sqrt

        mtx = (self._q.matrix()) / 2
        mtx_inv = mtx.inverse()
        boundary = [floor(sqrt(self._prec*abs(mtx_inv[i, i])))
                    for i in range(self._q.dim())]

    def compute_theta_series(self):
        """
        Compute the the Fourier coefficients of this Jacobi theta
        series. (See page 1 of [EZ1985]_ for a definition of Jacobi
        theta series.)

        INPUT:

        - This internal function should is called on an object where
          ``self._q`` is a positive-definite, unimodular quadratic
          form on a lattice and ``self._vec`` is a vector in the
          lattice.

        OUTPUT:

        - An object with the Jacobi form Fourier expansion of
        the Jacobi theta series associated to the quadratic form ._q
        and vector ._vec.

        EXAMPLES:

        We can create the Jacobi-theta series for the E8 lattice with
        the vector (0,1,0,0,0,0,0,0)::

            sage: ME8 = [[2, 0, 0, 0, 0, 0, 0, 1],\
                         [0, 2, 1, 1, 1, 1, 1, 1],\
                         [0, 1, 2, 1, 1, 1, 1, 1],\
                         [0, 1, 1, 2, 1, 1, 1, 1],\
                         [0, 1, 1, 1, 2, 1, 1, 1],\
                         [0, 1, 1, 1, 1, 2, 1, 1],\
                         [0, 1, 1, 1, 1, 1, 2, 0],\
                         [1, 1, 1, 1, 1, 1, 0, 2]]
            sage: QE8 = QuadraticForm(matrix(QQ,8,8,ME8))
            sage: Jacobi_Form(QE8,[0,1,0,0,0,0,0,0],4)
            Jacobi Form of weight 4 and index 1 with Fourier
            expansion: 1 + (Z^-2 + 56*Z^-1 + 126 + 56*Z + Z^2)*q +
            (126*Z^-2 + 576*Z^-1 + 756 + 576*Z + 126*Z^2)*q^2 +
            (56*Z^-3 + 756*Z^-2 + 1512*Z^-1 + 2072 + 1512*Z + 756*Z^2
            + 56*Z^3)*q^3 + (Z^-4 + 576*Z^-3 + 2072*Z^-2 + 4032*Z^-1 +
            4158 + 4032*Z + 2072*Z^2 + 576*Z^3 + Z^4)*q^4 + O(q^5)

        We can create the Jacobi-theta series for the D16^+ lattice
        with the vector (0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)::

            sage: MD16 = [[4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],\
                          [0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],\
                          [0, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],\
                          [0, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],\
                          [0, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],\
                          [0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],\
                          [0, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1],\
                          [0, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1],\
                          [0, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1],\
                          [0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1],\
                          [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1],\
                          [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1],\
                          [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1],\
                          [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1],\
                          [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 0],\
                          [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2]]
            sage: QD16 = QuadraticForm(matrix(QQ,16,16,MD16))
            sage: Jacobi_Form(QD16,[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],4)
            Jacobi Form of weight 8 and index 1 with Fourier
            expansion: 1 + (Z^-2 + 56*Z^-1 + 366 + 56*Z + Z^2)*q +
            (366*Z^-2 + 14016*Z^-1 + 33156 + 14016*Z + 366*Z^2)*q^2 +
            (56*Z^-3 + 33156*Z^-2 + 260712*Z^-1 + 462392 + 260712*Z +
            33156*Z^2 + 56*Z^3)*q^3 + (Z^-4 + 14016*Z^-3 + 462392*Z^-2
            + 1987392*Z^-1 + 2998638 + 1987392*Z + 462392*Z^2 +
            14016*Z^3 + Z^4)*q^4 + O(q^5)

        ALGORITHM:

        This algorithm is described in Section 5.2 of [deQ2010]_.  It
        involves iterating over all the vectors x in the lattice for
        which q(x)<= ``self._prec`` (the precision which we want for the
        Fourier expansion).  The vectors in the lattice which satisfy
        this property all lie within an n-dimensional ellipsoid,
        namely, q(x), where n is the dimension of the lattice.

        Notice if we intersect the ellipsoid defined by this equation
        with a hyperplane we obtain an (n-1)-dimensional ellipsoid. If
        intersect the resulting ellipsoid with another hyperplane we
        obtain an (n-2)-ellipsoid.  Repeating this process we
        eventually obtain a 0-dimensional point.

        If we iterate over all integral hyperplanes in one direction
        we get the (n-1)-ellipsoids that contain exactly those points
        (vectors x) such that q(x)<=._prec.  Following the recursive
        proceedure described in the last paragraph we can iterate over
        all of the integral points (vectors x) such that q(x)<=._prec.

        In order to perform this proceedure efficiently we will
        perform the following steps:

        #. Precomputations which allow us to effectively compute the
        changing position of the centre and boundaries of the
        ellipsoids as we interate through hyperplanes. The boundaries
        are computed using the method of Lagrange multipliers.

        #. Next we convert the python variables to C variables.

        #. As we iterate through the hyperplane, we iteratively
        compute the the values that change with each iteration (for
        example, the center and boundaries of the ellipse, and the
        values of the quadratic and bilinear forms restricted to the
        hyperplane). The computation is iterative in that it adjusts
        the value from the previous iteration rather the completely
        recomputing it.

        #. With each new hyperplane and resulting ellipsoid we recurse
        to the lower dimensions. We also iteratively compute the
        changing values.

        #. Finally we convert the number of vectors with each value of
        #q(x) and B(x,v) to coefficients for the Fourier expansion and
        #we return the object with the Fourier expansion.

        .. WARNING::

            Does not throw exceptions if given the wrong type of data.

        .. TODO::

            - use shortest basis and mention this in input
            - are the lattice condition in input the right ones?
            - general code cleanup
            - change tabs to space
            - throw exceptions
            - use relations between coefficients to make the algorithm far
              more efficient

        AUTHORS:

        - Victoria de Quehen
        - Andrew Fiori
        - David Roe
        """

        #### Variables that are precomputed
        # q = self._q = The quadratic form q(x) whose Jacobi theta
        # series we are computing.
        # self._v = The vector whose Jacobi theta series we are
        # computing.
        # dim = The dimension of the quadratic space.
        # self._prec = The highest exponent for q whose coefficient
        # we will compute.
        # error = The assumed bound on rounding errors from taking square
        # roots. This is used to round back to values known to be integral.
        # orig_radius = self._prec + error to ensure even with
        # rounding errors we still find all points <= self._prec.
        # z_prec = The highest exponent on zeta whose Fourier
        # coefficient we will compute. (This is a result of the bound
        # self._prec on q).
        # mtx = A, (where the quadratic form is q(x)=(1/2)(X^T)AX).
        # mtx_inv = The inverse matrix of MTX.
        # bil_mtx = The vector mtx*(self._v). This is used to compute
        # the bilinear form B(x,self._v).
        # mtx_final_col = The last column of mtx. This is used to
        # iteratively compute how the quadratic form changes as we
        # iterate.
        # mtxs[m] = Lower right (dim-m)x(dim-m)-submatrix of mtx.
        # mtxs_inv[m] = The inverse matrix of mtxs[m].
        # mtxs_inv_tip[m] = The upper left entry of mtx_inv[m]. This is
        # used to compute the boundary of iteration for mth-coordinate
        # of x.
        # adjust_center[m] = A vector calculated at the beginning that is
        # used to calculate how the center of the (dim-m)-dimensional
        # ellipsoid moves as we iterate.
        # adjust_radius[m] = A value used to compute how the bound for
        # the (dim-m) dimensional ellipse changes as we iterate.

        ##### Variables changed during iterations
        # x = The current lattice point on which we are iterating.
        # radius = The distance to the boundary for the current
        # hyperplane.
        # boundary[m] = The computed bounds for the current iteration
        # in the mth-dimension ellipsoid.
        # cen = The centre of the current ellipsoid.
        # coeff[n][r] = The number of vectors we have already found with
        # Q(x)=n and B(x,_v)=r. This is used at the end to calculate
        # the Fourier expansion of the Jacobi form.

        # num_hyperpl = The number of hyperplanes with which you have
        # currently intersected the ellipsoid. Thus the dimension of
        # the current ellipsoid is (dimension-num_hyperpl).

        # eval_bil = The current value of the bilinear form B(X,__V)
        # (computed iteratively).
        # eval_q = The current value of the quadratic form q
        # (computed iteratively).
        # adjust_eval_q = A variable used to iteratively compute
        # eval_q.
        # top = The upper bound on iteration for the lattice point x
        # in the current dimension.
        # diff = Variable to compute distance between boundary and
        # first lattice point on during an iteration.
        """
        This function computes the Fourier coeffiecients of this Jacobi
        theta series.

        In order to do this we must iterate over all the lattice points
        such that q(x) <= _prec.

        General Description of Algorithm:
        We iterate over all points of the lattice such that
        q(x) <= _prec by observing that if we intersect the ellipse
        defined by this equation with a hyperplane we get a lower dimensional
        ellipse. Thus, if we iterate over all integral hyperplanes in one
        direction and follow a recursive proceedure we can iterate over all
        of the integral points.
        In order to perform this proceedure efficiently we will perform
        the following steps:

        1. Precomputations which allow us to effectively compute the
        changing position of the centre and boundaries of the ellipses.

        2. As we iterate through the hyperplane, we iteratively compute
        the the values that change with each iteration (for example the
        the center of the ellipse, the values of the quadratic and
        bilinear forms). The computation is iterative in that it adjusts
        the value from the previous iteration rather the completely
        recomputing it.

        3. As we recurse to the lower dimensions, we also iteratively
        compute the changing values.

        .. TODO::

            - general code cleanup
            - change tabs to space
            - change variable names (consistent capitalization)
            - Improve variable descriptions
            - Most Arithmetic could be made integral, this would make things a little faster (LCM of the det of the upper principle minors for some stuff, 2 for mtx stuff)
        """
        #### Variables that are precomputed
        # q = The Quadratic Form.
        # self._vec = The vector whose Jacobi theta series we are computing.
        # dim = The dimension of the quadratic space.
        # self._prec = The highest exponent for q whose coefficient we will compute.
        # orig_radius = self._prec + epsilon to ensure even with rounding errors we still find all points <= self._prec.
        # z_prec = The highest exponent on zeta whose Fourier coefficient we will compute (this is a result of the bound on q).

        # mtx = 1/2A (where the quadratic form is q(x)=(1/2)(x^T)Ax).
        # mtx_inv = The inverse matrix of mtx.
        # bil_mtx = The vector A*(self._vec) (used to compute the bilinear form B(x,_vec)).
        # mtx_final_col = The last column of mtx, used to iteratively compute how the quadratic form changes as we iterate.
        # mtxs[m] = Lower right (dim-m)x(dim-m) submatrix of mtx.
        # mtxs_inv[m] = The inverse of mtxs[m]
        # mtxs_inv_tip[m] = The upper left entry of mtx_inv[m] which is used to compute the boundary of iteration for mth-coordinate of x.
        # adjust_center[m] = A vector calculated at the beginning which is used to calculate how the center of the (dim-m) dimensional ellipse moves as we iterate.
        # adjust_radius[m] = A value used to compute how the bound for the (dim-m) dimensional ellipse changes as we iterate.
        # error = Assumed bound on rounding errors from taking square roots, used to round back to values known to be integral.

        ##### Variables changed during iterations
        # x = The current lattice point on which we are iterating.
        # radius = The distance to the boundary for the current hyperplane.
        # boundary[m] = The computed bounds for the current iteration in the mth-dimension
        # cen = The centre of the current ellipse.
        # coeff[n][r] = The number of vectors we have already found with q(x)=n and B(x,_vec)=r.
        # num_hyperpl = The dimension we are currently iterating in.
        # eval_bil = The current value of B(x,_vec) (computed iteratively).
        # eval_q = The current value of q(x) (computed iteratively).
        # adjust_eval_q = A variable used to iteratively compute eval_q.
        # top = The upper bound on iteration for x in the current dimension.
        # diff = Variable to compute distance between boundary and first lattice point on during an iteration.

        #### Initialize Various Values ####
        cdef int z_prec = Integer(4*self._index*self._prec).isqrt() + 1
        q = self._q
        cdef int dim = q.dim()
        cdef int dim_minus_one = dim - 1
        cdef double error = .01
        cdef double radius = <double>(self._prec)+.01
        cdef double orig_radius = radius

        #### Precomputations ####
        mtx = (q.matrix())
        mtx_inv = mtx.inverse()
        mtxs = [copy(q.matrix())]
        mtxs_inv = [mtxs[0].inverse()]
        _adjust_center = [Matrix(QQ, [[0] for i in range(dim)])]
        _adjust_radius = [0]
        for m in range(1, dim):
             mtxs[m-1].subdivide(1, 1)
             M2 = mtxs[m-1].subdivision(1, 0)
             mtxs.append(copy(mtxs[m-1].subdivision(1,1)))
             mtxs_inv.append(mtxs[m].inverse())
             _adjust_center.append(-mtxs_inv[m]*M2)
             _adjust_radius.append(-(mtxs[m-1])[0,0]/2+(M2.transpose()*mtxs_inv[m]/2*M2)[0,0])
        _bil_mtx = mtx * self._vec
        LRing = LaurentSeriesRing(ZZ, 'Z')
        PRing = PowerSeriesRing(LRing, 'q')

        #### Convert Python Variables to C Variables ####
        cdef int i, k
        cdef int* bil_mtx = <int*>sage_malloc(dim*sizeof(int))
        for i from 0 <= i < dim:
            bil_mtx[i] = _bil_mtx[i]

        cdef int* mtx_final_col = <int*>sage_malloc(dim*sizeof(int))
        for i from 0 <= i < dim:
            mtx_final_col[i] = <int>(mtx[i, dim-1])

        cdef double* adjust_radius = <double*>sage_malloc(dim*sizeof(double))
        for i from 0 <= i < dim:
            adjust_radius[i] = _adjust_radius[i]

        cdef double* mtxs_inv_tip = <double*>sage_malloc(dim*sizeof(double))
        for i from 0 <= i < dim:
            mtxs_inv_tip[i] = mtxs_inv[i][0,0]

        cdef double** adjust_center = <double**>sage_malloc(dim*sizeof(double*))
        for k from 0 <= k < dim:
            adjust_center[k] = <double*>sage_malloc(dim*sizeof(double))
        for k from 0 <= k < dim:
            for i from 1 <= i <= k:
                adjust_center[i][k] = _adjust_center[i][k-i,0]

        cdef double* boundary = <double*>sage_malloc(dim*sizeof(double))
        for i from 0 <= i < dim:
            boundary[i] = 0
        boundary[0]=sqrt(<double>radius*(2)*(mtx_inv[0,0]))

        cdef int* x = <int*>sage_malloc(dim*sizeof(int))
        for i from 0 <= i < dim:
            x[i] = 0
        x[0] = <int>ceil(-boundary[0])

        cdef int** coeff = <int**>sage_malloc((self._prec+1)*sizeof(int*))
        for i from 0 <= i <= self._prec:
            coeff[i] = <int*>sage_malloc((z_prec*2)*sizeof(int))
            for k from 0 <= k < z_prec*2:
                coeff[i][k] = 0

        cdef double* cen = <double*>sage_malloc(dim*sizeof(double))
        for i from 0 <= i < dim:
            cen[i] = 0

        #### Initialize ####
        cdef int num_hyperpl = 0
        cdef int eval_bil = bil_mtx[0]*x[0] # eval_bil = <x, vec>, calculated iteratively.
        cdef int eval_q # eval_q = q(x), calculated iteratively.
        cdef int adjust_eval_q = mtx_final_col[0] * x[0]
        cdef int top
        cdef double diff

        #### Main Iteration ####
        while 1:
            top = <int>floor(boundary[num_hyperpl] + cen[num_hyperpl])
            while num_hyperpl < dim - 1 and x[num_hyperpl]<=top:
                num_hyperpl += 1

                ## Compute new the boundary while we advance in hyperplanes.
                radius += adjust_radius[num_hyperpl]*(x[num_hyperpl-1]-cen[num_hyperpl-1])**2
                if -error<radius and radius < 0:
                    # Correct any rounding errors.
                    radius = 0
                boundary[num_hyperpl]=sqrt(radius*mtxs_inv_tip[num_hyperpl]*2)
                cen[num_hyperpl] = 0
                ## Compute the new centre.
                for i from 1 <= i <= num_hyperpl:
                    cen[num_hyperpl]=cen[num_hyperpl]+adjust_center[i][num_hyperpl]*(x[i-1]-cen[i-1])
                ## Initialize the next iteration.
                x[num_hyperpl] = <int>ceil(-boundary[num_hyperpl]+cen[num_hyperpl])
                eval_bil += bil_mtx[num_hyperpl]*x[num_hyperpl]
                adjust_eval_q += x[num_hyperpl]*mtx_final_col[num_hyperpl]
                top=<int>floor( boundary[num_hyperpl]+ cen[num_hyperpl])

            ## Check the next vector.
            if x[num_hyperpl]<=top:
                ## In order to iteratively compute quadratic form we
                ## precompute its value at first integral point.
                diff = boundary[num_hyperpl]-cen[num_hyperpl]+x[num_hyperpl]
                eval_q= <int>(orig_radius + adjust_eval_q*diff - diff**2 * mtx_final_col[dim-1]/2 + error)
                # error so that it rounds to the nearest integer.
                ## Check the first one before iterating.
                coeff[eval_q][<int>(eval_bil+z_prec)] = coeff[eval_q][<int>(eval_bil+z_prec)]+ 1

            while x[num_hyperpl] < top:
                ## A tight iteration over final dimension.
                x[num_hyperpl]=x[num_hyperpl]+1
                eval_bil +=  bil_mtx[num_hyperpl]
                ## q(x+[0,..,1]) = q(x)+B(x,[0,...,0,1])+q([0,..,0,1])
                eval_q = <int>(eval_q + adjust_eval_q + mtx_final_col[dim-1]/2)
                ## B(x+1,1) = B(x,1) + 2q(1)
                adjust_eval_q += mtx_final_col[num_hyperpl]
                coeff[eval_q][<int>(eval_bil+z_prec)] = coeff[eval_q][<int>(eval_bil+z_prec)]+ 1

            while x[num_hyperpl] >=floor(boundary[num_hyperpl]+cen[num_hyperpl]) - 1 + error:
                ## Step back through to higher dimensional hyperplanes we have
                ## completed.
                radius = radius - adjust_radius[num_hyperpl]*(x[num_hyperpl-1]-cen[num_hyperpl-1])**2
                eval_bil -= bil_mtx[num_hyperpl]*x[num_hyperpl]
                adjust_eval_q -= mtx_final_col[num_hyperpl]*x[num_hyperpl]
                num_hyperpl=num_hyperpl-1
                if num_hyperpl < 0:
                    ## We are done now!
                    coeff_list = []
                    for i from 0 <= i <= self._prec:
                        coeff_list.append([])
                        for k from 0 <= k < 2*z_prec:
                            coeff_list[i].append(coeff[i][k])
                    tmp = [ LRing(coeff_list[i],-z_prec) for i in range(self._prec+1) ]
                    self._Fourier_expansion = PRing(tmp, self._prec+1)
                    #### Clean up the C variables ####
                    sage_free(bil_mtx)
                    sage_free(mtx_final_col)
                    sage_free(adjust_radius)
                    sage_free(mtxs_inv_tip)
                    for k from 0 <= k < dim:
                        sage_free(adjust_center[k])
                    sage_free(adjust_center)
                    sage_free(boundary)
                    sage_free(x)
                    for i from 0 <= i <= self._prec:
                        sage_free(coeff[i])
                    sage_free(coeff)
                    sage_free(cen)
                    return

            #### Advance the next iteration ####
            x[num_hyperpl] += 1
            eval_bil += bil_mtx[num_hyperpl]
            adjust_eval_q += mtx_final_col[num_hyperpl]

    def __init__(self, form=None, vec=None, prec=5, weight=None,
                 index=None, level=0, Fourier=None):
        """
        Calculate the Jacobi theta series for a given quadratic form and
        vector. See [EZ1985]_ for the definition of a Jacobi theta series.

        INPUT:

        2 or 3 objects in the following order:

        1. a positive-definite quadratic form with integral
           coefficients.

        2. an integral vector in the lattice associated to the
           quadratic form.

        3. (optional) the precision for the Jacobi form, that is,
           the highest power of q in the Fourier expansion.

        EXAMPLES::

            sage: D = DiagonalQuadraticForm(ZZ,[1,1,1,1])
            sage: Jacobi_Form(D,[0,1,0,0],5)
            Jacobi Form of weight 2 and index 1 with Fourier expansion:
            1 + (Z^-2 + 6 + Z^2)*q + (6*Z^-2 + 12 + 6*Z^2)*q^2 + (12*Z^-2
            + 8 + 12*Z^2)*q^3 + (Z^-4 + 8*Z^-2 + 6 + 8*Z^2 + Z^4)*q^4 +
            (6*Z^-4 + 6*Z^-2 + 24 + 6*Z^2 + 6*Z^4)*q^5 + O(q^6)

        ALGORITHMS:

        Jacobi forms are created using the function
        compute_theta_series. This algorithm is described in Section
        5.2 of [deQ2010]_. It uses the method of Lagrange multipliers
        to find the vectors x in the lattice for which the quadratic
        form q(x)<=the precision which we want for the Fourier
        expansion.

        .. TODO::

            - check weight, index, level, Fourier
            - input, output, example
            - input as Fourier series
            - check that the above really is a Jacobi form

        AUTHORS:

        - Victoria de Quehen
        - Andrew Fiori
        - David Roe
        """
        if not isinstance(form,QuadraticForm):
            raise TypeError("Form must be a quadratic form.")

        if form.base_ring() is not ZZ:
            try:
                form = QuadraticForm(form.matrix().change_ring(ZZ))
            except TypeError as err:
                raise TypeError("The matrix associated to the quadratic form must be integral (%s)."%err.msg)
        if not form.is_positive_definite():
            raise ValueError("The quadratic form must be positive-definite.")
        if form.dim() % 2 != 0:
            raise ValueError("The quadratic form must be even dimensional.")
        self._q = form
        self._weight = form.dim() // 2

        try:
            self._vec = vector(ZZ,vec)
        except (AttributeError, TypeError, ValueError):
            raise TypeError("vec must be a vector or list with integral entries")
        if len(vec) != form.dim():
            raise ValueError("The dimension of the quadratic form must equal the dimensional of the vector.")
        self._index = form(self._vec)

        prec = ZZ(prec)
        if prec < 0:
            raise ValueError("prec must be a positive integer")
        self._prec = prec

        self.compute_theta_series()

        # """
        #     Creates a Jacobi form object from the data without needing a basis.
        #     Note that this will not remember a quadratic form nor vector as it may not be a theta series.
        #     Input: arg1 is the Fourier series, arg2 is the level.
        # TODO - check the index, weight and level are the same
        #      - have minimum radius
        #      - check types
        # """
        # if not (form == None):
        #     self._q = form
        # if not (vec == None):
        #     self._vec = vec
        # self._prec =adjust_radius
        # self._index=index
        # self._weight=weight
        # self._level=level
        # self._Fourier_expansion=Fourier
        # return

    def __repr__(self):
        """
        Return the representation of self.
        
        .. TODO:: level
        """
        return "Jacobi Form of weight %i and index %i with Fourier expansion:\n %s"%(self._weight, self._index, self._Fourier_expansion)

    def _latex_(self):
        """
        Return the latex representation of self.

        .. TODO:: make pretty
        """
        return "Jacobi Form of weight %i and index %i with Fourier expansion:\n %s"%(self._weight, self._index, self._Fourier_expansion)

    def __getitem__(self, ij):
        """
        Return the c(i,j) coefficient of self.

        .. TODO::

            - Return 0 if 'j' is large/small
            - Return error if i > prec
            - do we want J[3] to return the coef of q^3 ?
        """
        i, j = ij
        i = int(i)
        j = int(j)
        # gets the coefficient of q^i
        L = self._Fourier_expansion.padded_list(i + 1)[i]
        # returns the coefficient of Z^j
        return L.list()[L.degree()+j]

    def __setitem__(self, ij, coeff):
        """
        TODO - should maybe disallow this...
        """
        raise NotImplementedError("Not allowed to set coeffs of Jacobi forms")

    def __eq__(self, right):
        """
        check equality with right

        TODO - NYI... be sure that this accounds for radius.
                      weight and index are implied by the above working (for non-constant)
                      level... not sure how to care.
        """
        if self._weight != right._weight:
            return False
        if self._index != right._index:
            return False
        return self._Fourier_expansion == right._Fourier_expansion

    def __add__(self, right):
        """
        Return sum of Jacobi forms.

        ToDo check the congruence subgroup and level and create a new Jacobi form

        ToDo We might be able to do more things in the event of theta series if they are
             in the same space and vectors in diff irreducible lattices, this is a theta series
        """
        if self._index != right._index :
            raise NotImplementedError("Oops, the index should be the same")
        if self._weight != right._weight :
            raise NotImplementedError("Oops, the weight should be the same")
        if self._level != right._level :
            raise NotImplementedError("Currently we do not know how to add modular forms of different levels")
        return Jacobi_Form(Fourier = self._Fourier_expansion+right._Fourier_expansion, adjust_radius = min(self._prec,right._prec), weight = self._weight, index = self._index)

    def __mul__(self, right):
        """
        Multiply Jacobi forms.

        Not Yet implemented

        The product of two jacobi theta series is a jacobi theta series
        ToDo we know the Quadratic form and vector to use, could track this
        """
        if self._level != right._level :
            raise NotImplementedError("Currently we do not know how to multiply modular forms of different levels")
        return Jacobi_Form(Fourier = self._Fourier_expansion*right._Fourier_expansion,adjust_radius = min(self._prec,right._prec) ,weight = self._weight+right._weight,index = self._index+right._index, level=self._level)

    def __call__(self, n1, n2):
        """
        Evaluate ... but on what?
        """
        return self._Fourier_expansion(n1)(n2)

    def weight(self):
        """
        Return the weight of the Jacobi form.

        EXAMPLES:
        """
        return self._weight

    def index(self):
        """
        Return the index of the Jacobi form.
        """
        return self._index

    def level(self):
        """
        Return the level of the Jacobi form.

        Not yet implemented.
        """
        raise NotImplementedError("Not yet implemented")
        return self._level

    def Fourier_expansion(self):
        """
        Return the Fourier expansion of the Jacobi form.
        """
        return self._Fourier_expansion

#    def compute_theta_series_precomp(self,mtx, mtx_inv, mtxs_inv,_adjust_center,_adjust_radius,vec,index):
        """
           This function computes the Fourier coeffiecients of this Jacobi
              theta series.
           In order to do this we must iterate over all the lattice points
              such that q(x) <= _prec.

           General Description of Algorithm:
               We iterate over all points of the lattice such that
           q(x) <= _prec by observing that if we intersect the ellipse
           defined by this equation with a hyperplane we get a lower dimensional
           ellipse. Thus, if we iterate over all integral hyperplanes in one
           direction and follow a recursive proceedure we can iterate over all
           of the integral points.
               In order to perform this proceedure efficiently we will perform
           the following steps:
               1. Precomputations which allow us to effectively compute the
               changing position of the centre and boundaries of the ellipses.
               2. As we iterate through the hyperplane, we iteratively compute
               the the values that change with each iteration (for example the
               the center of the ellipse, the values of the quadratic and
               bilinear forms). The computation is iterative in that it adjusts
               the value from the previous iteration rather the completely
               recomputing it.
               3. As we recurse to the lower dimensions, we also iteratively
               compute the changing values.

           TODO
                - general code cleanup
                - reduce boundary
                - change tabs to space
                - change variable names (consistend capitalization)
                - Improve variable descriptions
        """
        #### Variables that are precomputed
        # q = The Quadratic Form.
        # self._vec = The vector whose Jacobi theta series we are computing.
        # dim = The dimension of the quadratic space.
        # self._prec = The highest exponent for q whose coefficient we will compute.
        # orig_radius = self._prec + epsilon to ensure even with rounding errors we still find all points <= self._prec.
        # z_prec = The highest exponent on zeta whose Fourier coefficient we will compute (this is a result of the bound on q).

        # mtx = 1/2A (where the quadratic form is q(x)=(1/2)(x^T)Ax).
        # mtx_inv = The inverse matrix of mtx.
        # bil_mtx = The vector A*(self._vec) (used to compute the bilinear form B(x,_vec)).
        # mtx_final_col = The last column of mtx, used to iteratively compute how the quadratic form changes as we iterate.
        # mtxs[m] = Lower right (dim-m)x(dim-m) submatrix of mtx.
        # mtxs_inv[m] = The inverse of mtxs[m]
        # mtxs_inv_tip[m] = The upper left entry of mtx_inv[m] which is used to compute the boundary of iteration for mth-coordinate of x.
        # adjust_center[m] = A vector calculated at the beginning which is used to calculate how the center of the (dim-m) dimensional ellipse moves as we iterate.
             #### TODO - change to adjust_center ####
        # adjust_radius[m] = A value used to compute how the bound for the (dim-m) dimensional ellipse changes as we iterate.
             #### TODO - change variable name of adjust_radius to Adjustradius ####
        # error = Assumed bound on rounding errors from taking square roots, used to round back to values known to be integral.

        ##### Variables changed during iterations
        # x = The current lattice point on which we are iterating.
        # radius = The distance to the boundary for the current hyperplane.
             #### TODO - rename radius to radius ####
        # boundary[m] = The computed bounds for the current iteration in the mth-dimension
        # cen = The centre of the current ellipse.
        # coeff[n][r] = The number of vectors we have already found with q(x)=n and B(x,_vec)=r.
        # num_hyperpl = The dimension we are currently iterating in.
        # eval_bil = The current value of B(x,_vec) (computed iteratively).
        # eval_q = The current value of q(x) (computed iteratively).
        # adjust_eval_q = A variable used to iteratively compute eval_q.
        # top = The upper bound on iteration for x in the current dimension.
        # diff = Variable to compute distance between boundary and first lattice point on during an iteration.

        #### Initialize Various Values ####
        # cdef int z_prec = Integer(4*index*self._prec).isqrt() + 1
        # q = self._q
        # cdef int dim = q.dim()
        # cdef int dim_minus_one = dim - 1
        # cdef double error=.01
        # cdef double radius=<double>(self._prec)+.01
        # cdef double orig_radius = radius

        # #### Precomputations ####
        # _bil_mtx = 2*mtx*vec
        # LRing = LaurentSeriesRing(ZZ, 'Z')
        # PRing = PowerSeriesRing(LRing, 'q')
        # #### Convert Python Variables to C Variables ####
        # cdef int i, k
        # cdef int* bil_mtx = <int*>sage_malloc(dim*sizeof(int))
        # for i from 0 <= i < dim:
        #     bil_mtx[i] = _bil_mtx[i]

        # cdef double* mtx_final_col = <double*>sage_malloc(dim*sizeof(double))
        # for i from 0 <= i < dim:
        #     mtx_final_col[i] = mtx[i, dim-1]
        # cdef double* adjust_radius = <double*>sage_malloc(dim*sizeof(double))
        # for i from 0 <= i < dim:
        #     adjust_radius[i] = _adjust_radius[i]
        # cdef double* mtxs_inv_tip = <double*>sage_malloc(dim*sizeof(double))
        # for i from 0 <= i < dim:
        #     mtxs_inv_tip[i] = mtxs_inv[i][0,0]
        # cdef double** adjust_center = <double**>sage_malloc(dim*sizeof(double*))
        # for k from 0 <= k < dim:
        #     adjust_center[k] = <double*>sage_malloc(dim*sizeof(double))
        # for k from 0 <= k < dim:
        #     for i from 1 <= i <= k:
        #         adjust_center[i][k] = _adjust_center[i][k-i,0]
        # cdef double* boundary = <double*>sage_malloc(dim*sizeof(double))
        # for i from 0 <= i < dim:
        #     boundary[i] = 0
        # boundary[0]=sqrt(<double>radius*(mtx_inv[0,0]))
        # cdef int* x = <int*>sage_malloc(dim*sizeof(int))
        # for i from 0 <= i < dim:
        #     x[i] = 0
        # x[0] = <int>ceil(-boundary[0])
        # cdef int** coeff = <int**>sage_malloc((self._prec+1)*sizeof(int*))
        # for i from 0 <= i <= self._prec:
        #     coeff[i] = <int*>sage_malloc((z_prec*2)*sizeof(int))
        #     for k from 0 <= k < z_prec*2:
        #         coeff[i][k] = 0
        # cdef double* cen = <double*>sage_malloc(dim*sizeof(double))
        # for i from 0 <= i < dim:
        #     cen[i] = 0

        # #### Initialize ####
        # cdef int num_hyperpl=0
        # cdef double eval_bil = bil_mtx[0]*x[0] # eval_bil = <x, vec>, calculated iteratively
        # cdef int eval_q # eval_q = q(x), calculated iteratively
        # cdef int adjust_eval_q = <int>(2*mtx_final_col[0]*x[0])
        # cdef int top
        # cdef double diff

        # #### Main Iteration ####
        # while 1:
        #     top = <int>floor(boundary[num_hyperpl] + cen[num_hyperpl])
        #     while num_hyperpl < dim - 1 and x[num_hyperpl]<=top:
        #         num_hyperpl=num_hyperpl+1

        #         ## Compute new boundary while we advance in hyperplanes
        #         radius=radius+adjust_radius[num_hyperpl]*(x[num_hyperpl-1]-cen[num_hyperpl-1])**2
        #         if -error<radius and radius<0:
        #             # correct rounding errors
        #             radius=0
        #         boundary[num_hyperpl]=sqrt(radius*mtxs_inv_tip[num_hyperpl])
        #         cen[num_hyperpl]=0
        #         ## Compute the new centre
        #         for i from 1 <= i <= num_hyperpl:
        #             cen[num_hyperpl]=cen[num_hyperpl]+adjust_center[i][num_hyperpl]*(x[i-1]-cen[i-1])
        #         ## Initialize next iteration
        #         x[num_hyperpl] = <int>ceil(-boundary[num_hyperpl]+cen[num_hyperpl])
        #         eval_bil += bil_mtx[num_hyperpl]*x[num_hyperpl]
        #         adjust_eval_q += 2*x[num_hyperpl]*mtx_final_col[num_hyperpl]
        #         top=<int>floor( boundary[num_hyperpl]+ cen[num_hyperpl])

        #     ## Check the next vector
        #     if x[num_hyperpl]<=top:
        #         ## In order to iteratively compute quadratic form we
        #         ## precompute its value at first integral point
        #         diff = boundary[num_hyperpl]-cen[num_hyperpl]+x[num_hyperpl]
        #         eval_q= <int>(orig_radius + adjust_eval_q*diff - diff**2 * mtx_final_col[dim-1] + .5)
        #         # .5 so that it rounds to the nearest integer.  We shouldn't need nearly this much.
        #         ## Check the first one before iterating
        #         coeff[eval_q][<int>(eval_bil+z_prec)] = coeff[eval_q][<int>(eval_bil+z_prec)]+ 1

        #     while x[num_hyperpl] < top:
        #         ## Tight iteration over final dimension
        #         x[num_hyperpl]=x[num_hyperpl]+1
        #         eval_bil +=  bil_mtx[num_hyperpl]
        #         ## q(x+[0,..,1]) = q(x)+B(x,[0,...,0,1])+q([0,..,0,1])
        #         eval_q = <int>(eval_q + adjust_eval_q + mtx_final_col[dim-1])
        #                ## B(x+1,1) = B(x,1) + 2q(1)
        #         adjust_eval_q += 2*mtx_final_col[num_hyperpl]
        #         coeff[eval_q][<int>(eval_bil+z_prec)] = coeff[eval_q][<int>(eval_bil+z_prec)]+ 1

        #     while x[num_hyperpl] >=floor(boundary[num_hyperpl]+cen[num_hyperpl]) - 1 + error:
        #         ## Step back through to higher dimensional hyperplanes we have
        #         ## completed
        #         radius = radius - adjust_radius[num_hyperpl]*(x[num_hyperpl-1]-cen[num_hyperpl-1])**2
        #         eval_bil -= bil_mtx[num_hyperpl]*x[num_hyperpl]
        #         adjust_eval_q -= 2*mtx_final_col[num_hyperpl]*x[num_hyperpl]
        #         num_hyperpl=num_hyperpl-1
        #         if num_hyperpl < 0:
        #             ## We are done now!
        #             coeff_list = []
        #             for i from 0 <= i <= self._prec:
        #                 coeff_list.append([])
        #                 for k from 0 <= k < 2*z_prec:
        #                     coeff_list[i].append(coeff[i][k])
        #             tmp = [ LRing(coeff_list[i],-z_prec) for i in range(self._prec+1) ]
        #             returnresult = PRing(tmp, self._prec+1)
        #             #### Clean up C variables ####
        #             sage_free(bil_mtx)
        #             sage_free(mtx_final_col)
        #             sage_free(adjust_radius)
        #             sage_free(mtxs_inv_tip)
        #             for k from 0 <= k < dim:
        #                 sage_free(adjust_center[k])
        #             sage_free(adjust_center)
        #             sage_free(boundary)
        #             sage_free(x)
        #             for i from 0 <= i <= self._prec:
        #                 sage_free(coeff[i])
        #             sage_free(coeff)
        #             sage_free(cen)
        #             return Jacobi_Form(q, vec, self._prec,self._weight, index, self._level, returnresult)

        #     #### Advance the next iteration ####
        #     x[num_hyperpl]=x[num_hyperpl]+1
        #     eval_bil += bil_mtx[num_hyperpl]
        #     adjust_eval_q += 2*mtx_final_col[num_hyperpl]





#    def compute_all_theta_series(self,coefs, extra=0):
        """
            Code being messed with
        """
        #### Variables that are precomputed
        # q = The quadratic Form.
        # self._vec = The vector whose Jacobi theta series we are computing.
        # dim = The dimension of the quadratic space.
        # self._prec = The highest exponent for q whose coefficient we will compute.
        # orig_radius = self._prec + epsilon to ensure even with rounding errors we still find all points <= self._prec.
        # z_prec = The highest exponent on zeta whose Fourier coefficient we will compute (this is a result of the bound on q).

        # mtx = 1/2A (where the quadratic form is q(x)=(1/2)(x^T)Ax).
        # mtx_inv = The inverse matrix of mtx.
        # bil_mtx = The vector A*(self._vec) (used to compute the bilinear form B(x,_vec)).
        # mtx_final_col = The last column of mtx, used to iteratively compute how the quadratic form changes as we iterate.
        # mtxs[m] = Lower right (dim-m)x(dim-m) submatrix of mtx.
        # mtxs_inv[m] = The inverse of mtxs[m]
        # mtxs_inv_tip[m] = The upper left entry of mtx_inv[m] which is used to compute the boundary of iteration for mth-coordinate of x.
        # adjust_center[m] = A vector calculated at the beginning which is used to calculate how the center of the (dim-m) dimensional ellipse moves as we iterate.
             #### TODO - change to adjust_center ####
        # adjust_radius[m] = A value used to compute how the bound for the (dim-m) dimensional ellipse changes as we iterate.
             #### TODO - change variable name of adjust_radius to Adjustradius ####
        # error = Assumed bound on rounding errors from taking square roots, used to round back to values known to be integral.

        ##### Variables changed during iterations
        # x = The current lattice point on which we are iterating.
        # radius = The distance to the boundary for the current hyperplane.
             #### TODO - rename radius to radius ####
        # boundary[m] = The computed bounds for the current iteration in the mth-dimension
        # cen = The centre of the current ellipse.
        # coeff[n][r] = The number of vectors we have already found with q(x)=n and B(x,_vec)=r.
        # num_hyperpl = The dimension we are currently iterating in.
        # eval_bil = The current value of B(x,_vec) (computed iteratively).
        # eval_q = The current value of q(x) (computed iteratively).
        # adjust_eval_q = A variable used to iteratively compute eval_q.
        # top = The upper bound on iteration for x in the current dimension.
        # diff = Variable to compute distance between boundary and first lattice point on during an iteration.

        # tmpprec = self._prec
        # #### Initialize Various Values ####
        # cdef int z_prec = Integer(4*self._index*self._prec).isqrt() + 1
        # q = self._q
        # cdef int dim = q.dim()
        # cdef int dim_minus_one = dim - 1
        # cdef double error=.01
        # cdef double radius=<double>(self._prec)+.01
        # cdef double orig_radius = radius

        # #### Precomputations ####
        # mtx=(q.matrix())/2
        # mtx_inv=mtx.inverse()
        # mtxs = [copy(q.matrix())/2]
        # mtxs_inv=[mtxs[0].inverse()]
        # _adjust_center=[Matrix(QQ, [[0] for i in range(dim)])]
        # _adjust_radius=[0]
        # for m in range(1, dim):
        #      mtxs[m-1].subdivide(1,1)
        #      M2=mtxs[m-1].subdivision(1,0)
        #      mtxs.append(copy(mtxs[m-1].subdivision(1,1)))
        #      mtxs_inv.append(mtxs[m].inverse())
        #      _adjust_center.append(-mtxs_inv[m]*M2)
        #      _adjust_radius.append(-(mtxs[m-1])[0,0]+(M2.transpose()*mtxs_inv[m]*M2)[0,0])
        # _bil_mtx = 2*mtx*self._vec
        # LRing = LaurentSeriesRing(ZZ, 'Z')
        # PRing = PowerSeriesRing(LRing, 'q')

        # #### Convert Python Variables to C Variables ####
        # cdef int i, k
        # cdef int* bil_mtx = <int*>sage_malloc(dim*sizeof(int))
        # for i from 0 <= i < dim:
        #     bil_mtx[i] = _bil_mtx[i]
        # cdef double* mtx_final_col = <double*>sage_malloc(dim*sizeof(double))
        # for i from 0 <= i < dim:
        #     mtx_final_col[i] = mtx[i, dim-1]
        # cdef double* adjust_radius = <double*>sage_malloc(dim*sizeof(double))
        # for i from 0 <= i < dim:
        #     adjust_radius[i] = _adjust_radius[i]
        # cdef double* mtxs_inv_tip = <double*>sage_malloc(dim*sizeof(double))
        # for i from 0 <= i < dim:
        #     mtxs_inv_tip[i] = mtxs_inv[i][0,0]
        # cdef double** adjust_center = <double**>sage_malloc(dim*sizeof(double*))
        # for k from 0 <= k < dim:
        #     adjust_center[k] = <double*>sage_malloc(dim*sizeof(double))
        # for k from 0 <= k < dim:
        #     for i from 1 <= i <= k:
        #         adjust_center[i][k] = _adjust_center[i][k-i,0]
        # cdef double* boundary = <double*>sage_malloc(dim*sizeof(double))
        # for i from 0 <= i < dim:
        #     boundary[i] = 0
        # boundary[0]=sqrt(<double>radius*(mtx_inv[0,0]))
        # cdef int* x = <int*>sage_malloc(dim*sizeof(int))
        # for i from 0 <= i < dim:
        #     x[i] = 0
        # x[0] = <int>ceil(-boundary[0])
        # cdef double* cen = <double*>sage_malloc(dim*sizeof(double))
        # for i from 0 <= i < dim:
        #     cen[i] = 0

        # #### Initialize ####
        # cdef int num_hyperpl=0
        # cdef double eval_bil = bil_mtx[0]*x[0] # eval_bil = <x, vec>, calculated iteratively
        # cdef int eval_q # eval_q = q(x), calculated iteratively
        # cdef int adjust_eval_q = <int>(2*mtx_final_col[0]*x[0])
        # cdef int top
        # cdef double diff

        # LIST = []
        # VLIST = []
        # COUNT=0

        # self._prec = coefs

        # #### Main Iteration ####
        # while 1:
        #     top = <int>floor(boundary[num_hyperpl] + cen[num_hyperpl])
        #     while num_hyperpl < dim - 1 and x[num_hyperpl]<=top:
        #         num_hyperpl=num_hyperpl+1

        #         ## Compute new boundary while we advance in hyperplanes
        #         radius=radius+adjust_radius[num_hyperpl]*(x[num_hyperpl-1]-cen[num_hyperpl-1])**2
        #         if -error<radius and radius<0:
        #             # correct rounding errors
        #             radius=0
        #         boundary[num_hyperpl]=sqrt(radius*mtxs_inv_tip[num_hyperpl])
        #         cen[num_hyperpl]=0
        #         ## Compute the new centre
        #         for i from 1 <= i <= num_hyperpl:
        #             cen[num_hyperpl]=cen[num_hyperpl]+adjust_center[i][num_hyperpl]*(x[i-1]-cen[i-1])
        #         ## Initialize next iteration
        #         x[num_hyperpl] = <int>ceil(-boundary[num_hyperpl]+cen[num_hyperpl])
        #         eval_bil += bil_mtx[num_hyperpl]*x[num_hyperpl]
        #         adjust_eval_q += 2*x[num_hyperpl]*mtx_final_col[num_hyperpl]
        #         top=<int>floor( boundary[num_hyperpl]+ cen[num_hyperpl])

        #     ## Check the next vector
        #     if x[num_hyperpl]<=top:
        #         ## In order to iteratively compute quadratic form we
        #         ## precompute its value at first integral point
        #         diff = boundary[num_hyperpl]-cen[num_hyperpl]+x[num_hyperpl]
        #         eval_q= <int>(orig_radius + adjust_eval_q*diff - diff**2 * mtx_final_col[dim-1] + .5)
        #         # .5 so that it rounds to the nearest integer.  We shouldn't need nearly this much.
        #         ## Check the first one before iterating
        #         xx = vector([x[i] for i in range(dim)])
        #         J3 = compute_theta_series_precomp_nomaloc(self, mtx, mtx_inv, xx, eval_q, LRing, PRing, mtx_final_col, adjust_radius, mtxs_inv_tip, adjust_center)
        #         #J3 = self.compute_theta_series_precomp(mtx,mtx_inv,mtxs_inv,_adjust_center,_adjust_radius,xx,eval_q)
        #         #J3 = Jacobi_Form(q,xx,coefs)
        #         COUNT = COUNT+1
        #         if (COUNT % 100) == 0:
        #             print COUNT
        #         found = False
        #         for x in LIST:
        #             if J3 == x:
        #                 found = True
        #         if not found:
        #             LIST.append(J3)
        #             VLIST.append(xx)

        #     while x[num_hyperpl] < top:
        #         ## Tight iteration over final dimension
        #         x[num_hyperpl]=x[num_hyperpl]+1
        #         eval_bil +=  bil_mtx[num_hyperpl]
        #         ## q(x+[0,..,1]) = q(x)+B(x,[0,...,0,1])+q([0,..,0,1])
        #         eval_q = <int>(eval_q + adjust_eval_q + mtx_final_col[dim-1])
        #                ## B(x+1,1) = B(x,1) + 2q(1)
        #         adjust_eval_q += 2*mtx_final_col[num_hyperpl]
        #         xx = vector([x[i] for i in range(dim)])
        #         J3 = compute_theta_series_precomp_nomaloc(self, mtx, mtx_inv, xx, eval_q, LRing, PRing, mtx_final_col, adjust_radius, mtxs_inv_tip, adjust_center)
        #         #J3 = self.compute_theta_series_precomp(mtx,mtx_inv,mtxs_inv,_adjust_center,_adjust_radius,xx,eval_q)
        #         #J3 = Jacobi_Form(q,xx,coefs)
        #         COUNT = COUNT+1
        #         if (COUNT % 100) == 0:
        #             print COUNT
        #         found = False
        #         for x in LIST:
        #             if J3 == x:
        #                 found = True
        #         if not found:
        #             LIST.append(J3)
        #             VLIST.append(xx)

        #     while x[num_hyperpl] >=floor(boundary[num_hyperpl]+cen[num_hyperpl]) - 1 + error:
        #         ## Step back through to higher dimensional hyperplanes we have
        #         ## completed
        #         radius = radius - adjust_radius[num_hyperpl]*(x[num_hyperpl-1]-cen[num_hyperpl-1])**2
        #         eval_bil -= bil_mtx[num_hyperpl]*x[num_hyperpl]
        #         adjust_eval_q -= 2*mtx_final_col[num_hyperpl]*x[num_hyperpl]
        #         num_hyperpl=num_hyperpl-1
        #         if num_hyperpl < 0:
        #             ## We are done now!
        #             if extra > 0:
        #                 LIST=[]
        #                 self._prec = self._prec + extra
        #                 for xx in VLIST:
        #                     J3 = compute_theta_series_precomp_nomaloc(self, mtx, mtx_inv, xx, q(xx), LRing, PRing, mtx_final_col, adjust_radius, mtxs_inv_tip, adjust_center)
        #                     LIST.append(J3)

        #             #### Clean up C variables ####
        #             sage_free(bil_mtx)
        #             sage_free(mtx_final_col)
        #             sage_free(adjust_radius)
        #             sage_free(mtxs_inv_tip)
        #             for k from 0 <= k < dim:
        #                 sage_free(adjust_center[k])
        #             sage_free(adjust_center)
        #             sage_free(boundary)
        #             sage_free(x)
        #             sage_free(cen)
        #             self._prec = tmpprec
        #             return LIST

        #     #### Advance the next iteration ####
        #     x[num_hyperpl]=x[num_hyperpl]+1
        #     eval_bil += bil_mtx[num_hyperpl]
        #     adjust_eval_q += 2*mtx_final_col[num_hyperpl]
