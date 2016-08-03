from sage.matrix.constructor import matrix
from sage.matrix.matrix import is_Matrix
from sage.structure.element import is_Vector
from sage.sets.set import Set
from sage.matrix.constructor import column_matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.modules.free_module_element import vector
from sage.modules.free_module import VectorSpace
from sage.rings.integer import Integer
from encoder import Encoder
from decoder import Decoder, DecodingError
from sage.rings.integer_ring import ZZ
from sage.coding.relative_finite_field_extension import *
from sage.structure.sage_object import SageObject

class GabidulinCode(SageObject):
    
    _registered_encoders = {}
    _registered_decoders = {}
    
    def __init__(self, absolute_field, relative_field, length, dimension, \
            linearized_polynomial_ring, evaluation_points=None):
 
        FE = RelativeFiniteFieldExtension(absolute_field, relative_field)
        self._relative_finite_field_extension = FE
        m = absolute_field.degree()/relative_field.degree()
        self._m = m
        V = VectorSpace(relative_field, m)
        self._vector_space = V
        
        self._S = linearized_polynomial_ring #tests to be added

        if not length <= m or length not in ZZ or length < 1:
            raise ValueError("length of the code must be a positive integer less than or equal to \
                    the absolute_field_power which is %d" % m )
        self._length = length
        
        if not dimension <= length or dimension not in ZZ or dimension < 1:
            raise ValueError("dimension of the code must be a positive integer less than or equal \
                    to the length of the code which is %d" % length )
        self._dimension = dimension

        if evaluation_points is None:
            basis = FE.absolute_field_basis()
            evaluation_points = []
            for i in range(length):
                evaluation_points.append(basis[i])
        else:
            if not len(evaluation_points) == length:
                raise ValueError("The number of elements in the set should be equal to the length of the code")
            for i in range(length):
                if not evaluation_points[i] in absolute_field:
                    raise ValueError("Element does not belong to absolute field")
                basis = []
            for i in range(len(evaluation_points)):
                basis.append(FE.relative_field_representation(evaluation_points[i]))
            if V.linear_dependence(basis):
                raise ValueError("The elements provided are not linearly independent")
        self._evaluation_points = evaluation_points

    def _repr_(self):
        return "[%s, %s, %s] Linear Gabidulin Code over %s" \
                % (self.length(), self.dimension(),
                self.minimum_distance(), self.base_field())

    def _latex_(self):
        return "[%s, %s, %s] \\textnormal{ Linear Gabidulin Code over } %s"\
                % (self.length(), self.dimension() ,self.minimum_distance(),
                self.base_field()._latex_())

    def __eq__(self, other):
        return isinstance(other, GabidulinCode) \
                and self.relative_finite_field_extension().absolute_field() == other.relative_finite_field_extension().absolute_field() \
                and self.relative_finite_field_extension().relative_field() == other.relative_finite_field_extension().relative_field() \
                and self.length() == other.length() \
                and self.dimension() == other.dimension() \
                and self.evaluation_points() == other.evaluation_points() \
                
    def minimum_distance(self):
        return self.length() - self.dimension() + 1
    
    # parity_evaluation_points, dual_code and parity_check_matrix need to be fixed.
    def parity_evaluation_points(self):
        eval_pts = self.evaluation_points()
        n = self.length()
        k = self.dimension()
        Fqm = self.base_field()
        q = self.relative_finite_field_extension().relative_field().order()
        coefficient_matrix = matrix(Fqm, n - 1, n, lambda i,j: pow(eval_pts[j], pow(q, -n + k + 1 + i)))
        solution_space = coefficient_matrix.right_kernel()
        parity_eval_pts = solution_space.random_element()
        return parity_eval_pts
 
    def dual_code(self):
        parity_eval_pts = self.parity_evaluation_points()
        return GabidulinCode(self.relative_finite_field_extension().absolute_field(), 
                self.relative_finite_field_extension().relative_field(),
                self.length(), self.length() - self.dimension(), parity_eval_pts)

    def parity_check_matrix(self):
        E = GabidulinGeneratorMatrixEncoder(self.dual_code())
        return E.generator_matrix()

    def rank_weight(self, codeword):
        E = GabidulinPolynomialEvaluationEncoder(self)
        if is_Matrix(codeword):
            return codeword.rank()
        elif is_Vector(codeword):
            codeword_matrix = E.vector_to_matrix(codeword)
            return codeword_matrix.rank()

    def rank_distance(self, a, b):
        E = GabidulinPolynomialEvaluationEncoder(self)
        if is_Vector(a):
            a = E.vector_to_matrix(a)
        if is_Vector(b):
            b = E.vector_to_matrix(b)
        return (a - b).rank()

    def length(self):
        return self._length
    
    def dimension(self):
        return self._dimension
    
    def m(self): #rename
        return self._m

    def vector_space(self):
        return self._vector_space

    def relative_finite_field_extension(self):
        return self._relative_finite_field_extension

    def generator_matrix(self):
        E = GabidulinGeneratorMatrixEncoder(self)
        return E.generator_matrix()

    def base_field(self):
        return self.relative_finite_field_extension().absolute_field()

    def evaluation_points(self):
        return self._evaluation_points

    def message_space(self):
        return self._S


####################### encoders ###############################


####################### encoders ###############################


class GabidulinPolynomialEvaluationEncoder(Encoder):

    def __init__(self, code):
        super(GabidulinPolynomialEvaluationEncoder, self).__init__(code)

    def _repr_(self):
        return "Polynomial evaluation style encoder for %s" % self.code()

    def _latex_(self):
        return "\\textnormal{Polynomial evaluation style encoder for } %s" % self.code()._latex_()

    def __eq__(self, other):
        return isinstance(other, GabidulinPolynomialEvaluationEncoder) \
                and self.code() == other.code()

    def encode(self, p, form="vector"):
        C = self.code()
        M = C.message_space() 
        if p not in M:
            raise ValueError("The value to encode must be in %s" % M)
        if p.degree() >= C.dimension():
            raise ValueError("The polynomial to encode must have degree at most %s" % (C.dimension() - 1))
        eval_pts = C.evaluation_points()
        codeword = []
        for i in range(len(eval_pts)):
            codeword.append(p(eval_pts[i])) #call multi-point evaluation instead of these two lines.
        if form == "vector":
            return vector(codeword)
        elif form == "matrix":
            return self.vector_to_matrix(vector(codeword))
        else:
            return ValueError("The argument 'form' takes only either 'vector' or 'matrix' as valid input")

    def vector_to_matrix(self, codeword_vector):
        C = self.code()
        n = C.length()
        FE = C.relative_finite_field_extension()
        Fqm = FE.absolute_field()
        Fq = FE.relative_field()
        if not is_Vector(codeword_vector):
            raise TypeError("Input must be a vector")
        if not len(codeword_vector) == n:
            raise ValueError("Codeword vector must be of length %d" % n)   # Does not check if vector belongs to code
        codeword_matrix = matrix(Fq, C.m(), n, lambda i,j: FE.relative_field_representation(codeword_vector[j])[i])
        return codeword_matrix

    def matrix_to_vector(self, codeword_matrix):
        C = self.code()
        n = C.length()
        m = C.m()
        FE = self.relative_finite_field_extension()
        Fqm = FE.absolute_field()
        codeword_vector = []
        if not is_Matrix(codeword_matrix):
            raise TypeError("Input must be a matrix")
        if not codeword_matrix.ncols() == n:
            raise ValueError("Number of columns of matrix must be %d" % n)
        if not codeword_matrix.nrows() == m:
            raise ValueError("Number of columns of matrix must be %d" % m)  #Does not check if matrix belongs to code
        for i in range(codeword_matrix.ncols()):
            codeword_vector.append(FE.absolute_field_representation(codeword_matrix.column(i)))
        return vector(codeword_vector)

    def unencode_nocheck(self, c):
        C = self.code()
        eval_pts = C.evaluation_points()
        values = [c[i] for i in range(len(c))]
        p = C.message_space().interpolation_polynomial(eval_pts, values)
        return p

class GabidulinGeneratorMatrixEncoder(Encoder):

    def __init__(self, code):
        super(GabidulinGeneratorMatrixEncoder, self).__init__(code)

    def _repr_(self):
        return "Generator matrix style encoder for %s" % self.code()

    def _latex_(self):
        return "\\textnormal{Generator matrix style encoder for } %s" % self.code()._latex_()

    def __eq__(self, other):
        return isinstance(other, GabidulinGeneratorMatrixEncoder) \
                and self.code() == other.code()
    
    def generator_matrix(self):
        C = self.code()
        eval_pts = C.evaluation_points()
        k = C.dimension()
        sigma = C.message_space().twist_map()
        create_matrix_elements = lambda A,k,f: reduce(lambda L,x: [x] + \
                map(lambda l: map(f,l), L), [A]*k, [])
        return matrix(C.base_field(), C.dimension(), C.length(), \
                create_matrix_elements(eval_pts, C.dimension(), sigma))


####################### decoders ###############################


####################### decoders ###############################


class GabidulinGaoDecoder(Decoder):

    def __init__(self, code):
        if not isinstance(code, GabidulinCode):
            raise ValueError("code has to be a Gabidulin Code")
        super(GabidulinGaoDecoder, self).__init__(code, code.message_space(), "PolynomialEvaluation")

    def _repr_(self):
        return "Gao decoder for %s" % self.code()

    def _latex_(self):
        return "\\textnormal{Gao decoder for } %s" % self.code()._latex_()

    def __eq__(self, other):
        return isinstance(other, GabidulinGaoDecoder) \
            and self.code() == other.code()

    def _partial_xgcd(self, a, b, d_stop):
        C = self.code()
        S = C.message_space()
        if (a not in S) or (b not in S):
            raise ValueError("both the input polynomials must belong to %s" % S)
        if a.degree() < b.degree():
            raise ValueError("degree of first polynomial must be greater than or equal to degree of second polynomial")
        r_previous = a
        r_current = b
        u_previous = S.zero()
        u_current = S.gen()
        v_previous = u_current
        v_current = u_previous

        while r_current.degree() >= d_stop:
            r_hold = r_current
            q_current, r_current = r_current.right_quo_rem(r_previous)
            r_previous = r_hold
            u_hold = u_current
            u_current = u_previous - q_current*u_current
            u_previous = u_hold
            v_hold = v_current
            v_current = v_previous - q_current*v_current
            v_previous = v_hold

        return r_current, u_current, v_current

    def decode_to_message(self, r):
        C = self.code()
        length = len(r)
        if not length <= C.m() or length < 1:
            raise ValueError("length of the received code must be a positive integer \
                    less than or equal to the absolute_field_power which is %d" % m )
        eval_pts = C.evaluation_points()
        S = C.message_space()
        r_ = S.zero()
        for i in range(length):
            e_ = eval_pts[:i] + eval_pts[(i+1):]
            van_poly = S.minimal_vanishing_polynomial(e_) 
            r_ += r[i]*(van_poly/van_poly(eval_pts[i]))
        r_out, u_out, v_out = self._partial_xgcd(S.minimal_vanishing_polynomial(eval_pts), \
                r_, (C.length() + C.dimension())/2)
        quo, rem = r_out.left_quo_rem(u_out)
        print "rem", rem
        print "quo", quo
        if rem == S.zero():
            return quo
        else:
            raise ValueError("Decoding failure")

