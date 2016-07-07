from sage.matrix.constructor import matrix
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
    
    def __init__(self, absolute_field, relative_field, length, dimension, linearly_independent_elements=None):
        if not absolute_field.is_finite() or not absolute_field.is_field():
            raise ValueError("absolute_field must be a finite field")
        if not relative_field.is_finite() or not relative_field.is_field():
            raise ValueError("relative_field must be a finite field")
        p = relative_field.characteristic()
        q = relative_field.order()
        s = log(q, p)
        m = log(absolute_field.order(), q)
        sm = log(absolute_field.order(), p)
        if not s.divides(sm):
            raise ValueError("relative_field must be a subfield of the absolute_field")
        self._prime_field = relative_field.base_ring()
        self._relative_field = relative_field
        self._relative_field_power = s
        self._absolute_field = absolute_field
        self._absolute_field_power = sm
        if not length <= m or length not in ZZ or length < 1:
            raise ValueError("length of the code must be a positive integer less than or equal to the absolute_field_power which is %d" % m )
        self._length = length
        
        if not dimension <= length or dimension not in ZZ or dimension < 1:
            raise ValueError("dimension of the code must be a positive integer less than or equal to the length of the code which is %d" % length )
        self._dimension = dimension

#       V = VectorSpace(relative_field, m) is a vector space of dimension m over Fq. ``basis`` is the basis of Fqm over Fp.
#       Fqm over Fq is a vector space of dimension m. ``linearly_independent_elements`` is a set of n elements belonging to Fqm.
#       ``FE.relative_field_representation(linearly_independent_elements[i])`` gives elements belonging to Fq.
#       And any n such elements are linearly independent over V. Thus, we generate g_0, g_1,...., g_{n-1}
        V = VectorSpace(relative_field, m)
        self._vector_space = V
        FE = RelativeFiniteFieldExtension(absolute_field, relative_field)
        self._relative_finite_field_extension = FE
        
        if linearly_independent_elements is None:
            basis = FE.absolute_field_basis()
            linearly_independent_elements = []
            for i in range(length):
                linearly_independent_elements.append(basis[i])
        else:
            if not len(linearly_independent_elements) == length:
                raise ValueError("The number of elements in the set should be equal to the length of the code")
            for i in range(length):
                if not linearly_independent_elements[i] in absolute_field:
                    raise ValueError("Element does not belong to absolute field")
            basis = []
            for i in range(len(linearly_independent_elements)):
                basis.append(FE.relative_field_representation(linearly_independent_elements[i]))
            if V.linear_dependence(basis):
                raise ValueError("The elements provided are not linearly independent")
        self._evaluation_points = linearly_independent_elements

    def _repr_(self):
        return "[%s, %s, %s] Linear Gabidulin Code over %s" \
                % (self.length(), self.dimension(), \
                self.minimum_distance(), self.base_field())

    def _latex_(self):
        return "[%s, %s, %s] \\textnormal{ Linear Gabidulin Code over } %s"\
                % (self.length(), self.dimension() ,self.minimum_distance(),
                self.base_field()._latex_())

    def __eq__(self, other):
        return isinstance(other, GabidulinCode) \
                and self.absolute_field() == other.absolute_field() \
                and self.relative_field() == other.relative_field() \
                and self.length() == other.length() \
                and self.dimension() == other.dimension() \
                and self.basis() == other.basis() \
                and self.evaluation_points() == self.evaluation_points()
                
    def minimum_distance(self):
        return self.length() - self.dimension() + 1
    
    def parity_evaluation_points(self):
        eval_pts = self.evaluation_points()
        n = self.length()
        k = self.dimension()
        coefficient_matrix = matrix(self.base_field(), n - 1, n, \
                lambda i,j: pow(eval_pts[j], pow(self.relative_field().order(), -n + k + 1 + i)))
        solution_space = coefficient_matrix.right_kernel()
        parity_eval_pts = solution_space.random_element()
        return parity_eval_pts
 
    def dual_code(self):
        parity_eval_pts = self.parity_evaluation_points()
        return GabidulinCode(self.absolute_field(), self.relative_field(), \
                self.length(), self.length() - self.dimension(), parity_eval_pts)

    def parity_check_matrix(self):
        E = GabidulinGeneratorMatrixEncoder(self.dual_code())
        return E.generator_matrix()

    def length(self):
        return self._length
    
    def dimension(self):
        return self._dimension

    def vector_space(self):
        return self._vector_space

    def relative_finite_field_extension(self):
        return self._relative_finite_field_extension

    def base_field(self):
        return self._absolute_field

    def prime_field(self):
        return self._prime_field

    def absolute_field(self):
        return self._absolute_field

    def relative_field(self):
        return self._relative_field

    def absolute_field_power(self):
        return self._absolute_field_power

    def relative_field_power(self):
        return self._relative_field_power

    def evaluation_points(self):
        return self._evaluation_points


####################### encoders ###############################


####################### encoders ###############################


class GabidulinPolynomialEvaluationEncoder(Encoder):

    def __init__(self, code):
        super(GabidulinPolynomialEvaluationEncoder, self).__init__(code)
        Frob = code.base_field().frobenius_endomorphism()
        self._sigma = Frob
        self._R = code.base_field()['x', Frob]

    def _repr_(self):
        return "Polynomial evaluation style encoder for %s" % self.code()

    def _latex_(self):
        return "\\textnormal{Polynomial evaluation style encoder for } %s" % self.code()._latex_()

    def __eq__(self, other):
        return isinstance(other, GabidulinPolynomialEvaluationEncoder) \
                and self.code() == other.code()

    def encode(self, p):
        M = self.message_space()
        if p not in M:
            raise ValueError("The value to encode must be in %s" % M)
        C = self.code()
        if p.degree() >= C.dimension():
            raise ValueError("The polynomial to encode must have degree at most %s" % (C.dimension() - 1))
        eval_pts = C.evaluation_points()
        sigma = self.sigma()
        codeword = []
        elem = C.absolute_field().zero()
        for i in range(len(eval_pts)):
            for j in range(p.degree()):
                elem += p.coefficients()[j]*(pow(sigma(eval_pts[i]), j))
            codeword.append(elem)
            elem = C.absolute_field().zero()
        return vector(codeword)

    def sigma(self):
        return self._sigma

    def message_space(self):
        return self._R


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
        return matrix(C.base_field(), C.dimension(), C.length(), lambda i,j: pow(eval_pts[j], pow(C.relative_field().order(), i)))
