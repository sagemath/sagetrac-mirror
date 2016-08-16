from sage.coding.linear_code import AbstractLinearCode
from sage.rings.integer import Integer
from sage.coding.relative_finite_field_extension import *
from sage.matrix.constructor import matrix
from sage.matrix.matrix import is_Matrix
from sage.structure.element import is_Vector

def rank_weight(C, c):
    if is_Matrix(c):
        return c.rank()
    elif is_Vector(c):
        c = to_matrix_representation(C, c)
    return c.rank()

def rank_distance(C, a, b):
    if is_Vector(a):
        a = to_matrix_representation(C, a)
    if is_Vector(b):
        b = to_matrix_representation(C, b)
    return (a - b).rank()

def to_matrix_representation(C, v):
    if not is_Vector(v):
        raise TypeError("Input must be a vector")
    n = v.length()
    FE = C.field_extension()
    Fqm = FE.absolute_field()
    Fq = FE.relative_field()
    m = matrix(Fq, C.m(), n, lambda i,j: FE.relative_field_representation(v[j])[i])
    return m

def from_matrix_representation(C, m):
    if not is_Matrix(m):
        raise TypeError("Input must be a matrix")
    FE = C.field_extension()
    Fqm = FE.absolute_field()
    v = []
    for i in range(m.ncols()):
        v.append(FE.absolute_field_representation(m.column(i)))
    return vector(v)


class AbstractRankMetricCode(AbstractLinearCode):

    _registered_encoders = {} 
    _registered_decoders = {} 

#    def __init__(self, base_field, sub_field, length, dimension, field_extension=None, default_encoder_name, default_decoder_name):
    def __init__(self, base_field, sub_field, length, dimension, default_encoder_name, default_decoder_name, field_extension = None):

        super(AbstractRankMetricCode, self).__init__(base_field, length, default_encoder_name, default_decoder_name)
        if not isinstance(dimension, (int, Integer)):
            raise ValueError("dimension must be a Python int or a Sage Integer")
        if not sub_field.is_field():
            raise ValueError("'sub_field' must be a field (and {} is not one)".format(base_field))
        if not field_extension: #if field_extension is provided, then what? how to check?
            field_extension = RelativeFiniteFieldExtension(base_field, sub_field)
        self._base_field = base_field
        self._sub_field = sub_field
        self._length = length
        self._dimension = dimension
        self._field_extension = field_extension

    def base_field(self):
        return self._base_field

    def sub_field(self):
        return self._sub_field

    def field_extension(self):
        return self._field_extension

    @property
    def automorphism_group_gens(self):
        raise AttributeError("%s has no attribute 'automorphism_group_gens'" % self.__class__)
 
    @property
    def canonical_representative(self):
        raise AttributeError("%s has no attribute 'canonical_representative'" % self.__class__)

    @property
    def permutation_automorphism_group(self):
        raise AttributeError("%s has no attribute 'permutation_automorphism_group'" % self.__class__)

    @property
    def assmuss_mattson_designs(self):
        raise AttributeError("%s has no attribute 'assmus_mattson_designs'" % self.__class__)

    @property
    def binomial_moment(self):
        raise AttributeError("%s has no attribute 'binomial_moment'" % self.__class__)

  
