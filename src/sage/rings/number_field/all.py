from __future__ import absolute_import

from .number_field import (NumberField, NumberFieldTower, CyclotomicField, QuadraticField,
                           is_fundamental_discriminant, is_real_place)
from .number_field_element import NumberFieldElement

from .order import EquationOrder, GaussianIntegers, EisensteinIntegers

from .totallyreal import enumerate_totallyreal_fields_prim
from .totallyreal_data import hermite_constant
from .totallyreal_rel import enumerate_totallyreal_fields_all, enumerate_totallyreal_fields_rel

from .unit_group import UnitGroup

from sage.misc.lazy_import import lazy_import
lazy_import('sage.rings.number_field.arithgroup_nf', 'arithgroup_nf')
