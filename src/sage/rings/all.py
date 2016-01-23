"""
Rings
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Ring base classes
from ring import Ring, Field
from commutative_ring import CommutativeRing
from integral_domain import IntegralDomain
from dedekind_domain import DedekindDomain
from principal_ideal_domain import PrincipalIdealDomain
from euclidean_domain import EuclideanDomain

# Ring element base classes
from sage.structure.element import (CommutativeAlgebraElement,
        RingElement, CommutativeRingElement, IntegralDomainElement,
        DedekindDomainElement, PrincipalIdealDomainElement,
        EuclideanDomainElement, FieldElement)

# Ideals
from ideal import Ideal
ideal = Ideal

# Quotient
from quotient_ring import QuotientRing

# Infinities
from infinity import infinity, Infinity, InfinityRing, unsigned_infinity, UnsignedInfinityRing

# Rational integers.
from integer_ring import IntegerRing, ZZ, crt_basis
from integer import Integer

# Rational numbers
from rational_field import RationalField, QQ
from rational import Rational
Rationals = RationalField

# Integers modulo n.
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing, Zmod
from sage.rings.finite_rings.integer_mod import IntegerMod, Mod, mod
Integers = IntegerModRing

# Finite fields
from finite_rings.all import *

# Number field
from number_field.all import *

# Function field
from function_field.all import *

# Finite residue fields
from finite_rings.residue_field import ResidueField

# p-adic field
from padics.all import *
from padics.padic_printing import _printer_defaults as padic_printing

# Semirings
from semirings.all import *

# Real numbers
from real_mpfr import (RealField, RR,
                       create_RealNumber as RealNumber)   # this is used by the preparser to wrap real literals -- very important.
Reals = RealField

from real_double import RealDoubleField, RDF, RealDoubleElement

from real_lazy import RealLazyField, RLF, ComplexLazyField, CLF

from sage.rings.real_arb import RealBallField, RBF

# Polynomial Rings and Polynomial Quotient Rings
from polynomial.all import *


# Algebraic numbers
from qqbar import (AlgebraicRealField, AA,
                   AlgebraicReal,
                   AlgebraicField, QQbar,
                   AlgebraicNumber,
                   number_field_elements_from_algebraics)
from universal_cyclotomic_field import UniversalCyclotomicField, E

# Intervals
from real_mpfi import (RealIntervalField,
                       RIF,
                       RealInterval)

# Complex numbers
from complex_field import ComplexField
from complex_number import (create_ComplexNumber as ComplexNumber)
Complexes = ComplexField
from complex_interval_field import ComplexIntervalField
from complex_interval import (create_ComplexIntervalFieldElement as ComplexIntervalFieldElement)

from complex_double import ComplexDoubleField, ComplexDoubleElement, CDF

from complex_mpc import MPComplexField

from sage.rings.complex_arb import ComplexBallField, CBF

# Power series rings
from power_series_ring import PowerSeriesRing
from power_series_ring_element import PowerSeries

# Laurent series ring in one variable
from laurent_series_ring import LaurentSeriesRing
from laurent_series_ring_element import LaurentSeries

# Pseudo-ring of PARI objects.
from pari_ring import PariRing, Pari

# Big-oh notation
from big_oh import O

# Fraction field
from fraction_field import FractionField
Frac = FractionField

# c-finite sequences
from cfinite_sequence import CFiniteSequence, CFiniteSequences

from bernoulli_mod_p import bernoulli_mod_p, bernoulli_mod_p_single

from monomials import monomials

CC = ComplexField()
CIF = ComplexIntervalField()

from misc import composite_field

from sage.misc.lazy_import import lazy_import
lazy_import('sage.rings.invariant_theory', 'invariant_theory')
lazy_import('sage.arith.all', '*', deprecation=19879)

from fast_arith import prime_range

# continued fractions
from sage.rings.continued_fraction import (farey, convergents,
  continued_fraction, continued_fraction_list,
   Hirzebruch_Jung_continued_fraction_list)
# and deprecated continued fractions
from sage.rings.contfrac import (CFF, ContinuedFractionField)

# asymptotic ring
from asymptotic.all import *

# Register classes in numbers abc
import numbers_abc
