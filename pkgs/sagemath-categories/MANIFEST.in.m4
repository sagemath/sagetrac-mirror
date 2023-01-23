dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

global-include all__sagemath_categories.py
graft sage/categories
# Exclude what is already shipped in sagemath-objects
exclude sage/categories/action.*
exclude sage/categories/algebra_functor.*
exclude sage/categories/basic.*
exclude sage/categories/cartesian_product.*
exclude sage/categories/category*.*
exclude sage/categories/covariant_functorial_construction.*
exclude sage/categories/facade_sets.*
exclude sage/categories/functor.*
exclude sage/categories/homset.*
exclude sage/categories/homsets.*
exclude sage/categories/map.*
exclude sage/categories/morphism.*
exclude sage/categories/isomorphic_objects.*
exclude sage/categories/objects.*
exclude sage/categories/primer.*
exclude sage/categories/pushout.*
exclude sage/categories/quotients.*
exclude sage/categories/realizations.*
exclude sage/categories/sets_cat.*
exclude sage/categories/sets_with_partial_maps.*
exclude sage/categories/subobjects.*
exclude sage/categories/subquotients.*
exclude sage/categories/with_realizations.*

include sage/interfaces/abc.p*


include sage/rings/ideal.*
include sage/rings/ring.*
graft sage/typeset                     # dep of sage.categories.tensor

include sage/rings/quotient_ring*.p*
include sage/rings/homset.p*
include sage/rings/ideal*.p*
include sage/rings/localization.p*
include sage/rings/morphism.p*

include sage/rings/abc.*
include sage/rings/integer*.*
include sage/rings/rational*.*
include sage/rings/infinity.*

include sage/arith/*.*
include sage/misc/allocator.*
include sage/misc/functional.p*
include sage/misc/latex*.*
include sage/misc/html.p*
include sage/misc/table.p*
include sage/misc/map_threaded.p*
include sage/misc/binary_tree.p*
include sage/misc/derivative.p*
include sage/misc/mrange.p*

include sage/rings/finite_rings/element_base.*
include sage/rings/finite_rings/stdint.*
include sage/rings/finite_rings/finite_field_base.p*
include sage/rings/finite_rings/finite_field_constructor.py
include sage/rings/fast_arith.*
include sage/rings/finite_rings/integer_mod_limits.h
include sage/rings/finite_rings/integer_mod.pxd   # .pyx depends on pari

include sage/rings/number_field/number_field_base.*

include sage/rings/real_double.p*

include sage/rings/fraction_field.p*
include sage/rings/fraction_field_element.p*
include sage/rings/qqbar_decorators.p*

include sage/rings/padics/padic_generic.p*
include sage/rings/padics/local_generic.p*
include sage/rings/padics/local_generic_element.p*
include sage/rings/padics/precision_error.p*
include sage/rings/padics/misc.p*

include sage/rings/polynomial/polynomial_ring.p*
include sage/rings/polynomial/polynomial_ring_constructor.p*
include sage/rings/polynomial/polynomial_quotient_ring.p*
include sage/rings/polynomial/polynomial_singular_interface.p*
include sage/rings/polynomial/multi_polynomial_ring.p*
include sage/rings/polynomial/multi_polynomial_ring_base.p*
include sage/rings/polynomial/multi_polynomial_sequence.p*
include sage/rings/polynomial/multi_polynomial_ideal.p*
include sage/rings/polynomial/infinite_polynomial_ring.p*
include sage/rings/polynomial/commutative_polynomial.p*
include sage/rings/polynomial/polynomial_compiled.p*
include sage/rings/polynomial/polynomial_element.p*
include sage/rings/polynomial/polynomial_element_generic.p*
include sage/rings/polynomial/polynomial_fateman.p*
include sage/rings/polynomial/polynomial_quotient_ring_element.p*
include sage/rings/polynomial/multi_polynomial.p*
include sage/rings/polynomial/multi_polynomial_element.p*
include sage/rings/polynomial/polydict.p*
include sage/rings/polynomial/term_order.p*

global-exclude *.c
global-exclude *.cpp

global-exclude __pycache__
global-exclude *.py[co]
global-exclude *.bak
global-exclude *.so
global-exclude *~
prune .tox
prune build
prune dist
