dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.

dnl Include all from sagemath-objects (via m4 include)
include(`../sagemath_objects/src/MANIFEST.in')

# Extra in sagemath-categories:
global-include all__sagemath_categories.py
graft sage/categories
include sage/misc/prandom.*              # dep of sage/rings/ring
include sage/rings/ideal.*
include sage/rings/ring.*
graft sage/typeset                     # dep of sage.categories.tensor
# include sage/rings/integer*.*     # depends on cypari, flint - https://trac.sagemath.org/ticket/30022
# include sage/rings/rational*.*
# include sage/rings/infinity.*
