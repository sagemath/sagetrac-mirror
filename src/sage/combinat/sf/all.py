r"""
Symmetric function features that are imported by default in the interpreter namespace
"""
from sage.misc.lazy_import import lazy_import

# In the long run, this will be the single entry point
# Nothing else will be exported
from sf import SymmetricFunctions

# This is deprecated:
lazy_import('sage.combinat.sf.sfa', ['SymmetricFunctionAlgebra'])

# Advanced stuff:

lazy_import('sage.combinat.sf.kfpoly', 'KostkaFoulkesPolynomial')

from ns_macdonald import NonattackingFillings, AugmentedLatticeDiagramFilling, LatticeDiagram
