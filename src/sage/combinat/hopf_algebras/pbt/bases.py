from sage.combinat.hopf_algebras.pbt import PlanarBinaryTreeFunctions as PBT
from sage.combinat.hopf_algebras import register_as_realization

# P-basis
from fundamental_basis import Fundamental
register_as_realization(PBT, Fundamental, "P")

# Q-basis
from fundamental_dual_basis import FundamentalDual
register_as_realization(PBT, FundamentalDual, "Q")

# H-basis
from complete_basis import Complete
register_as_realization(PBT, Complete, "H")

# E-basis
from elementary_basis import Elementary
register_as_realization(PBT, Elementary, "E")