from sage.combinat.hopf_algebras import register_as_realization
from sage.combinat.hopf_algebras.fqsym import FreeQuasiSymmetricFunctions

# F-basis
from fundamental_basis import Fundamental
register_as_realization(FreeQuasiSymmetricFunctions, Fundamental, "F")

# G-basis
from fundamental_dual_basis import FundamentalDual
register_as_realization(FreeQuasiSymmetricFunctions, FundamentalDual, "G")

# E-basis
from elementary_basis import Elementary
register_as_realization(FreeQuasiSymmetricFunctions, Elementary, "E")

# n-basis
from elementary_dual_basis import ElementaryDual
register_as_realization(FreeQuasiSymmetricFunctions, ElementaryDual, "n")

# H-basis
from homogene_basis import Homogene
register_as_realization(FreeQuasiSymmetricFunctions, Homogene, "H")

# m-basis
from homogene_dual_basis import HomogeneDual
register_as_realization(FreeQuasiSymmetricFunctions, HomogeneDual, "m")
