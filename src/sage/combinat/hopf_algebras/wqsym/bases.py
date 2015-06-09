from sage.combinat.hopf_algebras import register_as_realization
from sage.combinat.hopf_algebras.wqsym import WordQuasiSymmetricFunctions as WQSym

# S-basis
from fundamental_basis import Fundamental
register_as_realization(WQSym, Fundamental, "S")

# M-basis
from fundamental_dual_basis import FundamentalDual
register_as_realization(WQSym, FundamentalDual, "M")

# Me-basis
from elementary_basis import Elementary
register_as_realization(WQSym, Elementary, "Me")

# Ms-basis
from homogene_basis import Homogene
register_as_realization(WQSym, Homogene, "Ms")
