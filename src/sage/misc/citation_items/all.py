from sage.misc.citation_items.axiom import Axiom_CitationItem
from sage.misc.citation_items.ecm import ECM_CitationItem
from sage.misc.citation_items.flint import FLINT_CitationItem
from sage.misc.citation_items.frobby import Frobby_CitationItem
from sage.misc.citation_items.gap import GAP_CitationItem
from sage.misc.citation_items.gmp import GMP_CitationItem
from sage.misc.citation_items.givaro import Givaro_CitationItem
from sage.misc.citation_items.kash import KASH_CitationItem
from sage.misc.citation_items.lie import LiE_CitationItem
from sage.misc.citation_items.linbox import Linbox_CitationItem
from sage.misc.citation_items.m4ri import M4RI_CitationItem
from sage.misc.citation_items.mpfi import MPFI_CitationItem
from sage.misc.citation_items.mpfr import MPFR_CitationItem
from sage.misc.citation_items.macaulay2 import Macaulay2_CitationItem
from sage.misc.citation_items.magma import Magma_CitationItem
from sage.misc.citation_items.maple import Maple_CitationItem
from sage.misc.citation_items.mathematica import Mathematica_CitationItem
from sage.misc.citation_items.maxima import Maxima_CitationItem
from sage.misc.citation_items.mupad import MuPAD_CitationItem
from sage.misc.citation_items.ntl import NTL_CitationItem
from sage.misc.citation_items.octave import Octave_CitationItem
from sage.misc.citation_items.pari import PARI_CitationItem
from sage.misc.citation_items.polybori import PolyBoRi_CitationItem
from sage.misc.citation_items.r import R_CitationItem
from sage.misc.citation_items.singular import Singular_CitationItem
from sage.misc.citation_items.symmetrica import Symmetrica_CitationItem
from sage.misc.citation_items.tachyon import Tachyon_CitationItem
from sage.misc.citation_items.gfan import gfan_CitationItem
from sage.misc.citation_items.ginac import ginac_CitationItem
from sage.misc.citation_items.matlab import matlab_CitationItem
from sage.misc.citation_items.mwrank import mwrank_CitationItem
from sage.misc.citation_items.numpy import numpy_CitationItem
from sage.misc.citation_items.povray import povray_CitationItem
from sage.misc.citation_items.qsieve import qsieve_CitationItem
from sage.misc.citation_items.scipy import scipy_CitationItem

citation_items = [
    Axiom_CitationItem(),
    ECM_CitationItem(),
    FLINT_CitationItem(),
    Frobby_CitationItem(),
    GAP_CitationItem(),
    gfan_CitationItem(),
    ginac_CitationItem(),
    GMP_CitationItem(),
    Givaro_CitationItem(),
    KASH_CitationItem(),
    LiE_CitationItem(),
    Linbox_CitationItem(),
    Macaulay2_CitationItem(),
    Magma_CitationItem(),
    Maple_CitationItem(),
    M4RI_CitationItem(),
    Mathematica_CitationItem(),
    matlab_CitationItem(),
    Maxima_CitationItem(),
    MPFI_CitationItem(),
    MPFR_CitationItem(),
    MuPAD_CitationItem(),
    mwrank_CitationItem(),
    NTL_CitationItem(),
    numpy_CitationItem(),
    Octave_CitationItem(),
    PARI_CitationItem(),
    PolyBoRi_CitationItem(),
    povray_CitationItem(),
    qsieve_CitationItem(),
    R_CitationItem(),
    scipy_CitationItem(),
    Singular_CitationItem(),
    Symmetrica_CitationItem(),
    Tachyon_CitationItem()
]

from citation_record import citation_record

__all__ = [citation_items, citation_record]
