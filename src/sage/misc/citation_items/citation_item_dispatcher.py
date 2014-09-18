r"""
Dispatching strings to citation items.

The citation framework, to speed things up, uses strings to register
ctiations.  In this module, we dispatch them to the corresponding
citation items.

AUTHOR:

- Martin Westerholt-Raum (martin@raum-brothers.eu)

EXAMPLES::

    sage: from sage.misc.citation.citation_item_dispatcher import get_citation_item
    sage: get_citation_item("singular")
    singular
"""

def get_citation_item(item_str):
    r"""
    Dispatch a string to a citation items.

    INPUT:

    - ``item_str`` -- A string.

    OUTPUT:

    A corresponding instance of class:`sage.misc.citation_items.citation_item.CitationItem`.

    EXAMPLES::

        sage: from sage.misc.citation_items.citation_item_dispatcher import get_citation_item
        sage: get_citation_item("ecm")
        ECM
    """
    return citation_items[item_str.lower_case()]()

import sage.misc.citation_items as citation_items_module
citation_items = {
    "axiom" : citation_items_module.axiom.Axiom_CitationItem,
    "ecm" : citation_items_module.ecm.ECM_CitationItem,
    "flint" : citation_items_module.flint.FLINT_CitationItem,
    "frobby" : citation_items_module.frobby.Frobby_CitationItem,
    "gap" : citation_items_module.gap.GAP_CitationItem,
    "gfan" : citation_items_module.gfan.gfan_CitationItem,
    "ginac" : citation_items_module.ginac.ginac_CitationItem,
    "givaro" : citation_items_module.givaro.Givaro_CitationItem,
    "gmp" : citation_items_module.gmp.GMP_CitationItem,
    "kash" : citation_items_module.kash.KASH_CitationItem,
    "lie" : citation_items_module.lie.LiE_CitationItem,
    "linbox" : citation_items_module.linbox.Linbox_CitationItem,
    "m4ri" : citation_items_module.m4ri.M4RI_CitationItem,
    "macaulay2" : citation_items_module.macaulay2.Macaulay2_CitationItem,
    "magma" : citation_items_module.magma.Magma_CitationItem,
    "maple" : citation_items_module.maple.Maple_CitationItem,
    "mathematica" : citation_items_module.mathematica.Mathematica_CitationItem,
    "matlab" : citation_items_module.matlab.matlab_CitationItem,
    "maxima" : citation_items_module.maxima.Maxima_CitationItem,
    "mpfi" : citation_items_module.mpfi.MPFI_CitationItem,
    "mpfr" : citation_items_module.mpfr.MPFR_CitationItem,
    "mupad" : citation_items_module.mupad.MuPAD_CitationItem,
    "mwrank" : citation_items_module.mwrank.mwrank_CitationItem,
    "ntl" : citation_items_module.ntl.NTL_CitationItem,
    "numpy" : citation_items_module.numpy.numpy_CitationItem,
    "octave" : citation_items_module.octave.Octave_CitationItem,
    "pari" : citation_items_module.pari.PARI_CitationItem,
    "polybori" : citation_items_module.polybori.PolyBoRi_CitationItem,
    "povray" : citation_items_module.povray.povray_CitationItem,
    "qsieve" : citation_items_module.qsieve.qsieve_CitationItem,
    "r" : citation_items_module.r.R_CitationItem,
    "scipy" : citation_items_module.scipy.scipy_CitationItem,
    "singular" : citation_items_module.singular.Singular_CitationItem,
    "symmetrica" : citation_items_module.symmetrica.Symmetrica_CitationItem,
    "tachyon" : citation_items_module.tachyon.Tachyon_CitationItem
}
