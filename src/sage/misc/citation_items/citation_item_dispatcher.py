r"""
Dispatching strings to citation items.

The citation framework, to speed things up, uses strings to register
ctiations.  In this module, we dispatch them to the corresponding
citation items.

AUTHOR:

- Martin Westerholt-Raum (martin@raum-brothers.eu)

EXAMPLES::

    sage: from sage.misc.citation_items.citation_item_dispatcher import get_citation_item
    sage: get_citation_item("singular")
    Singular
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
    return citation_items_by_string[item_str.lower()]()

import sage.misc.citation_items.all as citation_items
citation_items_by_string = {
    "axiom" : citation_items.Axiom_CitationItem,
    "ecm" : citation_items.ECM_CitationItem,
    "flint" : citation_items.FLINT_CitationItem,
    "frobby" : citation_items.Frobby_CitationItem,
    "gap" : citation_items.GAP_CitationItem,
    "gfan" : citation_items.gfan_CitationItem,
    "ginac" : citation_items.ginac_CitationItem,
    "givaro" : citation_items.Givaro_CitationItem,
    "gmp" : citation_items.GMP_CitationItem,
    "kash" : citation_items.KASH_CitationItem,
    "lie" : citation_items.LiE_CitationItem,
    "linbox" : citation_items.Linbox_CitationItem,
    "m4ri" : citation_items.M4RI_CitationItem,
    "macaulay2" : citation_items.Macaulay2_CitationItem,
    "magma" : citation_items.Magma_CitationItem,
    "maple" : citation_items.Maple_CitationItem,
    "mathematica" : citation_items.Mathematica_CitationItem,
    "matlab" : citation_items.matlab_CitationItem,
    "maxima" : citation_items.Maxima_CitationItem,
    "mpfi" : citation_items.MPFI_CitationItem,
    "mpfr" : citation_items.MPFR_CitationItem,
    "mupad" : citation_items.MuPAD_CitationItem,
    "mwrank" : citation_items.mwrank_CitationItem,
    "ntl" : citation_items.NTL_CitationItem,
    "numpy" : citation_items.numpy_CitationItem,
    "octave" : citation_items.Octave_CitationItem,
    "pari" : citation_items.PARI_CitationItem,
    "polybori" : citation_items.PolyBoRi_CitationItem,
    "povray" : citation_items.povray_CitationItem,
    "qsieve" : citation_items.qsieve_CitationItem,
    "r" : citation_items.R_CitationItem,
    "scipy" : citation_items.scipy_CitationItem,
    "singular" : citation_items.Singular_CitationItem,
    "symmetrica" : citation_items.Symmetrica_CitationItem,
    "tachyon" : citation_items.Tachyon_CitationItem
}
