from sage.misc.citation_items.citation_item import CitationItem

class MPFR_CitationItem( CitationItem ):
    _name = "MPFR"

    _re = [r"^sage.rings.real_mpfr",
           "^sage.rings.complex_number"]
