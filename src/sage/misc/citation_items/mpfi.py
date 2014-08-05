from sage.misc.citation_items.citation_item import CitationItem

class MPFI_CitationItem( CitationItem ):
    _name = "MPFI"

    _re = [r"^sage.rings.real_mpfi",
           "^sage.rings.complex_interval"]
