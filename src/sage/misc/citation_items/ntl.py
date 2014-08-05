from sage.misc.citation_items.citation_item import CitationItem

class NTL_CitationItem( CitationItem ):
    _name = "NTL"

    _re = [r"^sage.libs.ntl",
           "^sage.rings.finite_rings.element_ntl_gf2e"]
