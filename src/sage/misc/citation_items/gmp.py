from sage.misc.citation_items.citation_item import CitationItem

class GMP_CitationItem( CitationItem ):
    _name = "GMP"

    _re = [r"^sage.rings.integer.Integer"]
