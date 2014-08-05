from sage.misc.citation_items.citation_item import CitationItem

class MuPAD_CitationItem( CitationItem ):
    _name = "MuPAD"

    _re = [r"^sage.interfaces.mupad"]
