from sage.misc.citation_items.citation_item import CitationItem

class R_CitationItem( CitationItem ):
    _name = "R"

    _re = [r"^sage.interfaces.r"]
