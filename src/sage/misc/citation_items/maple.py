from sage.misc.citation_items.citation_item import CitationItem

class Maple_CitationItem( CitationItem ):
    _name = "Maple"

    _re = [r"^sage.interfaces.maple"]
