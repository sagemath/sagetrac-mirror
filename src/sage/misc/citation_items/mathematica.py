from sage.misc.citation_items.citation_item import CitationItem

class Mathematica_CitationItem( CitationItem ):
    _name = "Mathematica"

    _re = [r"^sage.interfaces.mathematica"]
