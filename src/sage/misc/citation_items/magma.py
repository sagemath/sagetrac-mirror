from sage.misc.citation_items.citation_item import CitationItem

class Magma_CitationItem( CitationItem ):
    _name = "Magma"

    _re = [r"^sage.interfaces.magma",
           "^sage.interfaces.magma_free"]
