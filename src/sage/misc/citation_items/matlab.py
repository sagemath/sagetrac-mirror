from sage.misc.citation_items.citation_item import CitationItem

class matlab_CitationItem( CitationItem ):
    _name = "matlab"

    _re = [r"^sage.interfaces.matlab"]
