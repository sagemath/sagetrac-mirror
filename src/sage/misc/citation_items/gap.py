from sage.misc.citation_items.citation_item import CitationItem

class GAP_CitationItem( CitationItem ):
    _name = "GAP"

    _re = [r"^sage.interfaces.gap"]
