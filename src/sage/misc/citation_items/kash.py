from sage.misc.citation_items.citation_item import CitationItem

class KASH_CitationItem( CitationItem ):
    _name = "KASH"

    _re = [r"^sage.interfaces.kash"]
