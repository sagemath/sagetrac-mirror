from sage.misc.citation_items.citation_item import CitationItem

class PARI_CitationItem( CitationItem ):
    _name = "PARI"

    _re = [r"^sage.libs.pari",
           "^sage.interfaces.gp"]
