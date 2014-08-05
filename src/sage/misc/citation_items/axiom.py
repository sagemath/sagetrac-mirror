from sage.misc.citation_items.citation_item import CitationItem

class Axiom_CitationItem( CitationItem ):
    _name = "Axiom"

    _re = [r"^sage.interfaces.axiom"]
