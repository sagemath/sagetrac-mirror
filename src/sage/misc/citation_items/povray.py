from sage.misc.citation_items.citation_item import CitationItem

class povray_CitationItem( CitationItem ):
    _name = "povray"

    _re = [r"^sage.interfaces.povray"]
