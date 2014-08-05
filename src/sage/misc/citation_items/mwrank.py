from sage.misc.citation_items.citation_item import CitationItem

class mwrank_CitationItem( CitationItem ):
    _name = "mwrank"

    _re = [r"^sage.interfaces.mwrank",
           "^sage.libs.mwrank"]
