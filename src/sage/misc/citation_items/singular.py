from sage.misc.citation_items.citation_item import CitationItem

class Singular_CitationItem( CitationItem ):
    _name = "Singular"

    _re = [r"^sage.interfaces.singular",
           "_libsingular",
           "^sage.libs.singular"]
