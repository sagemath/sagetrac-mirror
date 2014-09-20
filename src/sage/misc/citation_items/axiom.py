###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class Axiom_CitationItem( CitationItem ):
    _name = "Axiom"

    _bibtex = r"""
@Book{software-axiom,
    Author = {Jenks, Richard D. and Sutor, Robert S.},
    Title = {{Axiom. The scientific computation system.}},
    Pages = {xxiv + 742},
    Year = {1992},
    Publisher = {Springer-Verlag, Berlin},
}
    """
