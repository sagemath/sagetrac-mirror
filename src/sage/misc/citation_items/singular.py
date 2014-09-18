###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class Singular_CitationItem( CitationItem ):
    _name = "Singular"

    _bibtex = r"""
@misc {software-singular,
    title = {{Singular -- A computer algebra system for polynomial computations}},
    author = {Decker, Wolfram and Greuel, Gert-Martin and Pfister, Gerhard and Sch\"onemann, Hans},
    note = {{Version~3.1.6}},
    howpublished = {{\url{http://www.singular.uni-kl.de}}},
}
    """
