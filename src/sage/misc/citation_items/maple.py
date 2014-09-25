###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class Maple_CitationItem( CitationItem ):
    _name = "Maple"

    ## based on citation by ACM
    _bibtex = r"""
@misc{software-maple,
    title = {{Maple}},
    organization = {{Maplesoft}},
    howpublished = {{\url{http://www.maplesoft.com}}},
} 
    """
