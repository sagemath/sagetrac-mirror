###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class gfan_CitationItem( CitationItem ):
    _name = "gfan"

    ## provided by Gfan homepage
    _bibtex = r"""
@Misc{software-gfan,
    author = {Jensen, Anders N.},
    title = {{G}fan, a software system for {G}r{\"o}bner fans and tropical varieties},
    note = {Version~0.5},
    howpublished = {{\url{http://home.imf.au.dk/jensen/software/gfan/gfan.html}}}
} 
    """
