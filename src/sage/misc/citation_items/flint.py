###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class FLINT_CitationItem( CitationItem ):
    _name = "FLINT"

    ## based on FLINT homepage
    _bibtex = r"""
@misc{software-flint,
    author = {Hart, William and Johansson, Fredrik and Pancratz, Sebastian},
    title = {{FLINT: Fast Library for Number Theory}},
    year = {2014},
    note = {{Version~2.4.3}}
    howpublished = {{\url{http://flintlib.org}}}
}
    """
