###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class R_CitationItem( CitationItem ):
    _name = "R"

    ## based on FAQ at http://www.r-project.org/
    _bibtex = r"""
@Manual{software-r,
  title        = {{R: A Language and Environment for Statistical Computing}},
  author       = {{R Core Team}},
  organization = {{R Foundation for Statistical Computing}},
  address      = {Vienna, Austria},
  year         = {2014},
  note         = {{Version~3.1.1, \url{http://www.R-project.org}}},
}
    """
