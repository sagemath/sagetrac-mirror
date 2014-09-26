###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class matlab_CitationItem( CitationItem ):
    _name = "matlab"

    _bibtex = r"""
@misc{software-matlab,
    title = {{MATLAB}},
    organization = {{The MathWorks, Inc.}},
    howpublished = {{\url{http://www.mathworks.com/products/matlab}}},
}
    """
