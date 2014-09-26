###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class Octave_CitationItem( CitationItem ):
    _name = "Octave"

    _bibtex = r"""
 @misc{software-numpy,
    author = {{Eaton, John W. and others}}
    title        = {{GNU Octave}}
    howpublished = {{\url{http://www.gnu.org/software/octave}}},
}
    """
