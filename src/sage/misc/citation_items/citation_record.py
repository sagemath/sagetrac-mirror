r"""
A record of citation items.

This is a list of citation items, which allows the use to access LaTeX
cite statements and BibTeX code.


AUTHORS:

- Martin Westerholt-Raum (martin@raum-brothers.eu)

"""

###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.structure.sage_object import SageObject

class CitationRecord( SageObject, list ):

    def bibtex(self):
        return "\n\n".join(item.bibtex() for item in self)

    def print_bibtex(self):
        print self.bibtex()

    def latex_citation(self):
        return "\\cite{{{}}}".format(", ".join(item._latex_citation_label()
                                               for item in self))

    def print_latex_citation(self):
        print "LaTeX \cite:  (print BibTeX code by calling method print_bibtex())"
        print self.latex_citation()

    def _repr_(self):
        return repr(list(self))

    def _latex_(self):
        return latex(list(self))
