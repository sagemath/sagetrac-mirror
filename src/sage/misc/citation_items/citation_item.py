from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

import re

class CitationItem( SageObject, UniqueRepresentation ):
    def re(self):
        return [re.compile(s) for s in self._re]

    def _repr_(self):
        return self._name

    def _latex_(self):
        return r"\text{{{}}}".format(self._name)
