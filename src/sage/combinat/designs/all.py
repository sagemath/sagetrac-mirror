"""
Combinatorial design features that are imported by default in the interpreter namespace
"""
from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import

lazy_import("sage.combinat.designs.ext_rep", ['designs_from_XML',
                                              'designs_from_XML_url'],
            deprecation=27066)

from .block_design import BlockDesign

from .incidence_structures import IncidenceStructure

from .incidence_structures import IncidenceStructure as Hypergraph

from .covering_design import (CoveringDesign,
                             schonheim,
                             trivial_covering_design)

from . import design_catalog as designs

del absolute_import
