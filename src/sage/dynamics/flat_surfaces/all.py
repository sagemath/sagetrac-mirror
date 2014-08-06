from sage.misc.lazy_import import lazy_import
lazy_import("sage.dynamics.flat_surfaces.abelian_strata",
            ["AbelianStrata", "AbelianStratum"])
lazy_import("sage.dynamics.flat_surfaces.quadratic_strata",
            ["QuadraticStrata", "QuadraticStratum"])
lazy_import("sage.dynamics.flat_surfaces.homology",
            ["RibbonGraph","RibbonGraphWithAngles"])
lazy_import("sage.dynamics.flat_surfaces.separatrix_diagram",
            ["SeparatrixDiagram","CylinderDiagram"])

from origamis.all import *
