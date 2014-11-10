"""
Iwahori-Hecke algebras and their representations
"""
from sage.misc.lazy_import import lazy_import

lazy_import('sage.algebras.iwahori_hecke_algebras.iwahori_hecke_algebra', 'IwahoriHeckeAlgebra')

lazy_import('sage.algebras.iwahori_hecke_algebras.iwahori_hecke_algebra_representations',
            'LeftCellRepresentationOfHeckeAlgebra')
lazy_import('sage.algebras.iwahori_hecke_algebras.iwahori_hecke_algebra_representations',
            'RightCellRepresentationOfHeckeAlgebra')

lazy_import('sage.algebras.iwahori_hecke_algebras.iwahori_hecke_algebra_representations',
            'SpechtModuleWithMurphyBasis')
lazy_import('sage.algebras.iwahori_hecke_algebras.iwahori_hecke_algebra_representations',
            'SeminormalRepresentation')
lazy_import('sage.algebras.iwahori_hecke_algebras.iwahori_hecke_algebra_representations',
            'SeminormalRepresentation_Orthogonal')
lazy_import('sage.algebras.iwahori_hecke_algebras.iwahori_hecke_algebra_representations',
            'SeminormalRepresentation_Alternating')

