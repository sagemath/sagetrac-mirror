from permgroup_named import (SymmetricGroup, AlternatingGroup,
                       DihedralGroup, SplitMetacyclicGroup, SemidihedralGroup, CyclicPermutationGroup,
                       DiCyclicGroup, TransitiveGroup, PGL, PSL, PSp,PSU,PGU,
                       MathieuGroup, KleinFourGroup, QuaternionGroup,
                       PrimitiveGroup, PrimitiveGroups,
                       SuzukiGroup, TransitiveGroups, GeneralDihedralGroup)

from permgroup import  PermutationGroup, PermutationGroup_generic, PermutationGroup_subgroup, direct_product_permgroups

from permgroup_element import PermutationGroupElement

from permgroup_morphism import (PermutationGroupMorphism as PermutationGroupMap,
                                PermutationGroupMorphism_im_gens,
                                PermutationGroupMorphism_id)
PermutationGroupMorphism = PermutationGroupMorphism_im_gens

from cubegroup import CubeGroup, RubiksCube

from sage.misc.lazy_import import lazy_import
lazy_import("sage.groups.perm_gps", "permgroup", "pg", deprecation=18140)
