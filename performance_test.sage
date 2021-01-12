

# This file was *autogenerated* from the file performance_test.sage
from sage.all_cmdline import *   # import sage library

_sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_2 = Integer(2); _sage_const_4 = Integer(4); _sage_const_8 = Integer(8); _sage_const_31 = Integer(31); _sage_const_14 = Integer(14); _sage_const_20 = Integer(20); _sage_const_6 = Integer(6); _sage_const_17 = Integer(17); _sage_const_42 = Integer(42); _sage_const_71 = Integer(71); _sage_const_60 = Integer(60)
from sage.modules.fp_over_steenrod_algebra.fp_module import FP_Module

rels = [ [Sq(_sage_const_1 ),_sage_const_0 ,_sage_const_0 ,_sage_const_0 ], [Sq(_sage_const_2 ),_sage_const_0 ,_sage_const_0 ,_sage_const_0 ], [Sq(_sage_const_4 ),_sage_const_0 ,_sage_const_0 ,_sage_const_0 ], [Sq(_sage_const_8 ),_sage_const_0 ,_sage_const_0 ,_sage_const_0 ], [_sage_const_0 ,Sq(_sage_const_1 ),_sage_const_0 ,_sage_const_0 ], [_sage_const_0 ,Sq(_sage_const_2 ),_sage_const_0 ,_sage_const_0 ], [    _sage_const_0 ,Sq(_sage_const_4 ),_sage_const_0 ,_sage_const_0 ], [Sq(_sage_const_31 ),Sq(_sage_const_14 ),_sage_const_0 ,_sage_const_0 ], [_sage_const_0 ,Sq(_sage_const_20 ),_sage_const_0 ,_sage_const_0 ], [_sage_const_0 ,_sage_const_0 ,Sq(_sage_const_1 ),_sage_const_0 ], [_sage_const_0 ,_sage_const_0 ,Sq(_sage_const_2 ),_sage_const_0 ], [_sage_const_0 ,Sq(_sage_const_31 ),Sq(_sage_const_6 ),_sage_const_0 ],     [_sage_const_0 ,_sage_const_0 ,Sq(_sage_const_8 ),_sage_const_0 ], [_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,Sq(_sage_const_1 )], [_sage_const_0 ,_sage_const_0 ,Sq(_sage_const_31 ),Sq(_sage_const_2 )], [_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,Sq(_sage_const_4 )], [_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,Sq(_sage_const_8 )] ]

M = FP_Module([_sage_const_0 , _sage_const_17 , _sage_const_42 , _sage_const_71 ], SteenrodAlgebra(_sage_const_2 ), relations=rels)

# The performance will be measured when fp_morphism._resolve_kernel() in internal degrees >= 20:
M.resolution(_sage_const_2 , top_dim=_sage_const_60 )

