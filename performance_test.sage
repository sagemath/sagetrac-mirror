from sage.modules.fp_over_steenrod_algebra.fp_module import FP_Module

rels = [ [Sq(1),0,0,0], [Sq(2),0,0,0], [Sq(4),0,0,0], [Sq(8),0,0,0], [0,Sq(1),0,0], [0,Sq(2),0,0], [    0,Sq(4),0,0], [Sq(31),Sq(14),0,0], [0,Sq(20),0,0], [0,0,Sq(1),0], [0,0,Sq(2),0], [0,Sq(31),Sq(6),0],     [0,0,Sq(8),0], [0,0,0,Sq(1)], [0,0,Sq(31),Sq(2)], [0,0,0,Sq(4)], [0,0,0,Sq(8)] ]

M = FP_Module([0, 17, 42, 71], SteenrodAlgebra(2), relations=rels)

# The performance will be measured when fp_morphism._resolve_kernel() in internal degrees >= 20:
M.resolution(2, top_dim=60)
