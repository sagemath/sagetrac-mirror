from sage.misc.lazy_import import lazy_import
lazy_import('sage.algebras.vertex_algebras.vertex_algebra', 'VertexAlgebra')
lazy_import('sage.algebras.vertex_algebras.examples',
    ('VirasoroVertexAlgebra',
    'AffineVertexAlgebra',
    'FreeBosonsVertexAlgebra',
    'FreeFermionsVertexAlgebra',
    'WeylVertexAlgebra',
    'NeveuSchwarzVertexAlgebra',
    'N2VertexAlgebra',
    'AbelianVertexAlgebra'))
