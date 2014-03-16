"""
Internal test script

This is not meant to appear in the final Sage release.
"""

from sage.all import *

from diamond_cutting import *

import time

"""L = Lattice([[2, 0], [0, 3]])
g = plot_lattice(L)
g.show()

L = special_lattice('BodyCenteredCubic')
print L.embedded_basis()
plot_lattice(L).show(viewer='tachyon')
"""

L = special_lattice('SimpleCubic')
g = plot_lattice(L)
V = L.voronoi_cell()
g += V.plot()
g.show(viewer='tachyon')
g.save('lattice.png')

#L = Lattice([[GF(3)(1), 0], [0, 1]])
#print L
#print L.zero() + L.an_element()

#L_complex = Lattice([[i, 0], [0, 1]])

#L = Lattice([[2, 0], [0, 1]])
#print L([2, 1])

#L = Lattice([[6, 1], [9, 0]])
#print L

#L.voronoi_cell()

"""L = Lattice([[1, 2, 0], [4, 0, 6], [0, 8, 9]])
#V = L.voronoi_cell()
V = calculate_voronoi_cell(L.embedded_basis_matrix(), radius=None, debug=True)
print V
print V.Vrepresentation()
print L.embedded_basis_matrix()
print L.gram_matrix()
print L.discriminant()
print V.volume()
print V.lrs_volume()

L = Lattice(matrix([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]]))
#V = L.voronoi_cell(radius=100)
V = calculate_voronoi_cell(L.embedded_basis_matrix(), radius=847, debug=True)
print V.volume()
print V.lrs_volume()
print L.discriminant()
"""

cr63 = """
4.47    0.00    0.00    0.00    0.00    0.00    
-3.35    2.96    0.00    0.00    0.00    0.00    
1.34    -3.55    2.37    0.00    0.00    0.00    
-0.22    1.77    -3.55    2.05    0.00    0.00    
-0.22    -0.59    1.77    -3.76    1.53    0.00    
1.34    1.18    0.59    3.07    -2.29    1.32
"""
cr6 = """
 2.4494897427832E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00  
-1.6329931618555E+00   1.8257418583506E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00  
 4.0824829046386E-01  -1.8257418583506E+00   1.5811388300842E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00  
 0.0000000000000E+00   5.4772255750517E-01  -1.8973665961010E+00   1.4491376746189E+00   0.0000000000000E+00   0.0000000000000E+00  
 0.0000000000000E+00   0.0000000000000E+00   6.3245553203368E-01  -1.9321835661586E+00   1.3662601021279E+00   0.0000000000000E+00  
 4.0824829046386E-01   3.6514837167011E-01   3.1622776601684E-01   9.6609178307930E-01  -1.7078251276599E+00   1.3228756555323E+00  
"""
m2a = """
2 0
0 1
"""
m2b = """
-2 -2
2 -2
"""
es6 = """
1.41    0.00    0.00    0.00    0.00    0.00    
-0.71    1.22    0.00    0.00    0.00    0.00    
0.00    -0.82    1.15    0.00    0.00    0.00    
0.00    0.00    -0.87    1.12    0.00    0.00    
0.00    -0.82    -0.58    -0.45    0.37    0.00    
0.00    0.00    -0.87    -0.67    0.55    0.71
"""
lv6 = """
  1.4135770e+000  0.0000000e+000  0.0000000e+000  0.0000000e+000  0.0000000e+000  0.0000000e+000
  3.7281308e-001  1.3635287e+000  0.0000000e+000  0.0000000e+000  0.0000000e+000  0.0000000e+000
 -2.9499631e-001 -2.2516690e-001  1.4037368e+000  0.0000000e+000  0.0000000e+000  0.0000000e+000
 -3.7281308e-001 -2.8456344e-001 -8.7491668e-001  1.0063572e+000  0.0000000e+000  0.0000000e+000
  3.7281308e-001  2.8456344e-001 -1.7307180e-001 -4.5556135e-001  1.2412671e+000  0.0000000e+000
 -7.4569690e-001 -5.6918088e-001  3.4619130e-001 -5.5058489e-001 -6.4855953e-001  6.2011907e-001
 """
matas6 = """
 2.4494897427832E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00  
-4.0824829046386E-01   2.4152294576982E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00  
-4.0824829046386E-01  -4.8304589153965E-01   2.3664319132398E+00   0.0000000000000E+00   0.0000000000000E+00   0.0000000000000E+00  
-4.0824829046386E-01  -4.8304589153965E-01  -5.9160797830996E-01   2.2912878474779E+00   0.0000000000000E+00   0.0000000000000E+00  
-4.0824829046386E-01  -4.8304589153965E-01  -5.9160797830996E-01  -7.6376261582597E-01   2.1602468994693E+00   0.0000000000000E+00  
-4.0824829046386E-01  -4.8304589153965E-01  -5.9160797830996E-01  -7.6376261582597E-01  -1.0801234497346E+00   1.8708286933870E+00  
"""

def volume_simplex(S):
    n = len(S)
    V = matrix([vector(list(S[k] - S[0])) for k in range(1, n)]).transpose()
    return abs(V.determinant()) / factorial(n - 1)

def volume(P):
    from scipy.spatial import Delaunay
    if hasattr(P, 'vertices'):
        V = P.vertices()
    else:
        V = P
    D = Delaunay(V)
    vol = 0
    for simplex in D.vertices:
        vol += volume_simplex([D.points[index] for index in simplex])
    return vol

def test_lattice(data):
    M = []
    for line in data.splitlines():
        row = []
        for value in line.split(' '):
            value = value.strip()
            if value:
                value = QQ(float(value))
                row.append(value)
        if row:
            M.append(row)
    M = matrix(M)
    start = time.clock()
    V = calculate_voronoi_cell(M, radius=10, debug=True)
    stop = time.clock()
    print V.Vrepresentation()
    print V.Hrepresentation()
    print "Computed Voronoi cell in %f seconds" % (stop - start)
    print "Volume: %f" % volume(V)
    L = Lattice(M)
    print "Det:    %f" % L.determinant()
    
#test_lattice(cr63)

#test_lattice("""
#2 543 2 123
#432 425 123 5
#43 21 8 0
#45 7 0 1
#""")

m3 = """
1 2 0
0 3 4
3 1 0
"""

#test_lattice(m3)
    
#test_lattice(m2b)
#test_lattice(cr6)
#test_lattice(lv6)
#test_lattice(matas6)

#test_lattice("1 0 0 0 0\n0 2 1 0 0")

#L = Lattice([[1, 0], [2, 0], [0, 2]])
#print L.voronoi_cell().Vrepresentation()

#L = Lattice(quadratic_form=[[2,0], [0,2]])

#L = Lattice([[-1, -1, 0], [1, -1, 0], [0, 1, -1]])
#V = L.voronoi_cell()
#print V.lrs_volume()

"""L = Lattice(quadratic_form=[[1,0], [0,1]])
print L.closest_vector([0.4, 0.2])
print L.closest_vector([1.5, 2.5])
print L.closest_vector([1, 1])
print L.closest_vector([6, 0.6])

L = Lattice([[1,0,0], [0,1,0]])
print L.closest_vector([4, 2.6, 3.6])"""

"""print Lattice(basis=[[1, 0], [0, 1]])
print Lattice(basis=[[3, 2, 0], [1, 0, 1]])
print Lattice(basis=[[ZZ(1) / ZZ(2), 0], [0, 1]])
print Lattice(basis=[[0.5, 0], [0, 1]])
print Lattice(basis=[[1, 0], [0, 1]], base_ring=Zp(2))

L = Lattice([[6, 1], [9, 0]])
print L
print L.gram_matrix()
print L.determinant()
print L.discriminant()

L = Lattice([[2, 0], [0, 3]])
print L.shortest_vectors()

leech_lattice = [
    [8, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 2, 4, 2, 2, 2, 4, 2, 2, 2, 0, 0, 0, -3],
    [4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 1, 2, 1, 0, 0, -1],
    [4, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 2, 2, 1, 1, 1, 0, 0, -1],
    [4, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 0, 0, -1],
    [4, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, 0, -1],
    [4, 2, 2, 2, 2, 4, 2, 2, 2, 2, 2, 1, 2, 2, 1, 1, 2, 1, 2, 1, 0, 0, 0, -1],
    [4, 2, 2, 2, 2, 2, 4, 2, 2, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 0, 0, 0, -1],
    [2, 2, 2, 2, 2, 2, 2, 4, 1, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 0, 0, 1],
    [4, 2, 2, 2, 2, 2, 2, 1, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, -1],
    [4, 2, 2, 2, 2, 2, 2, 1, 2, 4, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 0, 1, 0, -1],
    [4, 2, 2, 2, 2, 2, 2, 1, 2, 2, 4, 2, 2, 1, 2, 1, 2, 1, 2, 1, 0, 0, 1, -1],
    [2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 4, 1, 2, 2, 2, 1, 2, 2, 2, 2, 1, 1, 1],
    [4, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 4, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, -1],
    [2, 2, 1, 1, 2, 2, 1, 2, 2, 2, 1, 2, 2, 4, 2, 2, 1, 2, 2, 2, 2, 2, 1, 1],
    [2, 1, 2, 1, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 4, 2, 1, 2, 2, 2, 2, 1, 2, 1],
    [2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 4, 1, 2, 2, 2, 2, 1, 1, 1],
    [4, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 4, 2, 2, 2, 1, 1, 1, -1],
    [2, 1, 2, 1, 2, 1, 1, 2, 2, 2, 1, 2, 1, 2, 2, 2, 2, 4, 2, 2, 2, 2, 1, 1],
    [2, 1, 1, 2, 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 4, 2, 2, 1, 2, 1],
    [2, 2, 1, 1, 2, 1, 2, 2, 2, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 4, 2, 1, 1, 1],
    [0, 1, 1, 1, 1, 0, 0, 2, 1, 0, 0, 2, 1, 2, 2, 2, 1, 2, 2, 2, 4, 2, 2, 2],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 4, 2, 2],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 2, 4, 2],
    [-3, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 2, 2, 2, 4]
]
L = Lattice(inner_product_matrix=leech_lattice)
print L.shortest_vectors_count()

L = Lattice([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
print L.shortest_vectors()
"""

"""
#GM = matrix([[0, 3], [3, -1]])
L = Lattice([[1, 0], [0, 1]])
GM = L.basis()
V = calculate_voronoi_cell(matrix(GM))
print V.Hrepresentation()
print V.Vrepresentation()
"""
