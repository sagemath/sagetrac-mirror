#Format
#Lattices = [[n,G],...,[n',G']]

def xid(x,n):
	return x*Matrix.identity(QQ,n)
	
Lattices = []

U = Matrix(QQ,2,2,[0,1,1,0])
A2 = Matrix(QQ,2,2,[2,-1,-1,2])
D4 = Matrix(QQ,4,4,[2,0,0,-1,0,2,0,-1,0,0,2,-1,-1,-1,-1,2])
E8 = Matrix(QQ,8,8,[2,0,0,-1,0,0,0,0,0,2,-1,0,0,0,0,0,0,-1,2,-1,0,0,0,0,-1,0,-1,2,-1,0,0,0,0,0,0,-1,2,-1,0,0,0,0,0,0,-1,2,-1,0,0,0,0,0,0,-1,2,-1,0,0,0,0
,0,0,-1,2])

# n = 1---------------
Lattices.append([1,Matrix.block_diagonal([U,U,U,-2*E8])])

# n = 2---------------
Lattices.append([2,Matrix.block_diagonal([U,3*U,3*U,-1*A2,-1*A2])])
# n = 3---------------
Lattices.append([3,Matrix.block_diagonal([U,2*U,2*U,-2*D4])])
# n = 4---------------
Lattices.append([4,Matrix.block_diagonal([U,4*U,4*U,xid(-2,2)])])
# n = 6---------------
Lattices.append([6,Matrix.block_diagonal([3*U,2*A2,-A2,-A2])])
# n = 9---------------
Lattices.append([9,Matrix.block_diagonal([2*U,2*U,2*U,xid(-4,2)])])
# n = 10---------------
Lattices.append([10,Matrix.block_diagonal([U,xid(4,2),xid(-4,3)])])
# n = 12---------------
M = Matrix(QQ,3,3,[6,2,2,2,6,-2,2,-2,6])
Lattices.append([12,Matrix.block_diagonal([M,xid(-2,2)])])
# n = 16---------------
Lattices.append([16,Matrix.block_diagonal([U,5*U,5*U])])
# n = 17---------------
Lattices.append([17,Matrix.block_diagonal([U,2*A2,-4*A2])])
# n = 18---------------
Lattices.append([18,Matrix.block_diagonal([U,6*U,6*U])])
# n = 21---------------
Lattices.append([21,Matrix.block_diagonal([2*U,2*U,2*U,xid(-8,1)])])
# n = 22---------------
Lattices.append([22,Matrix.block_diagonal([2*U,xid(4,2),xid(-4,2)])])
# n = 26---------------
Lattices.append([26,Matrix.block_diagonal([8*U,xid(2,1),xid(4,1)])])
# n = 30---------------
M = Matrix(QQ,2,2,[2,3,3,0])
Lattices.append([30,Matrix.block_diagonal([3*U,3*U,M])])
# n = 32---------------
M = Matrix(QQ,2,2,[4,2,2,6])
Lattices.append([32,Matrix.block_diagonal([5*U,M])])
# n = 33---------------
M = Matrix(QQ,2,2,[2,1,1,4])
Lattices.append([33,Matrix.block_diagonal([7*U,M])])
# n = 34---------------
Lattices.append([34,Matrix.block_diagonal([U,4*A2,xid(-12,1)])])
# n = 39---------------
Lattices.append([39,Matrix.block_diagonal([2*U,xid(4,1),xid(-4,1),xid(8,1)])])
# n = 40---------------
Lattices.append([40,Matrix.block_diagonal([xid(4,3),xid(-4,2)])])
# n = 46---------------
Lattices.append([46,Matrix.block_diagonal([A2,xid(6,1),xid(-18,1)])])
# n = 48---------------
Lattices.append([48,Matrix.block_diagonal([3*U,6*A2])])
# n = 49---------------
Lattices.append([49,Matrix.block_diagonal([2*U,2*A2,xid(-8,1)])])
# n = 51---------------
Lattices.append([51,Matrix.block_diagonal([2*U,xid(12,2)])])
# n = 54---------------
M = Matrix(QQ,3,3,[2,0,0,0,16,8,0,8,16])
Lattices.append([54,Matrix.block_diagonal([M])])
# n = 55---------------
Lattices.append([55,Matrix.block_diagonal([U,10*A2])])
# n = 56---------------
Lattices.append([56,Matrix.block_diagonal([xid(4,3),xid(-8,1)])])
# n = 61---------------
Lattices.append([61,Matrix.block_diagonal([3*U,4*A2])])
# n = 62---------------
M = Matrix(QQ,3,3,[6,0,3,0,6,3,3,3,12])
Lattices.append([62,Matrix.block_diagonal([M])])
# n = 63---------------
M = Matrix(QQ,3,3,[2,0,0,0,12,6,0,6,12])
Lattices.append([63,Matrix.block_diagonal([M])])
# n = 65---------------
Lattices.append([65,Matrix.block_diagonal([2*A2,xid(4,1),xid(-8,1)])])
# n = 70---------------
M1 = Matrix(QQ,3,3,[4,1,0,1,4,0,0,0,20])
M2 = Matrix(QQ,3,3,[4,2,2,2,6,1,2,1,16])
Lattices.append([70,Matrix.block_diagonal([M1]),Matrix.block_diagonal([M2])])
# n = 74---------------
M1 = Matrix(QQ,3,3,[2,1,0,1,4,0,0,0,28])
M2 = Matrix(QQ,3,3,[4,2,2,2,8,1,2,1,8])
Lattices.append([74,Matrix.block_diagonal([M1]),Matrix.block_diagonal([M2])])
# n = 75---------------
M = Matrix(QQ,4,4,[4,0,2,0,0,4,2,0,2,2,4,4,0,0,4,0])
Lattices.append([75,Matrix.block_diagonal([M])])
# n = 76---------------
M = Matrix(QQ,3,3,[4,0,0,0,8,0,0,0,12])
Lattices.append([76,M])
# n = 77---------------
M = Matrix(QQ,3,3,[4,0,0,0,8,4,0,4,8])
Lattices.append([77,Matrix.block_diagonal([M])])
# n = 78---------------
M = Matrix(QQ,3,3,[8,4,4,4,8,2,4,2,8])
Lattices.append([78,Matrix.block_diagonal([M])])
# n = 79---------------
M1 = Matrix(QQ,3,3,[2,1,0,1,8,0,0,0,12])
M2 = Matrix(QQ,3,3,[6,0,3,0,6,3,3,3,8])
Lattices.append([79,Matrix.block_diagonal([M1]),Matrix.block_diagonal([M2])])
# n = 80---------------
M = Matrix(QQ,3,3,[4,0,0,0,8,0,0,0,8])
Lattices.append([80,Matrix.block_diagonal([M])])
# n = 81---------------
M = Matrix(QQ,3,3,[4,0,2,0,4,2,2,2,12])
Lattices.append([81,Matrix.block_diagonal([M])])


# ------SANITY-CHECKS------
def checkDeterminante(L):
    for i in range(0,len(L)):
        g=L[i][1]
        G = Genus(g)
        print(L[i][0],g.det(),g.ncols())
        #print(G)

