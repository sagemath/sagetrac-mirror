import canon
# G = PermutationGroup([[(1,2,3,4,5,6,7,8)]])
# I = IntegerVectorsModPermutationGroup(G)

# 

G = PermutationGroup([[(1,8,14,12,3,7,13,9,2,5,16,11),(4,6,15,10)], [(1,13,10),(2,14,12,3,15,9,4,16,11),(5,6),(7,8)]]) 

I = IntegerVectorsModPermutationGroup(G)

SGS = G.strong_generating_system()

V = [2,1,0,0,2,0,0,0,1,1,0,0,1,0,0,0] 

def Perm2Vect(p):
    return canon.Perm([p(i)-1 for i in range(1, 9)])

def SGS2cpp(sgs):
    res = canon.PermListList()
    for trans in sgs:
        rest = canon.PermList()
        for x in trans:
            rest.append(Perm2Vect(x))
        res.append(rest)
    return res

SGScpp = SGS2cpp(SGS)

assert(canon.is_canonical(SGScpp, canon.Perm(V)) == I.is_canonical(V))

print "BEFORE :"
timeit("I.is_canonical(V)")
print "AFTER :"
timeit("canon.is_canonical(SGScpp, canon.Perm(V))")
