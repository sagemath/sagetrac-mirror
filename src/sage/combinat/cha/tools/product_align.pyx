"""
Fast computation of alignment sequence for WQSym, CQSym, PQSym

"""
from sage.rings.integer import Integer
from sage.rings.integer cimport Integer as Int

Z0 = Int(0)
Z1 = Int(1)
Z2 = Int(2)

cpdef Int bin(list st):
    """
    Construit un objet binaire a partir d'une liste de bit a 1::
    
        sage: bin([2])
        4
        sage: bin([0,2])
        5
        sage: bin([0,1])
        3
        sage: bin([0])
        1
    """
    cdef Int res = Z0
    cdef int i
    for i in st:
        res += Z1<<i
    return res
    
cpdef Int bin2(int n):
    """
    Construit un binaire de n 1 consécutifs::
    
        sage: bin2(1)
        1
        sage: bin2(2)
        3
        sage: bin2(3)
        7
    """
    cdef Int res = Z0
    cdef int i
    for i in range(n):
        res += Z1<<i
    return res

cpdef list set2(Int bn):
    """
    Retourne la liste des bit a un::
    
        sage: set2(4)
        [2]
        sage: set2(7)
        [0, 1, 2]
        sage: set2(9)
        [0, 3]
        sage: set2(bin([0,3]))
        [0, 3]
    """
    return [i for i, v in enumerate(bn.digits(2)) if v]

cpdef tuple park_func_to_blocks(list pf, int shift_pos=0):
    """
    créer des segment de bits::
    
        sage: park_func_to_blocks([1,3,1,1, 5, 6,6,7, 9])[1]
        [15, 16, 224, 256]
        sage: set2(15)
        [0, 1, 2, 3]
        sage: sage: set2(16)
        [4]
        sage: sage: set2(224)
        [5, 6, 7]
        sage: sage: set2(256)
        [8]

    le premier element retourne est un dico... a voir quoi en faire 
    """
    cdef int i, n, k
    cdef int maX = 0
    cdef int shift = 0
    cdef Int b = Z0
    import collections
    d = collections.defaultdict(list)
    n = len(pf)
    for i in range(n):
        d[pf[i] + shift].append(i + shift_pos)
    blocks = []
    for k in d.keys():
        n = len(d[k])
        if k > maX:
            blocks.append(b)
            b = Z0
            shift = maX
        maX += n
        b = (b<<(n)) + (bin2(n)<<shift)
    blocks.append(b)
    blocks.pop(0)
    return d, blocks
    
cpdef Int segment_f(list segment):
    """
    a partir d'un segment d'entier, on calcul l'entier associé.
    
    TESTS::
    
        sage: segment_f([1,2,4])
        7
    
    ATTENTION si le segment n'est pas disjoints::
    
        sage: segment_f([1,3]) # b(1)=1, b(3)=101
        3
    """
    cdef Int r = Z0
    for block in segment:
        r |= block
    return r
    
def lines_output(list l, list r):
     print set2(segment_f(l)), "::", set2(segment_f(r)) 
     s = "" 
     seet = set2(segment_f(l))
     for i in range( max(seet)+1 ):
         s += "__" if i in seet else "  "
     print s
     s = ""
     seet = set2(segment_f(r))
     for i in range( max(seet)+1 ):
         s += "__" if i in seet else "  "
     print s
    
def shuffle(list seg1, list seg2):
    """
    TESTS::
        
        sage: def ssh(s1,s2): 
        sage:     k = 0
        sage:     for (l,r) in shuffle(s1,s2):
        ....:         k+=1
        ....:         lines_output(l, r)
        ....:     print k
        sage: b1 = park_func_to_blocks([1])[1]
        sage: b2 = park_func_to_blocks([1])[1]
        sage: ssh(b1,b2)
        [0] :: [0]
        __
        __
        [1] :: [0]
          __
        __
        [0] :: [1]
        __
          __
        3
        sage: b1 = park_func_to_blocks([1, 1])[1]
        b1 sage: b2 = park_func_to_blocks([1])[1]
        sage: ssh(b1,b2)
        [0, 1] :: [0]
        ____
        __
        [1, 2] :: [0]
          ____
        __
        [0, 1] :: [1]
        ____
          __
        [0, 1] :: [2]
        ____
            __
        4
        sage: b1 = park_func_to_blocks([1,2])[1]
        sage: ssh(b1,b1)
        [0, 1] :: [0, 1]
        ____
        ____
        [0, 2] :: [0, 1]
        __  __
        ____
        [1, 2] :: [0, 1]
          ____
        ____
        [1, 2] :: [0, 2]
          ____
        __  __
        [1, 2] :: [0, 3]
          ____
        __    __
        [1, 3] :: [0, 2]
          __  __
        __  __
        [0, 1] :: [0, 2]
        ____
        __  __
        [0, 1] :: [1, 2]
        ____
          ____
        [0, 2] :: [1, 2]
        __  __
          ____
        [0, 3] :: [1, 2]
        __    __
          ____
        [0, 2] :: [1, 3]
        __  __
          __  __
        11
    """
    def corps(list seg1, list seg2):
        cdef int i,n
        cdef Int shift, seg12tmp
        cdef list segtmp
        cdef Int seg12 = segment_f(seg2) | segment_f(seg1)
        n = len(seg1)
        for i in range(n-1):
            shift = seg1[i]<<1
            if (shift ^ seg1[i+1]) == (shift | seg1[i+1]):
                segtmp = seg1[:i] + [shift] + seg1[i+1:]
                seg12tmp = segment_f(segtmp) | segment_f(seg2)
                if seg12tmp == seg12:
                    for tup in shuffle(segtmp, seg2): 
                        yield tup
        shift = seg1[-1]<<1
        segtmp = seg1[:-1] + [shift]
        seg12tmp = segment_f(segtmp) | segment_f(seg2)
        if seg12tmp == seg12 or seg12tmp == ((seg12<<1) | 1):
            for tup in shuffle(segtmp, seg2): 
                yield tup
    
    yield seg1, seg2
    for tup in corps(seg1, seg2):
        yield tup
    for s2,s1 in corps(seg2, seg1):
        yield s1, s2
            
            
            
            
            
            
            