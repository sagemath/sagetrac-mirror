from sage.combinat.parking_functions import ParkingFunction
## must be in parking_function
## use in pqsym
def to_park( word ):
    """
    TESTS::

        sage: from sage.combinat.some_hopf_algebra.parking_functions import to_park
        sage: to_park([1,1,1,5])
        [1, 1, 1, 4]
        sage: to_park([1,2,4,5])
        [1, 2, 3, 4]
        sage: to_park([3,1,2,4,5])
        [3, 1, 2, 4, 5]
    """
    if len( word ) == 0:
        return ParkingFunction( [] )
    import collections
    dic_pos = collections.defaultdict( list )
    for pos in range( len( word ) ):
        dic_pos[word[pos]].append( pos )
    i = 1; lt = min( word )
    w = {}
    for l in sorted( dic_pos.keys() ):
        i = min( i + abs( l - lt ), len( w.keys() ) + 1 )
        lt = l
        for pos in dic_pos[l]:
            w[pos] = i
    return ParkingFunction( w.values() )