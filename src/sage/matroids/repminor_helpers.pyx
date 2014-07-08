cpdef search_assgn(Ga,Gb,assgn,cand=None):
    """ Implement Ullman's algorithm for exact subgraph isomorphism 
    """
    uGa_S1 = set(Ga['S1']) - set(assgn.keys())
    uGb_S1= set(Gb['S1']) - set(assgn.values())
    uGa_S2 = set(Ga['S2']) - set(assgn.keys())
    uGb_S2= set(Gb['S2']) - set(assgn.values())
    uGa = uGa_S1 | uGa_S2
    uGb = uGb_S1 | uGb_S2
    print 'Ga',Ga
    print 'Gb',Gb
    if cand == None:
        # first call. populate cand
        cand = {}
        for va in set(Ga['S1']):
            cand[va]=[]
            for vb in set(Gb['S1']):
                print 'Ga:',va,'Gb:',vb,Ga[va][0],'?',Gb[vb][0]
                if Gb[vb][0] >= Ga[va][0]:
                    cand[va].append(vb)
        for va in set(Ga['S2']):
            cand[va]=[]
            for vb in set(Gb['S2']):
                print 'Ga:',va,'Gb:',vb,Ga[va][0],'?',Gb[vb][0]
                if Gb[vb][0] >= Ga[va][0]:
                    cand[va].append(vb)
    print cand
    # loop over uGa
    for v in uGa:
        # check that whenever vertex c of Gb is among candidates for vertex v of 
        # GA, then every neighbor of v has at least one candidate among 
        # neighbors of c.
        for c in cand[v]:
            valid_neighbours = set(Ga[1])
    return True

