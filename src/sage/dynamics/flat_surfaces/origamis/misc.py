r"""
Old code that might moved or be removed.
"""

def rauzy_move((ltop,lbot),l):
    r"""
    Start from a generalized permutation given as ((ltop,lbot),l)
    where ltop and lbot are lists with the top and bottom labels
    (the labels are assumed to pair up -- no check), and l is a
    dictionary of lengths for the labels used in ltop, lbot.

    rauzy_move will apply one Rauzy move at the rightmost end:

    - if rightmost top and bottom intervals have distinct lengths,
      do the usual thing

    - if rightmost top and bottom intervals have same label, error

    - if rightmost top and bottom intervals have same length but
      not same label, just remove one of the labels

    Note: this code was originally used to compute the connected component of an
    origami. This is no more ncessary but it might be a good idea to keep for
    double check
    """
    if len(ltop) < 2 and len(lbot) < 2:
        raise ValueError("cannot induce, too few intervals")
    ktop = ltop[-1]
    kbot = lbot[-1]
    lltop = ltop[:-1]
    llbot = lbot[:-1]
    ll = l.copy()
    if l[ktop] - l[kbot] > 0:
        # print "bottom interval shorter\n"

        # ltop[-1] is some symbol ktop, which appears twice: once as ltop[-1], and the other
        # time either as ltop[j] for j < len(ltop) - 1 or as lbot[j] for j < len(lbot) - 1;
        # as ltop[-1], ktop needs to be replaced by [ltop[-1],lbot[-1]],
        # the other occurrence of ktop needs to be replaced
        # - by [ltop[-1],lbot[-1]] if it is lbot[j]
        # - by [lbot[-1],ltop[-1]] if it is ltop[j]
        # the dictionary of lengths needs one update: ll[ktop] = l[ktop] - l[kbot]

        if ktop in lltop: lltop.insert(lltop.index(ktop),kbot)
        else: llbot.insert(llbot.index(ktop)+1,kbot)
        lltop.append(ktop)
        ll[ktop] = l[ktop] - l[kbot]

    elif l[ktop] - l[kbot] < 0:
        # print "top interval shorter\n"

        # lbot[-1] is some symbol kbot, which appears twice: once as lbot[-1], and the other
        # time either as lbot[j] for j < len(lbot) - 1 or as ltop[j] for j < len(ltop) - 1;
        # as lbot[-1], kbot needs to be replaced by [lbot[-1],ltop[-1]],
        # the other occurrence of ktop needs to be replaced
        # - by [lbot[-1],ltop[-1]] if it is ltop[j]
        # - by [ltop[-1],lbot[-1]] if it is lbot[j]
        # the dictionary of lenghts needs one update: ll[kbot] = l[kbot] - l[ktop]

        if kbot in llbot: llbot.insert(llbot.index(kbot),ktop)
        else: lltop.insert(lltop.index(kbot)+1,ktop)
        llbot.append(kbot)
        ll[kbot] = l[kbot] - l[ktop]

    elif ktop == kbot:
        # in this case the surface has a torus component... weird.
        raise ValueError("it seems the surface \n %s \n %s \n %s either is disconnected or has genus 1" %s(ltop,lbot,l))
    else:
        # ktop != kbot but the intervals labeled ktop and kbot have equal lengths
        # just remove ktop and kbot at the end of ltop and lbot, and match the
        # remaining ktop and kbot (keep only one of the names; also in the length dict)
        if ktop in lltop: lltop[lltop.index(ktop)] = kbot
        else: llbot[llbot.index(ktop)] = kbot
        ll.pop(ktop)

    return ((lltop,llbot),ll)

def iet_of_origami(o):
    r"""
    x, y are permutations which together represent an origami
    this origami is assumed to be a torus cover
    (ie corresponds to an abelian differential)

    First get a multi-iet by considering the first return to I, the union of antidiagonals
    of all squares, of the flow in direction (1,u), where u is the golden mean.
    The antidiagonals are parametrized by their horizontal coordinate, so that they appear
    to have length 1, therefore I has length n. Apply Rauzy moves until I has length 1.

    We could define an iet by using this permutation and this length vector
    p = Permutation([(lambda k: 2*y((k+1)//2) if k%2 else 2*x(k//2)-1)(j) for j in range(1,n+1)])
    l = [2-u,u-1] * n
    """
    from sage.dynamics.interval_exchanges.all import iet
    from sage.rings.number_field.all import NumberField
    from sage.rings.polynomial.all import PolynomialRing
    from sage.rings.all import RR

    sigma = o.r()
    tau = o.u()
    n = o.nb_squares()

    R = PolynomialRing(QQ,'x')
    x = R.gen()
    K = NumberField(x**2 - x - 1, 'a',embedding=RR(1.618))
    phi = K.gen()

    ltop = []
    for k in xrange(1,n+1):
        ltop.extend([(k,0),(k,1)])
    lbot = []
    for k in xrange(1,n+1):
        lbot.extend([(sigma.inverse()(k),1),(tau.inverse()(k),0)])
    l = {}
    for k in xrange(1,n+1):
        l[k,0] = 2 - RR(phi)
        l[k,1] = RR(phi) - 1

    while sum(l.values())>1:
        ((ltop, lbot), l) = rauzy_move((ltop,lbot),l)

    return iet.Permutation(ltop,lbot),l


