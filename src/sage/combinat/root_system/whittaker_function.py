def whittaker_function(weight, w, q1, q2):
    """
    Return the Whittaker function ?

    INPUT:

    - weight -- element of the weight lattice ?

    - w -- an element of the Weyl group ?

    - q1 -- a parameter ?

    - q2 -- a parameter ?

    EXAMPLES::
    
        sage: from sage.combinat.root_system.whittaker_function import whittaker_function
        sage: L = RootSystem(["A",2]).ambient_space()
        sage: K = QQ['v'].fraction_field()
        sage: v = K.gen()
        sage: omega = L.fundamental_weights()
        sage: W = L.weyl_group()

        sage: w = W.an_element()
        sage: whittaker_function(omega[1], w, -1, v)
        v^2*B[(-1, -1, 3)] + (v^2-v)*B[(-1, 1, 1)] + (v^2-v)*B[(-1, 0, 2)] + (v^2-v)*B[(0, -1, 2)] + (v^2-2*v+1)*B[(0, 0, 1)]

        sage: L = RootSystem(["A",1]).ambient_space()
        sage: K = QQ['v'].fraction_field()
        sage: v = K.gen()
        sage: omega = L.fundamental_weights()
        sage: W = L.weyl_group()

        sage: w = W.an_element()
        sage: whittaker_function(omega[1], w, -1, v)
        v^2*B[(-1, -1, 3)] + (v^2-v)*B[(-1, 1, 1)] + (v^2-v)*B[(-1, 0, 2)] + (v^2-v)*B[(0, -1, 2)] + (v^2-2*v+1)*B[(0, 0, 1)]
        sage: whittaker_function(omega[1], w, -1, v)
        (v^2-v)*B[(-1, 1)] + v^2*B[(-1, 0)] + (v^2-v)*B[(0, 0)]

        sage: whittaker_function(0*omega[1], W.one(), -1, v)
        B[(0, 0)]
        sage: whittaker_function(omega[1], W.one(), -1, v)
        B[(1, 0)]
    """
    from sage.categories.pushout import pushout
    K = pushout(q1.parent(), q2.parent())
    q1 = K(q1)
    q2 = K(q2)
    L = weight.parent()
    KL = L.algebra(K)
    T = KL.demazure_lusztig_operators(q1, q2, convention="haiman")
    rho = L.rho()
    result = T[w](KL.monomial(weight + rho))
    return result.map_support((-rho)._add_)
