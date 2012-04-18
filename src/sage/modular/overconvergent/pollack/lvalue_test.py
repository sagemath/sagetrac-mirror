from sage.all import *

def lvalue_using_OMS_split(A,p,D,r,prec1,prec2):
    """
    Computes special value of p-adic L-series in dim 2 when p splits

    INPUT:
        - A: modular symbols space
        - p: good ordinary prime that splits in Hecke eigenvalue field
        - D: discriminant of quadratic twist, currently needs to be a fundamental discriminant of a number field
        - r: rank of A
        - prec1: computation for 1st symbol done up to O(p^prec1)
        - prec2: computation for 2nd symbol done up to O(p^prec2)


    OUTPUT:
        - special value of p-adic L-series associated to A

    EXAMPLES:

    sage: M = ModularSymbols(93,2,1)
    sage: S = M.cuspidal_submodule()
    sage: Ne = S.new_subspace()
    sage: D = Ne.decomposition()
    sage: A = D[0]
    sage: from sage.modular.overconvergent.pollack.lvalue_test import lvalue_using_OMS_split
    sage: lvalue_using_OMS_split(A,11,1,2,4,5)

    
    sage: from sage.modular.overconvergent.pollack.lvalue_test import lvalue_using_OMS_split
    sage: N = 103
    sage: M = ModularSymbols(N,2,1)
    sage: S = M.cuspidal_submodule()
    sage: Ne = S.new_subspace()  
    sage: D = Ne.decomposition()
    sage: A = D[0]
    sage: lvalue_using_OMS_split(A,19,1,2,2,2)
    """
    coeff = ZZ(r)/ZZ(2)
    from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_decomposition    
    phi = form_modsym_from_decomposition(A)
    phis1 = phi.coerce_to_Qp(p,prec1)
    phis2 = phi.coerce_to_Qp(p,prec2)
    phi1 = phis1[0][0]
    phi2 = phis2[1][0]

    #these are ok
    #print "phi1 = ", phi1
    #print "phi2 = ", phi2
    psi1 = phis1[0][1]
    
    psi2 = phis2[1][1]
    f = A.q_eigenform(p+1,'alpha')
    ap = A.eigenvalue(p)
    
    phi1p = phi1.p_stabilize_ordinary(p,ZZ(psi1(ap)),prec1)
    #print "phi1p = ", phi1p
    phi2p = phi2.p_stabilize_ordinary(p,ZZ(psi2(ap)),prec2)
    #print "phi2p = ", phi2p
    Phi1 = phi1p.lift_to_OMS_eigen(p,prec1)
    #print "Phi1 = ", Phi1
    v1 = pLfunction_coef(Phi1,alpha1,coeff,D,p+1)
    print "v1 = ", v1
    R1 = Qp(p,prec1)['x']
    x = R1.gen()
    h1 = x**2-ZZ(psi1(ap))*x+p
    alpha1 = h1.roots()[0][0]
    Phi2 = phi2p.lift_to_OMS_eigen(p,prec2)
    print "Phi2 = ", Phi2
    R2 = Qp(p,prec2)['x']
    x = R2.gen()
    h2 = x**2-ZZ(psi2(ap))*x+p
    alpha2 = h2.roots()[0][0]
    
    v2 = pLfunction_coef(Phi2,alpha2,coeff,D,p+1)
    return v1*v2
