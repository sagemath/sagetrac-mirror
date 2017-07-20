from cypari2.paridecl cimport *
from cypari2.gen cimport Gen
from sage.libs.pari import pari
from cypari2.stack cimport new_gen
from cysignals.signals cimport sig_on

def pari_hyperellpadicfrobenius(Gen H,ulong p,prec):
    """
    Wrapper around pari's hyperellpadicfrobenius function.

    TESTS::

        sage: x = QQ['x'].gen()
        sage: H = HyperellipticCurve(x^5-5*x^3+7*x^2+4*x+1)
        sage: m = H.padic_frobenius_matrix(7,5) # indirect doctest
        sage: m.charpoly()
        (1 + O(7^5))*x^4 + (6 + O(7^5))*x^3 + (1 + 3*7 + O(7^5))*x^2 + (6*7 + O(7^5))*x + (7^2 + O(7^5))

    """
    sig_on()
    return new_gen(hyperellpadicfrobenius(RgX_to_FpX(H.g, utoi(p)),p,prec))

def pari_nfhyperellpadicfrobenius(Gen h,Gen t,ulong p,prec):
    """
    Wrapper around pari's nfhyperellpadicfrobenius function.

    TESTS::

        sage: x = QQ['x'].gen()
        sage: K.<a> = NumberField(x^3-5*x+1)
        sage: x = K['x'].gen()
        sage: H = HyperellipticCurve(x^7+x^5-3*(a+1)*x^2-4*(a+1)*x-1)
        sage: P = K.ideal(13).factor()[0][0]
        sage: m = H.padic_frobenius_matrix(P,5) # indirect doctest
        sage: m.charpoly()
        (1 + O(13^5))*x^6 + (5 + 13 + O(13^5))*x^5 + (3 + 5*13 + 13^2 + O(13^5))*x^4 + (4 + 3*13 + 13^3 + O(13^5))*x^3 + (3*13^2 + 5*13^3 + 13^4 + O(13^5))*x^2 + (5*13^4 + O(13^5))*x
    """
    cdef GEN pp = utoi(p)
    cdef GEN T
    cdef GEN H
    sig_on()
    T = RgX_to_FpX(t.g, pp)
    H = RgX_to_FpXQX(h.g, T, pp)
    return new_gen(nfhyperellpadicfrobenius(H,T,p,prec))
