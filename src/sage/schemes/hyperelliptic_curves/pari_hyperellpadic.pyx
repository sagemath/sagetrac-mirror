from cypari2.paridecl cimport *
from cypari2.gen cimport Gen
from sage.libs.pari import pari
from cypari2.stack cimport new_gen
from cysignals.signals cimport sig_on

def pari_hyperellpadicfrobenius(Gen H,ulong p,prec):
    sig_on()
    return new_gen(hyperellpadicfrobenius(RgX_to_FpX(H.g, utoi(p)),p,prec))

def pari_nfhyperellpadicfrobenius(Gen h,Gen t,ulong p,prec):
    cdef GEN pp = utoi(p)
    cdef GEN T
    cdef GEN H
    sig_on()
    T = RgX_to_FpX(t.g, pp)
    H = RgX_to_FpXQX(h.g, T, pp)
    return new_gen(nfhyperellpadicfrobenius(H,T,p,prec))
