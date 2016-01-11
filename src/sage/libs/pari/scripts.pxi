r"""
Declarations for functions defined in GP scripts compiled with gp2c.

Only the functions meant for public use should be declared here, not
the functions used privately within the script.

AUTHORS:

 - Jeroen Demeyer (2015-07-10): initial version (:trac:`15809`)

"""

# Denis Simon's scripts
cdef extern:
    long DEBUGLEVEL_qfsolve
    long DEBUGLEVEL_ell
    long LIM1
    long LIM3
    long LIMTRIV
    long ELLREDGENFLAG
    long COMPLETE
    long MAXPROB
    long LIMBIGPRIME
    long DEBUGLEVEL_res

    void init_simon()
    GEN ellrank(GEN ell, GEN help, long prec)
    GEN bnfellrank(GEN bnf, GEN ell, GEN help, GEN bigflag, GEN flag3, long prec)
