# coding=utf8
r"""
Beta-adic tools.

Beta-adic is a way to write numbers in the form

    - :math:`\sum_{i=0}^\infty \beta^i c_i`

where :math:`beta` is a element of a field (for example a complex number),
and the :math:`c_i` are varying in a finite set of digits.
The possible finite sequences of digits are given by a deterministic automaton.

AUTHORS:

- Paul Mercat (2013) -  I2M AMU Aix-Marseille Universite -initial version
- Dominique Benielli (2018) - Labex Archimede - I2M -
  AMU Aix-Marseille Universite - Integration in SageMath

EXAMPLES::

    # Tribonacci
    sage: m = BetaAdicSet(x^3-x^2-x-1, {0,1})
    sage: print(m)
    b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 1 state and 2 letters
    sage: ared = m.reduced_words_automaton()
    sage: print(ared)
    DetAutomaton with 4 states and an alphabet of 2 letters
"""
# *****************************************************************************
#  Copyright (C) 2013 Paul Mercat <mercatp@icloud.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
# *****************************************************************************
from sage.sets.set import Set
from sage.rings.qqbar import QQbar
from sage.rings.padics.factory import Qp
from libc.stdlib cimport malloc, free
from math import pi as pi_number
from cysignals.signals cimport sig_on, sig_off, sig_check
cimport sage.combinat.words.cautomata
from sage.combinat.words.cautomata cimport DetAutomaton, FreeAutomaton
from sage.combinat.words.cautomata_generators import DetAutomatonGenerators
from sage.rings.integer import Integer
from sage.combinat.words.morphism import WordMorphism
from sage.rings.number_field.number_field import NumberField

from cysignals.signals cimport sig_on, sig_off

cdef extern from "complex.h":
    cdef cppclass Complexe:
        double x, y

cdef extern from "Automaton.h":
    ctypedef char bool
    cdef cppclass State:
        int* f
        bool final

    cdef cppclass Automaton:
        State* e  # states
        int n   # number of states
        int na  # number of letters
        int i  # initial state

    cdef cppclass Transition:
        int l  # label
        int e  # arrival state

    cdef cppclass NState:
        Transition* a
        int n
        bool final
        bool initial

    cdef cppclass NAutomaton:
        NState* e  # states
        int n   # number of states
        int na  # number of letters

    Automaton CopyAutomaton (Automaton a, int nalloc, int naalloc)

cdef extern from "relations.h":
    cdef cppclass Element:
        int *c  # liste des n coeffs

    cdef cppclass PlaceArch:
        Complexe *c  # 1, b, b^2, ... for this place

    # class containing the informations needed to compute the relations automaton
    cdef cppclass InfoBetaAdic:
        int n         # degree
        Element bn    # expression of b^n as a polynome in b of degree < n
        Element b1    # expression of 1/b as a polynome in b of degree < n
        Element *c    # list of digits used for the calculation of relations automaton
        int nc        # number of digits
        int ncmax     # number of allocated digits
        PlaceArch *p  # list of na places
        double *cM    # square of the bound
        int na        # number of places

    Element NewElement(int n)
    void FreeElement(Element e)
    InfoBetaAdic allocInfoBetaAdic(int n, int na, int ncmax, int nhash, bool verb)
    void freeInfoBetaAdic(InfoBetaAdic *iba)
    Automaton RelationsAutomatonT(InfoBetaAdic *iba2, Element t, bool isvide, bool ext, bool verb)

cdef extern from "numpy/arrayobject.h":
    ctypedef int intp
    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
        cdef int flags

cimport numpy

from libc.stdint cimport uint8_t, uint32_t
from libc.math cimport log
from libc.math cimport ceil
from libc.math cimport floor
from libc.math cimport round
from libc.math cimport fabs

cdef extern from "draw.h":
    ctypedef unsigned char uint8
    cdef cppclass Color:
        uint8 r
        uint8 g
        uint8 b
        uint8 a
    cdef cppclass Surface:
        Color **pix
        int sx, sy
    cdef cppclass Complexe:
        double x
        double y
    cdef cppclass BetaAdic:
        Complexe b
        Complexe* t  # list of translations
        int n        # number of translations
        Automaton a
    cdef cppclass BetaAdic2:
        Complexe b
        Complexe* t  # list of translations
        int n        # number of translations
        Automaton* a
        int na
    ctypedef Color* ColorList

    # void *GetSDL_SurfaceFromNumpy (numpy.ndarray na)
    # void SDL_SurfaceToNumpy (void *ss, numpy.ndarray na)
    # void TestSDL()
    void SurfaceToNumpy (Surface *s, numpy.ndarray na)
    Surface NewSurface(int sx, int sy)
    void FreeSurface(Surface s)
    ColorList NewColorList(int n)
    void FreeColorList(ColorList l)
    Color randColor(int a)
    #    Automate NewAutomate (int n, int na)
    #    void FreeAutomate(Automate a)
    void FreeAutomatons(Automaton* a, int n)
    BetaAdic NewBetaAdic(int n)
    void FreeBetaAdic(BetaAdic b)
    BetaAdic2 NewBetaAdic2(int n, int na)
    void FreeBetaAdic2(BetaAdic2 b)
    int *DrawZoom(BetaAdic b, int sx, int sy, int n, int ajust, Color col, int nprec, double sp, int verb)
    Automaton UserDraw(BetaAdic b, int sx, int sy, int n, int ajust, Color col, double sp, int verb)
    #  void WordZone (BetaAdic b, int *word, int nmax)
    int *Draw(BetaAdic b, Surface s, int n, int ajust, Color col, int nprec, double sp, int verb)
    void Draw2(BetaAdic b, Surface s, int n, int ajust, Color col, double sp, int verb)
    void DrawList(BetaAdic2 b, Surface s, int n, int ajust, ColorList lc, double alpha, double sp, int nprec, int verb)
    void print_word(BetaAdic b, int n, int etat)


cdef uint32_t moy(uint32_t a, uint32_t b, float ratio):
    return <uint32_t><uint8_t>((a%256)*(1.-ratio) + (b%256)*ratio) | \
           (<uint32_t>(<uint8_t>(((a>>8)%256)*(1.-ratio) + ((b>>8)%256)*ratio)))<<8 | \
           (<uint32_t>(<uint8_t>(((a>>16)%256)*(1.-ratio) + ((b>>16)%256)*ratio)))<<16 | \
           (<uint32_t>(<uint8_t>((a>>24)*(1.-ratio) + (b>>24)*ratio)))<<24;


cdef double fmax(double a, double b):
    if a < b:
        return b
    return a

# plot the Rauzy fractal corresponding to the direction vector d,
# for the C-adic system given by the Cassaigne's algorithm
def plot_Cadic(numpy.ndarray dv, int sx=800, int sy=600,
               float mx=-2, float my=-2, float Mx=2, float My=2,
               int n=1000, int nptsmin=50000, int nptsmax=60000,
               bool verb=False, bool printl=True, bool get_ndarray=False):
    """
    Plot the Rauzy fractal corresponding to the direction vector ``dv``
    for the C-adic system given by the Cassaigne's algorithm

    INPUT:

        - ``dv``- ndarray array , direction vector

        - ``sx`` int -- (default 800) size of Image direction x

        - ``sy`` int -- (default 60) size of Image direction y

        - ``mx`` float -- (default -2)

        - ``my`` float -- (default -2)

        - ``Mx``  float -- (default 2)

        - ``My`` float  -- (default 2)

        - ``n`` int -- (default 1000)

        - ``nptsmin`` int -- (default 50000)

        - ``nptsmax`` int -- (default 60000)

        - ``verb`` bool -- (default ''False'')

        - ``printl`` bool -- (default ''False'')


        OUTPUT:

        Plot the Rauzy fractal corresponding to the direction vector dv.

        EXAMPLES::

            sage: from sage.arith.beta_adic import plot_Cadic
            sage: import numpy as np
            sage: plot_Cadic(np.array((random(), random(), random())))      # random

    """
    cdef numpy.ndarray l, d, im
    cdef int i, j, k, u, nA, i0, e, e0, npts, su, rsu
    cdef uint32_t x, y
    cdef uint32_t color
    cdef float fx, fy

    npts = 0
    color = 255 << 24
    colors = [255 | 255 << 24, 255 << 8 | 255 << 24, 255 << 16 | 255 << 24]

    import numpy as np
    d = np.empty(3, dtype=np.float)
    d[0] = <float>dv[0]
    d[1] = <float>dv[1]
    d[2] = <float>dv[2]
    from sage.combinat.words.morphism import WordMorphism
    s = WordMorphism('a->a,b->ac,c->b')
    t = WordMorphism('a->b,b->ac,c->c')
    auts = s.DumontThomas(proj=False)
    autt = t.DumontThomas(proj=False)
    aut = [auts, autt]
    A = [np.array(a) for a in auts.alphabet]
    nA = len(A)
    #if autt.alphabet != A:
    #    raise RuntimeError("The two Dumont-Thomas automata must have the same alphabet !")
    ms = s.incidence_matrix()
    mt = t.incidence_matrix()
    if verb:
        print("ms=%s"%ms)
        print("mt=%s"%mt)
    msi = (ms**(-1)).numpy()
    mti = (mt**(-1)).numpy()
    if verb:
        print("msi=%s"%msi)
        print("mti=%s"%mti)
    lm = [ms.numpy(), mt.numpy()]
    # compute an orthonormal basis
    v1 = np.array([1, -1, 0])
    v2 = np.array([1, 0, -1])
    v1 = v1 - v1.dot(d)/d.dot(d)*d
    v2 = v2 - v2.dot(d)/d.dot(d)*d
    from sage.functions.other import sqrt
    v1 = v1/sqrt(v1.dot(v1))
    v2 = v2/sqrt(v2.dot(v2))
    v2 = v2 - v1.dot(v2)*v1
    # Cassaigne's algorithm
    l = np.empty(n, dtype=np.int8)
    m = np.identity(3, dtype=np.int)
    v0 = np.zeros(3, dtype=np.int)
    v0[0] = 1
    su = 0
    for i in range(n):
        if d[0] > d[2]:
            d = msi.dot(d)
            l[i] = 0
        else:
            d = mti.dot(d)
            l[i] = 1
        m = m.dot(lm[l[i]])
        rsu = su
        su = sum(m.dot(v0))
        if rsu > nptsmin:
            n = i
            break
        if su > nptsmax:
            n = i
            break
        d = d/sum(d)
        if verb:
            print("d=%s" % d)
    if verb or printl:
        print("n=%s, l=%s"%(n, l[:n]))
    # Draw the Rauzy fractal
    im = np.empty([sy, sx], dtype=np.dtype(
        (np.uint32, {'r': (np.uint8, 0), 'g': (np.uint8, 1),
                     'b': (np.uint8, 2), 'a': (np.uint8, 3)})))
    # im.fill(0) #fill the image with transparent
    im.fill(255 | 255 << 8 | 255 << 16 | 255 << 24)  # fill with white

    if verb:
        print("A=%s" % A)
        print("nA=%s" % nA)

    p = [(np.zeros(3, dtype=np.int), 0, 0)]
    while len(p) > 0:
        k = len(p)-1
        u = l[n-k-1]
        # print("k=%s"%k)
        t, i, e = p[-1]
        # print("t=%s, i=%s, e=%s"%(t, i, e))
        if k == n:
            # we draw the point t
            # print(t)
            fx = (t.dot(v1) - mx)*sx/(Mx-mx)
            fy = (t.dot(v2) - my)*sy/(My-my)
            x = <uint32_t> fx
            y = <uint32_t> fy
            if verb:
                print(t)
                print(fx, fy)
                print(x, y)
                # print("")
            if x < sx and y < sy:
                if x+1 < sx and y+1 < sy:
                    im[y, x] = moy(im[y, x], colors[e], (1.-fx+x)*(1.-fy+y))
                    im[y, x+1] = moy(im[y, x+1], colors[e], (fx-x)*(1.-fy+y))
                    im[y+1, x] = moy(im[y+1, x], colors[e], (1.-fx+x)*(fy-y))
                    im[y+1, x+1] = moy(im[y+1, x+1], colors[e], (fx-x)*(fy-y))
                else:
                    im[y, x] = colors[e]
            npts += 1
            # increment
            # print("increment...")
            while True:
                t, i, e = p.pop()
                k = len(p)
                if k == 0:
                    break
                t0, i0, e0 = p[-1]
                u = l[n-k]
                # print("k=%s, u=%s, t=%s, i=%s, e=%s"%(k, u, t, i, e))
                while True:
                    i0 += 1
                    if i0 == nA or aut[u].succ(e0, i0) != -1:
                        break
                # print("i=%s"%i)
                if i0 != nA:
                    p[-1] = (t0, i0, e0)
                    p.append((lm[u].dot(t0)+A[i0], 0, aut[u].succ(e0, i0)))
                    break
        else:
            i = 0
            while i < nA and aut[u].succ(e, i) == -1:
                i += 1
            # print("starting i=%s k=%s u=%s t=%s e=%s"%(i, k, u, t, e))
            p[-1] = (t, i, e)
            p.append((lm[u].dot(t)+A[i], 0, aut[u].succ(e, i)))
        #for j2, (m2, t2, i2, e2) in enumerate(p):
            #print("%s : m=%s, t=%s, i=%s, e=%s"%(j2, m2, t2, i2, e2))
    if printl:
        print("%s pts computed."%npts)
    if get_ndarray:
        return im
    from PIL import Image
    return Image.fromarray(im, 'RGBA')


# plot the Rauzy fractal corresponding to the direction vector d,
# for the C-adic system given by the Cassaigne's algorithm
def plot_Cadic2(numpy.ndarray dv, int sx=800, int sy=600,
                float mx=-2, float my=-2, float Mx=2, float My=2,
                int n=40, bool verb=False, bool printl=True):
    cdef numpy.ndarray l, d, im
    cdef int i, j, k, u, nA, i0, e, e0, npts
    cdef uint32_t x, y
    cdef uint32_t color
    cdef float fx, fy

    npts = 0
    color = 255 << 24
    import numpy as np
    d = np.empty(3, dtype=np.float)
    d[0] = <float>dv[0]
    d[1] = <float>dv[1]
    d[2] = <float>dv[2]
    from sage.combinat.words.morphism import WordMorphism
    s = WordMorphism('a->a,b->ac,c->b')
    t = WordMorphism('a->b,b->ac,c->c')
    auts = s.DumontThomas(proj=False).mirror()
    autt = t.DumontThomas(proj=False).mirror()
    aut = [auts, autt]
    A = [np.array(a) for a in auts.alphabet]
    nA = len(A)
    # if autt.alphabet != A:
    #    raise RuntimeError("The two Dumont-Thomas automata must have the same alphabet !")
    ms = s.incidence_matrix()
    mt = t.incidence_matrix()
    if verb:
        print("ms=%s" % ms)
        print("mt=%s" % mt)
    msi = (ms**(-1)).numpy()
    mti = (mt**(-1)).numpy()
    if verb:
        print("msi=%s" % msi)
        print("mti=%s" % mti)
    lm = [ms.numpy(), mt.numpy()]
    # compute an orthonormal basis
    v1 = np.array([1,-1,0])
    v2 = np.array([1,0,-1])
    v1 = v1 - v1.dot(d)/d.dot(d)*d
    v2 = v2 - v2.dot(d)/d.dot(d)*d
    from sage.functions.other import sqrt
    v1 = v1/sqrt(v1.dot(v1))
    v2 = v2/sqrt(v2.dot(v2))
    v2 = v2 - v1.dot(v2)*v1
    # Cassaigne's algorithm
    l = np.empty(n, dtype=np.int8)
    for i in range(n):
        if d[0] > d[2]:
            d = msi.dot(d)
            l[i] = 0
        else:
            d = mti.dot(d)
            l[i] = 1
        d = d/sum(d)
        if verb:
            print("d=%s" % d)
    if verb or printl:
        print("l=%s" % l)
    # Draw the Rauzy fractal
    im = np.empty([sy, sx], dtype=np.dtype(
        (np.uint32, {'r': (np.uint8, 0), 'g': (np.uint8, 1),
                     'b': (np.uint8, 2), 'a': (np.uint8, 3)})))
    # im.fill(0) #fill the image with transparent
    im.fill(255 | 255 << 8 | 255 << 16 | 255 << 24)  # fill with white

    if verb:
        print("A=%s" % A)
        print("nA=%s" % nA)

    p = [(np.identity(3, dtype=np.int), np.zeros(3, dtype=np.int), 0, 0)]
    while len(p) > 0:
        k = len(p)-1
        u = l[k]
        #print("k=%s"%k)
        m, t, i, e = p[-1]
        #print("t=%s, i=%s, e=%s"%(t, i, e))
        if k == n-1:
            #we draw the point t
            #print(t)
            fx = (t.dot(v1) - mx)*sx/(Mx-mx)
            fy = (t.dot(v2) - my)*sy/(My-my)
            x = <uint32_t>fx
            y = <uint32_t>fy
            if verb:
                print(t)
                print(fx,fy)
                print(x,y)
                #print("")
            if x < sx and y < sy:
                #if x+1 < sx and y+1 < sy:
                #    im[y,x] = moy(im[y,x], color, (1.-fx+x)*(1.-fy+y))
                #    im[y,x+1] = moy(im[y,x+1], color, (fx-x)*(1.-fy+y))
                #    im[y+1,x] = moy(im[y+1,x], color, (1.-fx+x)*(fy-y))
                #    im[y+1,x+1] = moy(im[y+1,x+1], color, (fx-x)*(fy-y))
                #else:
                im[y,x] = color
            npts += 1
            #increment
            #print("increment...")
            while True:
                m, t, i, e = p.pop()
                k = len(p)
                if k == 0:
                    break
                m0, t0, i0, e0 = p[-1]
                u = l[k-1]
                nA = aut[u].n_succs(e0)
                # print("k=%s, u=%s, t=%s, i=%s, e=%s"%(k, u, t, i, e))
                while True:
                    i0 += 1
                    if i0 == nA or aut[u].succ(e0, i0) != -1:
                        break
                # print("i=%s"%i)
                if i0 != nA:
                    p[-1] = (m0, t0, i0, e0)
                    p.append((m, t0+m0.dot(A[i0]), 0, aut[u].succ(e0, i0)))
                    break
        else:
            i = 0
            nA = aut[u].n_succs(e)
            while i < nA and aut[u].succ(e, i) == -1:
                i += 1
            # print("starting i=%s k=%s u=%s t=%s e=%s"%(i, k, u, t, e))
            p[-1] = (m, t, i, e)
            p.append((m.dot(lm[u]), lm[u].dot(t)+A[i], 0, aut[u].succ(e, i)))
        # for j2, (m2, t2, i2, e2) in enumerate(p):
            # print("%s : m=%s, t=%s, i=%s, e=%s"%(j2, m2, t2, i2, e2))
    print("%s pts computed." % npts)
    from PIL import Image
    return Image.fromarray(im, 'RGBA')


# compute the p-adic absolute value
def absp(c, p, d):
    """
    Computation of the p-adic absolute value.

    INPUT:

    - ``c`` -- the algebraic number for which we compute the absolute value

    - ``p`` -- the prime number

    - ``d`` -- the degree

    OUTPUT:

    The p-adic absolute value.

    TESTS:

        sage: absp(1, 2, 3) # not implemented

    """
    return ((c.polynomial())(p).norm().abs())**(1/d)


cdef getElement(e, Element r, int n):
    cdef j
    p = e.lift()
    for j in range(n):
        r.c[j] = p[j]

cdef InfoBetaAdic initInfoBetaAdic(self,
                                   Ad=None, plus=True, nhash=1000003,
                                   verb=False) except *:
    b = self.b
    K = b.parent()
    A = self.a.alphabet
    if verb:
        print("A = %s" % A)

    if verb:
        print(K)

    # determine the places to consider
    parch = []
    for p in K.places():  # archimedian places
        if plus:
            if p(b).abs() > 1:
                parch += [p]
        else:
            if p(b).abs() < 1:
                parch += [p]
    pi = K.defining_polynomial()
    from sage.arith.misc import gcd
    # return the polynomial with integer coefficients and capacity 1
    pi = pi / gcd(pi.list())
    if verb:
        print("pi=%s" % pi)
    # list of concerned prime numbers
    lp = (Integer(pi.list()[0])).prime_divisors()
    if verb:
        print("lp=%s" % lp)
    # list of the considered ultrametric places
    pultra = []
    for p in lp:
        # find every places behind p in the field K
        k = Qp(p)
        Kp = k['a']
        a = Kp.gen()
        for f in pi(a).factor():
            kp = f[0].root_field('e')
            if kp == k:
                c = f[0].roots(kp)[0][0]
            else:
                c = kp.gen()
            if verb:
                print("c=%s (abs=%s)" % (c, (c.norm().abs())**(1/f[0].degree())))
            if plus:
                if (c.norm().abs())**(1/f[0].degree()) > 1:
                    pultra += [(c, f[0].degree())]
            else:
                if (c.norm().abs())**(1/f[0].degree()) < 1:
                    pultra += [(c, f[0].degree())]

    if verb:
        print("spaces: ")
        print(parch)
        print(pultra)

    if (len(pultra) > 0):
        raise ValueError("Not implemented for b algebraic non-integer.")
    # compute the max bound for each absolute value
    if Ad is None:
        Ad = Set([c-c2 for c in A for c2 in A])
    Ad = [K(c) for c in Ad]
    if verb:
        print("Ad = %s" % Ad)

    n = K.degree()
    na = len(parch)
    ncmax = len(Ad)
    cdef InfoBetaAdic i
    if verb:
        print("alloc...")
    sig_on()
    i = allocInfoBetaAdic(n, na, ncmax, nhash, verb)
    sig_off()
    cdef int j
    # initialize bn
    if verb:
        print("init bn...")
    getElement(b**n, i.bn, n)
    # initialize b1
    if verb:
        print("init b1...")
    getElement(1/b, i.b1, n)
    # initialize places
    if verb:
        print("init places...")
    for k in range(na):
        for j in range(n):
            i.p[k].c[j] = complex(parch[k](b**j))
    # initialize digits and bounds
    if verb:
        print("init digits...")
    initCdInfoBetaAdic(self, &i, Ad=Ad, parch=parch, verb=verb)
    return i

cdef initCdInfoBetaAdic(self, InfoBetaAdic *i, list Ad, list parch, verb=False):
    if verb:
        print("initCdInfoBetaAdic Ad = %s" % Ad)
    m = dict([])
    for p in parch:
        m[p] = max([p(c).abs() for c in Ad])/abs(1-p(self.b).abs())
    if verb:
        print("bounds : %s" % m)
    # conversion to C
    i.nc = len(Ad)
    if i.nc > i.ncmax:
        raise ValueError("Too much digits : %d > %d max (initialize BetaAdicSet with more digits)."%(i.nc, i.ncmax))
    for j, c in enumerate(Ad):
        getElement(c, i.c[j], i.n)
    for j, p in enumerate(parch):
        i.cM[j] = m[p]**2

cdef Complexe complex(c):
    cdef Complexe r
    r.x = c.real()
    r.y = c.imag()
    return r

cdef Color getColor(c):
    if len(c) < 4:
        raise ValueError("Colors must be defined by 4 float numbers between 0 and 1.")
    cdef Color r
    r.r = c[0] * 255
    r.g = c[1] * 255
    r.b = c[2] * 255
    r.a = c[3] * 255
    return r

cdef surface_to_img(Surface s):
    import numpy as np
    from PIL import Image
    # arr = np.empty([s.sy, s.sx], dtype=['uint8', 'uint8', 'uint8', 'uint8'])
    # arr = np.empty([s.sy, s.sx], dtype=[('r', 'uint8'), ('g', 'uint8'),('b', 'uint8'), ('a', 'uint8')])
    # arr = np.zeros([s.sy, s.sx], dtype=[('r', 'uint8'), ('g', 'uint8'),('b', 'uint8'), ('a', 'uint8')])
    arr = np.empty([s.sy, s.sx], dtype=np.dtype((np.uint32, {'r':(np.uint8,0), 'g':(np.uint8,1), 'b':(np.uint8,2), 'a':(np.uint8,3)})))

#    cdef int x, y
#    cdef Color c
#    for x in range(s.sx):
#        for y in range(s.sy):
#            c = s.pix[x][s.sy - y - 1]
#            #arr[y, x]['r'] = c.r
#            #arr[y, x]['g'] = c.g
#            #arr[y, x]['b'] = c.b
#            arr[y, x] = c.r | c.g << 8 | c.b << 16 | c.a<<24;
    sig_on()
    SurfaceToNumpy (&s, arr)
    sig_off()
    return Image.fromarray(arr, 'RGBA')

cdef Automaton getAutomaton(DetAutomaton a, list A, verb=False):
    cdef int i
    if verb:
        print("getAutomaton %s..." % a)
    cdef DetAutomaton fa
    cdef Automaton aut
    #if isinstance(a, DetAutomaton):
    if set(A).issubset(a.A):
        fa = a.permut(A, verb=verb)
    else:
        fa = a.bigger_alphabet(A)
    aut = fa.a[0]
    # free(fa.a)
    # fa.a = NULL
    aut = CopyAutomaton(aut, aut.n, aut.na);
    return aut
    # else:
    #    raise ValueError("DetAutomaton expected.")


def mahler(pi):
    from sage.rings.qqbar import AA
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.rational_field import RationalField
    QQ = RationalField()
    R = PolynomialRing(QQ, 'x')
    try:
        pi = R(pi)
    except:
        raise ValueError("The argument must be a polynomial over ZZ")
    pi.leading_coefficient()
    pi *= pi.denominator()
    rr = pi.roots(ring=QQbar)
    p = pi.leading_coefficient()
    for r in rr:
        if r[0] not in AA:
            rr.remove((r[0].conjugate(), r[1]))
        a = abs(r[0])
        if a > 1:
            p *= a
    return p


cdef BetaAdic getBetaAdic(m, prec=53, mirror=False, verb=False):
    from sage.rings.complex_field import ComplexField
    CC = ComplexField(prec)
    cdef BetaAdic b
    a = m.a.prune().minimize()
    if mirror:
        a = a.mirror().determinize().minimize()
    A = a.alphabet
    nA = a.n_letters

    b = NewBetaAdic(nA)
    b.b = complex(CC(m.b))
    for i, c in zip(range(b.n), A):
        b.t[i] = complex(CC(c))
    b.a = getAutomaton(a, A=A, verb=verb)
    return b

cdef BetaAdic2 getBetaAdic2(BetaAdicSet self, la=None,
                            prec=53, mirror=False, verb=False):
    if verb:
        print("getBetaAdic %s" % self)
    from sage.rings.complex_field import ComplexField
    CC = ComplexField(prec)
    cdef BetaAdic2 b
    if la is None:
        la = self.get_la(verb=verb)

    # check that every element of la is a DetAutomaton or convert it
    la = [getDetAutomaton(self, a) for a in la]
    
    # add the automaton of self as first element
    la = [self.a]+la

    # simplify each automaton
    la = [a.prune().minimize() for a in la]

    if mirror:
        la = [a.mirror().determinize().minimize() for a in la]

    if verb:
        print("la=%s" % la)

    A = set()
    for a in la:
        A.update(a.alphabet)
    A = list(A)
    if verb:
        print("A=%s" % A)

    b = NewBetaAdic2(len(A), len(la))
    b.b = complex(CC(self.b))
    d = {}
    for i, c in zip(range(b.n), A):
        b.t[i] = complex(CC(c))
        d[c] = i
    # automata
    for i in range(len(la)):
        b.a[i] = getAutomaton(getDetAutomaton(self, la[i]), A=A, verb=verb)
    return b

# used by substitution()
cdef fils(list tree, int e):
    """
    Return the list of leaves's sub-tree  starting at e.

    INPUT:

    - ``tree`` the tree.
    - ``e`` the starting node.

    OUTPUT:

    list of ``tree`` leaves's sub-tree  starting at e.

    """
    if tree[e] == []:  # e is a
        return [e]
    r = []
    for f in tree[e]:
        r += fils(tree, f)
    return r

# used by substitution()
# test if a is included in one of the pieces of l or not
cdef included(DetAutomaton a, list l, list lm):
    """
    Test if the piece described by the automaton a is included in one of the pieces of lm whose indices are in l.

    INPUT:

    - ``a`` DetAutomaton.

    - ``l`` list of indices of pieces to compare with a
    
    - ``lm`` list of pieces with corresponding translations

    OUTPUT:
    
    - ``True`` if a is empty
    
    - int, index of the piece of lm where a is included

    - ``None`` if a is not included in one of the pieces of l

    """
    # quickly test if we are included in a piece or not
    incl = False
    w = a.find_word()
    if w is None:
        if a.has_empty_language():
            print("Error: empty automata!")
        else:
            print("Error: empty word, but non-empty language!")
        return True
    lr = []
    for j in l:
        a2 = lm[j][0]
        if a2.rec_word(w):
            if a.included(a2):
                return j
            else:
                return None
    raise RuntimeError("******* Error: the word %s of a=%s is not in any pieces of the list l=%s ! *********" %(w, a, l))


# split a1 according to ba (returns a couple (a11, a12), with
# a11 in ba and a12 in the complementary)
def split_ba(i, tr, np, lm, m, aa, ap, verb=False):
    b = m.b
    a1 = lm[i][0]
    # TODO : use precomputed relations automata
    # for translations of the domain exchange
    at = m.Proj(aa, ap, t=(b**(-np))*tr)
    if at.intersect(a1):
        ar = at.intersection(a1)
        # determine if we are included in ba
        ar.zero_complete_op()
        if ar.has_same_language_as(a1):
            return (a1, None)
        else:
            # split into two pieces
            ar2 = a1.intersection(ar.complementary())
            ar2.zero_complete_op()
            return (ar, ar2)
    else:
        return (None, a1)


cdef class ImageIn:
    r"""
    This class permits to load an image and test
    if a point is in the image or outside (using transparency).

    INPUT:

    - file_name - The location of the image file.

    EXAMPLE::

        sage: from sage.arith.beta_adic import ImageIn
        sage: ImageIn("SomeImage.png")
        Traceback (most recent call last):
        ...
        IOError: [Errno 2] No such file or directory: 'SomeImage.png'

    """
    cdef numpy.ndarray img

    def __init__(self, file_name):
        import matplotlib.image as mpimg
        self.img = mpimg.imread(file_name)

    def __repr__(self):
        if self.img.ndim < 2:
            raise RuntimeError("the number of dimensions must be at least two")
        return "Image of size %sx%s" % (self.img.shape[1], self.img.shape[0])

    def __contains__(self, p):
        cdef int x,y,w,h
        if self.img.ndim < 2:
            raise RuntimeError("the number of dimensions must be at least two")
        h = self.img.shape[0]
        w = self.img.shape[1]
        from sage.rings.complex_field import ComplexField
        CC = ComplexField(53)
        try:
            p = CC(p)
            x = p.real
            y = p.imag
        except:
            x,y = p
        if x < 0 or y < 0 or x >= w or y >= h:
            return False
        return self.img[y][x][3] > .5

    @property
    def img(self):
        return self.img

    def height(self):
        if self.img.ndim < 1:
            raise RuntimeError("the number of dimensions must be at least one")
        return self.img.shape[0]

    def width(self):
        if self.img.ndim < 2:
            raise RuntimeError("the number of dimensions must be at least two")
        return self.img.shape[1]


def getDetAutomaton(self, a):
    if type(a) is BetaAdicSet:
        if self.b != a.b:
            raise ValueError("The two beta-adic sets must have the same" +
                             "b (here %s != %s).", self.b, a.b)
        a = a.a
    else:
        try:
            a = DetAutomaton(a)
        except Exception:
            raise ValueError("The argument a must be a BetaAdicSet or an automaton.")
    return a


cdef class BetaBase:
    r"""
    The purpose of this class is just to write more conveniently some computations.
    It is used in the computation of a substitution describing a BetaAdicSet.
    """

    def __init__(self, b):
        self.m = BetaAdicSet(b, DetAutomaton(None))
        self.b = self.m.b

    @property
    def b(self):
        return self.b

    def Proj(self, DetAutomaton a, DetAutomaton b, t=0, arel=None, bool only_aut=True):
        self.m.a = a
        return self.m.proj(b, t=-t, arel=arel, aut=only_aut)

    def relations_automaton(self, t=0, bool isvide=False, list Ad=None,
                            list A=None, list B=None,
                            bool couples=False, bool ext=False,
                            bool mirror=False,
                            bool prune=True, int nhash=1000003, int prec=53, int algo=3,
                            int coeff=1, bool verb=False):
        return self.m.relations_automaton(t=-t, isvide=isvide, Ad=Ad, A=A, B=B,
                                          couples=couples, ext=ext,
                                          mirror=mirror,
                                          prune=prune, nhash=nhash, prec=prec,
                                          algo=algo, coeff=coeff, verb=verb)

cdef getBetaAdicSet(BetaAdicSet self, a):
    if type(a) is BetaAdicSet:
        if self.b != a.b:
            raise ValueError("The two beta-adic sets must have the same b (here %s != %s).", self.b, a.b)
    elif type(a) is not DetAutomaton:
        try:
            a = DetAutomaton(a)
        except Exception:
            raise ValueError("The argument a must be a BetaAdicSet or an automaton.")
        a = BetaAdicSet(self.b, a)
    return a

cdef class BetaAdicSet:
    r"""
    Define a numeration in base b, i.e. set of numbers of the form

        :math:`\sum_{i=0}^\infty \beta^i c_i`

    where :math:`\beta` is an element of a field (for example a complex number),
    and the :math:`c_i` form a word recognized by a deterministic automaton ``a``.

    INPUT:

    - ``b`` -- number, base of the numeration.

    - ``a`` -- DetAutomaton, giving the allowed sequence of digits.


    EXAMPLES::

        sage: from sage.combinat.words.cautomata_generators import dag
        sage: m1 = BetaAdicSet(3, dag.AnyWord([0, 1, 3]))
        sage: print(m1)
        b-adic set with b root of x - 3, and an automaton of 1 state and 3 letters
        sage: m2 = BetaAdicSet((1 + sqrt(5)) / 2, dag.AnyWord([0, 1]))
        sage: print(m2)
        b-adic set with b root of x^2 - x - 1, and an automaton of 1 state and 2 letters
        sage: b = (x^3-x-1).roots(ring=QQbar)[0][0]
        sage: m3 = BetaAdicSet(b, dag.AnyWord([0, 1]))
        sage: print(m3)
        b-adic set with b root of x^3 - x - 1, and an automaton of 1 state and 2 letters

    """
    def __init__(self, b, a):
        r"""
        Construction of the b-adic with base ``b`` and automaton ``a``.

        EXAMPLES::

            sage: m1 = BetaAdicSet(3, dag.AnyWord([0, 1, 3]))
            sage: m1
            b-adic set with b root of x - 3, and an automaton of 1 state and 3 letters
            sage: c = Automaton({0:{1:'0',2:'1',3:'2'}, 2:{5:'1'}},initial_states=[0])
            sage: m3 = BetaAdicSet(m1.b, c)
            sage: m3
            b-adic set with b root of x - 3, and an automaton of 5 states and 3 letters
            sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
            sage: m
            b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 1 state and 2 letters
            sage: m1 = BetaAdicSet(3,[0,1])
            sage: m1
            b-adic set with b root of x - 3, and an automaton of 1 state and 2 letters

        """
        cdef int i, j
        from sage.rings.complex_field import ComplexField
        CC = ComplexField()
        if b not in CC:
            # raise ValueError("b must be a number.")
            from sage.rings.qqbar import QQ
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            K = PolynomialRing(QQ, 'x')
            try:
                pi = K(b)
                rr = [r[0] for r in pi.roots(ring=QQbar)]
                rrm = [r for r in rr if abs(r) < 1]
                if len(rrm) > 0:
                    b = rrm[0]
                else:
                    b = rr[0]
            except Exception:
                raise ValueError("b must be a number, or a polynomial over QQ")
        try:
            b = QQbar(b)
            pi = QQbar(b).minpoly()
            K = NumberField(pi, 'b', embedding=b)
            self.b = K.gen()
        except Exception:
            self.b = b

        if type(a) != DetAutomaton:
            try:
                a = DetAutomaton(a)
            except Exception:
                try:
                    a = list(a)
                except Exception:
                    raise ValueError("a must be an automaton or an iterable.")
                from sage.combinat.words.cautomata_generators import dag
                a = dag.AnyWord(a)
        self.a = a

        # test if letters of a are in K
        try:
            K = self.b.parent()
            self.a.A = [K(c) for c in self.a.A]
        except Exception:
            raise ValueError("Alphabet %s of the automaton is not in the field %s of b !"%(self.a.A, self.b.parent()))

    def __repr__(self):
        r"""
        Returns the string representation of the BetaAdicSet.

        EXAMPLES::

            sage: from sage.combinat.words.cautomata_generators import dag
            sage: BetaAdicSet((1+sqrt(5))/2, dag.AnyWord([0, 1]))
            b-adic set with b root of x^2 - x - 1, and an automaton of 1 state and 2 letters
            sage: BetaAdicSet(3, dag.AnyWord([0, 1, 3]))
            b-adic set with b root of x - 3, and an automaton of 1 state and 3 letters


        TESTS::

            sage: from sage.combinat.words.cautomata_generators import dag
            sage: m = BetaAdicSet(3/sqrt(2), dag.AnyWord([0, 1]))
            sage: repr(m)
            'b-adic set with b root of x^2 - 9/2, and an automaton of 1 state and 2 letters'

        """

        from sage.rings.qqbar import QQbar
        if self.b not in QQbar:
            str = "(%s)-adic set with an "%self.b
        else:
            K = self.b.parent()
            from sage.rings.rational_field import QQ
            if K.base_field() == QQ:
                str = "b-adic set with b root of %s, and an "%self.b.minpoly()
            else:
                if K.characteristic() != 0:
                    str = "b-adic set with b root of %s (in characteristic %s), and an "%(self.b.minpoly(), K.characteristic())
                else:
                    str = "b-adic set with b root of %s, and an "%K.modulus()
        str += "automaton of %s state"%self.a.a.n
        if self.a.a.n > 1:
            str += 's'
        str += " and %s letter" % (self.a.a.na)
        if self.a.a.na > 1:
            str += 's'
        return str

    def string(self):
        r"""
        Return a string that can be evaluated to recover the BetaAdicSet

        OUTPUT:
        Return a string to define a BetaAdicSet, this set can be obtained by the ``use_draw`` method

        EXAMPLES::

            sage: m1 = BetaAdicSet(3,[0,1])
            sage: m1.string()
            'BetaAdicSet((x - 3).roots(ring=QQbar)[0][0], DetAutomaton([[0], [(0, 0, 0), (0, 0, 1)]], A=[0, 1], i=0, final_states=[0]))'
            sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
            sage: m.string()
            'BetaAdicSet((x^3 - x^2 - x - 1).roots(ring=QQbar)[1][0], DetAutomaton([[0], [(0, 0, 0), (0, 0, 1)]], A=[0, 1], i=0, final_states=[0]))'

        """
        pi = self.b.minpoly()
        from sage.rings.qqbar import QQbar
        rr = pi.roots(ring=QQbar)
        for i, r in enumerate(rr):
            if r[0] == self.b:
                break
        from sage.rings.rational_field import QQ
        if len([c for c in self.a.A if c not in QQ]) == 0:
            return "BetaAdicSet((%s).roots(ring=QQbar)[%s][0], %s)"%(pi, i, self.a.string())
        else:
            return "m = BetaAdicSet((%s).roots(ring=QQbar)[%s][0], DetAutomaton(None))\nb=m.b\nBetaAdicSet(b, %s)"%(pi, i, self.a.string())

    @property
    def a(self):
        """
        Get the ``DetAutomaton`` ``a`` of the ``BetaAdicSet``

        OUTPUT:

        ``DetAutomaton`` ``a`` attribut

        EXAMPLES::

            sage: from sage.combinat.words.cautomata_generators import dag
            sage: m = BetaAdicSet((1+sqrt(5))/2, dag.AnyWord([0, 1]))
            sage: m.a
            DetAutomaton with 1 state and an alphabet of 2 letters

        """
        return self.a

    @property
    def b(self):
        """
        Get the number ``b`` of the ``BetaAdicSet``

        OUTPUT:

        number ``b`` attribut

        EXAMPLES::

            sage: m = BetaAdicSet((1+sqrt(5))/2, dag.AnyWord([0, 1]))
            sage: m.b
            b

        """
        return self.b

    def copy(self):
        """
        return a copy of  the ``BetaAdicSet``

        OUTPUT:

        a ``BetaAdicSet``

        EXAMPLES::

            sage: m = BetaAdicSet((1+sqrt(5))/2, [0, 1])
            sage: m.copy()
            b-adic set with b root of x^2 - x - 1, and an automaton of 1 state and 2 letters

        TESTS::
            
            sage: m = BetaAdicSet((1+sqrt(5))/2, dag.AnyLetter([0,1]))
            sage: m2 = m.copy()
            sage: m.a.set_final(0)
            sage: m.a == m2.a
            False

        """

        return BetaAdicSet(self.b, self.a.copy())

    def mirror(self):
        """
        Return the beta-adic set with the mirror automaton.

        OUTPUT:

        A ``BetaAdicSet`` with the mirror automaton as attribut ``a``

        EXAMPLES::

            sage: m = BetaAdicSet((1+sqrt(5))/2, [0, 1])
            sage: m.mirror()
            b-adic set with b root of x^2 - x - 1, and an automaton of 1 state and 2 letters

        """
        return BetaAdicSet(self.b, self.a.mirror())

    def is_included(self, a, verb=False):
        """
        Determine if the BetaAdicSet is included in the BetaAdicSet given by a.

        INPUT:

        - ``a`` - ``BetaAdicSet`` to compare
        - ``verb`` - Boolean (default: False) Display informations for debug.

        OUTPUT:

        ``True``  if the BetaAdicSet is included in the BetaAdicSet given
        by a  ``False`` otherwise


        EXAMPLES::

            sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
            sage: m1 = BetaAdicSet(x^3-x^2-x-1, [0,1,2])
            sage: m1.is_included(m)
            False
            sage: m.is_included(m1)
            True

        """
        a = getDetAutomaton(self, a)
        if verb:
            print("a=%s"%a)
        b = self.a.concat_zero_star()
        b.zero_complete_op()
        if verb:
            print("b=%s"%b)
        m = BetaAdicSet(self.b, a)
        ap = m.proj(b, aut=True)
        if verb:
            print("ap=%s" % ap)
        return b.included(ap)
        # return ap.has_same_language_as(b)

    def is_equal_to(self, a):
        """
        Determine if the ``BetaAdicSet`` is equal to the given ``BetaAdicSet``.

        INPUT:

        - ``a`` - ``BetaAdicSet`` to compare

        OUTPUT:

        ``True``  if the BetaAdicSet is equal in the BetaAdicSet given
        by a  ``False`` otherwise
 

        EXAMPLES::

            sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
            sage: m1 = BetaAdicSet(x^3-x^2-x-1, [0,1,2])
            sage: m1.is_equal_to(m)
            False


        """
        a = getBetaAdicSet(self, a)
        return self.is_included(a) and a.is_included(self)

    def is_empty(self):
        """
        Tell if the BetaAdicSet is empty.

        OUTPUT:

        ``True``  if the BetaAdicSet is empty in the BetaAdicSet given
        by a  ``False`` otherwise
 

        EXAMPLES::

            sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
            sage: m.is_empty()
            False
            sage: m = BetaAdicSet(3, [])
            sage: m.is_empty()
            True

        TESTS::

            sage: m = BetaAdicSet(1/(1+I), DetAutomaton([(0,0,0),(0,1,1)], i=0, final_states=[]))
            sage: m.is_empty()
            True
            sage: m = BetaAdicSet(1/(1+I), DetAutomaton([(0,0,0),(0,1,1)], i=0, final_states=[1]))
            sage: m.is_empty()
            False

        """
        return self.a.has_empty_language()

#    def _testSDL(self):
#        """
#        Open a window to test the SDL library used for graphical representation.
#
#        TESTS::
#
#            sage: m3 = BetaAdicSet(1/(1+I), dag.AnyWord([0, 1]))
#            sage: m3._testSDL()
#            Video Mode: 800x600 32 bits/pixel
#        """
#        sig_on()
#        TestSDL()
#        sig_off()

    def get_la(self, bool verb=False):
        """
        Return a list of automata corresponding to each final state of the automaton.
        For each state of self, give a copy of self but whose set of final states is this state.

        INPUT:

        -``verb`` -- bool (default ''False'') - set to ''True'' for verbose mode

        OUTPUT:
        Return a list of automata.

        EXAMPLES::

            sage: m=BetaAdicSet((1+sqrt(5))/2, dag.AnyWord([0, 1]))
            sage: m.get_la()
            [DetAutomaton with 1 state and an alphabet of 2 letters]
            
         #. plot a Rauzy fractal
            sage: m=WordMorphism('a->ab,b->ac,c->a').DumontThomas()
            sage: la = m.get_la()
            sage: la
            [DetAutomaton with 3 states and an alphabet of 2 letters,
             DetAutomaton with 3 states and an alphabet of 2 letters,
             DetAutomaton with 3 states and an alphabet of 2 letters]
            sage: m.plot_list(la)       #random
        """
        cdef DetAutomaton a = self.a.copy()
        # compute la
        la = []
        for v in range(a.a.n):
            a.set_final_states([v])
            la.append(a.copy())
        return la

#    def points_exact(self, n=None, i=None):
#        r"""
#        Returns a set of exacts values (in the number field of b)
#        corresponding to points of the b-adic set for words of length at most ``n``.
#
#        INPUT:
#
#        - ``n`` - integer (default: ``None``)
#          The number of iterations used to plot the fractal.
#          Default values: between ``5`` and ``16`` depending on the number
#          of generators.
#
#        - ``i`` - integer (default: ``None``)
#          State of the automaton of self taken as the initial state .
#
#        OUTPUT:
#
#            List of numbers, given with exact values.
#
#        EXAMPLES::
#
#            #. The dragon fractal::
#            sage: e = QQbar(1/(1+I))
#            sage: m=BetaAdicSet(e, dag.AnyWord([0, 1]))
#            sage: print(m)
#            b-adic set with b root of x^2 - x + 1/2, and an automaton of 1 state and 2 letters
#            sage: P = m.points_exact()
#            age: len(P)
#            65536
#            sage: P = m.points_exact(i=0)
#            sage: len(P)
#            65536
#        """
#        K = self.K
#        b = self.b
#        a = self.a
#        A = a.alphabet
#        ng = a.n_letters
#
#        if i is None:
#            i = a.initial_state
#
#        if n is None:
#            if ng == 2:
#                n = 16
#            elif ng == 3:
#                n = 9
#            else:
#                n = 5
#
#        if n == 0:
#            return [0]
#        else:
#            orbit_points = set()
#            V = set([v for c in A for v in [a.succ(i, c)] if v != -1])
#            orbit_points0 = dict()
#            for v in V:
#                orbit_points0[v] = self.points_exact(n=n-1, i=v)
#            for c in A:
#                v = a.succ(i, c)
#                if v is not None:
#                    orbit_points.update([b*p+c for p in orbit_points0[v]])
#        return orbit_points

    def user_draw(self, n=None,
                  int sx=800, int sy=600, bool ajust=True, int prec=53, color=(0, 0, 0, 255),
                  bool simplify=True, bool mirror=False, bool only_aut=False, bool verb=False):
        r"""
        Display a window where the user can draw a b-adic set based on the current b-adic set.
        Use keyboard p to reduce the size of the pen and the keyboard m to increse.
        Draw the figure with the the mouse and click to record the shape.

        INPUT:

        - ``n`` - integer (default: ``None``)
          The number of iterations used to plot the fractal.
          Default values: between ``5`` and ``16`` depending on the number
          of generators.

        - ``sx`` -- integer (default: ``800``) - width of the window

        - ``sy`` -- integer (default: ``600``) - height of the window

        - ``ajust``  -- boolean (default ``True``) - If True, change the zoom in order to fit the window.

        - ``prec`` -- integer (default: ``53``) - precision of computed values

        - ``color`` -- tuple (default: (0, 0, 0, 255)) - color in RGBA values

        - ``simplify`` -- (default: ``True``) - If True, minimize the result

        - ``only_aut`` -- (default: ``False``) - If True return a DetAutomaton, otherwise return a BetaAdicSet

        - ``verb`` -- (default ``False``) - set to ``True`` for verbose mod

        OUTPUT:

        A b-adic set, corresponding to what has been drawn by the user. Or only the automaton if only_aut was True.

        EXAMPLES::

            #. Draw a BetaAdicSet from the dragon fractal::

                sage: m = BetaAdicSet(1/(1+I), [0, 1])
                sage: P = m.user_draw()     # not tested (need the intervention of the user)
                sage: P.string()            # not tested

            #. Draw a BetaAdicSet from a Rauzy fractal::

                sage: m = WordMorphism('a->ab,b->ac,c->a').DumontThomas().mirror()
                sage: P = m.user_draw()     # not tested (need the intervention of the user)
                sage: P.plot()              # not tested

        """
        cdef BetaAdic b
        cdef Automaton a
        cdef DetAutomaton r
        b = getBetaAdic(self, prec=prec, mirror=mirror, verb=verb)
        # if verb:
        #    printAutomaton(b.a)
        # dessin
        cdef Color col
        col.r = color[0]
        col.g = color[1]
        col.b = color[2]
        col.a = color[3]
        if n is None:
            n = -1
        try:
            spr = self.a.prune().spectral_radius()
            sig_on()
            a = UserDraw(b, sx, sy, n, ajust, col, spr, verb)
            sig_off()
        except:
            raise NotImplementedError("This function does not exists in this version of Sage, because it uses the library SDL2. To use this function, you need to install the library SDL2 on your system and then recompile Sage from the sources, where you include the ticket https://trac.sagemath.org/ticket/21072.")
        r = DetAutomaton(None)
        r.a[0] = a
        r.A = self.a.A
        r.S = range(a.n)
        if simplify:
            r = r.minimize()
        if only_aut:
            return r
        else:
            return BetaAdicSet(self.b, r)

    def draw_zoom(self, n=None, int sx=800, int sy=600,
                  bool ajust=True, int prec=53, color=(0, 0, 0, 255),
                  int nprec=4, bool mirror=False, bool verb=False):
        r"""
        Display the BetaAdicSet in a window, with possibility for the user to zoom in.
        Use 'p' to zoom in, 'm' to zoom out, the arrows to translate the view, and 'Esc' to quit.
        You can also select a zone to zoom in with the mouse.

        INPUT:

        - ``n`` - integer (default: ``None``)
          The number of iterations used to plot the fractal.
          Default values: between ``5`` and ``16`` depending on the number
          of generators.

        - ``sx``  -- (default 800)

        - ``sy``  -- (default 600)

        - ``ajust``  -- (default ``True``) If ``True``, change the zoom in order to fit the window.

        - ``prec``  precision of computed values -- (default: ``53``)

        - ``color`` tuple of color in RGB values -- (default: (0, 0, 0, 255))

        - ``nprec`` int -- (default 4) - additional iterations for the drawing (if ``n`` is None).

        - ``mirror`` bool -- (default ``False) set to ``True`` to use the mirror of the automaton

        - ``verb`` -- (default ``False``) set to ``True`` for verbose mod

        OUTPUT:

        A word that corresponds to the place where we draw.

        EXAMPLES::

            #. The dragon fractal::

                sage: m = BetaAdicSet(1/(1+I), [0, 1])
                sage: w = m.draw_zoom()     # not tested (need the intervention of the user)


            #. Zoom in a complicated Rauzy fractal

                sage: s = WordMorphism('1->2,2->3,3->12')
                sage: m = s.DumontThomas().mirror(); m
                b-adic set with b root of x^3 - x - 1, and an automaton of 4 states and 2 letters
                sage: m.draw_zoom()         # not tested (need the intervention of the user)

        """
        cdef BetaAdic b
        b = getBetaAdic(self, prec=prec, mirror=mirror, verb=verb)
        # dessin
        cdef int *word
        cdef Color col
        cdef int i
        col.r = color[0]
        col.g = color[1]
        col.b = color[2]
        col.a = color[3]
        if n is None:
            n = -1
        try:
            spr = self.a.prune().spectral_radius()
            sig_on()
            word = DrawZoom(b, sx, sy, n, ajust, col, nprec, spr, verb)
            sig_off()
        except:
            raise NotImplementedError("This function does not exists in this version of Sage, because it uses the library SDL2. To use this function, you need to install the library SDL2 on your system and then recompile Sage from the sources where you include the ticket https://trac.sagemath.org/ticket/21072.")
        res = []
        if word is not NULL:
            for i in xrange(1024):
                if word[i] < 0:
                    break
                res.append(self.a.alphabet[word[i]])
            res.reverse()
        return res

    def plot(self, n=None, int sx=800, int sy=600,
             bool ajust=True, int prec=53, color=(0, 0, 0, 255),
             int nprec=4, bool mirror=False, bool verb=False):
        r"""
        Draw the beta-adic set. see processed examples on 

        INPUT:

        - ``n`` - integer (default: ``None``)
          The number of iterations used to plot the fractal.
          Default values: between ``5`` and ``16`` depending on the
          number of generators.

        - ``place`` - place of the number field of beta (default: ``None``)
          The place used to evaluate elements of the number field.

        - ``sx`` -- int (default: 800) - dimensions of the resulting in x dimension

        - ``sy`` -- int (default : 600) - dimensions of the resulting
          in y dimension image

        - ``ajust`` -- bool (default: ``True``) - adapt the drawing
          to fill all the image, with ratio 1 

        - ``prec`` - int (default: ``53``) - precision of returned values

        - ``color`` - list of four integers between 0
          and 255 (RGBA format, default: ``(0,0,0,255)``) Color of the drawing.

        - ``mirror`` bool -- (default ``False``) - set to ``True`` to use the mirror of the automaton

        - ``nprec`` int -- (default 4) - additionnal iterations (if n is ``None``)

        - ``verb`` - bool (default: ``False``)
          Print informations for debugging.

        OUTPUT:

            A Graphics object.

        EXAMPLES::

            #. The dragon fractal::

                sage: m = BetaAdicSet(1/(1+I), dag.AnyWord([0,1]))
                sage: m.plot()                                      # random

            #. Another dragon fractal::

                sage: m = BetaAdicSet(2*x^2+x+1, dag.AnyWord([0,1]))
                sage: m.plot()                                      # random

            #. The Rauzy fractal of the Tribonacci substitution::

                sage: s = WordMorphism('1->12,2->13,3->1')
                sage: m = s.DumontThomas().mirror()
                sage: m.plot()                                      # random

            #. The Rauzy fractal of the flipped Tribonacci substitution::

                sage: s = WordMorphism('1->12,2->31,3->1')
                sage: m = s.DumontThomas().mirror()
                sage: m.plot()                                      # random

            #. A non-Pisot Rauzy fractal::

                sage: s = WordMorphism({1:[3,2], 2:[3,3], 3:[4], 4:[1]})
                sage: m = s.DumontThomas().mirror()
                sage: m.plot()                                      # random
                sage: m = BetaAdicSet(1/m.b, m.a)
                sage: m.plot()                                      # random

            #. A part of the boundary of the dragon fractal::

                sage: m = BetaAdicSet(1/(1+I), dag.AnyWord([0,1]))
                sage: mi = m.intersection_words([0], [1])
                sage: mi.plot(nprec=6)                              # random

            #. A part of the boundary of the "Hokkaido" fractal::

                sage: s = WordMorphism('a->ab,b->c,c->d,d->e,e->a')
                sage: m = s.DumontThomas().mirror()
                sage: mi = m.intersection_words([0], [1])
                sage: mi.plot()                                     # random

            #. A limit set that look like a tiling but with holes::

                sage: m = BetaAdicSet(x^4 + x^3 - x + 1, [0,1])
                sage: m.plot()                                      # random

        """
        cdef Surface s
        cdef BetaAdic b
        cdef Automaton aut
        cdef int i, j
        sig_on()
        s = NewSurface(sx, sy)
        sig_off()
        sig_on()
        b = getBetaAdic(self, prec=prec, mirror=mirror, verb=verb)
        sig_off()
        if verb:
            print("b=%s+%s*I", b.b.x, b.b.y)
            print("n=%s" % b.n)
            for i in range(b.n):
                print("t[%s] = %s+%s*I" % (i, b.t[i].x, b.t[i].y))
            # print("a=%s"%b.a)
            for i in range(b.a.n):
                if b.a.e[i].final:
                    print("(%s) " % i)
                else:
                    print("%s " % i)
            aut = b.a;
            for i in range(aut.n):
                for j in range(aut.na):
                    print("%s -%s-> %s\n" % (i, j, aut.e[i].f[j]))
        cdef Color col
        col.r = color[0]
        col.g = color[1]
        col.b = color[2]
        col.a = color[3]
        if n is None:
            n = -1
        spr = self.a.prune().spectral_radius()
        sig_on()
        Draw(b, s, n, ajust, col, nprec, spr, verb)
        sig_off()
        sig_on()
        im = surface_to_img(s)
        sig_off()
        if verb:
            print("Free...")
        sig_on()
        FreeSurface(s)
        FreeBetaAdic(b)
        sig_off()
        return im

    def plot_list(self, list la=None, n=None,
                  int sx=800, int sy=600, bool ajust=True, int prec=53, colormap='hsv',
                  backcolor=None, float opacity=1., bool mirror=False,
                  int nprec=4, bool verb=False):
        r"""
        Draw the beta-adic set self, with color according to the list of automata or BetaAdicSets given.

        INPUT:

        - ``la``- list (default: ``None``)
          List of automata or BetaAdicSet to plot.

        - ``n`` - integer (default: ``None``)
          The number of iterations used to plot the fractal.

        - ``sx`` -- int (default: 800) - width of the result image

        - ``sy`` -- int (default : 600) - height of the result image

        - ``ajust`` -- bool (default: ``True``) - adapt the drawing to fill all the image, with
          ratio 1 (default: ``True``)

        - ``prec`` - precision of returned values (default: ``53``)

        - ``colormap`` - list of colors (default: ``hsv``)
          Colors of the drawing.

        - ``backcolor`` - (default: ``None``) list of four integers between 0
          and 255  .

        - ``opacity`` -- float (default: ``1.``)
          Transparency of the drawing coefficient.

        - ``mirror`` -- bool (default ``False) set to ``True`` to use the mirror of the automaton

        - ``nprec`` -- int (default 4) - additionnal iterations

        - ``verb`` -- bool (default: ``False``)
          Print informations for debugging.

        OUTPUT:

            A Graphics object.

        EXAMPLES::

            #. The Rauzy fractal of the Tribonacci substitution::

                sage: s = WordMorphism('1->12,2->13,3->1')
                sage: m = s.DumontThomas()
                sage: m.plot_list(mirror=True)  # random

            #. A non-Pisot Rauzy fractal::

                sage: s = WordMorphism({1:[3,2], 2:[3,3], 3:[4], 4:[1]})
                sage: m = s.DumontThomas()
                sage: m = BetaAdicSet(1/m.b, m.a)
                sage: m.plot_list(mirror=True)             # random
                sage: m = BetaAdicSet(m.b, m.a.mirror())
                sage: m.plot_list(mirror=True)             # random

            #. The dragon fractal and its boundary::

                sage: m = BetaAdicSet(1/(1+I), [0,1])
                sage: mi = m.intersection_words([0], [1])
                sage: m.plot_list([mi], n=19, colormap=[(.5,.5,.5,.5), (0,0,0,1.)])  # random

            #. The "Hokkaido" fractal and its boundary::

                sage: s = WordMorphism('a->ab,b->c,c->d,d->e,e->a')
                sage: m = s.DumontThomas().mirror()
                sage: mi = m.intersection_words([0], [1])                    # long time
                sage: m.plot_list([mi], colormap='gist_rainbow')             # not tested

            #. A limit set that look like a tiling::

                sage: m = BetaAdicSet(x^4 + x^3 - x + 1, [0,1])
                sage: m = m.reduced().mirror()
                sage: m.plot_list(mirror=True)                 # random

            #. Plot a domain exchange computed from a BetaAdicSet

                sage: s = WordMorphism('a->ab,b->c,c->d,d->e,e->a')
                sage: m = s.DumontThomas().mirror()
                sage: la = m.domain_exchange()              # long time
                sage: m.plot_list([a for t,a in la])        # not tested
        """
        cdef Surface s = NewSurface(sx, sy)
        cdef BetaAdic2 b
        sig_on()
        b = getBetaAdic2(self, la=la, prec=prec, mirror=mirror, verb=verb)
        sig_off()
        # dessin
        if n is None:
            n = -1

        # Manage colors
        if backcolor is None:
            backcolor = (.5, .5, .5, .5)
        cdef ColorList cl
        sig_on()
        cl = NewColorList(b.na)
        sig_off()
        if isinstance(colormap, list):
            # if b.na > len(colormap):
            #    raise ValueError("The list of color must contain at least %d elements."%b.na)
            for i in range(b.na):
                if i < len(colormap):
                    cl[i] = getColor(colormap[i])
                else:
                    cl[i] = randColor(255)
                sig_check()
        elif isinstance(colormap, str):
            from matplotlib import cm
            if not colormap in cm.datad.keys():
                raise ValueError("Color map %s not known (type 'from matplotlib import cm' and look at cm for valid names)" % colormap)
            colormap = cm.__dict__[colormap]
            cl[0] = getColor(backcolor)
            for i in range(b.na-1):
                cl[i+1] = getColor(colormap(float(i)/float(b.na-1)))
                sig_check()
        else:
            raise TypeError("Type of option colormap (=%s) must be list of colors or str" % colormap)
        spr = self.a.prune().spectral_radius()
        sig_on()
        DrawList(b, s, n, ajust, cl, opacity, spr, nprec, verb)
        sig_off()
        # enregistrement du résultat
        sig_on()
        im = surface_to_img(s)
        sig_off()
        if verb:
            print("Free...")
        sig_on()
        FreeSurface(s)
        if la is None:
            FreeAutomatons(b.a, b.na)
        FreeBetaAdic2(b)
        FreeColorList(cl)
        sig_off()
        return im

    def relations_automaton(self, t=0, bool isvide=False, list Ad=None, list A=None, list B=None,
                             bool couples=False, bool ext=False, bool mirror=False,
                             bool prune=True, int nhash=1000003, int prec=53, int algo=3, int coeff=1, bool verb=False):
        r"""
        Assume that beta is an algebraic integer.
        Compute the relation automaton of the beta-adic set
        (also called "zero automaton").
        It is the minimal deterministic automaton that recognizes
        the set of words a_0 a_1 ... a_n in Ad^* such that
        :math:`a_0 + beta*a_1 + ... + beta^n*a_n = 0`.
        If couples is True, then it describes the set of words over AxB
        :math:`(a_0, b_0) (a_1, b_1) ... (a_n, b_n)` such that
        :math:`a_0 + beta*a_1 + ... + beta^n*a_n = b_0 + beta*b_1 + ... + beta^n*b_n`.

        If ext is True, it describes the set of words that can be prolongated to an infinite relation
        in the contracting space (which is the product of copies of R, C and p-adic spaces corresponding to
        places of the number field for which beta has an absolute value less than one).

         INPUT:

        - ``t`` integer (default: 0) the translation of one of the side

        - ``isvide`` boolean - (default: ``False``) If isvide is True,
          it only checks if the automaton is trivial or not.

        - ``Ad`` - list (default: ``None``)
          Alphabet of differences A-B where A and B
          are the alphabets to compare.

        - ``A`` -  (default: ``None``) alphabet on one side
          (used if Ad is None)

        - ``B`` -  (default: ``None``) alphabet on the other side
          (used if Ad is None)

        - ``couples``  boolean - (default: ``False``) If ``True``, the alphabet of the resulting automaton is AxB. If ``False``, it is Ad (=A-B).

        - ``ext``  boolean - (default: ``False``)
          If ``True``, compute the automaton that describes infinite relations.

        - ``mirror``  boolean - (default: ``False``) If ``True``, return the mirror.

        - ``prune`` boolean - (default: ``True``) Prune the result or not.

        - ``nhash`` int (default: 1000003) Size of the hash table (only for algo 2).

        - ``prec`` int - (default:53)

        - ``algo`` int - (default: 3) Algorithm used (choose in the set {1,2,3}).

        - ``verb`` bool - (default: ``False``)
          Print informations for debugging.

        OUTPUT:

        A DetAutomaton whose language describe the set of relations.

        EXAMPLES::

            sage: m = BetaAdicSet(1/(1+I), [0,1,3])
            sage: m.relations_automaton()
            DetAutomaton with 49 states and an alphabet of 7 letters

            sage: m = BetaAdicSet(1/(1+I), [0,1])
            sage: m.relations_automaton()
            DetAutomaton with 1 state and an alphabet of 3 letters
            sage: m.relations_automaton(ext=True)
            DetAutomaton with 7 states and an alphabet of 3 letters
            sage: m.plot()          #random

        TESTS::

            sage: m = BetaAdicSet(x^3-x-1, [0,1])
            sage: a1 = m.relations_automaton(algo=1)
            sage: a2 = m.relations_automaton(algo=2)
            sage: a3 = m.relations_automaton(algo=3)
            sage: a1.has_same_language_as(a2)
            True
            sage: a2.has_same_language_as(a3)
            True

            sage: m = BetaAdicSet(1/pi, [0,1])
            sage: m.relations_automaton()
            Traceback (most recent call last):
            ...
            ValueError: b must live in a number field!


        """
        cdef InfoBetaAdic ib
        cdef Automaton a
        cdef Element e
        cdef DetAutomaton r
        cdef bool tb

        t0 = t
        if mirror is not None:
            try:
                tb = mirror
            except Exception:
                raise ValueError("mirror=%s must be a bool."%mirror)

        b = self.b
        K = b.parent()
        if not K.is_field():
            raise ValueError("b must live in a field!")
        if not K.is_exact() or not hasattr(K, 'abs_val'):
            raise ValueError("b must live in a number field!")
        pi = b.minpoly()
        pi = pi*pi.denominator()
        # alphabet
        if Ad is None:
            if A is None:
                A = self.a.A
            if B is None:
                B = self.a.A
            Ad = list(set([a1-b1 for a1 in A for b1 in B]))
        else:
            try:
                list(Ad[0])
                Ad = list(set([a1-b1 for a1, b1 in Ad]))
            except Exception:
                pass
        if verb:
            print("Ad=%s" % Ad)
        if ext:
            if algo == 1 and not mirror:
                if verb:
                    print("Algo 1 cannot be used with ext=True and mirror=False: change to algo 3.")
                algo = 3
            elif (algo == 3 or algo == 2) and mirror:
                if verb:
                    print("Algo 2 or 3 cannot be used with ext=True and mirror=True: change to algo 1.")
                algo = 1
        if verb:
            print("algo=%s, mirror=%s" % (algo, mirror))
        if algo == 1:
            if ext:
                b = 1/b
                pi = b.minpoly()
                pi = pi*pi.denominator()
                mirror = not mirror
            # find absolute values for which b is greater than one
            places = []
            narch = 0
            # archimedian places
            for p in K.places(prec=prec):
                if K.abs_val(p, b) > 1:
                    places.append(p)
                    narch+=1
            # ultra-metric places
            from sage.arith.misc import prime_divisors
            lc = pi.leading_coefficient()
            for p in prime_divisors(lc):
                for P in K.primes_above(p):
                    if K.abs_val(P, b, prec=prec) > 1:
                        places.append(P)
            if verb:
                print(places)
            # bounds
            bo = []
            for i, p in enumerate(places):
                if i < narch:
                    bo.append(
                        coeff*max(
                            [K.abs_val(p, x) for x in Ad])/(K.abs_val(p, b) - 1))
                else:
                    bo.append(
                        coeff*max(
                            [K.abs_val(p, x) for x in Ad])/K.abs_val(p, b))
            if verb:
                print("bounds=%s" % bo)
            # compute the automaton
            L = []
            S = [t]  # remaining state to look at
            d = dict()  # states already seen and their number
            d[t] = 0
            c = 1  # count the states seen
            while len(S) > 0:
                S2 = []
                for s in S:
                    for t in Ad:
                        ss = b*s + t
                        # test if we keep ss
                        keep = True
                        for p, m in zip(places, bo):
                            if K.abs_val(p, ss) > m + .00000001:
                                keep = False
                                break
                        if keep:
                            if not d.has_key(ss):
                                S.append(ss)
                                d[ss] = c
                                c += 1
                            L.append((d[s], d[ss], t))
                S = S2
            if d.has_key(0):
                r = DetAutomaton(L, A=Ad, i=d[0], final_states=[0])
            else:
                r = DetAutomaton(L, final_states=[0])
            if verb:
                print("before pruning: %s" % r)
            if not mirror:
                r = r.mirror_det()
            if prune:
                if verb:
                    print("prune...")
                if ext:
                    r = r.prune_inf()
                else:
                    r = r.prune()
            if ext:
                r.set_final_states(range(r.a.n))
        elif algo == 2:
            sig_on()
            ib = initInfoBetaAdic(self, Ad=Ad, plus=False, nhash=nhash, verb=verb)
            e = NewElement(ib.n)
            K = self.b.parent()
            t = K(t)
            getElement(t, e, ib.n)
            a = RelationsAutomatonT(&ib, e, isvide, ext, verb)
            r = DetAutomaton(None)
            r.a[0] = a
            if verb:
                print("a (%s etats)" % a.n)
                print("Free element...")
            FreeElement(e)
            r.A = Ad
            if verb:
                print("Free InfoBetaAdic...")
            freeInfoBetaAdic(&ib)
            sig_off()
            if isvide:
                return a.na != 0
            if prune:
                if verb:
                    print("prune...")
                if ext:
                    r = r.prune_inf()
                    r.set_final_states(r.states)
                else:
                    r = r.prune()
            if mirror:
                r = r.mirror_det()
        else:
            # find absolute values for which b is less than one
            places = []
            narch = 0
            # archimedian places
            for p in K.places(prec=prec):
                if K.abs_val(p, b) < 1:
                    places.append(p)
                    narch+=1
            # ultra-metric places
            from sage.arith.misc import prime_divisors
            for p in prime_divisors(pi(0)):
                for P in K.primes_above(p):
                    if K.abs_val(P, b, prec=prec) < 1:
                        places.append(P)
            if verb:
                print(places)
            # bounds
            bo = []
            for i, p in enumerate(places):
                if i < narch:
                    bo.append(
                        coeff*max(
                            [K.abs_val(p, x) for x in Ad])/(1 - K.abs_val(p, b)))
                else:
                    bo.append(
                        coeff*max([K.abs_val(p,x) for x in Ad]))
            if verb:
                print("bounds=%s" % bo)
            # compute the automaton
            L = []
            S = [t]  # remaining state to look at
            d = dict()  # states already seen and their number
            d[t] = 0
            c = 1  # count the states seen
            while len(S) > 0:
                S2 = []
                for s in S:
                    for t in Ad:
                        ss = (s - t)/b
                        # test if we keep ss
                        keep = True
                        for p, m in zip(places, bo):
                            if K.abs_val(p, ss) > m + .00000001:
                                if verb:
                                    print("|%s|=%s > %s"
                                          % (ss, K.abs_val(p, ss), m))
                                keep = False
                                break
                        if keep:
                            if not d.has_key(ss):
                                S.append(ss)
                                d[ss] = c
                                c += 1
                            L.append((d[s], d[ss], t))
                            # L.append((s, ss, t))
                S = S2
            if d.has_key(0):
                r = DetAutomaton(L, A=Ad, i=0, final_states=[d[0]])
            else:
                r = DetAutomaton(L, A=Ad, i=0, final_states=[])
            if verb:
                print("before pruning: %s" % r)
            if mirror:
                r = r.mirror_det()
            if prune:
                if verb:
                    print("prune...")
                if ext:
                    r = r.prune_inf()
                else:
                    r = r.prune()
            if ext:
                r.set_final_states(range(r.a.n))
        if couples:
            if A is None or B is None:
                raise ValueError("Alphabets A and B must be defined !")
            d = {}
            for c1 in A:
                for c2 in B:
                    if not d.has_key(c1-c2):
                        d[c1-c2] = []
                    d[c1-c2].append((c1, c2))
            if verb:
                print(d)
            r = r.duplicate(d, verb=verb)
        return r

#    def critical_exponent_aprox(self, niter=10, verb=False):
#        """
#        Return an approximation of the critical exponent.
#        This function is inefficient and returns a bad approximation.
#
#        INPUT:
#
#        - ``niter`` int (default: 10) number of iterations
#
#        - ``verb`` - bool (default: ``False``)
#          verbose mode
#
#        OUTPUT:
#        A 
#
#        EXAMPLES::
#
#        #.
#            sage: m = BetaAdicSet(1/(1+I), dag.AnyWord([0,1]))
#            sage: m.critical_exponent_aprox()
#            2.0
#
#        #.
#            sage: s = WordMorphism('1->12,2->13,3->1')
#            sage: m = s.DumontThomas()
#            sage: m.critical_exponent_aprox()
#            2.0994952521...
#
#        """
#        cdef set S, S2, S3
#        b = self.b
#        K = b.parent()
#        A = self.a.alphabet
#        S = set([K.zero()])
#        for i in range(niter):
#            S2 = set([])
#            for s in S:
#                for c in A:
#                    S2.add((s+c)/b)
#            # intervertit S et S2
#            S3 = S2
#            S2 = S
#            S = S3
#            if verb:
#                print(len(S))
#        #m = mahler((1/b).minpoly())
#        m = abs(b.n())
#        return (log(len(S)) / (niter * abs(log(m))))

    def complexity(self, list Ad=None, prec=None, bool verb=False):
        r"""
        Return a estimation of an upper bound of the number of states
        of the relations automaton.
        This estimation is obtained by computing the volume occupied by the lattice containing the BetaAdicSet, in the space product of completions of Q for every absolute archimedian value and p-adic absolute value for which beta as modulus different of one.

        INPUT:

         - ``Ad`` -- list (default: ``None``) - list of differences of digits

         - ``prec`` -- integer (default: ``None``) - precision used for the computation

         - ``verb`` - Boolean (default: ``False``) - Display informations for debug.

        OUTPUT:

        A positive integer.

        EXAMPLES::

            sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
            sage: m.complexity()
            109

            sage: m = BetaAdicSet(x^3-x-1, [0,1])
            sage: m.complexity()
            1125
        """
        b = self.b
        K = b.parent()
        pi = b.minpoly()
        pi = pi*pi.denominator()

        if verb:
            print(K)

        A = self.a.A
        if Ad is None:
            Ad = list(set([c1-c2 for c1 in A for c2 in A]))

        # archimedian places
        places = K.places(prec=prec)
        narch = len(places)
        # ultra-metric places
        from sage.arith.misc import prime_divisors
        lc = pi.leading_coefficient()*pi(0)
        for p in prime_divisors(lc):
            for P in K.primes_above(p):
                if K.abs_val(P, b, prec=prec) != 1:
                    places.append(P)
        if verb:
            print(places)
        # bounds
        bo = []
        vol = 1.
        for i, p in enumerate(places):
            if i < narch:
                bo.append(
                    max([K.abs_val(p, x) for x in Ad])/abs(1 - K.abs_val(p, b)))
                if verb:
                    print("bo = %s" % bo[-1])
                if p(b).imag() == 0:
                    vol *= 2*bo[-1]
                else:
                    vol *= pi_number*bo[-1]**2
            else:
                bo.append(max([K.abs_val(p, x) for x in Ad])/K.abs_val(p, b))
                vol *= bo[-1]
            if verb:
                print("vol = %s", vol)
        if verb:
            print("bounds=%s" % bo)
        # from sage.functions.other import ceil
        return <int>(ceil(vol))

    def intersection(self, BetaAdicSet m, t=0, bool ext=False, int algo=3, bool verb=False):
        r"""
        Assume that b is an algebraic number.
        Compute the intersection of two beta-adic sets.
        (This can also be done with proj() for ext=False.)

        INPUT:

        - ``m`` - the other beta-adic set

        - ``t`` - translate m by t

        - ``ext`` - bool (default: ``False``)
          If True, consider the adherences.

        - ``verb``- bool (default: ``False``)
          If True, verbose mode.

        OUTPUT:

        A BetaAdicSet.

        EXAMPLES::

            #. Compute the boundary of the dragon fractal (see intersection_words for an easier way) ::

                sage: m = BetaAdicSet(1/(1+I), dag.AnyWord([0,1]))
                sage: m1 = m.prefix([0])
                sage: m2 = m.prefix([1])
                sage: mi = m1.intersection(m2, ext=True)
                sage: mi
                b-adic set with b root of x^2 - x + 1/2, and an automaton of 21 states and 2 letters
                sage: mi.plot()     # random

            #. Compute the intersection of two Rauzy fractals (for the same beta)
                sage: m = WordMorphism("a->ab,b->ac,c->a").DumontThomas().mirror()      # Tribonnacci
                sage: m2 = WordMorphism("a->ab,b->ca,c->a").DumontThomas().mirror()     # flipped Tribonnacci
                sage: mi = m.intersection(m2)
                sage: mi
                b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 99 states and 2 letters
                sage: mi.plot()                          # random
                sage: WordMorphism(mi.substitution())    # long time (>4s)
                WordMorphism: a->c, b->ba, c->bd, d->bg, e->id, f->be, g->h, h->bddfb, i->bdj, j->bek, k->idfb

        TESTS::

            sage: m = WordMorphism("a->ab,b->a").DumontThomas().mirror()
            sage: m2 = WordMorphism("a->aab,b->a").DumontThomas().mirror()
            sage: mi = m.intersection(m2)
            Traceback (most recent call last):
            ...
            ValueError: The two beta-adic sets must have same beta.

        """
        cdef DetAutomaton a, ar, ai

        if self.b != m.b:
            raise ValueError("The two beta-adic sets must have same beta.")

        a = self.a.concat_zero_star().product(m.a.concat_zero_star()).prune().minimize()
        if verb:
            print("Product = %s" % a)

        ar = self.relations_automaton(ext=ext, t=t, algo=algo, couples=True, A=self.a.A, B=m.a.A, verb=verb)
        if verb:
            print("Arel = %s" % ar)

        ai = ar.intersection(a)
        if verb:
            print("ai = %s" % ai)

        ai = ai.proji(0)
        if verb:
            print("ai = %s" % ai)

        if ext:
            ai = ai.prune_inf()
        else:
            ai = ai.prune().minimize()
        ai.zero_complete_op()
        return BetaAdicSet(self.b, ai)

    def prefix(self, list w):
        """
        Return a BetaAdicSet like self but where we keep only words starting by w.

        INPUT:

        - ``w`` - list - word that we want as prefix

        OUTPUT:

        BetaAdicSet

        EXAMPLES::

            sage: m = BetaAdicSet(x^3-x^2-x-1, dag.AnyWord([0,1]))
            sage: mp = m.prefix([0, 1, 1, 1]); mp
            b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 5 states and 2 letters
            sage: m.plot_list([mp])     # random

        """
        return BetaAdicSet(self.b, self.a.prefix(w))

    def intersection_words(self, list w1, list w2, bool ext=True, bool verb=False):
        r"""
        Compute the intersection of the adherences of the two beta-adic sets
        corresponding to words with prefix w1 and prefix w2.

        INPUT:

        - ``w1``- word
          The first prefix.

        - ``w2``- word
          The second prefix.

        OUTPUT:

        A Automaton.

        EXAMPLES::

            #. Compute a part of the boundary of the dragon fractal::

                sage: m = BetaAdicSet(1/(1+I), [0,1])
                sage: m.intersection_words([0], [1])
                b-adic set with b root of x^2 - x + 1/2, and an automaton of 21 states and 2 letters

            #. Draw a part of the boundary of a Rauzy fractal::

                sage: m = WordMorphism('a->ab,b->ac,c->a').DumontThomas().mirror()
                sage: mi = m.intersection_words([0], [1])
                sage: mi.plot()                        # not tested
        """
        m1 = self.prefix(w1)
        m2 = self.prefix(w2)
        mi = m1.intersection(m2, ext=ext, verb=verb)
        return mi

    def reduced_words_automaton(self, bool full=False, int step=100,
                                bool mirror=False, int algo_rel=3, bool verb=False):
        r"""
        Compute the reduced words automaton for the alphabet of the automaton of self.
        See http://www.i2m.univ-amu.fr/perso/paul.mercat/Publis/
        Semi-groupes%20fortement%20automatiques.pdf
        for a definition of such automaton.
        The number beta is assumed algebraic.

        INPUT:

        - ``full`` - bool (default: False)
          If True, compute a reduced_words_automaton for the full set
          of words over the alphabet of the automaton of self.

        - ``mirror`` - bool (default: ``False``)
          If True, compute the mirror.

        - ``algo_rel`` - int (default ``3``)
          Algorithm used for the computation of the relations automaton.

        - ``verb`` - bool (default: ``False``)
          If True, print informations for debugging.

        - ``step`` - int (default: 100)
          number of steps done (used for debugging)

        OUTPUT:

        DetAutomaton.

        EXAMPLES::

            sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
            sage: ared = m.reduced_words_automaton()
            sage: ared
            DetAutomaton with 4 states and an alphabet of 2 letters

            sage: m = BetaAdicSet(x^3-x-1, [0,1])
            sage: ared = m.reduced_words_automaton()
            sage: ared
            DetAutomaton with 808 states and an alphabet of 2 letters

        TESTS::

            sage: m = BetaAdicSet(pi, [0,1])
            sage: m.reduced_words_automaton()
            Traceback (most recent call last):
            ...
            ValueError: b must live in a number field!

        """
        cdef list A
        cdef list Ad
        cdef list Adp
        cdef int nAd, nA
        cdef DetAutomaton arel
        cdef int ne, ei

        A = self.a.A
        nA = len(A)

        if full:
            # compute the relations automaton
            arel = self.relations_automaton(algo=algo_rel, mirror=mirror)
            if verb:
                print("arel = %s" % arel)
            if step == 1:
                return arel

            # add a new state
            ei = arel.a.i
            ne = arel.a.n  # new added state
            arel.add_state(True)
            arel.set_final(ei, final=False)  # the new state is final
            if step == 2:
                return arel

            Ad = arel.A
            nAd = len(Ad)

            # add edges from the new state (copy edges from the initial state)
            for j in range(nAd):
                arel.set_succ(ne, j, arel.succ(ei, j))
            if step == 3:
                return arel

            Adp = [i for i in range(
                nAd) if Ad[i] in [x-y for j, x in enumerate(A) for y in A[:j]]]

            # suppress some edges from the initial state
            for j in Adp:
                arel.set_succ(ei, j, -1)
            if step == 4:
                return arel

            # change edges that point to the initial state :
            # make them point to the new state
            for e in arel.states:
                if e != ei:
                    for j in range(nAd):
                        if arel.succ(e, j) == ei:
                            arel.set_succ(e, j, ne)
            if step == 5:
                return arel

            # project, determinise and take the complementary
            d = {}
            for a in A:
                for b in A:
                    if not d.has_key(a - b):
                        d[a-b] = []
                    d[a-b].append((a, b))
            if verb:
                print(d)
            arel = arel.duplicate(d)  # replace differences with couples
            d = {}
            for j in A:
                for i in A:
                    d[(i, j)] = i
            if verb:
                print(d)
                print(arel)
            arel = arel.determinize_proj(d, noempty=False, nof=True)  # , verb=True)
            # project on the first value of the couple, determinise and take the complementary
            if verb:
                print(arel)
            arel = arel.prune()
            if step == 10:
                return arel
            return arel.minimize()
        else:
            arel = self.relations_automaton(couples=True, ext=False)
            if verb:
                print("arel=%s" % arel)
            ap = self.a.product(self.a)
            if verb:
                print("ap=%s" % ap)
            ai = ap.intersection(arel)
            if verb:
                print("ai=%s" % ai)
            alex = DetAutomaton([(0, 0, (i, i)) for i in A]
                                + [(0, 1, (A[i], A[j]))
                                   for i in range(nA) for j in range(i)]
                                + [(1, 1, (i, j)) for i in A for j in A],
                                i=0, final_states=[1])
            if verb:
                print("alex=%s" % alex)
            ai = ai.intersection(alex)
            if verb:
                print("ai=%s" % ai)
            ai = ai.proji(0)
            if verb:
                print("ai=%s" % ai)
            ai.complementary_op()
            if verb:
                print("ai=%s" % ai)
            return ai.intersection(self.a)

    def reduced(self, bool mirror=False, int algo_rel=3, bool verb=False):
        r"""
        Compute a ``BetaAdicSet`` describing the same set, but with unicity (i.e. each point is described by an unique word).

        INPUT:

        - ``mirror`` bool -- (default ``False) set to ``True`` in order to compute the mirror of the reduced language

        - ``algo_rel`` - int (default ``2``)
          Algorithm used for the computation of the relations automaton.

        - ``verb`` - bool (default: ``False``)
          If True, print informations for debugging.


        OUTPUT:

        DetAutomaton.

        EXAMPLES::

            sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
            sage: ared = m.reduced()
            sage: ared
            b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 4 states and 2 letters
        """
        return BetaAdicSet(self.b, self.reduced_words_automaton(mirror=mirror,
                                                                algo_rel=algo_rel,
                                                                verb=verb))

    def critical_exponent_free(self, prec=None, bool verb=False):
        r"""
        Compute the critical exponent of the beta-adic set,
        assuming it is free (or reduced, i.e. there is no relation,
        i.e. every point is uniquely represented by a word).
        When the beta-adic set is moreover algebraic and conformal, then it is equal
        to the Hausdorff dimension of the limit set in the
        contracting space (R or C).

        Rk: beta-adic sets coming from WordMorphism.DumontThomas()
        are always free and algebraic.

        INPUT:

        - ``prec``- precision (default: ``None``)

        - ``verb``- bool (default: ``False``)
          If True, print informations for debugging.

        OUTPUT:

        A real number.

        EXAMPLES::

            #. Hausdorff dimension of limit set of 3-adic expansion with numerals set {0,1,3}::

                sage: m = BetaAdicSet(3, dag.AnyWord([0,1,3]))
                sage: mr = m.reduced()
                sage: mr.critical_exponent_free()
                log(y)/log(3) where y is the max root of x^2 - 3*x + 1, and 3 is root of x - 3.
                0.8760357589...

            #. Hausdorff dimension of limit set of phi-adic expansion with numerals set {0,1}::

                sage: m = BetaAdicSet((1+sqrt(5))/2, dag.AnyWord([0,1]))
                sage: m = m.reduced()
                sage: m.critical_exponent_free()
                log(y)/log(1.618033988749895?) where y is the max root of x^2 - x - 1, and 1.618033988749895? is root of x^2 - x - 1.
                1.0


            #. Hausdorff dimension of the boundary of the dragon fractal::

                sage: m = BetaAdicSet(1/(1+I), dag.AnyWord([0,1]))
                sage: mi = m.intersection_words(w1=[0], w2=[1])
                sage: mi.critical_exponent_free()
                log(y)/log(1.414213562373095?) where y is the max root of x^3 - x^2 - 2, and 1.414213562373095? is root of x^2 - 2.
                1.5236270862...


            #. Hausdorff dimension of the boundary of a Rauzy fractal::

                sage: s = WordMorphism('1->12,2->13,3->1')
                sage: m = s.DumontThomas()
                sage: mi = m.intersection_words(w1=[0], w2=[1])
                sage: mi.critical_exponent_free()
                log(y)/log(1.356203065626296?) where y is the max root of x^4 - 2*x - 1, and 1.356203065626296? is root of x^6 - x^4 - x^2 - 1.
                1.0933641642...

            #. Hausdorff dimension of a non-Pisot Rauzy fractal::

                sage: s = WordMorphism({1:[3,2], 2:[3,3], 3:[4], 4:[1]})
                sage: m = s.DumontThomas().mirror()
                sage: m.critical_exponent_free()
                log(y)/log(1.215716761013442?) where y is the max root of x^3 - x^2 + x - 2, and 1.215716761013442? is root of x^6 - x^4 + 2*x^2 - 4.
                1.5485260383...
        """
        M = self.a.adjacency_matrix()
        if verb:
            print("Eigen values...")
        e = M.eigenvalues()
        if verb:
            print("max...")
        y = max(e, key=abs)
        if verb:
            print("")
        #m = QQbar(mahler((1/self.b).minpoly()))
        m = QQbar(self.b)
        m = m*m.conjugate()
        m.simplify()
        m = m.sqrt()
        m.simplify()
        if m == 1:
            raise NotImplementedError("The computation of the critical exponent is not implemented for a number of absolute value 1.")
        if m < 1:
            m = 1/m
        print("log(y)/log(%s) where y is the max root of %s, and %s is root of %s." % (m, QQbar(y).minpoly(), m, m.minpoly()))
        y = y.n(prec)
        # from sage.functions.log import log
        m = m.n(prec)
        if verb:
            print("y=%s, m=%s" % (y, m))
        return log(y) / abs(log(m))

    def critical_exponent(self, prec=None, int algo_rel=3, bool verb=False):
        r"""
        Compute the critical exponent of the beta-adic set.
        If the beta-adic set is algebraic and conformal, then it is equal
        to the Hausdorff dimension of the limit set in the
        contracting space (R or C). If the beta-adic set is algebraic but not conformal,
        then this critical exponent is equal to the dimension of the limit set
        in the contracting space (product of R, C and p-adic spaces), for an appropriate notion of dimension.

        INPUT:

        - ``prec``- precision (default: ``None``)

        - ``algo_rel`` - int (default: ``2``)
          Algorithm used for the computation of the relations automaton.

        - ``verb``- bool (default: ``False``)
          If True, print informations for debugging.

        OUTPUT:

        A real number.

        EXAMPLES::

            #. Hausdorff dimension of limit set of 3-adic expansion with numerals set {0, 1, 3}::

                sage: m = BetaAdicSet(3, [0,1,3])
                sage: m.critical_exponent()
                log(y)/log(3) where y is the max root of x^2 - 3*x + 1, and 3 is root of x - 3.
                0.8760357589718848

            #. Hausdorff dimension of limit set of phi-adic expansion with numerals set {0, 1}::

                sage: m = BetaAdicSet((1+sqrt(5))/2, [0,1])
                sage: m.critical_exponent()
                log(y)/log(1.618033988749895?) where y is the max root of x^2 - x - 1, and 1.618033988749895? is root of x^2 - x - 1.
                1.0

            #. A non-conformal example::

                sage: P = x^7 - 2*x^6 + x^3 - 2*x^2 + 2*x - 1
                sage: b = P.roots(ring=QQbar)[3][0]
                sage: m = BetaAdicSet(b, [0,1])
                sage: m.critical_exponent()                    # long time
                log(y)/log(1.225816904767620?) where y is the max root of x^11 - 2*x^10 - 4*x^2 + 8*x + 2, and 1.225816904767620? is root of x^42 - 2*x^40 + 2*x^38 - 3*x^36 + 2*x^34 + x^32 - 8*x^30 - 3*x^28 + 6*x^26 + 10*x^24 + 4*x^22 + 4*x^20 + 14*x^18 + 6*x^16 - 11*x^14 - 21*x^12 + 20*x^10 - 9*x^8 + 2*x^6 + x^4 - 1.
                3.3994454205...

        .. SEEALSO::

            #. See more examples with :ref:'../../../thematic_tutorials/beta_adic_set.html'
            critical_exponent_free()

        """
        if verb:
            print("Computation of reduce words' automata")
        m = self.reduced(algo_rel=algo_rel, verb=verb)
        if verb:
            print("%s"%m.a)
        return m.critical_exponent_free(prec=prec, verb=verb)

    # complete the language of a
    def complete(self, list A=None,
                 bool ext=False, DetAutomaton arel=None, bool simplify=True, bool verb=False):
        r"""
        Return the language of all words over the alphabet A
        that describe points of the beta-adic set.
        If ``ext`` is True, it includes words that can be
        prolongated to infinite words that fall
        into the limit set.
        If ``A`` is None, take the alphabet of the automaton of self.

        INPUT:

        - ``A`` - list -- (default : ``None``) alphabet of the result.
          If None, takes the alphabet of the automaton of self.

        - ``ext`` - bool -- (default: ``False``)
          If ''ext'' is True, this also include words equal at infinity.

        - ``arel`` - DetAutomaton (default: ``None``)
            Automaton of relations (if already computed, this permits to
            avoid recomputing it).

        - ``simplify`` - bool (default: ``True``)
            Prune and minimize the result if True.

        - ``verb``- bool (default: ``False``)
          If True, print informations for debugging.

        OUTPUT:

        An automaton.

        EXAMPLES::

            sage: m = BetaAdicSet(3, [0,1,3])
            sage: m.complete([0,1,2])
            b-adic set with b root of x - 3, and an automaton of 2 states and 3 letters

            sage: m = WordMorphism('1->12,2->13,3->1').DumontThomas().mirror()
            sage: m.complete()
            b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 1 state and 2 letters
        """
        cdef DetAutomaton a
        if 0 not in self.a.A:
            a = self.a.bigger_alphabet([0]+self.a.A)
        else:
            a = self.a
        if A is None:
            A = self.a.A
        z = A.index(0)
        from sage.combinat.words.cautomata_generators import dag
        ap = a.concat_zero_star(z=z).product(dag.AnyWord(A))
        if arel is None:
            arel = self.relations_automaton(couples=True, ext=ext, A=a.A, B=A)
        ai = ap.intersection(arel)
        ai = ai.proji(1)
        if ext:
            ai = ai.prune_inf()
        ai.zero_complete_op()
        if simplify:
            ai = ai.prune().minimize()
        return BetaAdicSet(self.b, ai)

    # project the translation by t of self on the zero completion of a
    def proj(self, a, t=0, DetAutomaton arel=None, int algo=3, bool aut=False):
        r"""
        Project the translation by t of self on the zero completion of a.

        INPUT:

        - ``a`` - automaton

        - ``t`` - int (default : ``0``)
          The translation.

        - ``arel`` - DetAutomaton (default : ``None``)
          The relations automaton
          (if ``None`` compute it)

        - ``aut``  bool -- (default: ``False``)
          If True, returns only the DetAutomaton rather that the BetaAdicSet.

        OUTPUT:

        Return a DetAutomaton or a BetaAdicSet

        EXAMPLES::

            #. Use another alphabet to describe a part of the set

            sage: m = BetaAdicSet(3, [0,1,3])
            sage: m2 = BetaAdicSet(3, [0,1,2])
            sage: m.proj(m2)
            b-adic set with b root of x - 3, and an automaton of 2 states and 3 letters

            #. Intersection of the Rauzy fractal with a translated copy of itself

            sage: m = WordMorphism('a->ab,b->ac,c->a').DumontThomas().mirror()
            sage: m = m.proj(m, t=1+m.b)
            sage: m
            b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 12 states and 2 letters
            sage: m.plot()          # not tested

            #. Intersection of two Rauzy fractals (for the same beta)

            sage: m1 = WordMorphism('a->ab,b->ac,c->a').DumontThomas().mirror()
            sage: m2 = WordMorphism('a->ab,b->ca,c->a').DumontThomas().mirror()
            sage: m = m2.proj(m1)
            sage: m
            b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 99 states and 2 letters
            sage: m.plot()      # random

        """
        cdef DetAutomaton ai
        cdef DetAutomaton r

        a = getDetAutomaton(self, a)
        if arel is None:
            # compute the relations automaton with translation t
            arel = self.relations_automaton(t=t, couples=True, algo=algo,
                                            A=a.alphabet, B=self.a.alphabet)
        ai = arel.intersection(a.concat_zero_star().product(self.a.concat_zero_star()))
        r = ai.proji(0)
        r.zero_complete_op()
        if aut:
            return r
        else:
            return BetaAdicSet(self.b, r)

    def shift_op(self, w):
        """
        Shift the automaton of self by w ON PLACE.

        INPUT:

        - ``w`` list (word to shift) or letter

        OUTPUT:

        Return the shifted BetaAdicSet

        EXAMPLES::
            sage: m = BetaAdicSet(3, [0,1,3])
            sage: m.shift_op([0,1,0])
            sage: m
            b-adic set with b root of x - 3, and an automaton of 1 state and 3 letters

        """
        try:
            w = list(w)
            self.a.shift_list_op(w)
        except Exception:
            self.a.shift_op(w)

    def shift(self, w):
        """
        Shift the automaton of self by w.

        INPUT:

        - ``w`` list (word to shift), or letter


        OUTPUT:

        Return the shifted BetaAdicSet

        EXAMPLES::
            sage: m = BetaAdicSet(3, [0,1,3])
            sage: m.shift([0, 1, 0])
            b-adic set with b root of x - 3, and an automaton of 1 state and 3 letters

        """
        m = self.copy()
        m.shift_op(w)
        return m

    # used by Approx
    def _approx_rec(self, DetAutomaton a, test, f, x, int n, int n2):
        r"""
        used by approx

        INPUT:

        - ``a``  DetAutomaton
        - ``test``
        - ``f``
        - ``x``
        - ``n``  int
        - ``n2``  int


        OUTPUT:

        number of state or -1

        TESTS::

            sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
            sage: n = 13
            sage: pm = m.b.parent().places()[1]
            sage: test = lambda x: (pm(x).real())^2 + (pm(x).imag())^2 < .4
            sage: a = DetAutomaton(None, A=m.a.alphabet)
            sage: f = a.add_state(1)
            sage: e = m._approx_rec(a, test, f, 0, n, n)
            sage: e
            3537

        """
        if n == 0:
            if test(x):
                return f
            else:
                return -1
        else:
            e = dict()
            add = False
            for t in a.A:
                e[t] = self._approx_rec(a, test, f, x+t*self.b**(n2-n), n-1, n2)
                if e[t] != -1:
                    add = True
            if add:
                e3 = a.add_state(0)
                for t in self.a.A:
                    if e[t] != -1:
                        a.add_transition(e3, t, e[t])
                return e3
            return -1

    def approx(self, n, test, bool get_aut=False, bool simplify=True):
        """
        Gives a BetaAdicSet describing an approximation of a set defined by the
        characteritic function test, with the alphabet of the automaton of self.
        Rk: could be improved by drawing with the automaton of self
        .. see `thematic_tutorials  beta adic <../../../../thematic_tutorials/beta_adic_set.html>`_

        INPUT:

        - ``n`` -- int - number of iterations/depth of the approximation
        - ``test`` -- test function - function that associated
            to any element of the beta-adic-set, a boolean
        - ``get_aut``  bool -- (default ``False``)
          if ``True`` return only a DetAutomaton
        - ``simplify``  bool -- (default ``True``) set
          to ``True`` to minimize and prune the automaton of the result

        OUTPUT:

        Return a DetAutomaton or a BetaAdicSet

        EXAMPLES::

            #. BetaAdicSet approximating a disk
                sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
                sage: pm = m.b.parent().places()[1]
                sage: a = m.approx(13, lambda x: (pm(x).real())^2 + (pm(x).imag())^2 < .4 )
                sage: print(a)
                b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 201 states and 2 letters
                sage: a.plot()  # not tested

            #. BetaAdicSet approximating a square
                sage: m = WordMorphism('a->ab,b->ac,c->a').DumontThomas().mirror()
                sage: pm = m.b.parent().places()[1]
                sage: a = m.approx(14, lambda x: (pm(x).real())^2 < .3 and (pm(x).imag())^2 < .3 )
                sage: print(a)
                b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 236 states and 2 letters
                sage: m.plot_list([a])  # not tested

            #. Slide of the dragon fractal
                sage: m = BetaAdicSet(1/(1.+I), [0,1])
                sage: m2 = m.approx(12, lambda x: x.real()^2 < .1)
                sage: m2
                (0.500000000000000 - 0.500000000000000*I)-adic set with an automaton of 32 states and 2 letters
                sage: m.plot_list([m2])     # random

            #. BetaAdicSet approximating an image 
                sage: m = WordMorphism('1->12,2->13,3->1').DumontThomas().mirror()
                sage: from sage.arith.beta_adic import ImageIn
                sage: im = ImageIn("SomeImage.png")                                     # not tested
                sage: w = im.width()                                                    # not tested
                sage: h = im.height()                                                   # not tested
                sage: ma = max(w,h)                                                     # not tested
                sage: pm = m.b.parent().places()[1]                                     # not tested
                sage: m.approx(15, lambda x: (pm(x).conjugate()+.5*(1+I))*ma in im)     # not tested
        """
        cdef DetAutomaton a
        a = DetAutomaton(None, A=self.a.A)
        f = a.add_state(1)
        e = self._approx_rec(a, test, f, 0, n, n)
        for t in self.a.A:
            a.add_transition(f, t, f)
        a.a.i = e
        if simplify:
            a = a.minimize()
        if get_aut:
            return a
        else:
            return BetaAdicSet(self.b, a)

    def union(self, a):
        """
        Return the union of BetaAdicSet and automaton or BetaAdicSet a

        INPUT:

        - ``a`` - automaton or BetaAdicSet

        OUTPUT:

        Return the BetaAdicSet union of ``a`` and ``self.a``

        EXAMPLE::

            sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
            sage: a = dag.AnyWord([0, 1, 2, 4])
            sage: m.union(a)
            b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 1 state and 4 letters
            
            #. Disjoint union of two Rauzy fractals with same beta
                sage: s = WordMorphism('1->12,2->13,3->1')
                sage: t = WordMorphism('1->12,2->31,3->1')
                sage: a = s.DumontThomas().mirror().unshift([0,0])
                sage: b = t.DumontThomas().mirror().unshift([1,0,0,0,0])
                sage: m = a.union(b); m
                b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 13 states and 3 letters
                sage: m.plot()              # not tested
                sage: m.substitution()      # long time (>20s)
                {'a': ['h'],
                 'b': ['j'],
                 'c': ['c', 'd'],
                 'd': ['j', 'k', 'v'],
                 'e': ['c', 'd', 'j'],
                 'f': ['j', 'k'],
                 'g': ['c', 'd', 'j', 'k'],
                 'h': ['l', 'g'],
                 'i': ['z', 'c', 'f', 'z'],
                 'j': ['l', 'g', 'z'],
                 'k': ['c', 'd', 'j', 'l', 'h'],
                 'l': ['c', 'd', 'j', 'k', 'v'],
                 'm': ['z', 'c', 'f', 'z', 'e'],
                 'n': ['l', 'g', 'z', 'c', 'b'],
                 'o': ['l', 'j', 'k', 'v'],
                 'p': ['m', 'q', 'c', 'd'],
                 'q': ['x', 'i', 'o'],
                 'r': ['x', 'i', 'x'],
                 's': ['m', 'q', 'c', 'd', 'h'],
                 't': ['n', 'r', 'p', 'a', 'v'],
                 'u': ['l', 'j'],
                 'v': ['m', 'q'],
                 'w': ['n', 'r'],
                 'x': ['k', 'p', 't', 'w'],
                 'y': ['n', 'r', 's', 'v'],
                 'z': ['s', 'p', 'y', 'u']}
        """
        a = getDetAutomaton(self, a)
        return BetaAdicSet(self.b, self.a.union(a))

    def complementary(self, a):
        """
        Compute the complementary of the BetaAdicSet in the BetaAdicSet or automaton a.

        INPUT:

        - ``a`` -- :class:`BetaAdicSet`
          or :class:`sage.combinat.words.cautomata.DetAutomaton`
          in which we take the complementary

        OUTPUT:

        A BetaAdicSet.

        EXAMPLES::

            #. The Rauzy fractal with a hole

                sage: m = WordMorphism('a->ab,b->ac,c->a').DumontThomas().mirror()
                sage: m = m.unshift([1,0,0,0]).complementary(m); m
                b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 7 states and 2 letters
                sage: m.plot()      # not tested

            #. Complementary of a Rauzy fractal in another (for the same beta)

                sage: s = WordMorphism('1->12,2->13,3->1')
                sage: t = WordMorphism('1->12,2->31,3->1')
                sage: a = s.DumontThomas().mirror()
                sage: b = t.DumontThomas().mirror()
                sage: m = b.complementary(a)
                sage: m
                b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 99 states and 2 letters
                sage: m.plot()      # random
                sage: m = a.complementary(b)
                sage: m
                b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 42 states and 3 letters
                sage: m.plot()      # random

        """
        a = getDetAutomaton(self, a)
        return BetaAdicSet(self.b, self.proj(a, aut=True).complementary().intersection(a))

    def unshift(self, l):
        """
        Return a BetaAdicSet with a ``self.a`` unshifted by ``l``

        INPUT:

        - ``l``  list of indices of letters, or the index of a letter

        OUTPUT:

        Return a BetaAdicSet with an unshifted language

        EXAMPLE::

            sage: m = BetaAdicSet(x^3-x^2-x-1, dag.AnyWord([0,1]))
            sage: m.unshift(1)
            b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 2 states and 2 letters
            sage: m = BetaAdicSet(x^3-x^2-x-1, dag.AnyWord([0,1]))
            sage: m.unshift([0,1])
            b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 3 states and 2 letters
        """
        try:
            l = list(l)
            return BetaAdicSet(self.b, self.a.unshiftl(l))
        except:
            return BetaAdicSet(self.b, self.a.unshift(l))

    def diff(self, a):
        """
        Compute the difference of two beta-adic sets.
        Return a beta-adic set describing the set of differences of the two beta-adic sets.

        INPUT:

        - ``a`` - a BetaAdicSet or an automaton

        OUTPUT:

        Return the difference of the two beta-adic sets.


        EXAMPLES::

            sage: m = BetaAdicSet(x^3-x^2-x-1, dag.AnyWord([0,1]))
            sage: a = dag.AnyWord([0, 1, 2, 4])
            sage: m.diff(a)
            b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 1 state and 6 letters

            #. Difference of a Rauzy fractal with itself

            sage: m = WordMorphism('a->ab,b->ac,c->a').DumontThomas().mirror()
            sage: mdr = m.diff(m).reduced()
            sage: mdr
            b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 321 states and 3 letters
            sage: mdr.plot()        # random

            #. Covering of a Rauzy fractal by another

            sage: m = WordMorphism('a->ab,b->ac,c->a').DumontThomas().mirror().unshift([0,0,0,0,0,0,0,0,0,0,0,0,0])
            sage: m2 = WordMorphism('a->ab,b->ca,c->a').DumontThomas().mirror()
            sage: mdr = m.diff(m2).reduced()    # long time (8s)
            sage: mdr                           # long time
            b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 1879 states and 5 letters
            sage: mdr.plot()                    # random

        """
        a = getDetAutomaton(self, a)
        return BetaAdicSet(self.b, self.a.diff(a))

    def is_Pisot(self, bool verb=False):
        """
        Test if the number b is the conjugate of a Pisot number or not.

        INPUT:

        - ``verb`` bool -- (default : ``False``) set to ``True`` for verbose mode
          If true, explains why we return False when it happens.

        OUTPUT:

        Return ``True`` or ``False``

        EXAMPLES::

            sage: m = BetaAdicSet(x^2-x-1, [0,1])
            sage: m.is_Pisot()
            True

            sage: m = BetaAdicSet(x^4-2*x^3+x^2-2*x+1, [0,1])
            sage: m.is_Pisot()
            False

            sage: m = BetaAdicSet(1+I, [0,1])
            sage: m.is_Pisot(verb=True)
            There is a conjugate of modulus greater than one which is not real.
            False

        """
        try:
            if not self.b.is_integral():
                if verb:
                    print("b is not an algebraic integer.")
                return False
            pi = self.b.minpoly()
            rr = [r[0] for r in pi.roots(ring=QQbar)]
            np = 0
            for r in rr:
                if abs(r) > 1:
                    if np != 0:
                        if verb:
                            print("There are more than one conjugate of modulus > 1.")
                        return False
                    from sage.rings.qqbar import AA
                    if r not in AA:
                        if verb:
                            print("There is a conjugate of modulus greater than one which is not real.")
                        return False
                    if r < 0:
                        if verb:
                            print("There is conjugate of modulus greater than one which is negative.")
                        return False
                    np = 1
                elif abs(r) == 1:
                    if verb:
                        print("There is a conjugate of modulus one.")
                    return False
            if np == 0:
                if verb:
                    print("There is no conjugate of modulus > 1.")
                return False
            return True
        except Exception:
            if verb:
                print("b is not an algebraic number.")
            return False

    def points(self, int n=1000, int npts=10000):
        """
        Compute points (in the number field of b) corresponding to words of length k recognized by the automaton,
        where k is at most n, and the total number of points is approximatively npts.
        Return (k, list of couples (state, point))

        INPUT:

        - ``n`` - integer (default: 1000)
          The maximum number of iterations.

        - ``npts`` - integer (default: 10000 )
          Approximation of a bound on the number of points computed.

        OUTPUT:

        Return (k, list of couples (state, point)),
        where k is the number of iterations computed.

        EXAMPLES::

            #. The dragon fractal::
                sage: m = BetaAdicSet(1/(1+I), [0, 1])
                sage: print(m)
                b-adic set with b root of x^2 - x + 1/2, and an automaton of 1 state and 2 letters
                sage: P = m.points()
                sage: P[0]
                13
                sage: len(P[1])
                8192
                sage: points([x.n() for i,x in P[1]], aspect_ratio=1)   # long time
                Graphics object consisting of 1 graphics primitive
        """
        cdef int i, j, k, f, nA
        nA = self.a.a.na
        l = self.a.prune().spectral_radius()
        n = min(n, <int>(log(<double>npts)/log(<double>l)))
        r = [(self.a.a.i, 0)]
        bn = 1
        for i in range(n):
            rr = []
            for j, t in r:
                for k in range(nA):
                    f = self.a.a.e[j].f[k]
                    if f != -1:
                        rr.append((f, t + bn*self.a.A[k]))
            bn = bn*self.b
            r = rr
        return (n, r)

    def zero_ball(self, p, int npts=1000):
        """
        Compute the radius of a ball centered at 0 and that covers the BetaAdicSet for the place p.
        We assume that abs(p(self.b)) < 1.

        INPUT:

        - ``p`` - archimedian place

        - ``npts`` - integer (default: 10000 )
            Approximation of the number of points computed to find the bound.

        """
        pts = self.points(npts=npts)
        M = abs(p(self.b**pts[0]))*max([abs(p(c))
                                        for c in self.a.A])/abs(1-abs(self.b))
        return max([abs(p(c[1]))+M for c in pts[1]])

    def diameter(self, p, int n=10, bool verb=False):
        """
        Compute an upper bound of the diameter of the BetaAdicSet for the place p.
        The error has order p(self.b)^n.
        (The algorithm used here is not optimal.)

        INPUT:

        - ``p`` - archimedian place used to compute the diameter

        - ``n`` - integer (default: 10) - number of iterations

        - ``verb`` bool -- (default : ``False``) set to ``True`` for verbose mode

        OUTPUT:

        double

        EXAMPLES::

            sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
            sage: p = m.b.parent().places()[1]
            sage: m.diameter(p)
            2.93122910465427
        """
        cdef int i, j, k, f, f2, nrr, nA
        cdef double d, dmm, dm2
        nA = self.a.a.na
        r = [(self.a.a.i, 0)]
        bn = 1
        M = max([abs(p(c)) for c in self.a.A])/abs(1-abs(self.b.n()))
        import numpy as np
        for i in range(n):
            rr = []
            for j, t in r:
                for k in range(nA):
                    f = self.a.a.e[j].f[k]
                    if f != -1:
                        rr.append((f, t + bn*self.a.A[k]))
            bn = bn*self.b
            if verb:
                print("rr : %s elements" % len(rr))
            r = []
            # compute the diameter of the set rr (this could be improved)
            dmm = 0
            dm = np.zeros(len(rr), dtype=np.float)
            v = np.empty(len(rr), dtype=np.complex)
            for f, (j, t) in enumerate(rr):
                v[f] = p(t)
            nrr = len(rr)
            for f in range(nrr):
                dm2 = 0
                for f2 in range(nrr):
                    d = abs(v[f] - v[f2])
                    if d > dm2:
                        dm2 = d
                dmm = fmax(dmm, dm2)
                dm[f] = dm2
            if verb:
                print("dmm = %s" % dmm)
            M2 = 2*abs(p(bn))*M
            if i == n-1:
                return dmm+M2
            for f, (j, t) in enumerate(rr):
                if dm[f]+M2 >= dmm:
                    r.append((j, t))
            if verb:
                print("r : %s elements" % len(r))

    def translations_iterator(self, bool test_Pisot=True, int ndiam=20, bool verb=False):
        """
        Compute a list of numbers containing the positive
        part of the BetaAdicSet, ordered in the expanding direction.
        Assume that self.b is the conjugate of a Pisot number.

        INPUT:

        - ``test_Pisot``  bool -- (default : ``True``) test if b is the conjugate of a Pisot number as needed

        - ``verb`` bool -- (default : ``False``) set to ``True`` for verbose mode

        - ``ndiam`` int  -- (default : 20): number of iterations
          used for the estimation of the diameter

        OUTPUT:
        Return an iterator.


        EXAMPLES::

            sage: m = BetaAdicSet((x^3-x^2-x-1).roots(ring=QQbar)[1][0], dag.AnyWord([0,1]))
            sage: m.translations_iterator().next()
            -b + 2

        """
        cdef int n, i, j
        if test_Pisot:
            if not self.is_Pisot():
                raise ValueError("b must be the conjugate of a Pisot number")
        # take a basis of the lattice
        d = self.b.minpoly().degree()
        B = [self.b**i for i in range(d)]
        # compute the min of the differences for every place
        Bd = set([a-b for a in B for b in B if a != b])
        K = self.b.parent()
        n = -2147483648
        # from sage.functions.other import ceil
        # from sage.functions.log import log
        P = [p for p in K.places() if abs(p(self.b)) < 1]
        M = [self.diameter(p, n=ndiam) for p in P]

        for i, p in enumerate(P):
            m = min([abs(p(b)) for b in Bd])
            if verb:
                print("p=%s, m=%s, M=%s" % (p,m,M))
                print("%s" % (log(m/(2*M[i]))/log(abs(p(self.b)))))
            n = max(n, 1+<int>floor(log(m/(2*M[i]))/log(abs(p(self.b)))))
        if verb:
            print("n=%s" % n)
        # multiply the bound by this power of b
        bn = self.b**n
        M = [M[i]*abs(p(bn)) for i, p in enumerate(P)]
        # compute the matrix corresponding to the multiplication by M to the left
        from sage.matrix.constructor import identity_matrix
        I = identity_matrix(d)
        pi = self.b.minpoly()
        pi /= pi.leading_coefficient()
        from sage.matrix.constructor import matrix

        m = matrix(
            [I[i] for i in range(1, d)] +
            [[-c for c in pi.list()[:d]]]).transpose()

        if verb:
            print("m=%s" % m)
        # compute the Perron-Frobenius eigenvector
        from sage.modules.free_module_element import vector
        v = vector(max(
            [r[1][0] for r in m.right_eigenvectors()], key=lambda x: x[0]))
        v /= sum(v)
        vB = vector(B)
        if verb:
            print("v=%s" % v)
        r = []
        from itertools import count
        for j in count(start=1):
            vi = vector([<int> round(j * x) for x in v])
            t = vi * vB
            if t == 0:
                continue
            if verb:
                print("j=%s, t=%s" % (j, t))
            # test if t is in the domain
            keep = True
            for i, p in enumerate(P):
                if abs(p(t)) > M[i]:
                    keep = False
                    break
            if keep:
                yield t/bn

    def translations_diff_iterator(self, bool test_Pisot=True, 
                                   int ndiam=20, bool verb=False):
        """
        Compute a list that contains the set of positive differences of points of the BetaAdicSet.
        The list is increasing for the expanding place.
        Assume that self.b is a Pisot number.

        INPUT:

        - ``test_Pisot``  bool -- (default : ``True``) : test if b is 
          the conjugate of a Pisot number as needed
          B : basis of a lattice containing the BetaAdicSet

        - ``ndiam`` int -- (default : 20): number of iterations used
          for the computation of the estimation of the diameter

        - ``verb`` bool -- (default : ``False``) set to ``True`` for verbose mode

        OUTPUT:
        Return an iterator.

        EXAMPLES::

            sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
            sage: it = m.translations_diff_iterator()
            sage: it.next()
            -b + 2
            sage: it.next()
            -b^2 + 2*b

        """
        cdef int n, i, j
        if test_Pisot:
            if not self.is_Pisot():
                raise ValueError("b must be the conjugate of a Pisot number")
        # take a basis of the lattice
        d = self.b.minpoly().degree()
        B = [self.b**i for i in range(d)]
        # compute the min of the differences for every place
        Bd = set([a-b for a in B for b in B if a != b])
        K = self.b.parent()
        n = -2147483648
        # from sage.functions.other import ceil
        # from sage.functions.log import log
        P = [p for p in K.places() if abs(p(self.b)) < 1]
        M = [self.diameter(p, ndiam) for p in P]
        for i, p in enumerate(P):
            m = min([abs(p(b)) for b in Bd])
            if verb:
                print("p=%s, m=%s, M=%s" % (p, m,M))
                print("%s" % (log(m/(2*M[i]))/log(abs(p(self.b)))))
            n = max(n, 1+<int> floor(log(m/(2*M[i])) / log(abs(p(self.b)))))
        if verb:
            print("n=%s" % n)
        # multiply the bound by this power of b
        bn = self.b**n
        M = [M[i]*abs(p(bn)) for i, p in enumerate(P)]
        # compute the matrix corresponding to the multiplication by M to the left
        from sage.matrix.constructor import identity_matrix
        I = identity_matrix(d)
        pi = self.b.minpoly()
        pi /= pi.leading_coefficient()
        if verb:
            print("pi=%s" % pi)
        from sage.matrix.constructor import matrix

        m = matrix(
            [I[i] for i in range(1,d)] +
            [[-c for c in pi.list()[:d]]]).transpose()

        if verb:
            print("m=%s" % m)
        # compute the Perron-Frobenius eigenvector
        from sage.modules.free_module_element import vector
        v = vector(max([r for r in m.right_eigenvectors()],
                        key=lambda x: abs(x[0]))[1][0])
        v /= sum(v)
        vB = vector(B)
        if verb:
            print("v=%s" % v)
        r = []
        from itertools import count
        for j in count(start=1):
            vi = vector([<int>round(j * x) for x in v])
            t = vi*vB
            if t == 0:
                continue
            if verb:
                print("j=%s, t=%s"%(j,t))
            # test if t is in the domain
            keep = True
            for i, p in enumerate(P):
                if abs(p(t)) > M[i]:
                    if verb:
                        print("%s > %s"%(abs(p(t)), M[i]))
                    keep = False
                    break
            if keep:
                yield t/bn

# TO BE ADDED LATER
#    def interior(self, verb=False):
#        r"""
#        We assume that self.b is a Pisot number.
#        Compute a BetaAdicSet describing the interior for the topology for which open sets are
#        sets of points of that are in the BetaAdicSet with same beta but with a language that recognize every words over the alphabet of self.a and that projects to an open set of P, for the natural projection on the contracting space (which is a product of copies of R, C and p-adic fields).
#        """
#        def S(a, verb=False):
#            F = []
#            for e in a.states:
#                ok = True
#                for l in range(len(a.alphabet)):
#                    if a.succ(e, l) != e:
#                        ok = False
#                        break
#                    if ok:
#                        F.append(e)
#                        a2 = a.copy()
#                        a2.set_final_states(F)
#                        if verb:
#                            print "F =",F
#                        return a2
#        arel = self.relations_automaton(t=0,couples=True)
#        if verb:
#            print "arel =",arel
#            #arel.plot()
#        ap = dag.AnyWord(self.a.alphabet).product(self.a.concat_zero_star())
#        ai = ap.intersection(arel)
#        if verb:
#            print "ai =",ai
#        aip = ai.proji(0)
#        if verb:
#            print "aip =", aip
#        aip.zero_complete_op()
#        af = S(aip.minimize())
#        af.zero_complete_op()
#        if verb:
#            print "af =",af
#        af2 = af.minimize().intersection(self.a)
#        af2 = af2.prune()
#        if verb:
#            print "af2 =",af2
#            print af2.equals_langages(self.a)
#        return BetaAdicSet(self.b, af2) 
#

    def domain_exchange(self, n=None, int algo=1, int algo_rel=3, bool test_Pisot=True,
                        int ndiam=30, bool verb=False):

        """
        Compute the domain exchange describing the BetaAdicSet.
        Assume that self.b is a Pisot number.
        Return a list of (translation, BetaAdicSet).

        INPUT:

        - ``n`` - int -- (default: ``None``)

        - ``algo`` - int -- (default: 1)
            algorithm used to compute the list of translations

        - ``algo_rel`` - int -- (default: 3)
            algorithm used to compute the relations automaton

        - ``test_Pisot``  bool -- (default: ``True``)
            test if b is the conjugate of a Pisot number as needed

        - ``ndiam`` int -- (default: 30) : number of iterations used for
           the estimation of the diameter

        - ``verb`` bool -- (default: ``False``) set to ``True`` 
          for verbose mode


        OUTPUT:
        List of tuple ``BetaAdicSet``

        EXAMPLES::

            #Domain exchange of the Tribonnacci substitution
            sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
            sage: l = m.domain_exchange(); l
            [(b^2 - b - 1,
              b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 4 states and 2 letters),
             (b - 1,
              b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 4 states and 2 letters),
             (1,
              b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 4 states and 2 letters)]
            sage: m.plot_list([a for t,a in l])             # not tested
            sage: m.plot_list([a.proj(m, t) for t,a in l])  # random
            <PIL.Image.Image image mode=RGBA size=800x600 at 0x7F57DFF3BC10>

            # A more complicated domain exchange
            sage: m = BetaAdicSet(x^3 - x^2 - x - 1, DetAutomaton([[0, 1], [(0, 17, 0), (0, 4, 1), (1, 16, 0), (2, 17, 0), (2, 4, 1), (3, 17, 0), (4, 17, 0), (5, 7, 0), (5, 0, 1), (6, 5, 0), (6, 0, 1), (7, 6, 0), (8, 10, 0), (9, 8, 0), (9, 0, 1), (10, 9, 0), (11, 15, 0), (11, 1, 1), (12, 14, 0), (12, 11, 1), (13, 8, 0), (13, 2, 1), (14, 13, 0), (14, 18, 1), (15, 5, 0), (15, 2, 1), (16, 17, 0), (16, 0, 1), (17, 17, 0), (17, 0, 1), (18, 16, 0), (18, 3, 1)]], i=12, final_states=[0, 1, 2, 3, 4, 16, 17, 18]))
            sage: l = m.domain_exchange(); l
            [(b^2 - b - 1,
              b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 5 states and 2 letters),
             (b - 1,
              b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 9 states and 2 letters),
             (1,
              b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 8 states and 2 letters),
             (2,
              b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 6 states and 2 letters),
             (2*b - 1,
              b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 7 states and 2 letters)]
            sage: m.plot_list([a for t,a in l])             # random
            <PIL.Image.Image image mode=RGBA size=800x600 at 0x7F57DFE50450>
            sage: m.plot_list([a.proj(m, t) for t,a in l])  # random
            <PIL.Image.Image image mode=RGBA size=800x600 at 0x7F57DFE50BD0>
            
            # See the thematic tutorial for more examples

        """
        if algo == 1:
            if verb:
                print("compute translations...")
            it = self.translations_diff_iterator(test_Pisot=test_Pisot,
                                                 ndiam=ndiam, verb=verb)
        else:
            if verb:
                print("diff...")
            md = self.diff(test_Pisot=test_Pisot)
            if verb:
                print("compute translations...")
            it = md.translations_iterator(verb=verb, ndiam=ndiam)
        m = self.copy()
        # from sage.combinat.words.cautomata_generators import dag
        # a = self.a.intersection(dag.AnyWord([0], A2=self.a.A).complementary())
        a = self.a.copy()
        r = []
        if n is None:
            n = -1
        for t in it:
            if not t.is_integral():
                if verb:
                    print("t=%s not integral", t)
                continue
            if verb:
                print("t=%s" % t)
            mi = m.intersection(m, -t, algo=algo_rel)
            mia = mi.a.intersection(a)
            if not mia.has_empty_language():
                if verb:
                    print("not empty ! mia=%s" % mia)
                mi = BetaAdicSet(m.b, mia)
                r.append((t, mi))
                a = a.intersection(mi.a.complementary())
                if a.has_empty_language():
                    return r
            n -= 1
            if n == 0:
                return r

    def substitution(self, DetAutomaton ap=None,
                              np=None, list lt=None, bool need_included=True,
                              bool get_aut=False, step=None, bool verb=False):
        r"""
        Assume that b is a conjugate of a Pisot number.
        Compute a substitution whose discrete line is this BetaAdicSet.

        Return a substitution given as a dictionnary. If get_aut is True, return also a list of (translation, automaton) describing each piece of the Rauzy fractal.

        INPUT:

        - ``ap``- DetAutomaton (default: ``None``)
            Language used to do the computations: we project everything on it. If ap is None, we use the automaton of self.

        - ``np``- int (default: ``None``)
            Power of beta used for the computing. The BetaAdicSet
            must be b^np invariant.
            If np is None, take the smallest possible positive integer.

        - ``lt``- list (default: None)
            List of (DetAutomaton, translations) that describes the pieces exchange,
            where the translation is an element of the integer ring.
            If None, compute it by calling self.domain_exchange().

        - ``get_aut``- Bool (default: ``False``)
            If True, gives also the list of automata.

        - ``verb``- bool (default: ``True``)
          If True, print informations about the computing.

        OUTPUT:

        A word morphism given by a dictionnary.

        EXAMPLES::

            #. Tribonnacci::

                sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
                sage: m.substitution()
                {1: [3], 2: [3, 1], 3: [3, 2]}
            
            #. Example with infinitely many connected components and where zero is not an inner point
            
                sage: m = BetaAdicSet(x^3-x^2-x-1, dag.AnyWord([0]).concat(dag.Word([1,0,0,0])).concat(dag.AnyWord([0,1])))
                sage: WordMorphism(m.substitution())
                WordMorphism: a->c, b->ba, c->d, d->h, e->gi, f->jma, g->fma, h->be, i->l, j->bma, k->ga, l->fe, m->gka

            #. Substitution whose Rauzy fractal approximate a disk

                sage: m = BetaAdicSet(x^3-x^2-x-1, [0,1])
                sage: pm = m.b.parent().places()[1]
                sage: a = m.approx(13, lambda x: (pm(x).real())^2 + (pm(x).imag())^2 < .4 )
                sage: s = WordMorphism(a.substitution())    # long time (>30s)
                sage: s.rauzy_fractal_plot()                # not tested

            #. Find a substitution whose Rauzy fractal is what the user draw

                sage: m = WordMorphism('a->ab,b->ac,c->a').DumontThomas().mirror()
                sage: m = m.user_draw()                     # not tested (need the intervention of the user)
                sage: s = WordMorphism(m.substitution())    # not tested
                sage: s.rauzy_fractal_plot()                # not tested
           
            #. The Tribonnacci Rauzy fractal with a hole

                sage: m = WordMorphism('a->ab,b->ac,c->a').DumontThomas().mirror()
                sage: m = m.unshift([1,0,0,0]).complementary(m); m
                b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 7 states and 2 letters
                sage: m.substitution()
                {1: [5], 2: [5, 1], 3: [3, 5, 1], 4: [3], 5: [6, 2], 6: [4, 2]}

            #. Disjoint union of two Rauzy fractals with same beta

                sage: s = WordMorphism('1->12,2->13,3->1')
                sage: t = WordMorphism('1->12,2->31,3->1')
                sage: a = s.DumontThomas().mirror().unshift([0,0])
                sage: b = t.DumontThomas().mirror().unshift([1,0,0,0,0])
                sage: m = a.union(b)
                sage: WordMorphism(m.substitution())        # long time (>15s)
                WordMorphism: a->h, b->j, c->cd, d->jkv, e->cdj, f->jk, g->cdjk, h->lg, i->zcfz, j->lgz, k->cdjlh, l->cdjkv, m->zcfze, n->lgzcb, o->ljkv, p->mqcd, q->xio, r->xix, s->mqcdh, t->nrpav, u->lj, v->mq, w->nr, x->kptw, y->nrsv, z->spyu
                
             #. A substitution whose Rauzy fractal is the union of a Cantor set and an interval
             
                sage: a = DetAutomaton([(0,1,0), (1,1,0), (1,1,1), (1,1,2), (0,2,2), (2,2,0), (2,2,2)], i=0)
                sage: m = BetaAdicSet(1-sqrt(2), a)
                sage: m.substitution()
                {1: [5],
                2: [3, 1, 4],
                3: [3, 1, 7],
                4: [9, 2],
                5: [3, 1],
                6: [6, 6, 8],
                7: [9, 2, 6, 8],
                8: [6],
                9: [7]}

             #. See more examples in the thematic tutorial

        TESTS::

            sage: m = BetaAdicSet(1/(1+I), [0,1])
            sage: m.substitution()
            Traceback (most recent call last):
            ...
            ValueError: The number b of the BetaAdicSet must be for the conjugate of a Pisot number.

        """
        cdef DetAutomaton a
        
        #test if b is a Pisot number
        if not self.is_Pisot():
            raise ValueError("The number b of the BetaAdicSet must be for the conjugate of a Pisot number.")
        #ensure that the alphabet of a contains 0
        a = self.a.concat_zero_star()
        a.zero_complete_op()
        A = a.A
        try:
            l0 = A.index(0)
        except:
            A = [0]+A
            a = a.bigger_alphabet([0]+A)
            l0 = 0
        #compute the domain exchange if necessary
        if lt is None:
            if verb:
                print("Compute the domain exchange...")
            l = self.domain_exchange(test_Pisot=False)
            lt = [(m.a,t) for t,m in l]
        if verb:
            print("Domain exchange with %s pieces."%len(lt))
        m = BetaBase(self.b)
        if ap is None:
            ap = a
        if verb:
            print("ap=%s" % ap)
        cdef DetAutomaton aa
        if not a.included(ap):
            aa = a.copy()
            aa.zero_complete_op()
            # check that Qap contains Qaa
            if not m.Proj(ap, aa).has_same_language_as(aa) and need_included:
                raise ValueError("The g-beta-set described by a is not included in the one described by ap.")
            # project aa on ap
            aa = m.Proj(aa, ap)
        else:
            aa = a
        if verb:
            print("aa=%s" % aa)
        A = aa.A
        # test if np is big enough
        if np is None:
            for i in range(1, 300):
                ba = aa.unshift(0, i)
                ba.zero_complete_op()
                if m.Proj(aa, ba).has_same_language_as(ba):
                    ba = m.Proj(ba, aa)
                    np = i
                    break
            if np is None:
                raise ValueError('The g-beta-expansion must be b^np invariant for some natural integer np.')
        else:
            ba = aa.unshift(0, np)
            ba.zero_complete_op()
            if not m.Proj(aa, ba).has_same_language_as(ba):
                raise ValueError('The g-beta-expansion must be b^np invariant (here np=%s).' % np)
            ba = m.Proj(ba, aa)
        if verb:
            print("np = %s" % np)
        # compute the pieces exchange
        if lt is None:
            raise NotImplementedError("You have to compute the domain exchange yourself for the moment !")
        lt = [(m.Proj(a, ap), t) for a, t in lt]
        if verb:
            print("Exchange of %s pieces" % len(lt))
        # compute the induction from the list of (piece, translation)
        # precomputation
        if verb:
            print("Pre-computation...")
        arel = dict()
        for a, t in lt:
            arel[t] = m.relations_automaton(t=-t, A=aa.A,
                                             B=ap.A, couples=True)
            if verb:
                print("arel[%s]=%s" % (t, arel[t]))
        if verb:
            print("ba : %s" % ba)
        # tree of subdivision of the pieces
        tree = [range(1, len(lt) + 1)] + [[] for i in range(len(lt))] 
        if verb:
            print("initial tree: %s" % tree)
        lm = [(aa, 0)] + lt  # list of pieces, translations
        if verb:
            print("lm = %s" % lm)
        # browse each piece (given by the list of pieces)
        d = [[] for i in range(len(lm))]
        if verb:
            print("d = %s" % d)
        lf = range(1, len(lm))  # list of leaf
        if verb:
            print("lf = %s" % lf)

        from copy import copy

        if verb:
            print("\n**********************\n   Step 1   \n**********************")

        # étape 1 : complétion des mots
        for i, (a1, t1) in enumerate(lm):
            if tree[i] != []:
                continue  # this piece is not a leaf
            if verb:
                print("\nCompute the piece %s/%s (%s, %s)..." % (i, len(lm), a1, t1))
                # print "lf = %s"%lf
                # print "d = %s"%d
                # print "tree = %s"%tree
            tr = 0  # total translation
            if d[i] != []:
                if d[i][-1] == -1:
                    continue  # the computation for this piece was already finished
                # go to the end of the word
                for j in d[i]:
                    if j < 0:
                        break
                    tr += lm[j][1]
            # compute b^np*a + tr
            a = a1.unshift(0, np).prune().minimize()
            if tr != 0:
                if verb:
                    print("Translation of %s..." % tr)
            # m.move2(t=-tr, a=a)
            # TODO : do not recompute this automaton already computed
            a = m.Proj(a, ap, t=-tr)
            if a.has_empty_language():
                raise RuntimeError("Empty language when projecting a on ap for t=%s"%(-tr))
            while True:
                # split according to other pieces
                j = included(a, lf, lm)
                if j is None:
                    # find the pieces that intersect a
                    l = []
                    for j in lf:
                        if lm[j][0].intersect(a):
                            l.append(j)
                    if len(l) < 2:
                        print("Error : intersection with %s piece but not included !!!"%len(l))
                    if verb:
                        print("Subdivision on %s pieces..." % len(l))
                    # compute intersections (split a1)
                    for j in l:
                        a2 = lm[j][0]
                        # cut a according to a2
                        # m.move2(t=tr, a=a2) #translate a2 de -tr
                        a = m.Proj(a2, ap.unshift(0, np), t=tr)
                        a.shift_op(0, np)  # multiply by b^(-np)
                        a = a.prune().minimize()
                        k = len(lm)  # index of the new piece
                        lf.append(k)  # new leaf
                        tree[i].append(k)
                        tree.append([])
                        from copy import copy
                        # print copy
                        d.append(copy(d[i]))
                        d[k].append(j)  # add the next translation
                        # add the new piece to the list
                        lm.append((a.intersection(a1), t1))
                        # split according to ba
                        (ab, abc) = split_ba(k, tr+lm[j][1], np,
                                             lm, m, aa, ap, verb)
                        if ab is None:
                            if verb:
                                print("k=%s, tr=%s+%s : computation to continue" % (k, tr, lm[j][1]))
                        else:
                            if abc is None:
                                if verb:
                                    print("tr=%s : computation finished" % tr)
                                d[k].append(-1)  # indicate that the computation of this piece is terminated
                            else:
                                if verb:
                                    print("tr=%s : subdivision of %s according to ba (new %s)..." % (tr, i, len(lm)))
                                lf.append(len(lm))  # new leaf
                                tree[k].append(len(lm))
                                tree.append([])
                                d.append(copy(d[k]))
                                # indicate that the computation is terminated for this piece (for the first step)
                                d[len(lm)].append(-1)
                                lm.append((ab, t1))
                                lf.append(len(lm))  # new leaf
                                tree[k].append(len(lm))
                                tree.append([])
                                d.append(copy(d[k]))
                                lm.append((abc, t1))
                                lf.remove(k)  # the piece k is no more a leaf
                    lf.remove(i)  # the piece i is no more a leaf
                    # computation ended for this piece since it is no more a leaf
                    break
                else:
                    # add the piece to the list and translate
                    d[i].append(j)
                    # if verb: print "Translation by %s..."%lm[j][1]
                    # m.move2(t=-lm[j][1], a=a, ar=arel[lm[j][1]])
                    a = m.Proj(a, ap, t=-lm[j][1], arel=arel[lm[j][1]])  
                    tr += lm[j][1]
                # split according to ba
                (ab, abc) = split_ba(i, tr, np, lm, m, aa, ap, verb)
                if ab is None:
                    pass
                    # if verb: print "tr=%s : computation to continue"%tr
                else:
                    if abc is None:
                        if verb:
                            print("tr=%s : end of computation" % tr)
                    else:
                        if verb:
                            print("tr=%s : subdivision of %s according to ba (new %s)..." % (tr, i, len(lm)))
                        lf.append(len(lm))  # new leaf
                        tree[i].append(len(lm))
                        tree.append([])
                        d.append(copy(d[i]))
                        # indicate that the computation is finished for this piece (for the first step)
                        d[len(lm)].append(-1)
                        lm.append((ab, t1))
                        lf.append(len(lm))  # new leaf
                        tree[i].append(len(lm))
                        tree.append([])
                        d.append(copy(d[i]))
                        lm.append((abc, t1))
                        lf.remove(i)  # the piece i is no more a leaf
                    break  # computation terminated for this piece (for the first step)

        if verb:
            print("\n*************\n   Step 2   \n*************")

        # second step : replacement of the letters that are not leaves
        while True:
            end = True
            for i in lf:
                a1, t1 = lm[i]
                if verb:
                    print("\nPiece %s/%s..." % (i, len(lm)))
                tr = 0  # total translation
                if d[i] == []:
                    print("Error : empty leaf !!!!")
                # got to the end of the word
                for ij, j in enumerate(d[i]):
                    if j < 0:
                        break
                    if tree[j] != []:  # we have to recompute this letter
                        # compute b^np*a + tr
                        a = a1.unshift(0, np).prune().minimize()
                        if tr != 0:
                            if verb:
                                print("Translation of %s..." % tr)
                        # m.move2(t=-tr, a=a)
                        # TODO : do not recompute this automaton already computed
                        a = m.Proj(a, ap, t=-tr)
                        # split according to the other pieces
                        f = fils(tree, j)
                        if verb:
                            print("Split by %s pieces" % len(f))
                        k = included(a, f, lm)
                        if k is None:
                            end = False
                            # find pieces that intersect a
                            l = []
                            for k in lf:
                                if lm[k][0].intersect(a):
                                    l.append(k)
                            if len(l) < 2:
                                print("Error : intersection with %s pieces but not included !!!" % len(l))
                            if verb:
                                print("Subdivision of %s pieces..." % len(l))
                            # compute intersections (split a1)
                            for j2 in l:
                                a2 = lm[j2][0]
                                # cut according to a2
                                # a = m.move2(t=tr, a=a2) #translate a2 de -tr
                                # a.zero_complete_op()
                                a = m.Proj(a2, ap.unshift(0, np), t=tr)
                                a.shift_op(0, np)  # multiply by b^(-np)
                                a = a.prune().minimize()
                                k = len(lm)  # index of the new piece
                                lf.append(k)  # new leaf
                                tree[i].append(k)
                                tree.append([])
                                d.append(copy(d[i]))
                                d[k][ij] = j2  # replace the letter
                                # add the new piece to the list
                                lm.append((a.intersection(a1), t1))
                            lf.remove(i)  # i is not more a leaf
                            if verb:
                                print("break...")
                            break  # the piece is no more a leaf
                        else:
                            # replace the letter
                            d[i][ij] = k
                    tr += lm[j][1]
            if end:
                break
        lf = [i for i in range(len(tree)) if tree[i] == []]
        # compute the substitution
        s = dict()
        for i in lf:
            if d[i][-1] < 0:
                d[i].pop()
            s[i] = d[i]
        # recode the substitution
        l = s.keys()
        dl = dict()  # inverse of l
        if len(l) > 9 and len(l) < 27:
            for i, k in enumerate(l):
                dl[k] = chr(i+ord('a'))
        else:
            for i, k in enumerate(l):
                dl[k] = i+1
        d = dict()
        for i in l:
            d[dl[i]] = [dl[j] for j in s[i]]
        if get_aut:
            return d, [(a, t) for i, (a, t) in enumerate(lm) if tree[i] == []]
        else:
            return d

