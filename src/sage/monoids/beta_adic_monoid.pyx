# coding=utf8
r"""
Beta-adic Monoids tools.

Beta-adic monoids are finitely generated monoids with generators of the form

    - :math:`x -> \beta  x + c`

where beta is a element of a field (for example a complex number),
and c is varying in a finite set of numerals.
It permits to describe beta-adic expansions, that is writing of numbers of the form

    - :math:`x = c_{0} + c_{1} \beta + c_{2} \beta^{2} +  ...`
      for :math:`c_{i}` 's in a finite set of numerals.

AUTHORS:

- Paul Mercat (2013) -  I2M AMU Aix-Marseille Universite -initial version
- Dominique Benielli (2018) - Labex Archimede - I2M -
  AMU Aix-Marseille Universite - Integration in SageMath

EXAMPLES::

    sage: pi = x^3-x^2-x-1 # Tribonacci
    sage: b = pi.roots(ring=QQbar)[1][0]
    sage: m = BetaAdicMonoid(b, {0,1})
    sage: print(m)
    Monoid of b-adic expansion with b root of x^3 - x^2 - x - 1 and numerals set {0, 1}
    sage: pp = m.b.parent().places()[0]
    sage: print(pp)
    Ring morphism:
      From: Number Field in b with defining polynomial x^3 - x^2 - x - 1
      To:   Real Field with 106 bits of precision
      Defn: b |--> 1.839286755214161132551852564671
    sage: pm = m.b.parent().places()[1]
    sage: print(pm)
    Ring morphism:
      From: Number Field in b with defining polynomial x^3 - x^2 - x - 1
      To:   Complex Field with 53 bits of precision
      Defn: b |--> -0.419643377607080 + 0.606290729207199*I
    sage: ared = m.reduced_words_automaton2()
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
from monoid import Monoid_class
from sage.rings.qqbar import QQbar
from sage.rings.padics.all import *
from libc.stdlib cimport malloc, free

#from sage.structure.factory import UniqueFactory
#from sage.misc.cachefunc import cached_method
from cysignals.signals cimport sig_on, sig_off
cimport sage.combinat.words.cautomata
from sage.combinat.words.cautomata cimport DetAutomaton, FreeAutomaton
from sage.combinat.words.cautomata_generators import DetAutomatonGenerators
from sage.rings.integer import Integer
from sage.rings.number_field.all import *
# from sage.structure.parent_gens import normalize_names
# from free_monoid_element import FreeMonoidElement


# from sage.combinat.words.automata import Automaton
# from sage.combinat.words.cautomata import DetAutomaton

# from sage.structure.factory import UniqueFactory
# from sage.misc.cachefunc import cached_method
from cysignals.signals cimport sig_on, sig_off


# garde la composante fortement connexe de 0
# def prune(a, K):
#     """
#     Return the strongly connex component
#
#     INPUT:
#
#     - ``a`` a tree
#     - ``K``
#
#     OUTPUT:
#
#     the the strongly  connex component of a which
#     correspond the the K zeros
#
#     EXAMPLES::
#
#         sage:
#         sage: prune()
#     """
#     for s in a.strongly_connected_components_subgraphs():
#         if K.zero() in s:
#             return s

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

cdef extern from "relations.h":
    cdef cppclass Element:
        int *c  # liste des n coeffs

    cdef cppclass PlaceArch:
        Complexe *c  # 1, b, b^2, ... pour cette place

    # structure contenant les infos nécessaires pour calculer l'automate des relations
    cdef cppclass InfoBetaAdic:
        int n         # degre
        Element bn    # expression of b^n as a polynome in b of degree < n
        Element b1    # expression of 1/b as a polynome in b of degree < n
        Element *c    # list of figures used for the calculation of  relations' automaton
        int nc        # number of figures
        int ncmax     # number of allocated figures
        PlaceArch *p  # list of na places
        double *cM    # square of valeurs absolues max
        int na        # number of va

    Element NewElement(int n)
    void FreeElement(Element e)
    InfoBetaAdic allocInfoBetaAdic(int n, int na, int ncmax, bool verb)
    void freeInfoBetaAdic(InfoBetaAdic iba)
    Automaton RelationsAutomaton(InfoBetaAdic iba2, bool isvide, bool ext, bool verb)
    Automaton RelationsAutomatonT(InfoBetaAdic iba2, Element t, bool isvide, bool ext, bool verb)

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
        Complexe* t  # liste des translations
        int n        # nombre de translations
        Automaton a
    cdef cppclass BetaAdic2:
        Complexe b
        Complexe* t  # liste des translations
        int n        # nombre de translations
        Automaton* a
        int na
    ctypedef Color* ColorList
    #    cdef cppclass SDLImage:
    #        void *img

    void* OpenImage(const char *file_name)
    bool InImage(void* img, int x, int y)
    int ImageWidth(void *img)
    int ImageHeight(void *img)
    void CloseImage(void* img)
    void TestSDL()
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
    void DrawZoom(BetaAdic b, int sx, int sy, int n, int ajust, Color col, double coeff, int verb)
    Automaton UserDraw(BetaAdic b, int sx, int sy, int n, int ajust, Color col, int only_pos, int verb)
    #  void WordZone (BetaAdic b, int *word, int nmax)
    int *WordDrawn()
    void Draw(BetaAdic b, Surface s, int n, int ajust, Color col, double coeff, int verb)
    void Draw2(BetaAdic b, Surface s, int n, int ajust, Color col, int verb)
    void DrawList(BetaAdic2 b, Surface s, int n, int ajust, ColorList lc, double alpha, int verb)
    void print_word(BetaAdic b, int n, int etat)

# calcul de la valeur absolue p-adique (car non encore implémenté autrement)
def absp(c, p, d):
    """
    Computation of the p-adic absolute value (not yet implemented)

    INPUT:

    - ``c`` first argument

    - ``p`` second argument

    - ``d`` third argument

    OUTPUT:

    The p-adic absolute value

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
                                   Cd=None, plus=True,
                                   verb=False) except *:
    # compute all the data in sage
    # K = NumberField((1/self.b).minpoly(), 'b', embedding=QQbar(1/self.b))
    b = self.b
    K = b.parent()
    #    b = K.gen()
    #    C = [c.lift()(1/b) for c in self.C]
    C = self.C
    if verb:
        print("C = %s" % C)

    if verb:
        print(K)

    # détermine les places qu'il faut considérer
    parch = []
    for p in K.places():  # places archimédiennes
        if plus:
            if p(b).abs() > 1:
                parch += [p]
        else:
            if p(b).abs() < 1:
                parch += [p]
    pi = K.defining_polynomial()
    from sage.arith.misc import gcd
    # rend le polynôme à coefficients entiers et de contenu 1
    pi = pi / gcd(pi.list())
    if verb:
        print("pi=%s" % pi)
    # liste des nombres premiers concernés
    lp = (Integer(pi.list()[0])).prime_divisors()
    if verb:
        print("lp=%s" % lp)
    # liste des place  s ultramétriques considérées
    pultra = []
    for p in lp:
        # détermine toutes les places au dessus de p dans le corps de nombres K
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

    self.parch = parch

    if (len(pultra) > 0):
        raise ValueError("Not implemented for b algebraic non-integer.")
    # calcule les bornes max pour chaque valeur absolue
    if Cd is None:
        Cd = Set([c-c2 for c in C for c2 in C])
    else:
        #        Cd = [K(c).lift()(1/b) for c in Cd]
        Cd = [K(c) for c in Cd]
    if verb:
        print("Cd = %s" % Cd)

    #    m = dict([])
    #    for p in parch:
    #        m[p] = max([p(c).abs() for c in Cd])/abs(1-p(b).abs())
    #    for p, d in pultra:
    #        m[p] = max([absp(c, p, d) for c in Cd])
    #    if verb: print "bornes : %s"%m
    # convert the data in C
    n = K.degree()
    na = len(parch)
    ncmax = len(set([c-c2 for c in self.C for c2 in self.C]))
    cdef InfoBetaAdic i
    if verb:
        print("alloc...")
    sig_on()
    i = allocInfoBetaAdic(n, na, ncmax, verb)
    sig_off()
    cdef int j
    # initialise bn
    if verb:
        print("init bn...")
    getElement(b**n, i.bn, n)
    # initialise b1
    if verb:
        print("init b1...")
    getElement(1/b, i.b1, n)
    # initialise les places
    if verb:
        print("init places...")
    for k in range(na):
        for j in range(n):
            i.p[k].c[j] = complex(parch[k](b**j))
    # initialise les chiffres et bornes
    if verb:
        print("init chiffres...")
    initCdInfoBetaAdic(self, &i, Cd=Cd, verb=verb)
    return i

cdef initCdInfoBetaAdic(self, InfoBetaAdic *i, Cd, verb=False):
    # recalcule les bornes max pour chaque valeur absolue
    #    if Cd is None:
    #        Cd = Set([c-c2 for c in self.C for c2 in self.C])
    Cd = list(Cd)
    if verb:
        print("Cd = %s" % Cd)
    m = dict([])
    for p in self.parch:
        m[p] = max([p(c).abs() for c in Cd])/abs(1-p(self.b).abs())
    if verb:
        print("bornes : %s" % m)
    #    for p, d in self.pultra:
    #        m[p] = max([absp(c, p, d) for c in Cd])
    # conversion en C
    i.nc = len(Cd)
    if i.nc > i.ncmax:
        raise ValueError("Trop de chiffres : %d > %d max (initialiser le BetaAdicMonoid avec plus de chiffre)."%(i.nc, i.ncmax))
    for j, c in enumerate(Cd):
        getElement(c, i.c[j], i.n)
    for j, p in enumerate(self.parch):
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
    arr = np.zeros([s.sy, s.sx], dtype=[('r', 'uint8'), ('g', 'uint8'),
                                        ('b', 'uint8'), ('a', 'uint8')])
    cdef int x, y
    cdef Color c
    for x in range(s.sx):
        for y in range(s.sy):
            c = s.pix[x][s.sy - y - 1]
            arr[y, x][0] = c.r
            arr[y, x][1] = c.g
            arr[y, x][2] = c.b
            arr[y, x][3] = c.a
    return Image.fromarray(arr, 'RGBA')
    # img.save("/Users/mercat/Desktop/output.png")
    # img.save(file)

cdef Automaton getAutomate(a, d, list C, iss=None, verb=False):
    cdef int i
    if verb:
        print("getAutomate %s..." % a)
    cdef DetAutomaton fa
    if isinstance(a, DetAutomaton):
        fa = a
        fa.permut_op(C, verb=verb)
        # fa = fa.permut(C, verb=verb)
        return fa.a[0]
    else:
        raise ValueError("DetAutomaton expected.")

#    #assume in the following that a is a Automaton
#    lv = a.vertices()
#    if hasattr(a, 'F'):
#        F = a.F
#    else:
#        F = lv
#    #alloue l'automate
#    cdef Automate r
#    r = NewAutomaton(a.num_verts(), len(a.alphabet))
#    #réindice les sommets
#    dv = {}
#    for u,i in zip(lv, range(len(lv))):
#        dv[u] = i
#        if u in F:
#            r.e[i].final = 1
#    if verb:
#        print len(lv)
#    #copie l'automate en C
#    le = a.edges()
#    if verb:
#        print "len(le)=%s"%len(le)
#    for u,v,l in le:
#        if d.has_key(l):
#            #if dv.has_key(u) and dv.has_key(v):
#            r.e[dv[u]].f[d[l]] = dv[v]
#            #else:
#            #   print "Erreur : pas de clef %s ou %s !"%(u,v)
#        else:
#            print "Erreur : pas de clef %s !"%l
#    if verb:
#        print "I..."
#    if iss is not None:
#        r.i = iss
#    else:
#        if hasattr(a, 'I') and len(a.I) > 0:
#            r.i = dv[list(a.I)[0]]
#        else:
#            r.i = -1
#            #raise ValueError("The initial state must be defined !")
#    if verb:
#        print "...getAutomate"
#    return r

cdef BetaAdic getBetaAdic(input_a, prec=53, ss=None, tss=None, iss=None,
                          transpose=True, add_letters=True, verb=False):
    from sage.rings.complex_field import ComplexField
    CC = ComplexField(prec)
    cdef BetaAdic b
    if ss is None:
        if hasattr(input_a, 'ss'):
            ss = input_a.ss
        else:
            ss = input_a.default_ss()
    else:
        if transpose and tss is None:
            if verb:
                print("Transpose computation...\n")  # "Calcul de la transposée...\n"
            tss = ss.mirror().determinize()
            self.tss = tss
            if verb:
                print(tss)
    if transpose:
        if tss is None:
            if hasattr(input_a, 'tss'):
                tss = input_a.tss
            else:
                if verb:
                    print("Transpose computation...\n")  # "Calcul de la transposée...\n"
                tss = ss.mirror().determinize()
                input_a.tss = tss
                if verb:
                    print(tss)
        a = tss
    else:
        a = ss

    if add_letters:
        C = set(input_a.C)
        C.update(a.alphabet)
    else:
        C = a.alphabet
    C = list(C)

    b = NewBetaAdic(len(C))
    b.b = complex(CC(input_a.b))
    d = {}
    for i, c in zip(range(b.n), C):
        b.t[i] = complex(CC(c))
        d[c] = i
    # automaton
    # if isinstance(a, DetAutomaton):
    #    a = a.permut(C, verb=verb)
    b.a = getAutomate(a, d, C=C, iss=iss, verb=verb)
    return b

cdef BetaAdic2 getBetaAdic2(input_a, la=None, ss=None, tss=None,
                            prec=53, add_letters=True, verb=False):
    if verb:
        print("getBetaAdic %s" % input_a)
    from sage.rings.complex_field import ComplexField
    #from sage.combinat.words.automata import Automaton
    CC = ComplexField(prec)
    cdef BetaAdic2 b
    if la is None:
        la = input_a.get_la(ss=ss, tss=tss, verb=verb)

    if add_letters:
        C = set(input_a.C)
        for a in la:
            C.update(a.alphabet)
    C = list(C)

    b = NewBetaAdic2(len(C), len(la))
    b.b = complex(CC(input_a.b))
    d = {}
    for i, c in zip(range(b.n), C):
        b.t[i] = complex(CC(c))
        d[c] = i
    # automata
    for i in range(len(la)):
        # if isinstance(la[i], DetAutomaton):
        #    la[i] = la[i].permut(C, verb=verb)
        b.a[i] = getAutomate(la[i], d, C=C, iss=None, verb=verb)
    return b

def PrintWord(m, n):
    """
    Print of beta adic 

    INPUT:

    - ``m`` first word argument

    - ``n`` second word argument


    OUTPUT:

    Print the word

    TESTS:

        sage:import sage.monoids.beta_adic_monoid as mn

    """
    b = getBetaAdic(m, prec=53, ss=None, tss=None, iss=None, transpose=False,
                    add_letters=True, verb=False)
    print_word(b, n, b.a.i)


# ##used by compute_substitution()
# donne la liste des feuilles du sous-arbre partant de e
def fils(tree, e):
    """
    Return the list of sheet's sub-tree  starting on e.

    INPUT:

    - ``tree`` the tree.
    - ``e`` the starting sheet.

    OUTPUT:

    list of ``tree`` sheet's sub-tree  starting on e.

    """
    if tree[e] == []:  # e is a
        return [e]
    r = []
    for f in tree[e]:
        r += fils(tree, f)
    return r


# teste si a est inclus dans un des morceaux de l ou pas
def included(a, l, lm):
    """
    Return the index of

    INPUT:

    - ``a`` word to find in ``l``.

    - ``l`` including word to test

    OUTPUT:

    - the word in ``l`` if a is present

    - ``True`` if the automata is empty

    - ``None``
    
    

    """
    # teste vite fait si l'on est inclus dans un morceau ou pas
    incl = False
    w = a.find_word()
    if w is None:
        print("Error : empty automata!")
        return True
    lr = []
    for j in l:
        a2 = lm[j][0]
        if a2.rec_word(w):
            if a.included(a2):
                return j
            else:
                return None
    print("******* Error : word %s is conconize by any automata ! *********" % w)
    return None


# split a1 selon ba (rend un couple (a11, a12) avec a11 la partie
# dans ba et a12 celle disjointe de ba)
def split_ba(i, tr, np, lm, m, aa, ap, verb=False):
    b = m.b
    a1 = lm[i][0]
    # teste l'intersection avec ba
    # at = m.move2(t=(b**(-np))*tr, a=aoc)
    # TODO : utiliser les automates des relations précalculés
    # pour les translations de l'échange
    at = m.Proj(aa, ap, t=(b**(-np))*tr)
    if at.intersect(a1):
        ar = at.intersection(a1)
        # détermine si l'on est inclus dans ba
        ar.zero_completeOP()
        if ar.equals_langages(a1):
            return (a1, None)
        else:
            # on subdivise en deux
            ar2 = a1.intersection(ar.complementary())
            ar2.zero_completeOP()
            return (ar, ar2)
    else:
        return (None, a1)


# split a1 selon baoc (rend un couple (a11, a12) avec a11 la partie dans baoc
# et a12 celle disjointe de baoc)
def split_baoc(i, tr, np, lm, m, aoc, verb=False):
    b = m.b
    a1 = lm[i][0]
    # teste l'intersection avec baoc
    at = m.move2(t=(b**(-np))*tr, a=aoc)
    # TODO : utiliser les automates des relations précalculés
    # pour les translations de l'échange
    if at.intersect(a1):
        ar = at.intersection(a1)
        # détermine si l'on est inclus dans baoc
        ar.zero_completeOP()
        if ar.equals_langages(a1):
            return (a1, None)
        else:
            # on subdivise en deux
            return (ar, a1.intersection(ar.complementary()))
    else:
        return (None, a1)


cdef class ImageIn:
    cdef void** s

    def __cinit__(self):
        self.s = <void **>malloc(sizeof(void*))

    def __init__(self, file_name):
        self.s[0] = OpenImage(file_name)

    def __dealloc__(self):
        CloseImage(self.s[0])
        free(self.s)

    def __repr__(self):
        w = ImageWidth(self.s[0])
        h = ImageHeight(self.s[0])
        return "Image of size %sx%s" % (w, h)

    def __contains__(self, p):
        from sage.rings.complex_field import ComplexField
        CC = ComplexField(53)
        if p in CC:
            return InImage(self.s[0], p.real(), p.imag())
        else:
            return InImage(self.s[0], p[0], p[1])

    def width(self):
        return ImageWidth(self.s[0])

    def height(self):
        return ImageHeight(self.s[0])


class BetaAdicMonoid(Monoid_class):
    r"""
    It is the beta-adic monoid generated by the set of affine transformations
    `:math:`{x -> b*x + c | c \in \mathbb C}`.

    INPUT:

    - ``b`` -- coefficient ''b'' of affine transformation
      :math:`{x -> b*x + c | c \in \mathbb C}`

    - ``C` -- coefficient ''C''of affine transformation  
      :math:`{x -> b*x + c | c \in \mathbb C}`


    EXAMPLES::

        sage: m1 = BetaAdicMonoid(3, {0,1,3})
        sage: print(m1)
        Monoid of 3-adic expansion with numerals set {0, 1, 3}
        sage: m2 = BetaAdicMonoid((1 + sqrt(5)) / 2, {0,1})
        sage: print(m2)
        Monoid of b-adic expansion with b root of x^2 - x - 1 and numerals set {0, 1}
        sage: b = (x^3-x-1).roots(ring=QQbar)[0][0]
        sage: m3 = BetaAdicMonoid(b, {0,1})
        sage: print(m3)
        Monoid of b-adic expansion with b root of x^3 - x - 1 and numerals set {0, 1}

    """
    def __init__(self, b, C):
        r"""
        Construction of the b-adic monoid generated by the set of
        affine transformations ``{x -> b*x + c | c in C}``.

        EXAMPLES::

            sage: m1 = BetaAdicMonoid(3, {0,1,3})
            sage: m1
            Monoid of 3-adic expansion with numerals set {0, 1, 3}
            sage: m2 = BetaAdicMonoid((1+sqrt(5))/2, {0,1})
            sage: m2
            Monoid of b-adic expansion with b root of x^2 - x - 1 and numerals set {0, 1}
            sage: b = (x^3-x-1).roots(ring=QQbar)[0][0]
            sage: m3 = BetaAdicMonoid(b, {0,1})
            sage: m3
            Monoid of b-adic expansion with b root of x^3 - x - 1 and numerals set {0, 1}

        """
        # print "init BAM with (%s,%s)"%(b,C)
        if b in QQbar:
            #            print b
            pi = QQbar(b).minpoly()
            K = NumberField(pi, 'b', embedding=QQbar(b))
        else:
            K = b.parent()
            try:
                K.places()
            except:
                print("b=%s must be a algebraic number, ring %s not accepted." % (b, K))

#        print K
        self.b = K.gen()  # beta (element of an NumberField)
#        print "b="; print self.b
        self.C = Set([K(c) for c in C])  # set of numerals
#        print "C="; print self.C

    def gen(self, i):
        r"""
        Return the element of C of index i.

        INPUT:

        - ``i`` -- index to return

        OUTPUT:

        Return the element of C of index ``i``

        EXAMPLES::

            sage: m1 = BetaAdicMonoid(3, {0,1,3})
            sage: m1.gen(2)
            3
            sage: m1.gen(3)
            IndexError                                Traceback (most recent call last)
            ...
            IndexError: list index out of range

        """
        if i >= len(self.C):
            raise IndexError()

        return self.C[i]

    def ngens(self):
        r"""
        Return the number of elements of C.

        OUTPUT:

        Return the number of elements of C.

        EXAMPLES::

            sage: m1 = BetaAdicMonoid(3, {0,1,3})
            sage: m1.ngens()
            3

        """
        return len(self.C)

    def _repr_(self):
        r"""
        Returns the string representation of the beta-adic monoid.

        EXAMPLES::

            sage: BetaAdicMonoid((1+sqrt(5))/2, {0,1})
            Monoid of b-adic expansion with b root of x^2 - x - 1 and numerals set {0, 1}
            sage: BetaAdicMonoid(3, {0,1,3})
            Monoid of 3-adic expansion with numerals set {0, 1, 3}

        TESTS::

            sage: m=BetaAdicMonoid(3/sqrt(2), {0,1})
            sage: repr(m)
            'Monoid of b-adic expansion with b root of x^2 - 9/2 and numerals set {0, 1}'
        """

        str = ""
        if hasattr(self, 'ss'):
            if self.ss is not None:
                str = " with subshift of %s states" % self.ss.num_verts()

        from sage.rings.rational_field import QQ
        if self.b in QQ:
            return "Monoid of %s-adic expansion with numerals set %s" % (self.b, self.C) + str
        else:
            K = self.b.parent()
            if K.base_field() == QQ:
                return "Monoid of b-adic expansion with b root of %s and numerals set %s" % (self.b.minpoly(), self.C) + str
            else:
                if K.characteristic() != 0:
                    return "Monoid of b-adic expansion with b root of %s and numerals set %s, in characteristic %s"%(self.b.minpoly(), self.C, K.characteristic()) + str
                else:
                    return "Monoid of b-adic expansion with b root of %s and numerals set %s"%(K.modulus(),self.C) + str

    def _testSDL(self):
        """
        Open mode video for graphical representation

        TESTS::

        sage: b = (x^3-x-1).roots(ring=QQbar)[0][0]
        sage: m3 = BetaAdicMonoid(b, {0,1})
        sage: m3._testSDL()
        Video Mode: 800x600 32 bits/pixel
        """
        sig_on()
        TestSDL()
        sig_off()

#     def default_ss(self, C=None):
#         r"""
#         Returns the full subshift (given by an Automaton) corresponding
#         to the beta-adic monoid.
# 
#         EXAMPLES::
# 
#             sage: m=BetaAdicMonoid((1+sqrt(5))/2, {0,1})
#             sage: m.default_ss()
#             Finite automaton with 1 states
#         """
#         if C is None:
#             C = self.C
#         ss = Automaton()
#         ss.allow_multiple_edges(True)
#         ss.allow_loops(True)
#         ss.add_vertex(0)
#         for c in C:
#             ss.add_edge(0, 0, c)
#         ss.I = [0]
#         ss.F = [0]
#         ss.A = C
#         return ss

    def default_ss(self, C=None):
        r"""
        Returns the full subshift (given by an Automaton) corresponding
        to the beta-adic monoid.

        INPUT:

        - ``C`` -- (default -- None)  alphabet

        OUTPUT:
        Returns the full subshift (given by an DetAutomaton) corresponding
        to the beta-adic monoid.

        EXAMPLES::

            sage: m=BetaAdicMonoid((1+sqrt(5))/2, {0,1})
            sage: m.default_ss()
            DetAutomaton with 1 states and an alphabet of 2 letters

        """
        if C is None:
            C = self.C
        #ss = Automaton()
        ss = DetAutomatonGenerators().AnyWord(C)
        #ss.allow_multiple_edges(True)
        #ss.allow_loops(True)
        #ss.add_vertex(0)
        #for c in C:
        #    ss.add_edge(0, 0, c)
        #ss.I = [0]
        #ss.F = [0]
        #ss.A = C
        return ss

    # liste des automates donnant le coloriage de l'ensemble limite
    def get_la(self, ss=None, tss=None, verb=False):
        """
        Retun the ''la of betaAdic
        INPUT:

        - ``ss`` -- (default ''None'') ``DetAutomaton``

        -``tss`` -- (default ''None'') transition

        -``verb`` -- (default ''False'') set to ''True'' for verbose mode

        OUTPUT:
        Return the list of automaton giving color of limit set

        EXAMPLES::

            sage: m=BetaAdicMonoid((1+sqrt(5))/2, {0,1})
            sage: m.get_la()
            [DetAutomaton with 1 states and an alphabet of 2 letters,
             DetAutomaton with 1 states and an alphabet of 2 letters]
        """

        if hasattr(self, 'la'):
            return self.la
        if tss is None:
            if hasattr(self, 'tss'):
                tss = self.tss
                if not isinstance(tss, DetAutomaton):
                    tss = DetAutomaton(tss)
        if ss is None:
            if hasattr(self, 'ss'):
                ss = self.ss
                if not isinstance(ss, DetAutomaton):
                    ss = DetAutomaton(ss)
        if ss is None:
            if tss is None:
                tss = DetAutomaton(None).full(list(self.C))
                # raise ValueError("la, ss, or tss must be defined !")
            if verb:
                print("Compute the transposition of tss=%s..." % tss)
            ss = tss.mirror().determinize()  # ze()
            if verb:
                print(ss)
                print("simplify...")
            ss = ss.prune_inf().prune()  # prune0_simplify()
            if verb:
                print(ss)
        if tss is None:
            if ss is None:
                ss = DetAutomaton(None).full(list(self.C))  # self.default_ss()
                # raise ValueError("la, ss, or tss must be defined !")
            if verb:
                print("Compute the transposition of ss=%s..." % ss)
            tss = ss.mirror().determinize()  # ze()
            if verb:
                print(tss)
                print("simplify...")
            tss = tss.prune_inf().prune()  # 0_simplify()
            if verb:
                print(tss)

        # compute la
        a = {}
        for v in ss.states:  # ss.vertices():
            a[v] = ss.copy()  # Automaton(ss) ##
            a[v].set_final_states([v])
            if verb:
                print("Compute the transposition...")
            a[v] = a[v].mirror().determinize()  # ze()
            if verb:
                print(a[v])
                print("simplify...")
            a[v] = a[v].prune_inf().prune()
            if verb:
                print(a[v])
        return [tss]+a.values()

    def points_exact(self, n=None, ss=None, iss=None):
        r"""
        Returns a set of exacts values (in the number field of beta)
        corresponding to the drawing of the limit set of the beta-adic monoid.

        INPUT:

        - ``n`` - integer (default: ``None``)
          The number of iterations used to plot the fractal.
          Default values: between ``5`` and ``16`` depending on the number
          of generators.

        - ``ss`` - DetAutomaton (default: ``None``)
          The subshift to associate to the beta-adic monoid for this drawing.

        - ``iss`` - set of initial states of the automaton
          ss (default: ``None``)

        OUTPUT:

            list of exact values

        EXAMPLES::

            #. The dragon fractal::
            sage: e = QQbar(1/(1+I))
            sage: m=BetaAdicMonoid(e, {0,1})
            sage: print(m)
            Monoid of b-adic expansion with b root of x^2 - x + 1/2 and numerals set {0, 1}
            sage: P = m.points_exact()
            Give the begin state iss of the automaton ss !
            sage: P = m.points_exact(iss=0)
            sage: len(P)
            65536
        """
        # global co

        C = self.C
        K = self.b.parent()
        b = self.b
        ng = C.cardinality()

        if ss is None:
            if hasattr(self, 'ss'):
                ss = self.ss
            else:
                ss = self.default_ss()

        if iss is None:
            if hasattr(ss, 'I'):
                iss = [i for i in ss.I][0]
            if iss is None:
                print("Give the begin state iss of the automaton ss !")
                return

        if n is None:
            if ng == 2:
                n = 16
            elif ng == 3:
                n = 9
            else:
                n = 5
        if n == 0:
            # donne un point au hasard dans l'ensemble limite
            # co = co+1
            return [0]
        else:
            orbit_points = set()
            V = set([v for c in C for v in [ss.succ(iss, c)] if v is not None])
            orbit_points0 = dict()
            for v in V:
                orbit_points0[v] = self.points_exact(n=n-1, ss=ss, iss=v)
            for c in C:
                v = ss.succ(iss, c)
                if v is not None:
                    orbit_points.update([b*p+c for p in orbit_points0[v]])
                    # orbit_points = orbit_points.union([b*p+c for p in self.points_exact(n=n-1, ss=ss, iss=v)])
            # orbit_points0 = self.points_exact(n-1)
            # orbit_points = Set([])
            # for c in C:
            #    v = self.succ(i, c)
            #    if v is not None:
            #       orbit_points = orbit_points.union(Set([b*p+c for p in orbit_points0]))
        # print "no=%s"%orbit_points.cardinality()
        return orbit_points

    def points(self, n=None, place=None, ss=None, iss=None, prec=53):
        r"""
        Returns a set of values (real or complex) corresponding to the drawing
        of the limit set of the beta-adic monoid.

        INPUT:

        - ``n`` - integer (default: ``None``)
          The number of iterations used to plot the fractal.
          Default values: between ``5`` and ``16`` depending on the number
          of generators.

        - ``place`` - place of the number field of beta (default: ``None``)
          The place we should use to evaluate elements of the number field
          given by points_exact()

        - ``ss`` - DetAutomaton (default: ``None``)
          The subshift to associate to the beta-adic monoid for this drawing.

        - ``iss`` - set of initial states of the automaton
          ss (default: ``None``)

        - ``prec`` - precision of returned values (default: ``53``)

        OUTPUT:

            list of real or complex numbers

        EXAMPLES::

            #. The dragon fractal::

                sage: e = QQbar(1/(1+I))
                sage: m = BetaAdicMonoid(e, {0,1})
                sage: P = m.points(iss=0)     # long time (360 s)
                sage: len(P)
                32768
        """

        C = self.C
        K = self.b.parent()
        b = self.b
        ng = C.cardinality()

        if n is None:
            if ng == 2:
                n = 18
            elif ng == 3:
                n = 14
            elif ng == 4:
                n = 10
            elif ng == 5:
                n = 7
            else:
                n = 5
            from sage.functions.log import log
            n = int(5.2/-log(abs(self.b.n(prec=prec))))

        from sage.rings.complex_field import ComplexField
        CC = ComplexField(prec)
        if place is None:
            if abs(b) < 1:
                # garde la place courante
                # place = lambda x: CC(x.n())
                return [CC(c).n(prec) for c in self.points_exact(n=n, ss=ss, iss=iss)]
            else:
                # choisis une place
                places = K.places()
                place = places[0]
                for p in places:
                    if abs(p(b)) < 1:
                        place = p
                        # break

        # from sage.rings.qqbar import QQbar
        # from sage.rings.qqbar import QQbar, AA
        # if QQbar(self.b) not in AA:
        #    #print "not in AA !"
        #    return [(place(c).conjugate().N().real(), place(c).conjugate().N().imag()) for c in self.points_exact(n=n, ss=ss, iss=iss)]
        # else:
        #    #print "in AA !"
        #    return [place(c).conjugate().N() for c in self.points_exact(n=n, ss=ss, iss=iss)]
        return [place(c).n(prec) for c in self.points_exact(n=n, ss=ss, iss=iss)]

#          if n == 0:
#             #donne un point au hasard dans l'ensemble limite
#             return [0]
#         else:
#             orbit_points0 = self.points(n-1)
#             orbit_points = Set([])
#             for c in C:
#                 orbit_points = orbit_points.union(Set([place(b)*p+place(c) for p in orbit_points0]))
#         return orbit_points

    def user_draw(self, n=None, tss=None, ss=None, iss=None,
                  sx=800, sy=600, ajust=True, prec=53, color=(0, 0, 0, 255),
                  method=0, add_letters=True, only_pos=False, verb=False):
        r"""
        Returns a set of values (real or complex) corresponding to the drawing
        of the limit set of the beta-adic monoid.

        INPUT:

        - ``n`` - integer (default: ``None``)
          The number of iterations used to plot the fractal.
          Default values: between ``5`` and ``16`` depending on the number
          of generators.

        -``tss`` -- (default ''None'') transition

        - ``ss`` - DetAutomaton (default: ``None``)
          The subshift to associate to the beta-adic monoid for this drawing.

        - ``iss`` - set of initial states of the automaton
          ss (default: ``None``)

        - ``sx``  -- (default 800)

        - ``sy``  -- (default 600)

        - ``ajust``  -- (default ``True``)

        - ``prec``  precision of returned values -- (default: ``53``)

        - ``color`` tuple of color in RGB values -- (default: (0, 0, 0, 255))

        - ``method`` -- (default 0)

        - ``add_letters`` -- (default ``True``)

        - ``only_pos`` -- (default ``False``)

        - ``verb`` -- (default ``False``) set ti ``True`` for verbose mod

        OUTPUT:

        list of real or complex numbers

        EXAMPLES::

            #. The dragon fractal::

                sage: e = QQbar(1/(1+I))
                sage: m = BetaAdicMonoid(e, {0,1})
                sage: P = m.user_draw()     # long time (360 s)

        """
        if tss is None:
            tss = self.reduced_words_automaton2()
        cdef BetaAdic b
        cdef Automaton a
        cdef DetAutomaton r
        sig_on()
        b = getBetaAdic(self, prec=prec, tss=tss, ss=ss, iss=iss,
                        add_letters=add_letters, transpose=True, verb=verb)
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
        if method == 0:
            a = UserDraw(b, sx, sy, n, ajust, col, only_pos, verb)
        elif method == 1:
            print("Not implemented !")
            return
        sig_off()
        r = DetAutomaton(None)
        r.a[0] = a
        r.A = list(self.C)
        return r

    def draw_zoom(self, n=None, tss=None, ss=None, iss=None,
                  sx=800, sy=600, ajust=True, prec=53, color=(0, 0, 0, 255),
                  method=0, add_letters=True, coeff=8., verb=False):
        r"""
        Returns a set of values (real or complex) corresponding to the drawing
        of the limit set of the beta-adic monoid.

        INPUT:

        - ``n`` - integer (default: ``None``)
          The number of iterations used to plot the fractal.
          Default values: between ``5`` and ``16`` depending on the number
          of generators.

        -``tss`` -- (default ''None'') transition

        - ``ss`` - DetAutomaton (default: ``None``)
          The subshift to associate to the beta-adic monoid for this drawing.

        - ``iss`` - set of initial states of the automaton
          ss (default: ``None``)

        - ``sx``  -- (default 800)

        - ``sy``  -- (default 600)

        - ``ajust``  -- (default ``True``)

        - ``prec``  precision of returned values -- (default: ``53``)

        - ``color`` tuple of color in RGB values -- (default: (0, 0, 0, 255))

        - ``method`` -- (default 0)

        - ``add_letters`` -- (default ``True``)

        - ``coeff`` -- (default 8.)

        - ``verb`` -- (default ``False``) set ti ``True`` for verbose mod

        OUTPUT:

        list of real or complex numbers

        EXAMPLES::

            #. The dragon fractal::

                sage: e = QQbar(1/(1+I))
                sage: m = BetaAdicMonoid(e, {0,1})
                sage: P = m.draw_zoom()     # long time (360 s)

        """
        if tss is None:
            tss = self.reduced_words_automaton2()

        cdef BetaAdic b
        b = getBetaAdic(self, prec=prec, tss=tss, ss=ss, iss=iss,
                        add_letters=add_letters, transpose=True, verb=verb)
        # if verb:
        #     printAutomaton(b.a)
        # dessin
        cdef Color col
        col.r = color[0]
        col.g = color[1]
        col.b = color[2]
        col.a = color[3]
        if n is None:
            n = -1
        if method == 0:
            sig_on()
            DrawZoom(b, sx, sy, n, ajust, col, coeff, verb)
            sig_off()
        elif method == 1:
            print("Not implemented !")
            return


    # give a word corresponding to one of the previously drawn points
    def word_drawn(self):
        r"""
        Draw word, return the list of words

        OUTPUT:

        return list of words

        EXAMPLES::

            #. The dragon fractal::

                sage: e = QQbar(1/(1+I))
                sage: m = BetaAdicMonoid(e, {0,1})
                sage: P = m.draw_zoom()     # long time (360 s)

        """
        sig_on()
        cdef int *word = WordDrawn()
        sig_off()
        cdef int i
        res = []
        for i in xrange(1024):
            if word[i] < 0:
                break
            res.append(word[i])
        res.reverse()
        return res

    def plot2(self, n=None, tss=None, ss=None, iss=None, sx=800, sy=600,
              ajust=True, prec=53, color=(0, 0, 0, 255), method=0, add_letters=True,
              coeff=8., verb=False):
        r"""
        Draw the limit set of the beta-adic monoid (with or without subshift).

        INPUT:

        - ``n`` - integer (default: ``None``)
          The number of iterations used to plot the fractal.
          Default values: between ``5`` and ``16`` depending on the number of generators.

        - ``place`` - place of the number field of beta (default: ``None``)
          The place used to evaluate elements of the number field.

        - ``tss`` (default: ``None``) transition

        - ``ss`` - Automaton (default: ``None``)
          The subshift to associate to the beta-adic monoid for this drawing.

        - ``iss`` - set of initial states of the
          automaton ss (default: ``None``)

        - ``sx`` -- (default: 800) dimensions of the resulting in x dimension

        - ``sy`` -- (default : 600) dimensions of the resulting
          in y dimension image

        - ``ajust`` - adapt the drawing to fill all the image,
          with ratio 1 (default: ``True``)

        - ``prec`` - precision of returned values -- (default: ``53``)

        - ``color`` - list of three integer between 0
          and 255 (default: ``(0,0,255,255)``) Color of the drawing.

        - ``method`` -- (default : 0)

        - ``add_letters`` -- (default ``True``)

        - ``coeff`` -- (default 8.)

        - ``verb`` - bool (default: ``False``)
          Print informations for debugging.

        OUTPUT:

            A Graphics object.

        EXAMPLES::

        #. The dragon fractal::

            sage: e = QQbar(1/(1+I))
            sage: m = BetaAdicMonoid(e, {0,1})
            sage: m.plot2()     # long time
            ValueError                                Traceback (most recent call last)
            ValueError: DetAutomaton expected.

            Exception ValueError: ValueError('DetAutomaton expected.',) in 'sage.monoids.beta_adic_monoid.getAutomate' ignored

        #. The Rauzy fractal of the Tribonacci substitution::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: m = s.rauzy_fractal_beta_adic_monoid()
            sage: m.plot2()     # long time

        #. A non-Pisot Rauzy fractal::

            sage: s = WordMorphism({1:[3,2], 2:[3,3], 3:[4], 4:[1]})
            sage: m = s.rauzy_fractal_beta_adic_monoid()
            sage: m.b = 1/m.b
            sage: m.plot2(tss=m.ss)     # long time

        #. The dragon fractal and its boundary::

            sage: e = QQbar(1/(1+I))
            sage: m = BetaAdicMonoid(e, {0,1})
            sage: ssi = m.intersection_words(w1=[0], w2=[1])     # long time
            sage: m.plot2(tss=ssi)                               # long time
            sage: m.plot2()                                      # long time

        #. The "Hokkaido" fractal and its boundary::

            sage: s = WordMorphism('a->ab,b->c,c->d,d->e,e->a')
            sage: m = s.rauzy_fractal_beta_adic_monoid()
            sage: ssi = m.intersection_words(w1=[0], w2=[1])                  # long time
            sage: m.plot2(la=[la[0], ssi]+la[1:], colormap='gist_rainbow')    # long time

        #. A limit set that look like a tiling::

            sage: P=x^4 + x^3 - x + 1
            sage: b = P.roots(ring=QQbar)[2][0]
            sage: m = BetaAdicMonoid(b, {0,1})
            sage: m.plot2(19)                                   # long time

        """
        # cdef char *file
        # from sage.misc.temporary_file import tmp_filename
        # file_name = tmp_filename()
        # file = file_name
        # if verb:
        #    print("file=%s"%file_name)
        sig_on()
        cdef Surface s = NewSurface(sx, sy)
        cdef BetaAdic b
        if tss is not None:
            tss = tss.prune()
        b = getBetaAdic(self, prec=prec, tss=tss, ss=ss, iss=iss,
                        add_letters=add_letters, transpose=True, verb=verb)
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
        if method == 0:
            Draw(b, s, n, ajust, col, coeff, verb)
        elif method == 1:
            print("Not implemented !")
            return
            # lv = s.rauzy_fractal_projection_exact().values()
            # for i,v in zip(range(len(lv)),lv):
            #        b.t[i] = complex(CC(v))
            # Draw2(b, s, n, ajust, col, verb)
        # enregistrement du résultat
        im = surface_to_img(s)
        if verb:
            print("Free...")
        FreeSurface(s)
        #  if not isinstance(tss, DetAutomaton):
        #  FreeAutomaton(&b.a)
        FreeBetaAdic(b)
        sig_off()
        # from PIL import Image
        return im  # Image.open(file_name+'.png')

    def plot3(self, n=None, la=None, ss=None, tss=None,
              sx=800, sy=600, ajust=True, prec=53, colormap='hsv',
              backcolor=None, opacity=1., add_letters=True, verb=False):
        r"""
        Draw the limit set of the beta-adic monoid with colors.

        INPUT:

        - ``n`` - integer (default: ``None``)
          The number of iterations used to plot the fractal.
          Default values: between ``5`` and ``16`` depending on the number of generators.

        - ``place`` - place of the number field of beta (default: ``None``)
          The place we should use to evaluate elements of the number field.

        - ``ss`` - Automaton (default: ``None``)
          The subshift to associate to the beta-adic monoid for this drawing.

        - ``iss`` - set of initial states of the
          automaton ss (default: ``None``)

        - ``sx, sy`` - dimensions of the resulting
          image (default : ``800, 600``)

        - ``ajust`` - adapt the drawing to fill all the image, with
          ratio 1 (default: ``True``)

        - ``prec`` - precision of returned values (default: ``53``)

        - ``colormap`` - list of colors (default: ``hsv``)
          Colors of the drawing.

        - ``opacity``- float (default: ``1.``)
          Transparency of the drawing.

        - ``verb`` - bool (default: ``False``)
          Print informations for debugging.

        OUTPUT:

            A Graphics object.

        EXAMPLES::

        #. The dragon fractal::
            sage: e = QQbar(1/(1+I))
            sage: m=BetaAdicMonoid(e, {0,1})
            sage: m.plot3()     # long time

        #. The Rauzy fractal of the Tribonacci substitution::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: m = s.rauzy_fractal_beta_adic_monoid()
            sage: m.plot3()     # long time

        #. A non-Pisot Rauzy fractal::

            sage: s = WordMorphism({1:[3,2], 2:[3,3], 3:[4], 4:[1]})
            sage: m = s.rauzy_fractal_beta_adic_monoid()
            sage: m.b = 1/m.b
            sage: m.plot3(tss=m.ss)     # long time

        #. The dragon fractal and its boundary::

            sage: e = QQbar(1/(1+I))
            sage: m = BetaAdicMonoid(e, {0,1})
            sage: ssi = m.intersection_words(w1=[0], w2=[1])                               # long time
            sage: m.plot2(tss=ssi) #plot the boundary                                      # long time
            sage: m.plot3(la=[m.default_ss(), ssi], colormap=[(.5,.5,.5,.5), (0,0,0,1.)])  # long time

        #. The "Hokkaido" fractal and its boundary::

            sage: s = WordMorphism('a->ab,b->c,c->d,d->e,e->a')
            sage: m = s.rauzy_fractal_beta_adic_monoid()
            sage: ssi = m.intersection_words(w1=[0], w2=[1])                # long time

            sage: m.plot3(la=[la[0], ssi]+la[1:], colormap='gist_rainbow')  # long time

        #. A limit set that look like a tiling::

            sage: P = x^4 + x^3 - x + 1
            sage: b = P.roots(ring=QQbar)[2][0]
            sage: m = BetaAdicMonoid(b, {0,1})
            sage: m.plot2(19)                                   # long time

        """
        if tss is not None:
            tss = tss.prune()
        cdef Surface s = NewSurface(sx, sy)
        cdef BetaAdic2 b
        sig_on()
        b = getBetaAdic2(self, la=la, ss=ss, tss=tss, prec=prec, add_letters=add_letters, verb=verb)
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
        if isinstance(colormap, list):
            # if b.na > len(colormap):
            #    raise ValueError("The list of color must contain at least %d elements."%b.na)
            for i in range(b.na):
                if i < len(colormap):
                    cl[i] = getColor(colormap[i])
                else:
                    cl[i] = randColor(255)
        elif isinstance(colormap, str):
            from matplotlib import cm
            if not colormap in cm.datad.keys():
                raise ValueError("Color map %s not known (type 'from matplotlib import cm' and look at cm for valid names)" % colormap)
            colormap = cm.__dict__[colormap]
            cl[0] = getColor(backcolor)
            for i in range(b.na-1):
                cl[i+1] = getColor(colormap(float(i)/float(b.na-1)))
        else:
            raise TypeError("Type of option colormap (=%s) must be list of colors or str" % colormap)
        DrawList(b, s, n, ajust, cl, opacity, verb)
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
        else:
            la[0] = la[0].prune()
#        else:
#            for i,a in enumerate(la):
#                if not isinstance(a, DetAutomaton):
#                    FreeAutomaton(&b.a[i])
        FreeBetaAdic2(b)
        FreeColorList(cl)
        sig_off()
        return im

    def plot(self, n=None, place=None, ss=None, iss=None,
             prec=53, point_size=None, color='blue', verb=False):
        r"""
        Draw the limit set of the beta-adic monoid (with or without subshift).

        INPUT:

        - ``n`` - integer (default: ``None``)
          The number of iterations used to plot the fractal.
          Default values: between ``5`` and ``16`` depending on the number of generators.

        - ``place`` - place of the number field of beta (default: ``None``)
          The place we should use to evaluate elements of the number field given by points_exact()

        - ``ss`` - Automaton (default: ``None``)
          The subshift to associate to the beta-adic monoid for this drawing.

        - ``iss`` - set of initial states of the automaton ss (default: ``None``)

        - ``prec`` - precision of returned values (default: ``53``)

        - ``point_size`` - real (default: ``None``)
          Size of the plotted points.

        - ``verb`` - bool (default: ``False``)
          Print informations for debugging.

        OUTPUT:

            A Graphics object.

        EXAMPLES::

        #. The dragon fractal::
            sage: e = QQbar(1/(1+I))
            sage: m=BetaAdicMonoid(e, {0,1})
            sage: m.plot()     # long time

        #. The Rauzy fractal of the Tribonacci substitution::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: m = s.rauzy_fractal_beta_adic_monoid()
            sage: m.plot()     # long time
            Graphics object consisting of 1 graphics primitive

        #. A non-Pisot Rauzy fractal::

            sage: s = WordMorphism({1:[3,2], 2:[3,3], 3:[4], 4:[1]})
            sage: m = s.rauzy_fractal_beta_adic_monoid()
            sage: m.b = 1/m.b
            sage: m.ss = m.ss.mirror().determinize().minimize()
            sage: m.plot()     # long time
            Graphics object consisting of 1 graphics primitive

        #. The dragon fractal and its boundary::

            sage: e = QQbar(1/(1+I))
            sage: m = BetaAdicMonoid(e, {0,1})
            sage: p1 = m.plot()      # long time
            Graphics object consisting of 1 graphics primitive
            sage: ssi = m.intersection_words(w1=[0], w2=[1])     # long time
            sage: p2 = m.plot(ss = ssi, n=18)                    # long time
            sage: p1+p2                                          # long time

        #. The "Hokkaido" fractal and its boundary::

            sage: s = WordMorphism('a->ab,b->c,c->d,d->e,e->a')
            sage: m = s.rauzy_fractal_beta_adic_monoid()
            sage: p1 = m.plot()                                     # long time
            Graphics object consisting of 1 graphics primitive
            sage: ssi = m.intersection_words(w1=[0], w2=[1])        # long time
            sage: p2 = m.plot(ss=ssi, n=40)                         # long time
            sage: p1+p2                                             # long time

        #. A limit set that look like a tiling::

            sage: P = x^4 + x^3 - x + 1
            sage: b = P.roots(ring=QQbar)[2][0]
            sage: m = BetaAdicMonoid(b, {0,1})
            sage: m.plot(18)                                    # long time
            Graphics object consisting of 1 graphics primitive

        """

        global co

        co = 0
        orbit_points = self.points(n=n, place=place, ss=ss, iss=iss, prec=prec)
        if verb:
            print("co=%s" % co)

        # Plot points size
        if point_size is None:
            point_size = 1

        # Make graphics
        from sage.plot.plot import Graphics
        G = Graphics()

        # dim = self.b.minpoly().degree()

        from sage.rings.qqbar import QQbar, AA
        if QQbar(self.b) not in AA:  # 2D plots
            from sage.all import points
            G = points(orbit_points, size=point_size, color=color)
        else:  # 1D plots
            from sage.all import plot
            G += plot(orbit_points, thickness=point_size, color=color)
#            if plotbasis:
#                from matplotlib import cm
#                from sage.plot.arrow import arrow
#                canonicalbasis_proj = self.rauzy_fractal_projection(eig=eig, prec=prec)
#                for i,a in enumerate(alphabet):
#                    x = canonicalbasis_proj[a]
#                    G += arrow((-1.1,0), (-1.1,x[0]), color=cm.__dict__["gist_gray"](0.75*float(i)/float(size_alphabet))[:3])
#        else:
#            print "dimension too large !"
        G.set_aspect_ratio(1)

        return G

    def _relations_automaton_rec(self, current_state, di, parch,
                                 pultra, m, Cd, ext, verb=False, niter=-1):
        r"""
        Used by relations_automaton()
        """

        if niter == 0:
            return di

        global count
        if verb and count % 10000 == 0:
            print(count)
        if count == 0:
            return di
        count -= 1

        b = self.b
        if not di.has_key(current_state):
            di[current_state] = dict([])
        # parcours les transitions partant de current_state
        for c in Cd:
            # e = b*current_state + c
            # calcule l'état obtenu en suivant la transition c
            e = (current_state + c)/b
            # calcule l'état obtenu en suivant la transition c
            # if verb: print "b=%s, e=%s, cur=%s, c=%s, di=%s"%(b, e, current_state, c, di)
            if not di.has_key(e):  # détermine si l'état est déjà dans le dictionnaire
                ok = True
                # calcule les valeurs abolues pour déterminer si l'état n'est pas trop grand
                for p in parch:
                    if not ext:
                        if p(e).abs() >= m[p]:
                            ok = False
                            break
                    else:
                        if p(e).abs() > m[p]+.000000001:
                            ok = False
                            break
                if not ok:
                    continue  # cesse de considérer cette transition
                for p, d in pultra:
                    if absp(e, p, d) > m[p]:
                        # if verb: print "abs(%s)=%s trop grand !"%(e, absp(e, p, d))
                        ok = False
                        break
                if ok:
                    # on ajoute l'état et la transition à l'automate
                    di[current_state][e] = c
                    di = self._relations_automaton_rec(current_state=e, di=di, parch=parch, pultra=pultra, m=m, Cd=Cd, ext=ext, verb=verb, niter=niter-1)
            else:
                # ajoute la transition
                di[current_state][e] = c
        return di

#     def relations_automaton(self, ext=False, ss=None, noss=False, Cd=None,
#                             verb=False, step=100, limit=None, niter=None):
#         r"""
#         Compute the relations automaton of the beta-adic monoid
#         (with or without subshift).
#         See http://www.latp.univ-mrs.fr/~paul.mercat/Publis/
#         Semi-groupes%20fortement%20automatiques.pdf for a definition
#         of such automaton (without subshift).
# 
#         INPUT:
# 
#         - ``ext`` - bool (default: ``False``)
#           If True, compute the extended relations automaton (which permit to describe infinite words in the monoid).
# 
#         - ``ss`` - Automaton (default: ``None``)
#           The subshift to associate to the beta-adic monoid for this operation.
# 
#         - ``noss`` - bool (default: ``False``)
# 
# 
#         - ``verb`` - bool (default: ``False``)
#           If True, print informations for debugging.
# 
#         - ``step`` - int (default: ``100``)
#           Stop to an intermediate state of the computing to verify that all is right.
# 
#         - ``limit``- int (default: None)
#           Stop the computing after a number of states limit.
# 
#         OUTPUT:
# 
#         A Automaton.
# 
#         EXAMPLES:::
# 
#             sage: m = BetaAdicMonoid(3, {0,1,3})
#             sage: m.relations_automaton()
#             Finite automaton with 3 states
# 
#             sage: b = (x^3-x-1).roots(ring=QQbar)[0][0]
#             sage: m = BetaAdicMonoid(b, {0,1})
#             sage: m.relations_automaton()
#             Finite automaton with 179 states
# 
#         REFERENCES:
# 
#         ..  [Me13] Mercat P.
#             Bull. SMF 141, fascicule 3, 2013.
# 
#         """
#         if not noss:
#             a = self.relations_automaton(ext=ext, ss=None,
#                                          noss=True, Cd=Cd, verb=verb,
#                                          step=step, limit=limit)
#             if not step:
#                 return a
#             step = step-1
#             if ss is None:
#                 if hasattr(self, 'ss'):
#                     ss = self.ss
#                 else:
#                     return a  # pas de sous-shift
#             if not step:
#                 return ss
#             step = step-1
#             d = dict()
#             for u in self.C:
#                 for v in self.C:
#                     if not d.has_key(u - v):
#                         d[u-v] = []
#                     d[u - v] += [(u, v)]
#             if not step:
#                 return d
#             step = step - 1
#             ss = ss.prune0_simplify()
#             P = ss.product(A=ss)
#             # P = P.prune0_simplify()
#             if not step:
#                 return P
#             step = step-1
#             a.relabel2(d)
#             if not step:
#                 return a
#             step = step-1
#             a = a.intersection(A=P)
#             if not step:
#                 return a
#             step = step-1
#             a = a.prune0_simplify()
#             if not step:
#                 return a
#             step = step-1
#             if not ext:
#                 a.pruneF()
#                 if not step:
#                     return a
#                 step = step - 1
#             # a = a.determinize(A=a.A, noempty=True)
#             # if not step:
#             #    return a
#             # step = step-1
#             # return a
#             return a.minimize()
# 
#         K = self.C[0].parent()
#         b = self.b
# 
#         if verb:
#             print(K)
# 
#         # détermine les places qu'il faut considérer
#         parch = []
#         for p in K.places():  # places archimédiennes
#             if p(b).abs() < 1:
#                 parch += [p]
#         pi = K.defining_polynomial()
#         from sage.arith.misc import gcd
#         # rend le polynôme à coefficients entiers et de contenu 1
#         pi = pi / gcd(pi.list())
#         # den = pi.constant_coefficient().denominator()
#         # lp = (pi.list()[pi.degree()].numerator()*den).prime_divisors()
#         # liste des nombres premiers concernés
#         lp = (Integer(pi.list()[0])).prime_divisors()
#         pultra = []  # liste des places ultramétriques considérées
#         for p in lp:
#             # détermine toutes les places au dessus
#             # de p dans le corps de nombres K
#             k = Qp(p)
#             Kp = k['a']
#             a = Kp.gen()
#             for f in pi(a).factor():
#                 kp = f[0].root_field('e')
# #                c = kp.gen()
#                 if kp == k:
#                     c = f[0].roots(kp)[0][0]
#                 else:
#                     c = kp.gen()
#                 if verb:
#                     print("c=%s (abs=%s)" % (c, (c.norm().abs())**(1/f[0].degree())))
#                 if (c.norm().abs())**(1/f[0].degree()) < 1:  # absp(c, c, f[0].degree()) > 1:
#                     pultra += [(c, f[0].degree())]
# 
#         if verb:
#             print("position: ")
#             print(parch)
#             print(pultra)
# 
#         # calcule les bornes max pour chaque valeur absolue
#         if Cd is None:
#             Cd = Set([c-c2 for c in self.C for c2 in self.C])
#         if verb:
#             print("Cd = %s" % Cd)
#         m = dict([])
#         for p in parch:
#             m[p] = max([p(c).abs() for c in Cd])/abs(1-p(p.domain().gen()).abs())
#         for p, d in pultra:
#             m[p] = max([absp(c, p, d) for c in Cd])
# 
#         if verb:
#             print("m = %s" % m)
# 
#         if verb:
#             print(K.zero().parent())
# 
#         global count
#         # print limit
#         if limit is None:
#             count = -1
#         else:
#             count = limit
#         if niter is None:
#             niter = -1
#         # print count
#         if verb:
#             print("Way...")
#         di = self._relations_automaton_rec(current_state=K.zero(), di=dict([]),
#                                            parch=parch, pultra=pultra, m=m,
#                                            Cd=Cd, ext=ext, verb=verb, niter=niter)
# 
#         if count == 0:
#             print("Maximum number of reached states.")
#         else:
#             if verb:
#                 if limit is None:
#                     print("%s covered states." % (-1 - count))
#                 else:
#                     print("%s covered states." % (limit - count))
# 
#         # a = Automaton([K.zero()], [K.zero()], di)
# 
#         # if verb: print "di = %s"%di
# 
#         res = Automaton(di, loops=True)  # , multiedges=True)
# 
#         if verb:
#             print("Avant emondation : %s" % res)
# 
#         res.I = [K.zero()]
#         res.A = Cd  # Set([c-c2 for c in self.C for c2 in self.C])
#         if verb:
#             print("Emondation...")
#         if not ext:
#             res.F = [K.zero()]
#             res.prune()
#         else:
#             # res = res.prune0_simplify() #pour retirer les états puits
#             res.prune0()
#             res.F = res.vertices()
#         return res

#     def relations_automaton2(self, verb=False, step=100,
#                              limit=None, niter=None):
#         r"""
# 
#         Do the same as relations_automaton, but avoid recursivity in order
#         to avoid the crash of sage.
# 
#         INPUT:
# 
#         - ``verb`` - bool (default: ``False``)
#           If True, print informations for debugging.
# 
#         - ``step`` - int (default: ``100``)
#           Stop to an intermediate state of the computing to verify
#           that all is right.
# 
#         - ``limit``- int (default: None)
#           Stop the computing after a number of states limit.
# 
#         OUTPUT:
# 
#         A Automaton.
#         """
# 
#         K = self.C[0].parent()
#         b = self.b
# 
#         if verb:
#             print(K)
# 
#         # détermine les places qu'il faut considérer
#         parch = []
#         for p in K.places():  # places archimédiennes
#             if p(b).abs() < 1:
#                 parch += [p]
#         pi = K.defining_polynomial()
#         from sage.arith.misc import gcd
#         # rend le polynôme à coefficients entiers et de contenu 1
#         pi = pi/gcd(pi.list())  
#         # liste des nombres premiers concernés
#         lp = (Integer(pi.list()[0])).prime_divisors()
#         pultra = []  # liste des places ultramétriques considérées
#         for p in lp:
#             # détermine toutes les places au dessus de p dans le corps de nombres K
#             k = Qp(p)
#             Kp = k['a']
#             a = Kp.gen()
#             for f in pi(a).factor():
#                 kp = f[0].root_field('e')
#                 if kp == k:
#                     c = f[0].roots(kp)[0][0]
#                 else:
#                     c = kp.gen()
#                 if verb:
#                     print("c=%s (abs=%s)" % (c, (c.norm().abs())**(1/f[0].degree())))
#                 if (c.norm().abs())**(1/f[0].degree()) < 1:  # absp(c, c, f[0].degree()) > 1:
#                     pultra += [(c, f[0].degree())]
# 
#         if verb:
#             print("positions: ")
#             print(parch)
#             print(pultra)
# 
#         # calcule les bornes max pour chaque valeur absolue
#         Cd = Set([c-c2 for c in self.C for c2 in self.C])
#         if verb:
#             print("Cd = %s" % Cd)
#         m = dict([])
#         for p in parch:
#             m[p] = max([p(c).abs() for c in Cd])/abs(1-p(p.domain().gen()).abs())
#         for p, d in pultra:
#             m[p] = max([absp(c, p, d) for c in Cd])
# 
#         if verb:
#             print("m = %s" % m)
# 
#         if verb:
#             print(K.zero().parent())
# 
#         if limit is None:
#             count = -1
#         else:
#             count = limit
#         if niter is None:
#             niter = -1
# 
#         if verb:
#             print("Way...")
# 
#         di = dict([])
#         S = [K.zero()]  # set of states to look at
#         iter = 0
#         while len(S) != 0:
#             if iter == niter:
#                 break
#             for s in S:
#                 S.remove(s)
#                 if not di.has_key(s):
#                     di[s] = dict([])
#                     if count % 10000 == 0:
#                         print(count)
#                     count -= 1
#                     if count == 0:
#                         iter = niter-1  # to break the main loop
#                         break
#                 for c in Cd:  # parcours les transitions partant de current_state
#                     e = (s + c)/b  # calcule l'état obtenu en suivant la transition c
#                     # if verb:
#                     # print "b=%s, e=%s, cur=%s, c=%s, di=%s"%(b, e, current_state, c, di)
#                     # détermine si l'état est déjà dans le dictionnaire
#                     if not di.has_key(e):
#                         ok = True
#                         # calcule les valeurs abolues pour déterminer si l'état n'est pas trop grand
#                         for p in parch:
#                             if p(e).abs() >= m[p]:
#                                 ok = False
#                                 break
#                         if not ok:
#                             continue  # cesse de considérer cette transition
#                         for p, d in pultra:
#                             if absp(e, p, d) > m[p]:
#                                 # if verb: print "abs(%s)=%s trop grand !"%(e, absp(e, p, d))
#                                 ok = False
#                                 break
#                         if ok:
#                             # on ajoute l'état et la transition à l'automate
#                             di[s][e] = c
#                             S.append(e)
#                         else:
#                             # ajoute la transition
#                             di[s][e] = c
#                             iter += 1
# 
#                 if count == 0:
#                     print("Maximum number of reached states.")
#                     return
#                 else:
#                     if verb:
#                         if limit is None:
#                             print("%s états parcourus." % (-1 - count))
#                         else:
#                             print("%s covered states." % (limit - count))
# 
#         res = Automaton(di, loops=True) # , multiedges=True)
# 
#         if verb:
#             print("Before emondation : %s" % res)
# 
#         res.I = [K.zero()]
#         res.A = Set([c-c2 for c in self.C for c2 in self.C])
#         res.F = [K.zero()]
#         if verb:
#             print("Emondation...")
#         res.prune()
#         return res

    def relations_automaton3(self, t=0, isvide=False, Cd=None,
                             ext=False, verb=False):
        r"""
        Compute the relation automaton of the beta-adic monoid.
        For beta algebraic integer only.
        If isvide is True, it only checks if the automaton is trivial or not.
        Cd is the set of differences A-B where A and B
        are the alphabets to compare.
        t is the translation of one of the side
        (initial state of the automaton).
        ext : automate des relations à l'infini ou pas.


        TESTS::

            sage: pi = x^3-x^2-x-1
            sage: b = pi.roots(ring=QQbar)[1][0]
            sage: m = BetaAdicMonoid(b, [0,1])
            sage: m .relations_automaton3()
            DetAutomaton with 7 states and an alphabet of 3 letters
        """
        if Cd is None:
            Cd = Set([c-c2 for c in self.C for c2 in self.C])
        Cd = list(Cd)
        cdef InfoBetaAdic ib
        sig_on()
        ib = initInfoBetaAdic(self, Cd=Cd, plus=True, verb=verb)
        cdef Automaton a
        a = RelationsAutomaton(ib, isvide, ext, verb)
        r = DetAutomaton(None)
        r.a[0] = a
        r.A = list(Cd)
        freeInfoBetaAdic(ib)
        sig_off()
        if isvide:
            return a.na != 0
        if ext:
            r2 = r.prune_inf()
            r2.set_final_states(r2.states())
        else:
            r2 = r.prune()
        return r2.mirror_det()

    def relations_automaton4(self, t=0, isvide=False, Cd=None, A=None, B=None,
                             couples=False, ext=False, transp=False,
                             prune=True, verb=False):
        r"""
        Compute the relation automaton of the beta-adic monoid.
        For beta algebraic integer only.
        If isvide is True, it only checks if the automaton is trivial or not.
        Cd is the set of differences A-B where A and B
        are the alphabets to compare.
        t is the translation of one of the side
        (initial state of the automaton).
        ext : automate des relations à l'infini ou pas.
        """
        if Cd is None:
            if A is None or B is None:
                Cd = Set([c-c2 for c in self.C for c2 in self.C])
            else:
                Cd = Set([a1-b1 for a1 in A for b1 in B])
        Cd = list(Cd)
        sig_on()
        cdef InfoBetaAdic ib
        ib = initInfoBetaAdic(self, Cd=Cd, plus=False, verb=verb)
        cdef Automaton a
        cdef Element e
        e = NewElement(ib.n)
        K = self.b.parent()
        t = K(t)
        getElement(t, e, ib.n)
        a = RelationsAutomatonT(ib, e, isvide, ext, verb)
        r = DetAutomaton(None)
        r.a[0] = a
        if verb:
            print("a (%s etats)" % a.n)
            print("Free element...")
        FreeElement(e)
        r.A = Cd
        if verb:
            print("Free InfoBetaAdic...")
        freeInfoBetaAdic(ib)
        sig_off()
        if isvide:
            return a.na != 0
        cdef DetAutomaton r2
        if prune:
            if verb:
                print("prune...")
            if ext:
                r2 = r.prune_inf()
                r2.set_final_states(r2.states())
            else:
                r2 = r.prune()
        else:
            r2 = r
        if transp:
            r2 = r2.mirror_det()
        if couples:
            if A is None or B is None:
                raise ValueError("Alphabets A and B must be defined !")
            d={}
            for c1 in A:
                for c2 in B:
                    if not d.has_key(c1-c2):
                        d[c1-c2] = []
                    d[c1-c2].append((c1, c2))
            if verb:
                print(d)
            r2 = r2.duplicate(d, verb=verb)
        if verb:
            print("return...")
        return r2

    def critical_exponent_aprox(self, niter=10, verb=False):
        b = self.b
        K = b.parent()
        C = self.C
        S = set([K.zero()])
        for i in range(niter):
            S2 = set([])
            for s in S:
                for c in C:
                    S2.add((s+c)/b)
            # intervertit S et S2
            S3 = S2
            S2 = S
            S = S3
            if verb:
                print(len(S))
        from sage.functions.log import log
        print("%s" % (log(len(S)).n() / (niter * abs(log(abs(b.n()))))))

    def complexity(self, verb=False):
        r"""
        Return a estimation of an upper bound of the number of states
        of the relations automaton.

        INPUT:

         - ``verb`` - Boolean (default: False) Display informations for debug.

        OUTPUT:

        A positive real number.

        EXAMPLES::

            sage: b = (x^3-x^2-x-1).roots(ring=QQbar)[1][0]
            sage: m = BetaAdicMonoid(b, [0,1])
            sage: m.complexity()
            108.523461211014
        """
        K = self.C[0].parent()
        b = self.b

        if verb:
            print(K)

        # détermine les places qu'il faut considérer
        parch = K.places()
        r = len(parch)
        pi = K.defining_polynomial()
        from sage.arith.misc import gcd
        pi = pi/gcd(pi.list())
        # rend le polynôme à coefficients entiers et de contenu 1
        # den = pi.constant_coefficient().denominator()
        # lp = (pi.list()[pi.degree()].numerator()*den).prime_divisors()
        # liste des nombres premiers concernés
        lp = (Integer(pi.list()[0])).prime_divisors()
        # liste des nombres premiers concernés
        pultra = []
        # liste des places ultramétriques considérées
        for p in lp:
            # détermine toutes les places au dessus de p dans le corps de nombres K
            k = Qp(p)
            Kp = k['a']
            a = Kp.gen()
            for f in pi(a).factor():
                kp = f[0].root_field('e')
#                c = kp.gen()
                if kp == k:
                    c = f[0].roots(kp)[0][0]
                else:
                    c = kp.gen()
                if verb:
                    print("c=%s (abs=%s)" % (c, (c.norm().abs())**(1 / f[0].degree())))
                if c.norm().abs() != 1:  # absp(c, c, f[0].degree()) > 1:
                    pultra += [(c, f[0].degree())]

        if verb:
            print("places: ")
            print(parch)
            print(pultra)

        # calcule les bornes max pour chaque valeur absolue
        Cd = Set([c-c2 for c in self.C for c2 in self.C])
        if verb:
            print("Cd = %s" % Cd)
        vol = 1.
        # from sage.rings.real_mpfr import RR
        for p in parch:
            if (p(b)).imag() == 0:
                vol *= 2 * max([p(c).abs() for c in Cd])/abs(1-p(b).abs())
                if verb:
                    print("real place  %s" % p)
            else:
                vol *= 3.1415926535*(max([p(c).abs() for c in Cd])/abs(1-p(b).abs()))**2
                if verb:
                    print("complex place %s" % p)
            # vol *= max([p(c).abs() for c in Cd])/abs(1-p(p.domain().gen()).abs())
            # vol *= max(1, max([p(c).abs() for c in Cd])/abs(1-p(p.domain().gen()).abs()))
        for p, d in pultra:
            vol *= max([(c.polynomial())(p).norm().abs() for c in Cd])
            # vol *= max([absp(c, p, d) for c in Cd])
            # vol *= max(1, max([absp(c, p, d) for c in Cd]))
        # from sage.functions.other import sqrt
        # return vol/(K.regulator()*(sqrt(r+1).n()))
        return vol

    # def infinite_relations_automaton (self, verb=False):
    #    a = self.relations_automaton (ext=True, verb=verb)
    #    #retire la composante fortement connexe de l'etat initial
    #    K = self.b.parent()
    #    for s in a.strongly_connected_components_subgraphs():
    #        if K.zero() in s:
    #            a.delete_vertices(a.strongly_connected_components_subgraphs()[0].vertices())
    #            return a
    # calcule le sous-shift correspondant à l'intersection des deux
    # monoïdes avec sous-shifts
    def intersection(self, ss, ss2=None, Iss=None,
                     Iss2=None, ext=True, verb=False):
        r"""
        Compute the intersection of two beta-adic monoid with subshifts

        INPUT:

        - ``ss``- Automaton (default: ``None``)
          The first subshift to associate to the beta-adic monoid for this operation.

        - ``ss2``- Automaton (default: ``None``)
          The second subshift to associate to the beta-adic monoid for this operation.

        - ``Iss``- set of states of ss (default: ``None``)

        - ``Iss2``- set of states of ss2 (default: ``None``)

        - ``ext`` - bool (default: ``True``)
          If True, compute the extended relations automaton (which permit to describe infinite words in the monoid).  

        - ``verb``- bool (default: ``False``)
          If True, print informations for debugging.

        OUTPUT:

        A Automaton.

        EXAMPLES::

            #. Compute the boundary of the dragon fractal (see intersection_words for a easier way) ::

                sage: e = QQbar(1/(1+I))
                sage: m = BetaAdicMonoid(e, {0,1})
                sage: ss=m.default_ss()
                sage: iss=ss.I[0]
                sage: ss0 = ss.prefix(w=[0], i=iss)
                sage: ss1 = ss.prefix(w=[1], i=iss)
                sage: ssi = m.intersection(ss=ss0, ss2=ss1)
                sage: ssd = ssi.determinize(A=m.C, noempty=True)
                sage: ssd = ssd.prune0_simplify()
                sage: m.plot(ss = ssd, n=19)     # long time
        """

        m = None

        if ss2 is None:
            if hasattr(self, 'ss'):
                ss2 = self.ss
            else:
                raise ValueError("Only one sub-shift given, I need two !")
            if type(ss) == type(BetaAdicMonoid(2, {0, 1})):
                m = ss
                ss = m.ss
                if m.b != self.b:
                    raise ValueError("Cannot compute the intersection of two beta-adic monoids with differents beta.")
                m.C = m.C.union(self.C)
                self.C = self.C.union(m.C)
                if hasattr(m, 'ss'):
                    m.ss.A = m.C
                else:
                    raise ValueError("Only one sub-shift given, I need two !")
                self.ss.A = self.C

        if Iss is None:
            if hasattr(ss, 'I'):
                Iss = ss.I
            if Iss is None:
                Iss = [ss.vertices()[0]]
        if Iss2 is None:
            if hasattr(ss2, 'I'):
                Iss2 = ss2.I
            if Iss2 is None:
                Iss2 = [ss2.vertices()[0]]
            if verb:
                print("Iss = %s, Iss2 = %s" % (Iss, Iss2))

        a = ss.product(ss2)
        if verb:
            print("Product = %s" % a)

        ar = self.relations_automaton(ext=ext, noss=True)
        if verb:
            print("Arel = %s" % ar)

        # maps actual edges to the list of corresponding couple
        m = dict([])
        for c in self.C:
            for c2 in self.C:
                if m.has_key(c - c2):
                    m[c-c2] += [(c, c2)]
                else:
                    m[c-c2] = [(c, c2)]
        if verb:
            print("m = %s" % m)

        L = a.alphabet  # a.edge_labels()
        LA = ar.alphabet  # ar.edge_labels()
        d = dict([])
        for u, v in L:
            for ka in LA:
                for u2, v2 in m[ka]:
                    if u == u2 and v == v2:
                        d[((u,v), ka)] = u
                        break
                    else:
                        d[((u,v), ka)] = None
        if verb:
            print("d = %s" % d)
        p = a.product(A=ar, d=d)
        # I = [((i,i2),self.b.parent().zero()) for i in Iss for i2 in Iss2]
        # if verb: print "I = %s"%I
        # p.pruneI(I=I)
        # if verb: print "%s"%p
        # p.prune0(I=I)
        p.I = [((i, i2), self.b.parent().zero()) for i, i2 in zip(Iss, Iss2)]
        p = p.prune0_simplify()
        if m is not None:
            ssd = p.determinize2(A=self.C, noempty=True)
            ssd = ssd.prune0_simplify()
            return ssd
        return p

    # calcule le sous-shift correspondant à l'intersection
    # des deux monoïdes avec sous-shifts, utilise des DetAutomaton
    def intersection2(self, ss1, ss2, verb=False):
        r"""
        Compute the intersection of two beta-adic monoid with subshifts given by DetAutomaton

        INPUT:

        - ``ss``- DetAutomaton (default: ``None``)
          The first subshift to associate to the beta-adic monoid for this operation.

        - ``ss2``- DetAutomaton (default: ``None``)
          The second subshift to associate to the beta-adic monoid for this operation.

        - ``verb``- bool (default: ``False``)
          If True, print informations for debugging.

        OUTPUT:

        A DetAutomaton.

        EXAMPLES::

            #. Compute the boundary of the dragon fractal (see intersection_words for a easier way) ::

                sage: e = QQbar(1/(1+I)) 
                sage: m = BetaAdicMonoid(e, {0,1})
                sage: import sage.combinat.words.cautomata
                sage: from sage.combinat.words.cautomata import DetAutomaton
                sage: ss0 = DetAutomaton([(0,1,0)]+[(1,1,l) for l in m.C], i=0, final_states=[1])
                sage: ss1 = DetAutomaton([(0,1,1)]+[(1,1,l) for l in m.C], i=0, final_states=[1])
                sage: ssi = m.intersection2(ss0, ss1)
                sage: m.plot2(tss = ssi)     # long time
        """
        a = self.relations_automaton3()
        a = a.prune_inf()
        a.set_final_states(a.states())
        ssp = ss1.product(ss2)
        ssp = ssp.prune()
        d = {}
        for (la1, la2) in ssp.alphabet:
            for lb in a.alphabet:
                if lb == la1-la2:
                    d[((la1, la2), lb)] = (la1, la2)
        ssi = ssp.product(a, d)
        ssi = ssi.prune_inf()
        ssi = ssi.prune()
        d = {}
        for (l1, l2) in ssi.alphabet:
            d[(l1, l2)] = l1
        ssi = ssi.determinize_proj(d)
        ssi = ssi.prune_inf()
        ssi = ssi.prune()
        return ssi.minimize()

    def intersection_words(self, w1, w2, ss=None, iss=None):
        r"""
        Compute the intersection of two beta-adic monoid with
        subshifts corresponding to two prefix

        INPUT:

        - ``w1``- word
          The first prefix.

        - ``w2``- word
          The second prefix.

        - ``ss``- Automaton (default: ``None``)
          The subshift to associate to the beta-adic monoid for this operation.

        - ``iss``- initial state of ss (default: ``None``)

        OUTPUT:

        A Automaton.

        EXAMPLES::

            #. Compute the boundary of the dragon fractal::

                sage: e = QQbar(1/(1+I))
                sage: m = BetaAdicMonoid(e, {0,1})
                sage: m.intersection_words(w1=[0], w2=[1])
                Finite automaton with 21 states

            #. Draw the intersection of two sub-sets of a limit set::

                sage: m = BetaAdicMonoid(1/(1+I), {0,1})
                sage: ssi = m.intersection_words(w1=[0], w2=[1])
                sage: m.plot(n=10, ss=ssi)                        # long time
        """
        if ss is None:
            if hasattr(self, 'ss'):
                ss = self.ss
            else:
                ss = self.default_ss()
        if iss is None:
            if hasattr(ss, 'I'):
                iss = ss.I[0]
            if iss is None:
                iss = ss.vertices()[0]
        ss1 = ss.prefix(w=w1, i=iss)
        ss2 = ss.prefix(w=w2, i=iss)
        ssi = self.intersection(ss=ss1, ss2=ss2)
        ssd = ssi.determinize2(A=self.C, noempty=True)
        ssd = ssd.prune0_simplify()
        return ssd

    #     - ``aut`` - DetAutomaton (default: ``None``, full language)
    #       Automaton describing the language in which we live.
    def reduced_words_automaton2(self, step=100,
                                 verb=False, transpose=False):  # , DetAutomaton aut=None):
        r"""
        Compute the reduced words automaton of the beta-adic monoid
        (without considering the automaton of authorized words).
        See http://www.latp.univ-mrs.fr/~paul.mercat/Publis/
        Semi-groupes%20fortement%20automatiques.pdf for a definition of such automaton.
        Use fast C-functions but works only for algebraic integer.
        (Consider using reduced_words_automaton() if you're not in this case.)

        INPUT:

        - ``verb`` - bool (default: ``False``)
          If True, print informations for debugging.

        OUTPUT:

        DetAutomaton.

        EXAMPLES::

            sage: pi = x^3-x^2-x-1
            sage: b = pi.roots(ring=QQbar)[1][0]
            sage: m = BetaAdicMonoid(b, {0,1})
            sage: ared = m.reduced_words_automaton2()
            sage: ared
            DetAutomaton with 4 states and an alphabet of 2 letters

        """

        # compute the relations automaton
        Cd = list(set([c-c2 for c in self.C for c2 in self.C]))
        Cdp = [k for k in range(len(Cd)) if Cd[k] in [self.C[j]-self.C[i] for i in range(len(self.C)) for j in range(i)]] #indices des chiffres strictements négatifs dans Cd
        arel = self.relations_automaton3(Cd=Cd, ext=False)
        arel = arel.prune()
        if transpose:
            arel = arel.mirror_det()
        if verb:
            print("arel = %s" % arel)
        if step == 1:
            return arel

        # add a new state
        cdef int ne, ei
        ei = arel.initial_state
        ne = arel.n_states  # new added state
        arel.add_state(True)
        arel.set_final_state(ei, final=False)  # it is the new final state
        if step == 2:
            return arel

        # add edges from the new state (copy edges from the initial state)
        for j in range(len(Cd)):
            arel.set_succ(ne, j, arel.succ(ei, j))
        if step == 3:
            return arel

        # suppress some edges from the initial state
        for j in Cdp:
            arel.set_succ(ei, j, -1)
        if step == 4:
            return arel

        # change edges that point to the initial state :
        # make them point to the new state
        for e in arel.states:
            if e != ei:
                for j in range(len(Cd)):
                    if arel.succ(e, j) == ei:
                        arel.set_succ(e, j, ne)
        if step == 5:
            return arel

        # project, determinise and take the complementary
        d = {}
        for a in self.C:
            for b in self.C:
                if not d.has_key(a - b):
                    d[a-b] = []
                d[a-b].append((a, b))
        if verb:
            print(d)
        arel = arel.duplicate(d)  # replace differences with couples
        d = {}
        for j in self.C:
            for i in self.C:
                d[(i, j)] = i
        if verb:
            print(d)
            print(arel)
        arel = arel.determinize_proj(d, noempty=False, nof=True)  # , verb=True)
        # project on the first value of the couple, determinise and take the complementary
        if verb:
            print(arel)
        arel = arel.prune()
        arel = arel.prune_inf()
        if step == 10:
            return arels
        return arel.minimize()

#     def reduced_words_automaton(self, ss=None, Iss=None, ext=False,
#                                 verb=False, step=None, arel=None):
#         r"""
#         Compute the reduced words automaton of the beta-adic monoid (with or without subshift).
#         See http://www.latp.univ-mrs.fr/~paul.mercat/Publis/Semi-groupes%20fortement%20automatiques.pdf for a definition of such automaton (without subshift).
# 
#         WARNING: It seems there is a bug : result may be incorrect if ss is not None.
# 
#         INPUT:
# 
#         - ``ss``- Automaton (default: ``None``)
#           The first subshift to associate to the beta-adic monoid for this operation.
# 
#         - ``Iss``- set of states of ss (default: ``None``)
# 
#         - ``ext`` - bool (default: ``True``)
#           If True, compute the extended relations automaton (which permit to describe infinite words in the monoid).  
# 
#         - ``verb`` - bool (default: ``False``)
#           If True, print informations for debugging.
# 
#         - ``step`` - int (default: ``None``)
#           Stop to a intermediate state of the computing to make verifications.
# 
#         - ``arel`` - Automaton (default: ``None``)
#           Automaton of relations.
# 
#         OUTPUT:
# 
#         A Automaton.
# 
#         EXAMPLES::
# 
#             #. 3-adic expansion with numerals set {0,1,3}::
# 
#                 sage: m = BetaAdicMonoid(3, {0,1,3})
#                 sage: m.reduced_words_automaton()
#                 Finite automaton with 2 states
# 
#             #. phi-adic expansion with numerals set {0,1}::
# 
#                 sage: m = BetaAdicMonoid((1+sqrt(5))/2, {0,1})
#                 sage: m.reduced_words_automaton()
#                 Finite automaton with 3 states
# 
#             #. beta-adic expansion with numerals set {0,1} where beta is the plastic number::
#                 sage: b = (x^3-x-1).roots(ring=QQbar)[0][0]
#                 sage: m = BetaAdicMonoid(b, {0,1})
#                 sage: m.reduced_words_automaton()        # long time
#                 Finite automaton with 5321 states
#         """
#         if ss is None:
#             if hasattr(self, 'ss'):
#                 ss = self.ss
#                 if hasattr(self.ss, 'I'):
#                     Iss = self.ss.I
# 
#         if step is None:
#             step = 1000
# 
#         K = self.C[0].parent()
# 
#         if verb:
#             print("Computation of relations's automata")
#             # "Calcul de l'automate des relations..."
# 
#         if arel is None:
#             a = self.relations_automaton(noss=True)
#         else:
#             a = arel
# 
#         if verb:
#             print(" -> %s" % a)
# 
#         if step == 1:
#             return ("relations's automata", a)
# 
#         # add a state copy of K.0 (it will be the new initial state)
#         a.add_vertex('O')
# 
#         #        #add transitions to K.0 to 'O'
#         #        for f, d, l in a.incoming_edges(K.zero(), labels=True):
#         #            if f == K.zero():
#         #                a.add_edge('O', 'O', l)
#         #            else:
#         #                a.add_edge(f, 'O', l)
# 
#         # subset of positives labels
#         Cdp = []
#         for i in range(self.C.cardinality()):
#             for j in range(i):
#                 Cdp += [self.C[i] - self.C[j]]
# 
#         # redirect positives transitions from K.0
#         for f, d, l in a.outgoing_edges(K.zero(), labels=True):
#             if l in Cdp:
#             #                a.delete_edge(K.zero(), d, l)
#                 # add the edge
#                 a.add_edge('O', d, l)
# 
#         a.add_edge('O', 'O', a.edge_label(K.zero(), K.zero()))
# 
#         if verb:
#             print(a.incoming_edges(K.zero(), labels=True))
# 
#         # remove outgoing edges from K.0 (except from K.0 to K.0)
#         for f, d, l in a.outgoing_edges(K.zero(), labels=True):
#             if f != d:
#                 a.delete_edge(f, d, l)
# 
#         if step == 2:
#             a.I = ['O']
#             a.F = Set([K.zero()])
#             return ("automaton of ordoned relations", a)
#         a.pruneI(I=['O'])
# 
#         if step == 3:
#             return ("pruned automaton of ordoned relations", a)
# 
#         if ss is not None:  # not full sub-shift
#             if Iss is None:
#                 Iss = [ss.vertices()[0]]
#             # maps actual edges to the list of corresponding couple
#             m = dict([])
#             for c in self.C:
#                 for c2 in self.C:
#                     if m.has_key(c - c2):
#                         m[c-c2] += [(c, c2)]
#                     else:
#                         m[c-c2] = [(c, c2)]
#             # if verb: print "m=%s"%m
# 
#             # calculate the 'product to the right' of a with ss
#             d = dict([])
#             La = a.edge_labels()
#             Lss = ss.edge_labels()
#             for ka in La:
#                 for kss in Lss:
#                     d[(ka, kss)] = None
#                     for k in m[ka]:
#                         if k[1] == kss:
#                             d[(ka, kss)] = k[0]
#                             break
# 
#             # if verb: print "d=%s"%d
#             if verb:
#                 # "avant produit : a=%s (%s etats)"%(a, a.num_verts())
#                 print(" before product : a=%s (%s states)" % (a, a.num_verts()))
#             a = a.product(A=ss, d=d)
#             if verb:
#                 print(" after product : a=%s" % a)
#             if step == 4:
#                 return ("non reduce general words automata", a) #"automate des mots généraux non réduits", a)
# 
#             I = [('O', i) for i in Iss]
#             nof = Set([(K.zero(), i) for i in ss.vertices()])
# 
#             # if verb: print "I=%s, F=%s"%(I, nof)
# 
#             if ext:
#                 # a.emondeI(I=I)
#                 # a.emonde0(I=I) #pour retirer les états puits
#                 a = a.emonde0_simplify(I=I)
#             else:
#                 a = a.emonde0_simplify(I=I)
#                 a.emonde(I=I, F=nof)
#             # a.emondeI(I=I)
#             # a.emondeF(F=nof)
#             # if step == 4:
#             #    return ("automate des mots généraux non réduits, émondé", a)
#             # a.emondeF(F=nof)
# 
#             if verb:
#                 print("After emondation : a=%s" % a)
#             if step == 5:
#                 return("emonded automaton of non reducted general words", a)
# 
#             # return a
#         else:
#             # maps actual edges to element of self.C (the left part when writted c-c2)
#             m = dict([])
#             for c in self.C:
#                 for c2 in self.C:
#                     if m.has_key(c-c2):
#                         m[c-c2] += [c]
#                     else:
#                         m[c-c2] = [c]
#             # if verb: print "m=%s"%m
# 
#             a.allow_multiple_edges(True)
#             # replace each label by its mapping
#             for f, d, l in a.edges():
#                 a.delete_edge(f, d, l)
#                 for l2 in m[l]:
#                     a.add_edge(f, d, l2)
# 
#             I = ['O']
#             nof = Set([K.zero()])
# 
#         a.I = I
#         a.F = nof
#         a.C = self.C
# 
#         if verb:
#             print("Before determinisation : a=%s" % a)
#         if step == 6:
#             return ("emonded automaton of non reducted general words", a)
# 
#         # rend l'automate plus simple
#         a = a.emonde0_simplify()
# 
#         if verb:
#             print("simplification : a=%s" % a)
# 
#         if verb:
#             print("Determinization...")
#         # determinize
#         ad = a.determinize2(nof=a.F)
#         # ad = a.determinize(nof=a.F, verb=False)
#         # ad = a.determinize(I, self.C, nof, verb=verb)
# 
#         if verb:
#             print(" -> %s" % ad)
#         if step == 7:
#             return ("automate des mots généraux réduits", ad)
# 
#         if ss is not None:  # not full sub-shift
#             # calculate the intersection with ss
#             ad = ad.emonde0_simplify()
#             ad = ad.intersection(ss)
#             if verb:
#                 print("after intersection : a=%s" % ad)
# 
#         if step == 8:
#             return ("automaton of reduces words", ad)
# 
#         # F2=[e for e in a.vertices() nof in e[0]]
#         # if verb: print "I2=%s"%I2 #, F2=%s"%(I2,F2)
#         ad.A = self.C
#         # ad.emondeI(I=I2) #, F=F2)
#         ad = ad.emonde0_simplify()
#         ad.F = ad.vertices()
# 
#         if verb:
#             print("after emondation : a=%s" % ad)
#         if step == 9:
#             return ("emonded automaton of reduced words", ad)
# 
#         return ad

    def critical_exponent_free(self, ss=None, prec=None, verb=False):
        r"""
        Compute the critical exponent of the beta-adic monoid
        (with or without subshift), assuming it is free.
        See http://www.latp.univ-mrs.fr/~paul.mercat/
        Publis/Semi-groupes%20fortement%20automatiques.pdf
        for a definition (without subshift).

        INPUT:

        - ``ss``- Automaton (default: ``None``)
          The first subshift to associate to the beta-adic monoid
          for this operation.

        - ``prec``- precision (default: ``None``)

        - ``verb``- bool (default: ``False``)
          If True, print informations for debugging.

        OUTPUT:

        A real number.

        EXAMPLES::

            #. Hausdorff dimension of limit set of 3-adic expansion with numerals set {0,1,3}::

                sage: m = BetaAdicMonoid(3, {0,1,3})
                sage: m.critical_exponent_free(m.reduced_words_automaton())
                log(y)/log(|3|) where y is the max root of x^2 - 3*x + 1
                0.8760357589...

            #. Hausdorff dimension of limit set of phi-adic expansion with numerals set {0,1}::

                sage: m = BetaAdicMonoid((1+sqrt(5))/2, {0,1})
                sage: m.critical_exponent_free(m.reduced_words_automaton())
                log(y)/log(|b|) where y is the max root of x^2 - x - 1
                1.0000000000...

            #. Hausdorff dimension of the boundary of the dragon fractal::

                sage: e = QQbar(1/(1+I))
                sage: m = BetaAdicMonoid(e, {0,1})
                sage: ssi = m.intersection_words(w1=[0], w2=[1])
                sage: m.critical_exponent_free(ss=ssi)
                log(y)/log(|b|) where y is the max root of x^3 - x^2 - 2
                1.5236270862...

            #. Hausdorff dimension of the boundary of a Rauzy fractal::

                sage: s = WordMorphism('1->12,2->13,3->1')
                sage: m = s.rauzy_fractal_beta_adic_monoid()
                sage: ssi = m.intersection_words(w1=[0], w2=[1])
                sage: m.critical_exponent_free(ss=ssi)
                log(y)/log(|b|) where y is the max root of x^4 - 2*x - 1
                1.0933641642...

            #. Hausdorff dimension of a non-Pisot Rauzy fractal::

                sage: s = WordMorphism({1:[3,2], 2:[3,3], 3:[4], 4:[1]})
                sage: m = s.rauzy_fractal_beta_adic_monoid()
                sage: m.b = 1/m.b
                sage: m.ss = m.ss.mirror().determinize().minimize()
                sage: m.critical_exponent_free()
                log(y)/log(|1/2*b^2 - 1/2*b + 1/2|) where y is the max root of x^3 - x^2 + x - 2
                1.5485260383...
        """

        if ss is None:
            if hasattr(self, 'ss'):
                ss = self.ss
        if ss is None:
            y = self.C.cardinality()
            print("log(%s)/log(|%s|)" % (y, self.b))
        else:
            if verb:
                print("")
            M = ss.adjacency_matrix()
            if verb:
                print("Eigen values...")
            e = M.eigenvalues()
            if verb:
                print("max...")
            y = max(e, key=abs)
            if verb:
                print("")
            print("log(y)/log(|%s|) where y is the max root of %s" % (self.b, QQbar(y).minpoly()))
            y = y.n(prec)
        from sage.functions.log import log
        b = self.b.n(prec)
        if verb:
            print("y=%s, b=%s" % (y, b))
        return abs(log(y) / log(abs(b))).n()

    def critical_exponent(self, ss=None, prec=None, verb=False):
        r"""
        Compute the critical exponent of the beta-adic monoid
        (with or without subshift).
        See http://www.latp.univ-mrs.fr/~paul.mercat/Publis/
        Semi-groupes%20fortement%20automatiques.pdf for a definition (without subshift).

        INPUT:

        - ``ss``- Automaton (default: ``None``)
          The first subshift to associate to the beta-adic monoid
          for this operation.

        - ``prec``- precision (default: ``None``)

        - ``verb``- bool (default: ``False``)
          If True, print informations for debugging.

        OUTPUT:

        A real number.

        EXAMPLES::

            #. Hausdorff dimension of limit set of 3-adic expansion with numerals set {0, 1, 3}::

                sage: m = BetaAdicMonoid(3, {0,1,3})
                sage: m.critical_exponent()
                log(y)/log(|3|) where y is the max root of x^2 - 3*x + 1
                0.8760357589...

            #. Hausdorff dimension of limit set of phi-adic expansion with numerals set {0, 1}::

                sage: m = BetaAdicMonoid((1+sqrt(5))/2, {0,1})
                sage: m.critical_exponent()
                log(y)/log(|b|) where y is the max root of x^2 - x - 1
                1.0000000000...

            #. Critical exponent that is not the Hausdorff dimension of the limit set::

                sage: P = x^7 - 2*x^6 + x^3 - 2*x^2 + 2*x - 1
                sage: b = P.roots(ring=QQbar)[3][0]
                sage: m = BetaAdicMonoid(b, {0,1})
                sage: m.critical_exponent()                    # long time
                log(y)/log(|b|) where y is the max root of x^11 - 2*x^10 - 4*x^2 + 8*x + 2
                3.3994454205...

            #. See more examples with doc critical_exponent_free()

        """
        #     #. Critical exponent that is not the Hausdorff dimension of the limit set::
        #
        #     sage:

        if ss is None:
            if hasattr(self, 'ss'):
                ss = self.ss
        if verb:
            # Calcul de l'automate des mots réduits...\n")
            print("Computation of reduce words' automata")
        a = self.reduced_words_automaton(ss=ss, verb=verb)
        return self.critical_exponent_free(prec=prec, ss=a, verb=verb)

    # test if 0 is an inner point of the limit set
    def ZeroInner(self, verb=False):

        if not hasattr(self, 'ss'):
            self.ss = self.default_ss()

        if verb:
            print("relations automaton...")

        ar = self.relations_automaton(ext=True)
        ar.complementary()

        if verb:
            print("complementary : %s" % ar)
        # return ar

        a = self.default_ss().product(self.ss)

        if verb:
            print("a = %s" % a)
        # return a

        # a = ar.intersection(a)

        if verb:
            print("product...")

        L = a.edge_labels()
        Lr = ar.edge_labels()
        d = dict([])
        for k in L:
            for kr in Lr:
                if k == kr and k[0] == 0:
                    d[(k, kr)] = 0
                else:
                    d[(k, kr)] = None
        a2 = a.product(A=ar, d=d)

        if verb:
            print("product = %s" % a2)
        # return a2

        # test if there is a cycle in the graph
        if a2.is_directed_acyclic():
            print("0 is not an inner point.")
        else:
            ss = self.ss
            self.ss = None
            print("Zero is an inner point iff the %s has non-empty interior." % self)
            self.ss = ss

    # complete the langage of a
    def complete(self, DetAutomaton a, C=None,
                 ext=False, arel=None, verb=False):
        r"""
        Return an automaton that recognize the language of all words over
        C that represent elements recognized by a.
        If ext is True, this also include words equal at infinity.

        INPUT:

        - ``a`` - A DetAutomaton.

        - ``C`` - list of digits (default : ``None``).
        """
        if C is None:
            C = self.C
        ap = DetAutomaton(None).full(list(C)).product(a)
        if ext:
            ap = ap.prefix_closure()
        if verb:
            if ap.n_states < 100:
                ap.plot()
            print("ap=%s" % ap)
        d = dict()
        Cd = [c - c2 for c in C for c2 in a.A]
        for c in Cd:
            d[c] = []
        for c in C:
            for c2 in a.A:
                d[c - c2].append((c, c2))
        if arel is None:
            arel = self.relations_automaton3(Cd=Cd, ext=ext).duplicate(d)
        else:
            arel = arel.duplicate(d)
        if verb:
            if arel.n_states < 100:
                arel.plot()
            print("arel=%s" % arel)
        ai = ap.intersection(arel)
        if ext:
            ai = ai.prune_inf()
        ai = ai.prune()
        ai = ai.minimize()
        if verb:
            if ai.n_states < 100:
                ai.plot()
            print("ai=%s" % ai)
        d = dict()
        for c in C:
            for c2 in a.A:
                d[(c, c2)] = c
        ac = ai.determinize_proj(d)
        if ext:
            ac = ac.prune_inf()
        ac = ac.prune()
        ac = ac.minimize()
        return ac

    # donne l'automate décrivant l'adhérence de l'ensemble limite avec un nouvel alphabet C
    def adherence(self, tss=None, C=None, C2=None,
                  ext=False, verb=False, step=None):
        if tss is None:
            if hasattr(self, 'tss'):
                tss = self.tss
            else:
                tss = self.default_ss()
        if C is None:
            C = list(set(self.C))
        if C2 is None:
            C2 = list(set(tss.alphabet))
        if verb:
            print("Calcul de l'automate des relations...")
        Cd = list(set([c1 - c2 for c1 in C2 for c2 in C]))
        if verb:
            print("Cd=%s" % Cd)
        a = self.relations_automaton3(Cd=Cd, ext=ext)
        if verb:
            print(" -> %s" % a)
        if step == 1:
            return a
        if ext:
            a = a.prune_inf()
        a = a.prune()
        if verb:
            print(" Après émondation : %s" % a)
        if step == 2:
            return a
        d = {}
        for c1 in C2:
            for c2 in C:
                if not d.has_key(c1 - c2):
                    d[c1-c2] = []
                d[c1-c2].append((c1, c2))
        if verb:
            print(d)
        a2 = a.duplicate(d, verb=verb)
        if verb:
            print(a2.alphabet)
            print(a2)
        if step == 3:
            return a2
        ap = tss.product(DetAutomaton(self.default_ss(C=C)), verb=verb)
        if ext:
            ap = ap.prefix_closure()
        if step == 4:
            return ap
        a2 = ap.intersection(a2)
        if step == 5:
            return a2
        if ext:
            a2 = a2.prune_inf()
        a2 = a2.prune()
        if step == 6:
            return a2
        a2 = a2.minimize()
        if step == 7:
            return a2
        if verb:
            print("determine...")
        d = {}
        for c1, c2 in a2.alphabet:
            d[(c1, c2)] = c2
        a2 = a2.determinisze_proj(d, verb=verb)
        if step == 8:
            return a2
        if verb:
            print(" -> %s" % a2)
        if ext:
            a2 = a2.prune_inf()
        a2 = a2.prune()
        if step == 9:
            return a2
        if verb:
            print("After simplification : %s" % a2)
        return a2.minimize()

    # obsolete
    # donne l'automate décrivant le translaté de +t, avec les chiffres C
    # obsolete use move2
    def move(self, t, DetAutomaton tss=None, list C=None, step=None):
        if tss is None:
            if hasattr(self, 'tss'):
                if isinstance(self.tss, DetAutomaton):
                    tss = self.tss
                else:
                    tss = DetAutomaton(self.tss)
            else:
                tss = DetAutomaton(self.default_ss())
        if C is None:
            C = list(set(self.C))

        A = tss.alphabet
        k = self.b.parent()
        nA = list(set([k(a+t2) for a in A for t2 in [0,t]]))
        a = tss.bigger_alphabet(nA)

        # add a new state
        cdef int ne, ei
        ei = a.initial_state
        ne = a.n_states  # new added state
        a.add_state(a.is_final(ei))
        a.set_initial_state(ne)  # it is the new initial state
        if step == 2:
            return a

        # add edges from the new state (copy edges
        # from the initial state and move them)
        cdef s
        for j in range(len(A)):
            a.set_succ(ne, nA.index(A[j] + t), tss.succ(ei, j))
        if step == 3:
            return a

        # compute the adherence of the new automaton
        return self.adherence(tss=a, C=C, C2=nA)

    # project a translated by t on b (do the same as move2)
    def Proj(self, a, b, t=0, arel=None):
        m2 = BetaAdicMonoid(self.b, set(a.alphabet+b.alphabet))
        if arel is None:
            # compute the relations automaton with translation t
            arel = m2.relations_automaton4(t=t, couples=True,
                                           A=a.alphabet, B=b.alphabet)
        ai = arel.intersection(a.zero_complete2().product(b))
        d = {}
        for i in a.alphabet:
            for j in b.alphabet:
                d[(i, j)] = j
        r = ai.proj(d)
        r.zero_completeOP()
        return r


    # donne l'automate décrivant le translaté de +t de a,
    # avec les chiffres A au départ et B à l'arrivée, le tout dans l'ensemble décrit par b
    def move2(self, t, DetAutomaton a=None, DetAutomaton b=None,
              list A=None, list B=None, ar=None, verb=False):
        if a is None:
            if hasattr(self, 'tss'):
                if isinstance(self.tss, DetAutomaton):
                    a = self.tss
                else:
                    a = DetAutomaton(self.tss)
            else:
                a = DetAutomaton(self.default_ss())
        if b is None:
            b = DetAutomaton(self.default_ss())
            #            if hasattr(self, 'tss'):
            #                if isinstance(self.tss, DetAutomaton):
            #                    b = self.tss
            #                else:
            #                    b = DetAutomaton(self.tss)
            #            else:
            #                b = DetAutomaton(self.default_ss())
        if A is None:
            A = list(set(a.A))
        if B is None:
            B = list(set(b.A))
        if ar is None:
            # compute the relations automaton with translation t
            ar = self.relations_automaton4(t=t, A=A, B=B,
                                           couples=True, verb=False)
        # compute the product of a and b
        if verb:
            print("product...")
        ap = a.zero_complete2().product(b.zero_complete2())
        if verb:
            print("ap = %s" % ap)
        # compute the intersections
        if verb:
            print("intersection...")
        ai = ar.intersection(ap)
        if verb:
            print("ai = %s" % ai)
        if verb:
            print("min...")
        ai = ai.minimize()
        if verb:
            print("ai = %s" % ai)
        # project on one side
        d = {}
        for c1 in A:
            for c2 in B:
                d[(c1, c2)] = c2
        if verb:
            print("d=%s" % d)
        if verb:
            print("determinize...")
        ai = ai.determinize_proj(d, verb=verb)
        if verb:
            print("ai=%s" % ai)
        if verb:
            print("min")
        ai.zero_completeOP()
        return ai.prune().minimize()

    # return the automaton recognizing the division by beta and translation t
    def shift(self, DetAutomaton aa,
              DetAutomaton bb, t=0, verb=False):
        # détermine la liste des successeurs de l'état initial
        cdef int i
        cdef int e
        l = set()
        for j in range(aa.a.na):
            e = aa.a.e[aa.a.i].f[j]
            if e != -1:
                l.add(j)
        l = list(l)
        if verb:
            print("états à considérer : %s" % l)
        # calcule l'union des automates translatés
        a = DetAutomaton(None)
        A = aa.alphabet
        a.setAlphabet(A)
        for i in range(0, len(l)):
            a2 = aa.copy()
            a2.set_initial_state(aa.a.e[aa.a.i].f[l[i]])
            a2 = a2.prune().minimize()
            a2 = self.move2(t=-(A[l[i]]-t)/self.b, a=a2, b=bb)
            a = a.union(a2)
        return a

    # calcule l'intersection des ensembles limites
    def intersection3(self, DetAutomaton a, DetAutomaton b, ext=True):
        a2 = self.complete(a, ext=ext)
        b2 = self.complete(b, ext=ext)
        return a2.intersection(b2).prune()

    # determine if the limit sets intersect
    def intersect(self, DetAutomaton a, DetAutomaton b, ext=True, verb=False):
        a2 = self.complete(a, ext=ext)
        if verb:
            print("a2=%s" % a2)
        b2 = self.complete(b, ext=ext)
        if verb:
            print("b2=%s" % b2)
        return not a2.intersection(b2).is_empty(ext)

    # Dit si toutes les tuiles du pavages autosimilaire sont connexes ou pas
    def is_all_connected(self, DetAutomaton a=None, ext=True, verb=False):
        if a is None:
            if hasattr(self, 'tss'):
                if isinstance(self.tss, DetAutomaton):
                    a = m.tss
                else:
                    a = DetAutomaton(self.tss)
            else:
                a = DetAutomaton(None).full(list(self.C))

        from sage.graphs.graph import Graph
        n = a.n_states
        na = len(a.alphabet)
        d = dict([])  # dictionnaire des automates complétés
        if verb:
            print("Automaton of relations...")
        arel = self.relations_automaton3(ext=ext)
        for i in range(n):
            if verb:
                print("piece %s" % i)
            g = Graph({j: {} for j in range(na) if a.succ(i, j) != -1})
            # compute the neighboorhood graph of the piece i
            for j in g.vertices():
                for k in g.vertices():
                    if k >= j:
                        continue
                    if verb:
                        print(" intersection %s et %s..." % (j, k))
                    la = []
                    la.append(a.piece(j, e=i).minimize())
                    la.append(a.piece(k, e=i).minimize())
                    for l in range(2):
                        a2 = d.get(la[l])
                        # récupère le complété dans le dictionnaire s'il y est
                        if a2 is None:
                            if verb:
                                print("  complete %s..." % la[l])
                            a2 = self.complete(la[l], ext=ext)
                            d[la[l]] = a2
                            la[l] = a2
                        else:
                            la[l] = a2
                            if verb:
                                print("  already calculated !")
                    if verb:
                        print("  intersect...")
                    # if self.intersect(a.piece(j, e=i), a.piece(k, e=i), ext=ext, verb=verb):
                    if la[0].intersect(la[1], ext=ext):
                        if verb:
                            print("  yes !")
                        g.add_edge(j, k)
            if not g.is_connected():
                return False
        return True

    # Dit si l'ensemble limite est connexe ou pas
    def is_connected(self, DetAutomaton a=None):
        if a is None:
            if hasattr(self, 'tss'):
                a = DetAutomaton(self.tss)
            else:
                a = DetAutomaton(None).full(self.C)

        n = a.n_states
        na = len(a.alphabet)
        rules = [[[l] for l in a.alphabet] for i in range(n)]
        gvois = [Graph(na) for i in range(n)]  # graphe des morceaux voisins
        gnvois = [Graph(na) for i in range(n)]  # graphe des morceaux non voisins

        # liste des morceaux dont il faut tester
        # la connexité du graphe de voisinage
        m = [a.initial_state]

        while len(m) > 0:
            i = m.pop()
            gvois[i].connected_component_containing_vertex(i)

        raise ValueError("Not implemented !")

    # Dit si l'ensemble limite est simplement connexe ou pas
    def is_simply_connected(self, DetAutomaton a=None):
        if a is None:
            if hasattr(self, 'tss'):
                a = DetAutomaton(self.tss)
            else:
                a = DetAutomaton(None).full(self.C)
        # TODO!
        else:
            raise ValueError("Not implemented !")

    # used by Approx
    def Approx_rec(self, DetAutomaton a, test, f, x, int n, int n2):
        """
        EXAMPLE::

            sage: m = BetaAdicMonoid((x^3-x^2-x-1).roots(ring=QQbar)[1][0], {0,1})
            sage: pm = m.b.parent().places()[1]
            sage: a = m.Approx(13, lambda x: (pm(x).real())^2 + (pm(x).imag())^2 < .4 )
            sage: print(a)
            DetAutomaton with 3538 states and an alphabet of 2 letters
            sage: aoc = m.zero_complete(a) #get the canonical representation of the g-b-set
            sage: aoc
            DetAutomaton with 182 states and an alphabet of 2 letters
        """
        if n == 0:
            if test(x):
                return f
            else:
                return -1
        else:
            e = dict()
            add = False
            for t in self.C:
                e[t] = self.Approx_rec(a, test, f, x+t*self.b**(n2-n), n-1, n2)
                if e[t] != -1:
                    add = True
            if add:
                e3 = a.add_state(0)
                for t in self.C:
                    if e[t] != -1:
                        a.add_edge(e3, t, e[t])
                return e3
            return -1

    # gives a automaton describing a approximation of a set defined by
    # the characteritic function test
    # rk : can be improve using a reduced words automaton
    def Approx(self, n, test):  # , ared=None):
        #        if ared is None:
        #            ared = m.reduced_words_automaton2()
        a = DetAutomaton(None)
        a.setAlphabet(list(self.C))
        f = a.add_state(1)
        e = self.Approx_rec(a, test, f, 0, n, n)
        for t in self.C:
            a.add_edge(f, t, f)
        a.set_initial_state(e)
        return a

    def zero_complete(self, a, verb=False, test=False):
        if verb:
            print("a = %s" % a)
        a2 = a.zero_complete2()
        a2.zero_completeOP()
        if verb:
            print("a2 = %s (après zero-complétion)" % a2)
        aoc = self.move2(0, a2)
        if verb:
            print("aoc = %s" % aoc)
        aoc.zero_completeOP()
        aoc = aoc.zero_complete2()
        return aoc

    # calcule la liste triée (par rapport à la place >1) des
    # premiers points dans omega-omega
    #
    # THIS FUNCTION IS UNCORRECT AND VERY INEFFICIENT ! SHOULD BE IMPROVED.
    #
    def compute_translations(self, DetAutomaton aoc, imax=None, verb=False):
        b = self.b
        p = b.parent().places()
        if verb:
            print(p)
        if abs(p[0](b)) > 1:
            pp = p[0]
            pm = p[1]
        else:
            pp = p[1]
            pm = p[0]
        l = []
        if imax is None:
            imax = 20
        bound = max([abs(pm(x)) for x in self.C])/(1-abs(pm(b)))
        if verb:
            print("bound = %s" % bound)
        if b.minpoly().degree() == 3:
            for i in range(-imax, imax):
                for j in range(-imax, imax):
                    for k in range(-imax, imax):
                        x = i + b*j + b*b*k
                        if abs(pm(x)) <= bound and pp(x) > 0:
                        # if abs(pm(x)-z) < r + .5 and pp(x) > 0:
                            l.append(x)
        elif b.minpoly().degree() == 2:
            for i in range(-imax, imax):
                for j in range(-imax, imax):
                    x = i + b*j
                    if abs(pm(x)) <= bound and pp(x) > 0:
                    # if abs(pm(x)-z) < r + .5 and pp(x) > 0:
                        l.append(x)
        else:
            raise NotImplemented
        l.sort(key=pp)
        return l

    # décrit les mots de a de longueur n partant de e (utilisé par compute_translation2)
    def Parcours(self, A, a, e, t, n, bn):
        # print "Parcours e=%s t=%s n=%s bn=%s"%(e,t,n,bn)
        if n == 0:
            if a.is_final(e):
                return [t]
            else:
                return []
        else:
            l = []
            for i in range(len(A)):
                f = a.succ(e, i)
                if f != -1:
                    l += self.Parcours(A, a, f, t + bn * A[i], n - 1, bn*self.b)
            if a.is_final(e):
                l.append(t)
            return l

    def compute_translations2(self, DetAutomaton aoc, imax=None, verb=False):
        b = self.b
        A = self.C
        p = b.parent().places()
        if verb:
            print(p)
        if abs(p[0](b)) > 1:
            pp = p[0]
            pm = p[1]
        else:
            pp = p[1]
            pm = p[0]
        if imax is None:
            imax = 8
        # compute a reduced version of aoc
        # ared = self.reduced_words_automaton2()
        # aoc = aoc.intersection(ared)
        # compute aoc-aoc
        d = dict()
        for t1 in self.C:
            for t2 in self.C:
                d[(t1, t2)] = t1 - t2
        ad = aoc.product(aoc)
        ad = ad.determinize_proj(d)
        if verb:
            print("ad = %s" % ad)
        # compute the reduced words automaton with the difference alphabet
        # m2 = BetaAdicMonoid(b, set([t1 - t2 for t1 in A for t2 in A]))
        # if verb: print "m2 = %s"%m2
        # ar2 = m2.reduced_words_automaton2()
        # if verb: print "ar2 = %s"%ar2
        # intersect
        # adr = ad.intersection(ar2)
        # if verb: print "adr = %s"%adr
        adr = ad
        # compute the list of points
        if verb:
            print("Ways %s..." % imax)
        l = self.Parcours(adr.alphabet, adr, adr.initial_state, 0, imax, 1)
        if verb:
            print("%s computed points" % len(l))
        # sort
        if verb:
            print("sort...")
        l = list(set(l)) #avoid repetitions
        l.sort(key=pp)
        return l[l.index(0)+1:]

    # verify that a piece exchange is correct
    # w e assume that pieces are disjoints
    # WARNING : IF DOES NOT CONTAIN 0 THE PIECE EXCHANGE
    # IS NOT NECESSARLY INVERSIBLE
    def correct_morceaux(self, DetAutomaton aoc, list lm, bool verb=False):
        u = DetAutomaton(None)
        u.setAlphabet(aoc.alphabet)
        # verify that the union of the pieces is equal to aoc
        for a, t in lm:
            u = u.union(a)
        if not u.equals_langages(aoc):
            if verb: 
                print("Les morceaux ne recouvrent pas !")
            return False
        # verify that the union of the translated pieces is equal to aoc
        u = DetAutomaton(None)
        u.setAlphabet(aoc.alphabet)
        u.add_state(True)
        u.set_initial_state(0)
        # reconnait 0 (le seul élément à ne pas avoir forcément d'inverse)
        u.add_edge(0, 0, 0)  
        for a, t in lm:
            u = u.union(self.move2(-t, a))
        if not u.equals_langages(aoc):
            if verb:
                print("Translated pieces do not overlap !")
            return False
        return True

    def compute_morceaux(self, DetAutomaton aoc, lt=None, method=1, imax=None,
                         verb=False, stop=-1):
        r"""
        Compute the domain exchange describing the
        g-beta-expansion given by the automaton aoc.

        INPUT:

        - ``aoc``- DetAutomaton
            Automaton of the g-beta-expansion.

        - ``verb``- bool (default: ``False``)
          If True, print informations about the computing.

        OUTPUT:

        A list of (DetAutomaton, translation).

        EXAMPLES::

            #. Full Tribonnacci::

                sage: m = BetaAdicMonoid((x^3-x^2-x-1).roots(ring=QQbar)[1][0], {0, 1})
                sage: pm = m.b.parent().places()[1]
                sage: a = m.Approx(13, lambda x: (pm(x).real())^2 + (pm(x).imag())^2 < .4 )
                sage: aoc = m.zero_complete(a)
                sage: m.compute_morceaux(aoc) # long time
                [(DetAutomaton with 85 states and an alphabet of 2 letters, 1),
                 (DetAutomaton with 71 states and an alphabet of 2 letters, b^2 - b),
                 (DetAutomaton with 173 states and an alphabet of 2 letters, b^2),
                 (DetAutomaton with 69 states and an alphabet of 2 letters, b^2 + 1),
                 (DetAutomaton with 59 states and an alphabet of 2 letters, b^2 + b +
                1),
                 (DetAutomaton with 38 states and an alphabet of 2 letters, b^2 + b),
                 (DetAutomaton with 119 states and an alphabet of 2 letters, b),
                 (DetAutomaton with 80 states and an alphabet of 2 letters, b + 1)]
        """
        m = self
        A = aoc.A
        if lt is None:
            if verb:
                print("Compute the list of translations...")
            if method == 1:
                lt = self.compute_translations(aoc, imax, verb)
            else:
                lt = self.compute_translations2(aoc, imax, verb)
        if verb:
            print("list of %s points." % len(lt))
        if verb:
            print("Compute the pieces...")
        u = DetAutomaton(None)
        u.setAlphabet(list(A))
        uc = u.complementary()
        at = dict()
        for i, t in enumerate(lt):
            if i == stop:
                break
            if verb:
                print("t=%s" % t)
            at[t] = m.move2(t=t, a=aoc, b=aoc)  # , verb=verb)
            # if verb: print "intersection..."
            at[t].zero_completeOP()
            # if verb: print "intersection..."
            at[t] = at[t].intersection(uc)
            at[t].zero_completeOP()
            # if verb: print "union..."
            u = u.union(at[t]).prune().minimize()
            # if verb: print "prune..."
            at[t] = at[t].prune().minimize()
            # teste si c'est fini
            # if verb: print "compl..."
            uc = u.complementary()
            ai = uc.intersection(aoc)
            # if verb: print "empty..."
            if ai.is_empty():
                break  # c'est fini !
            if at[t].is_empty():
                if verb:
                    print("Empty Automaton !")
                del at[t]
            else:
                if verb:
                    print(at[t])
        if stop == -1 and not ai.is_empty():
            raise ValueError("Incorrect list of translations for" +
                             "the computation of piece exchange." +
                             "Try to increase imax.")

    # décrit les éléments de a de longueur n (utilisé par compute_morceaux2)
    #
    #  NE FONCTIONNE PAS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #
    def ParcoursUnic(self, a, pp, n, verb=False):
        r = dict()
        A = a.alphabet
        b = self.b
        p = [(a.initial_state, 0, 1)]  # pile des éléments à traiter
        r[0] = None
        res = []
        bnm = b**n
        while len(p) != 0:
            e, t, bn = p.pop()
            r[t] = None
            if verb:
                print("%s, %s, %s" % (e, t, bn))
            # print("p = %s"%p)
            for l in A:
                f = a.succ(e, l)
                if f != -1:
                    t2 = t+bn*l
                    # print "l=%s, t2 = %s, pp(bn)=%s, pp(bnm)=%s"%(l, t2, pp(bn), pp(bnm))
                    if pp(bn) <= pp(bnm):
                        p.append((f, t2, bn*b))
                        if a.is_final(f):
                            if not r.has_key(t2) and pp(t2) > 0:
                                res.append(t2)
        return res

    def compute_morceaux2(self, DetAutomaton aoc, iplus=1,
                          reduit=False, verb=False, step=None):
        r"""
        Compute the domain exchange describing the g-beta-expansion given by the automaton aoc.

        INPUT:

        - ``aoc``- DetAutomaton
            Automaton of the g-beta-expansion.

        - ``verb``- bool (default: ``False``)
          If True, print informations about the computing.

        OUTPUT:

        A list of (DetAutomaton, translation).

        EXAMPLES::

            #. Full Tribonnacci::

                sage: m = BetaAdicMonoid((x^3-x^2-x-1).roots(ring=QQbar)[1][0], {0, 1})
                sage: pm = m.b.parent().places()[1]
                sage: a = m.Approx(13, lambda x: (pm(x).real())^2 + (pm(x).imag())^2 < .4 )
                sage: aoc = m.zero_complete(a)
                sage: m.compute_morceaux2(aoc) # long time
                [(DetAutomaton with 85 states and an alphabet of 2 letters, 1),
                 (DetAutomaton with 71 states and an alphabet of 2 letters, b^2 - b),
                 (DetAutomaton with 173 states and an alphabet of 2 letters, b^2),
                 (DetAutomaton with 69 states and an alphabet of 2 letters, b^2 + 1),
                 (DetAutomaton with 59 states and an alphabet of 2 letters, b^2 + b +
                1),
                 (DetAutomaton with 38 states and an alphabet of 2 letters, b^2 + b),
                 (DetAutomaton with 119 states and an alphabet of 2 letters, b),
                 (DetAutomaton with 80 states and an alphabet of 2 letters, b + 1)]       
        """
        m = self
        A = aoc.A
        b = self.b
        p = b.parent().places()
        d = dict()
        for t1 in A:
            for t2 in A:
                d[(t1, t2)] = t2 - t1
        if verb:
            print(p)
        if abs(p[0](b)) > 1:
            pp = p[0]
            # pm = p[1]
        else:
            pp = p[1]
            # pm = p[0]
        if reduit:
            if verb:
                print("Compute ared...")
            ared = self.reduced_words_automaton2()
            if verb:
                print("Compute aocr...")
            aocr = aoc.intersection(ared)
        if verb:
            print("Compute the pieces...")
        u = DetAutomaton(None)
        u.setAlphabet(list(A))
        uc = u.complementary().intersection(aoc)
        at = dict()
        while True:
            #############
            # compute uc-aoc
            if verb:
                print("product...")
            if reduit:
                ucr = uc.intersection(ared)
                ad = ucr.product(aocr)
            else:
                ad = uc.product(aoc)
            if verb:
                print("determinize...")
            ad = ad.determinize_proj(d)
            if verb:
                print("minimise...")
            ad = ad.prune().minimize()
            ad.zero_completeOP()
            if verb:
                print("ad = %s" % ad)
            # compute the list of points
            imax = 3
            fin = False
            while True:
                if imax >= 1000:
                    return "imax too big!", ad, uc, [(at[t], t) for t in at.keys()]
                if verb:
                    print("Way %s..." % imax)
                l = self.Parcours(ad.alphabet, ad, ad.initial_state, 0, imax, 1)
                if len(l) == 0:
                    imax += 1
                    continue
                if verb:
                    print("%s computed points" % len(l))
                # sort
                if verb:
                    print("sort...")
                l = list(set(l)) #avoid repetitions
                l.sort(key=pp)
                if 0 in l:
                    l = l[l.index(0) + 1:]
                else:
                    return "0 non présent !", ad
                    raise ValueError("Erreur : 0 non présent dans la différence !!!")
                if fin:
                    break
                if l != []:
                    # if iplus == 0:
                    #    break
                    fin = True
                    imax += iplus
                imax += 1
            t = l[0]
            if verb:
                print(l[:10])
            #############
            if verb:
                print("t=%s" % t)
            if at.has_key(t):
                return "tr deja vu", [(at[t], t) for t in at.keys()], uc
                raise ValueError("Error : computed translation already seen!!!")
            at[t] = m.move2(t=t, a=aoc)  # , verb=verb)
            # if verb: print "intersection..."
            at[t].zero_completeOP()
            # if verb: print "intersection..."
            if not uc.included(aoc):
                raise ValueError("Error : uc is not included in aoc !!!")
            at[t] = at[t].intersection(uc)
            at[t].zero_completeOP()
            # if verb: print "prune..."
            at[t] = at[t].prune().minimize()
            # if verb: print "union..."
            u = u.union(at[t]).prune().minimize()
            at[t] = at[t].prune()
            # teste si c'est fini
            # if verb: print "compl..."
            uc = u.complementary()
            uc = uc.intersection(aoc)
            # if verb: print "empty..."
            if at[t].is_empty():
                return "tr incorrecte", [(at[t], t) for t in at.keys()], uc
                raise ValueError("Error : incorrect computed translation!!!")
            else:
                if verb:
                    print("at[%s]=%s, uc=%s" % (t, at[t], uc))
            if uc.is_empty():
                break  # c'est fini !
            if step is not None:
                if step == 0:
                    return "Arret", [(at[t], t) for t in at.keys()], uc
        if not uc.is_empty():
            raise ValueError("Error: piece exchange doesn't pave!!!")
        return [(at[t], t) for t in at.keys()]

    def compute_morceaux3(self, DetAutomaton a, DetAutomaton ap=None,
                          bound=100, iplus=1, getad=False, verb=False,
                          need_included=True, step=None):
        r"""
        Compute the domain exchange describing the g-beta-expansion given by the automaton aoc.

        INPUT:

        - ``a``- DetAutomaton
          Automaton of the g-beta-expansion.
        -  ``ap``- DetAutomaton (default: ``None``)
           Langage used for the computations. Everything is projected on it.
        - ``verb``- bool (default: ``False``)
          If True, print informations about the computing.

        OUTPUT:

        A list of (DetAutomaton, translation).

        EXAMPLES::

            # Full Tribonnacci::

                sage: m = BetaAdicMonoid((x^3-x^2-x-1).roots(ring=QQbar)[1][0], {0,1})
                sage: pm = m.b.parent().places()[1]
                sage: a = m.Approx(13, lambda x: (pm(x).real())^2 + (pm(x).imag())^2 < .4 )
                sage: aoc = m.zero_complete(a)
                sage: m.compute_morceaux3(aoc) # long time
                [(DetAutomaton with 86 states and an alphabet of 2 letters, 1),
                 (DetAutomaton with 164 states and an alphabet of 2 letters, b^2),
                 (DetAutomaton with 70 states and an alphabet of 2 letters, b^2 + 1),
                 (DetAutomaton with 60 states and an alphabet of 2 letters, b^2 + b +
                1),
                 (DetAutomaton with 39 states and an alphabet of 2 letters, b^2 + b),
                 (DetAutomaton with 120 states and an alphabet of 2 letters, b),
                 (DetAutomaton with 81 states and an alphabet of 2 letters, b + 1)]
        """
        m = self
        if ap is None:
            if hasattr(self, 'tss'):
                ap = DetAutomaton(self.tss)
            else:
                ap = DetAutomaton(None).full(list(self.C))
        cdef DetAutomaton aa
        if not a.included(ap):
            if verb:
                print("Project a on ap...")
            aa = a.copy()
            # check that Qap contains Qaa
            if not m.Proj(ap, aa).equals_langages(aa) and need_included:
                raise ValueError("The g-beta-set described by a is not included in the one described by ap.")
            # project aa on ap
            aa = m.Proj(aa, ap)
        else:
            aa = a
        A = ap.A
        if verb:
            print("A=%s, aa=%s, ap=%s" % (A, aa, ap))
        b = self.b
        p = b.parent().places()
        d = dict()
        for t1 in A:
            for t2 in A:
                d[(t1, t2)] = t2 - t1
        if abs(p[0](b)) > 1:
            pp = p[0]
            # pm = p[1]
        else:
            pp = p[1]
            # pm = p[0]
        if verb:
            print("pp=%s" % pp)
        # check that the alphabet of ap is non-negative
        for i in ap.alphabet:
            if pp(i) < 0:
                raise ValueError("The alphabet of ap must be non-negative in the expanding direction !")
        if verb:
            print("Compute the pieces...")
        u = DetAutomaton(None)
        u.setAlphabet(list(A))
        uc = u.complementary().intersection(aa)
        at = dict()
        while True:
            #############
            # compute uc-aa
            if verb:
                print("product...")
            ad = uc.product(aa)
            if verb:
                print("determinize...")
            ad = ad.determinize_proj(d)
            if verb:
                print("minimise...")
            ad = ad.prune().minimize()
            if verb:
                print("ad = %s" % ad)
            # project on ap
            if verb:
                print("proj on ap...")
            ad = m.Proj(ad, ap)
            if verb:
                print("ad = %s" % ad)
            if getad:
                return ad
            #############
            # compute the list of points
            R = set()
            while len(R) == 0:
                if verb:
                    print("bound=%s" % bound)
                S = set()
                R = set()
                TS = [(ad.initial_state, 0, 1)]
                while len(TS) > 0:
                    TS2 = []
                    # print TS
                    for s, t, bn in TS:
                        for i in ad.alphabet:
                            ss = ad.succ(s, i)
                            if ss != -1:
                                tt = t + i*bn
                                if pp(tt) < bound:
                                    if (ss, tt) not in S:
                                        S.add((ss, tt))
                                        TS2.append((ss, tt, bn * m.b))
                                        if ad.is_final(ss) and pp(tt) > 0:
                                            R.add(tt)
                    TS = TS2
                bound = bound * 2
            l = sorted(R, key=pp)
            if verb:
                print(l)
            #############
            # compute the pieces for these translations
            for t in l:
                if verb:
                    print("t=%s" % t)
                if at.has_key(t):
                    return "tr deja vu",at, uc
                    raise ValueError("Erreur : translation calculée déjà vue !!!")
                # pourrait valoir le coup de projeter sur aa pour l'éfficacité des calculs ?
                at[t] = m.Proj(aa, ap, t=t)
                # if verb: print "intersection..."
                at[t] = at[t].intersection(uc)
                if at[t].is_empty():
                    at.pop(t)
                    continue
                # if verb: print "union..."
                u = u.union(at[t])
                # if verb: print "compl..."
                uc = u.complementary()
                uc = uc.intersection(aa)
                if uc.is_empty():
                    break  # c'est fini !
                else:
                    if verb:
                        print("at[%s]=%s, uc=%s" % (t, at[t], uc))
            if uc.is_empty():
                break
        return [(at[t], t) for t in at.keys()]

    def compute_substitution(self, DetAutomaton a=None,
                             np=None, lt=None, method=2,
                             method_tr=1, iplus=2, imax=None,
                             get_aut=False, verb=True):
        r"""Compute a substitution whose fixed point is the g-beta-expansion given
        by the beta-adic monoid with automaton a.

        INPUT:

        - ``a``- DetAutomaton (default: ``None``)
            Automaton of the g-beta-expansion.

        - ``np``- int (default: ``None``)
            Power of beta for the computing. The g-beta-expansion must be b^np invariant.

        - ``lt``- list of elements of the integer rings
            List of translations to compute the pieces exchange.

        - ``get_aut``- Bool (default: ``False``)
            If True, gives also the list of automata.

        - ``verb``- bool (default: ``True``)
          If True, print informations about the computing.

        OUTPUT:

        A word morphism given by a dictionnary.

        EXAMPLES::

            #. Full Tribonnacci::

                sage: m = BetaAdicMonoid((x^3-x^2-x-1).roots(ring=QQbar)[1][0], {0, 1})
                sage: m.compute_substitution(verb=False)      # long time
                {1: [1, 3], 2: [1], 3: [1, 2]}

        """
        m = self
        if a is None:
            if hasattr(self, 'tss'):
                a = DetAutomaton(self.tss)
            else:
                a = DetAutomaton(None).full(list(self.C))
        a.zero_completeOP()
        A = a.A
        # complete a
        aoc = m.move2(t=0, a=a)
        aoc.zero_completeOP()
        aoc = aoc.prune().minimize()
        if verb:
            print("aoc = %s" % aoc)
        # test if np is big enough
        if np is None:
            for i in range(1, 300):
                baoc = aoc.unshift(0, i)
                if baoc.included(aoc):
                    np = i
                    break
            if np is None:
                raise ValueError('The g-beta-expansion must be b^np invariant for some natural integer np.')
        else:
            baoc = aoc.unshift(0, np)
            if not baoc.included(aoc):
                raise ValueError('The g-beta-expansion must be b^np invariant (here np=%s).' % np)
        if verb:
            print("np = %s" % np)
        # compute the pieces exchange
        if lt is None:
            if verb:
                print("compute the pieces exchange...")
            if method == 1:
                lt = self.compute_morceaux(aoc, method=method_tr,
                                           imax=imax, verb=verb)
            else:
                lt = self.compute_morceaux2(aoc, iplus, verb=verb)
        if verb:
            print("Exchange of %s pieces" % len(lt))
        # calcule l'induction à partir de la liste de (morceau, translation)
        # précalculs
        if verb:
            print("Pre- computation...")
        arel = dict()
        for a, t in lt:
            arel[t] = m.relations_automaton4(t=-t, A=a.A, B=m.C, couples=True)
            if verb:
                print("arel[%s]=%s" % (t, arel[t]))
        baoc = aoc.unshift(0, np)  # multiplie par b
        baoc = m.move2(t=0, a=baoc, b=aoc)  # complete
        baoc.zero_completeOP()
        if verb:
            print("baoc : %s" % baoc)
        # arbre de subdivision des morceaux
        arbre = [range(1, len(lt) + 1)] + [[] for i in range(len(lt))]
        if verb:
            print("initial tree: %s" % arbre)
        lm = [(aoc, 0)] + lt  # liste des morceaux, translations
        if verb:
            print("lm = %s" % lm)
        # parcours de chaque morceau (donné par une liste de morceaux)
        d = [[] for i in range(len(lm))]
        if verb:
            print("d = %s" % d)
        lf = range(1, len(lm))  # liste des feuilles
        if verb:
            print("lf = %s" % lf)

        from copy import copy

        if verb:
            print("\n********************\n   Step 1   \n********************")

        # étape 1 : complétion des mots
        for i, (a1, t1) in enumerate(lm):
            if arbre[i] != []:
                continue  # ce morceau n'est pas une feuille
            if verb:
                print("\nComputation of piece %s/%s (%s, %s)..." % (i, len(lm), a1, t1))
                # print "lf = %s"%lf
                # print "d = %s"%d
                # print "arbre = %s"%arbre
            tr = 0  # translation totale
            if d[i] != []:
                if d[i][-1] == -1:
                    continue  # le morceau était déjà fini de calculé
                # va à la fin du mot
                for j in d[i]:
                    if j < 0:
                        break
                    tr += lm[j][1]
                # calcule b^np*a + tr
            a = a1.unshift(0, np).prune().minimize()
            if tr != 0:
                if verb:
                    print("Translation of %s..." % tr)
                # TODO : ne pas recalculer cet automate déjà calculé
                a = m.move2(t=-tr, a=a)
            while True:
                # split selon les autres morceaux
                j = included(a, lf, lm)
                if j is None:
                    # détermine les morceaux qui intersectent a
                    l = []
                    for j in lf:
                        if lm[j][0].intersect(a):
                            l.append(j)
                    if len(l) < 2:
                        print("Error : intersection with %s pieces not included !!!" % len(l))
                    if verb:
                        print("Subdivision on %s pieces..." % len(l))
                    # calcule les intersections (découpe en morceaux a1)
                    for j in l:
                        a2 = lm[j][0]
                        # découpe a selon a2
                        a = m.move2(t=tr, a=a2)  # translate a2 de -tr
                        a.zero_completeOP()
                        a.shiftOP(0, np)  # multiplie par b^(-np)
                        a = a.prune().minimize()
                        k = len(lm)  # indice du nouveau morceau
                        lf.append(k)  # nouvelle feuille
                        arbre[i].append(k)
                        arbre.append([])
                        from copy import copy
                        # print copy
                        d.append(copy(d[i]))
                        d[k].append(j)  # ajoute la translation suivante
                        # ajoute le nouveau morceau à la liste
                        lm.append((a.intersection(a1), t1))
                        # split selon baoc
                        (ab, abc) = split_baoc(k, tr + lm[j][1], np,
                                               lm, m, aoc, verb)
                        if ab is None:
                            if verb:
                                print("k=%s, tr=%s+%s : computation to continu" % (k, tr, lm[j][1]))
                        else:
                            if abc is None:
                                if verb:
                                    print("tr=%s : end of computation" % tr)
                                # indique que le calcul de ce morceau est terminé
                                d[k].append(-1)
                            else:
                                if verb:
                                    print("tr=%s : subdivision of %s by baoc (new %s)..." % (tr, i, len(lm)))
                                lf.append(len(lm))  # nouvelle feuille
                                arbre[k].append(len(lm))
                                arbre.append([])
                                d.append(copy(d[k]))
                                # indique que le calcul est terminé pour ce morceau (pour l'étape 1)
                                d[len(lm)].append(-1)
                                lm.append((ab, t1))
                                lf.append(len(lm))  # nouvelle feuille
                                arbre[k].append(len(lm))
                                arbre.append([])
                                d.append(copy(d[k]))
                                lm.append((abc, t1))
                                # le morceau k n'est plus une feuille
                                lf.remove(k)
                    lf.remove(i)  # le morceau i n'est plus une feuille
                    # calcul fini pour ce morceau puisque ce n'est plus une feuille
                    break
                else:
                    # ajoute le morceau à la liste et translate
                    d[i].append(j)
                    # if verb: print "Translation de %s..."%lm[j][1]
                    a = m.move2(t=-lm[j][1], a=a, ar=arel[lm[j][1]])
                    tr += lm[j][1]
                # split selon baoc
                (ab, abc) = split_baoc(i, tr, np, lm, m, aoc, verb)
                if ab is None:
                    # if verb: print "tr=%s : calcul à continuer"%tr
                    pass
                else:
                    if abc is None:
                        if verb:
                            print("tr=%s : end of computation" % tr)
                    else:
                        if verb:
                            print("tr=%s : subdivision of %s by baoc (new %s)..." % (tr, i, len(lm)))
                        lf.append(len(lm))  # nouvelle feuille
                        arbre[i].append(len(lm))
                        arbre.append([])
                        d.append(copy(d[i]))
                        # indique que le calcul est terminé pour ce morceau (pour l'étape 1)
                        d[len(lm)].append(-1)
                        lm.append((ab, t1))
                        lf.append(len(lm))  # nouvelle feuille
                        arbre[i].append(len(lm))
                        arbre.append([])
                        d.append(copy(d[i]))
                        lm.append((abc, t1))
                        lf.remove(i)  # le morceau i n'est plus une feuille
                    break  # calcul fini pour ce morceau (pour cette étape)
        if verb:
            print("\n*************\n   Step 2   \n*************")

        # étape 2 : remplacement des lettres qui ne sont pas des feuilles
        while True:
            end = True
            for i in lf:
                a1, t1 = lm[i]
                if verb:
                    print("\nPiece %s/%s..." % (i, len(lm)))
                tr = 0  # translation totale
                if d[i] == []:
                    print("Error : empty sheet !!!!")
                # va à la fin du mot
                for ij, j in enumerate(d[i]):
                    if j < 0:
                        break
                    if arbre[j] != []:  # il faut recalculer cette lettre
                        # calcule b^np*a + tr
                        a = a1.unshift(0, np).prune().minimize()
                        if tr != 0:
                            if verb:
                                print("Translation of %s..." % tr)
                            # TODO : ne pas recalculer cet automate déjà calculé
                            a = m.move2(t=-tr, a=a)
                        # split selon les autres morceaux
                        f = fils(arbre, j)
                        if verb:
                            print("Split selon %s morceaux" % len(f))
                        k = included(a, f, lm)
                        if k is None:
                            end = False
                            # détermine les morceaux qui intersectent a
                            l = []
                            for k in lf:
                                if lm[k][0].intersect(a):
                                    l.append(k)
                            if len(l) < 2:
                                print("Error : intersection with %s pieces but not included !!!" % len(l))
                            if verb:
                                print("Subdivision en %s morceaux..." % len(l))
                            # calcule les intersections (découpe en morceaux a1)
                            for j2 in l:
                                a2 = lm[j2][0]
                                # découpe selon a2
                                a = m.move2(t=tr, a=a2)  # translate a2 de -tr
                                a.zero_completeOP()
                                a.shiftOP(0, np)  # multiplie par b^(-np)
                                a = a.prune().minimize()
                                k = len(lm)  # indice du nouveau morceau
                                lf.append(k)  # nouvelle feuille
                                arbre[i].append(k)
                                arbre.append([])
                                d.append(copy(d[i]))
                                d[k][ij] = j2  # remplace la lettre
                                # ajoute le nouveau morceau à la liste
                                lm.append((a.intersection(a1), t1))
                            lf.remove(i)  # i n'est plus une feuille
                            if verb:
                                print("break...")
                            break  # le morceau n'est plus une feuille
                        else:
                            # remplace la lettre
                            d[i][ij] = k
                    tr += lm[j][1]
            if end:
                break
        # compute the substitution
        s = dict()
        for i in lf:
            if d[i][-1] < 0:
                d[i].pop()
            s[i] = d[i]
        # recode the substitution
        l = s.keys()
        dl = dict()  # inverse of l
        for i, k in enumerate(l):
            dl[k] = i+1
        d = dict()
        for i in s:
            d[dl[i]] = [dl[j] for j in s[i]]
        if get_aut:
            return d, [(a, t) for i, (a, t) in enumerate(lm) if arbre[i] == []]
        else:
            return d

    def compute_substitution2(self, DetAutomaton a, DetAutomaton ap=None,
                              np=None, lt=None, method=2, method_tr=1,
                              iplus=2, imax=None, need_included=True,
                              get_aut=False, verb=True):
        r"""
        Compute a substitution whose fixed point is the g-beta-expansion given
        by the beta-adic monoid with automaton a.

        INPUT:

        - ``a``- DetAutomaton
            Automaton of the g-beta-expansion.

        - ``ap``- DetAutomaton (default: ``None``)
            Language used to do the computations : we project everything on it.

        - ``np``- int (default: ``None``)
            Power of beta for the computing. The g-beta-expansion
            must be b^np invariant.

        - ``lt``- list of elements of the integer rings
            List of translations to compute the pieces exchange.

        - ``get_aut``- Bool (default: ``False``)
            If True, gives also the list of automata.

        - ``verb``- bool (default: ``True``)
          If True, print informations about the computing.

        OUTPUT:

        A word morphism given by a dictionnary.

        EXAMPLES::

            #. Full Tribonnacci::

                sage: m = BetaAdicMonoid((x^3-x^2-x-1).roots(ring=QQbar)[1][0], {0, 1})
                sage: m.compute_substitution(verb=False)          # long time
                {1: [1, 3], 2: [1], 3: [1, 2]}

        """
        m = self
        if ap is None:
            if hasattr(self, 'tss'):
                ap = DetAutomaton(self.tss)
            else:
                ap = DetAutomaton(None).full(list(self.C))
        if verb:
            print("ap=%s" % ap)
        cdef DetAutomaton aa
        if not a.included(ap):
            aa = a.copy()
            aa.zero_completeOP()
            # check that Qap contains Qaa
            if not m.Proj(ap, aa).equals_langages(aa) and need_included:
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
                ba.zero_completeOP()
                if m.Proj(aa, ba).equals_langages(ba):
                    ba = m.Proj(ba, aa)
                    np = i
                    break
            if np is None:
                raise ValueError('The g-beta-expansion must be b^np invariant for some natural integer np.')
        else:
            ba = aa.unshift(0, np)
            ba.zero_completeOP()
            if not m.Proj(aa, ba).equals_langages(ba):
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
        # calcule l'induction à partir de la liste de (morceau, translation)
        # précalculs
        if verb:
            print("Pre-computation...")
        arel = dict()
        for a, t in lt:
            arel[t] = m.relations_automaton4(t=-t, A=aa.A,
                                             B=ap.A, couples=True)
            if verb:
                print("arel[%s]=%s" % (t, arel[t]))
        if verb:
            print("ba : %s" % ba)
        # arbre de subdivision des morceaux
        arbre = [range(1, len(lt) + 1)] + [[] for i in range(len(lt))] 
        if verb:
            print("initial tree: %s" % arbre)
        lm = [(aa, 0)] + lt  # liste des morceaux, translations
        if verb:
            print("lm = %s" % lm)
        # parcours de chaque morceau (donné par une liste de morceaux)
        d = [[] for i in range(len(lm))]
        if verb:
            print("d = %s" % d)
        lf = range(1, len(lm))  # liste des feuilles
        if verb:
            print("lf = %s" % lf)

        from copy import copy

        if verb:
            print("\n**********************\n   Step 1   \n**********************")

        # étape 1 : complétion des mots
        for i, (a1, t1) in enumerate(lm):
            if arbre[i] != []:
                continue  # ce morceau n'est pas une feuille
            if verb:
                print("\nCalcul du morceau %s/%s (%s, %s)..." % (i, len(lm), a1, t1))
                # print "lf = %s"%lf
                # print "d = %s"%d
                # print "arbre = %s"%arbre
            tr = 0  # translation totale
            if d[i] != []:
                if d[i][-1] == -1:
                    continue  # le morceau était déjà fini de calculé
                # va à la fin du mot
                for j in d[i]:
                    if j < 0:
                        break
                    tr += lm[j][1]
            # calcule b^np*a + tr
            a = a1.unshift(0, np).prune().minimize()
            if tr != 0:
                if verb:
                    print("Translation de %s..." % tr)
            # m.move2(t=-tr, a=a)
            # TODO : ne pas recalculer cet automate déjà calculé
            a = m.Proj(a, ap, t=-tr)
            while True:
                # split selon les autres morceaux
                j = included(a, lf, lm)
                if j is None:
                    # détermine les morceaux qui intersectent a
                    l = []
                    for j in lf:
                        if lm[j][0].intersect(a):
                            l.append(j)
                    if len(l) < 2:
                        print("Error : intersection with %s piece but not included !!!"%len(l))
                    if verb:
                        print("Subdivision on %s pieces..." % len(l))
                    # calcule les intersections (découpe en morceaux a1)
                    for j in l:
                        a2 = lm[j][0]
                        # découpe a selon a2
                        # m.move2(t=tr, a=a2) #translate a2 de -tr
                        a = m.Proj(a2, ap.unshift(0, np), t=tr)
                        a.shiftOP(0, np)  # multiplie par b^(-np)
                        a = a.prune().minimize()
                        k = len(lm)  # indice du nouveau morceau
                        lf.append(k)  # nouvelle feuille
                        arbre[i].append(k)
                        arbre.append([])
                        from copy import copy
                        # print copy
                        d.append(copy(d[i]))
                        d[k].append(j)  # ajoute la translation suivante
                        # ajoute le nouveau morceau à la liste
                        lm.append((a.intersection(a1), t1))
                        # split selon ba
                        (ab, abc) = split_ba(k, tr+lm[j][1], np,
                                             lm, m, aa, ap, verb)
                        if ab is None:
                            if verb:
                                print("k=%s, tr=%s+%s : calcul à continuer" % (k, tr, lm[j][1]))
                        else:
                            if abc is None:
                                if verb:
                                    print("tr=%s : calcul terminé" % tr)
                                d[k].append(-1)  # indique que le calcul de ce morceau est terminé
                            else:
                                if verb:
                                    print("tr=%s : subdivision of %s by ba (new %s)..." % (tr, i, len(lm)))
                                lf.append(len(lm))  # nouvelle feuille
                                arbre[k].append(len(lm))
                                arbre.append([])
                                d.append(copy(d[k]))
                                # indique que le calcul est terminé pour ce morceau (pour l'étape 1)
                                d[len(lm)].append(-1)
                                lm.append((ab, t1))
                                lf.append(len(lm))  # nouvelle feuille
                                arbre[k].append(len(lm))
                                arbre.append([])
                                d.append(copy(d[k]))
                                lm.append((abc, t1))
                                lf.remove(k)  # le morceau k n'est plus une feuille
                    lf.remove(i)  # le morceau i n'est plus une feuille
                    # calcul fini pour ce morceau puisque ce n'est plus une feuille
                    break
                else:
                    # ajoute le morceau à la liste et translate
                    d[i].append(j)
                    # if verb: print "Translation de %s..."%lm[j][1]
                    # m.move2(t=-lm[j][1], a=a, ar=arel[lm[j][1]])
                    a = m.Proj(a, ap, t=-lm[j][1], arel=arel[lm[j][1]])  
                    tr += lm[j][1]
                # split selon ba
                (ab, abc) = split_ba(i, tr, np, lm, m, aa, ap, verb)
                if ab is None:
                    pass
                    # if verb: print "tr=%s : calcul à continuer"%tr
                else:
                    if abc is None:
                        if verb:
                            print("tr=%s : end of computation" % tr)
                    else:
                        if verb:
                            print("tr=%s : subdivision of %s by ba (new %s)..." % (tr, i, len(lm)))
                        lf.append(len(lm))  # nouvelle feuille
                        arbre[i].append(len(lm))
                        arbre.append([])
                        d.append(copy(d[i]))
                        # indique que le calcul est terminé pour ce morceau (pour l'étape 1)
                        d[len(lm)].append(-1)
                        lm.append((ab, t1))
                        lf.append(len(lm))  # nouvelle feuille
                        arbre[i].append(len(lm))
                        arbre.append([])
                        d.append(copy(d[i]))
                        lm.append((abc, t1))
                        lf.remove(i)  # le morceau i n'est plus une feuille
                    break  # calcul fini pour ce morceau (pour cette étape)

        if verb:
            print("\n*************\n   Step 2   \n*************")

        # étape 2 : remplacement des lettres qui ne sont pas des feuilles
        while True:
            end = True
            for i in lf:
                a1, t1 = lm[i]
                if verb:
                    print("\nPiece %s/%s..." % (i, len(lm)))
                tr = 0  # translation totale
                if d[i] == []:
                    print("Error : empty sheet !!!!")
                # va à la fin du mot
                for ij, j in enumerate(d[i]):
                    if j < 0:
                        break
                    if arbre[j] != []:  # il faut recalculer cette lettre
                        # calcule b^np*a + tr
                        a = a1.unshift(0, np).prune().minimize()
                        if tr != 0:
                            if verb:
                                print("Translation of %s..." % tr)
                        # m.move2(t=-tr, a=a)
                        # TODO : ne pas recalculer cet automate déjà calculé
                        a = m.Proj(a, ap, t=-tr)
                        # split selon les autres morceaux
                        f = fils(arbre, j)
                        if verb:
                            print("Split by %s pieces" % len(f))
                        k = included(a, f, lm)
                        if k is None:
                            end = False
                            # détermine les morceaux qui intersectent a
                            l = []
                            for k in lf:
                                if lm[k][0].intersect(a):
                                    l.append(k)
                            if len(l) < 2:
                                print("Error : intersection with %s pieces but not included !!!" % len(l))
                            if verb:
                                print("Subdivision of %s pieces..." % len(l))
                            # calcule les intersections (découpe en morceaux a1)
                            for j2 in l:
                                a2 = lm[j2][0]
                                # découpe selon a2
                                # a = m.move2(t=tr, a=a2) #translate a2 de -tr
                                # a.zero_completeOP()
                                a = m.Proj(a2, ap.unshift(0, np), t=tr)
                                a.shiftOP(0, np)  # multiplie par b^(-np)
                                a = a.prune().minimize()
                                k = len(lm)  # indice du nouveau morceau
                                lf.append(k)  # nouvelle feuille
                                arbre[i].append(k)
                                arbre.append([])
                                d.append(copy(d[i]))
                                d[k][ij] = j2  # remplace la lettre
                                # ajoute le nouveau morceau à la liste
                                lm.append((a.intersection(a1), t1))
                            lf.remove(i)  # i n'est plus une feuille
                            if verb:
                                print("break...")
                            break  # le morceau n'est plus une feuille
                        else:
                            # remplace la lettre
                            d[i][ij] = k
                    tr += lm[j][1]
            if end:
                break
        # compute the substitution
        s = dict()
        for i in lf:
            if d[i][-1] < 0:
                d[i].pop()
            s[i] = d[i]
        # recode the substitution
        l = s.keys()
        dl = dict()  # inverse of l
        for i, k in enumerate(l):
            dl[k] = i+1
        d = dict()
        for i in s:
            d[dl[i]] = [dl[j] for j in s[i]]
        if get_aut:
            return d, [(a, t) for i, (a, t) in enumerate(lm) if arbre[i] == []]
        else:
            return d

