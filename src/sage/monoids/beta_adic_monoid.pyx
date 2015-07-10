# coding=utf8
r"""
Beta-adic Monoids

AUTHORS:

- Paul Mercat (2013)

Beta-adic monoids are finitely generated monoids with generators of the form
	x -> beta*x + c
where beta is a element of a field (for example a complex number),
and c is varying in a finite set of numerals.
It permits to describe beta-adic expansions, that is writing of numbers of the form
	x = c_0 + c_1*beta + c_2*beta^2 + ...
for c_i's in a finite set of numerals.
"""

#*****************************************************************************
#  Copyright (C) 2013 Paul Mercat <mercatp@icloud.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#	This code is distributed in the hope that it will be useful, 
#	but WITHOUT ANY WARRANTY; without even the implied warranty 
#	of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# 
#  See the GNU General Public License for more details; the full text 
#  is available at:
#
#				  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer import Integer
from sage.rings.number_field.all import *
#from sage.structure.parent_gens import normalize_names
#from free_monoid_element import FreeMonoidElement

from sage.sets.set import Set
from monoid import Monoid_class
from sage.rings.qqbar import QQbar
from sage.rings.padics.all import *
from sage.combinat.words.automata import Automaton
#from sage.combinat.words.cautomata import FastAutomaton

#from sage.structure.factory import UniqueFactory
#from sage.misc.cachefunc import cached_method
include "sage/ext/interrupt.pxi"

#calcul de la valeur absolue p-adique (car non encore implémenté autrement)
def absp (c, p, d):
	return ((c.polynomial())(p).norm().abs())**(1/d)

#garde la composante fortement connexe de 0
def emonde (a, K):
	for s in a.strongly_connected_components_subgraphs():
		if K.zero() in s:
			return s

cdef extern from "complex.h":
	cdef cppclass Complexe:
		double x,y

cimport sage.combinat.words.cautomata
from sage.combinat.words.cautomata cimport Automate, FastAutomaton
from sage.combinat.words.cautomata import FastAutomaton

#cdef extern from "Automaton.h":
#	cdef cppclass Etat:
#		int* f
#		int final
#	cdef cppclass Automate:
#		int n
#		int na
#		Etat* e
#		int i

cdef extern from "Automaton.h":
	ctypedef char bool

cdef extern from "relations.h":
	cdef cppclass Element:
		int *c	#liste des n coeffs
	
	cdef cppclass PlaceArch:
		Complexe *c	#1, b, b^2, ... pour cette place

	#structure contenant les infos nécessaires pour calculer l'automate des relations
	cdef cppclass InfoBetaAdic:
		int n		#degré
		Element bn   #expression de b^n comme un polynome en b de degré < n
		Element *c   #liste des chiffres utilisés pour le calcul de l'automate des relations
		int nc	   #nombre de chiffre
		int ncmax	#nombre de chiffres alloués
		PlaceArch *p #liste des na places
		double *cM   #carré des valeurs absolues max
		int na		 #nombre de va
	
	InfoBetaAdic allocInfoBetaAdic (int n, int na, int ncmax, bool verb)
	void freeInfoBetaAdic (InfoBetaAdic iba)
	Automate RelationsAutomaton (InfoBetaAdic iba2, bool isvide, bool ext, bool verb)

cdef getElement (e, Element r, int n):
	cdef j
	p = e.lift()
	for j in range(n):
		r.c[j] = p[j]

cdef InfoBetaAdic initInfoBetaAdic (self, Cd=None, verb=False) except *:
	#compute all the data in sage
	K = NumberField((1/self.b).minpoly(), 'b', embedding=QQbar(1/self.b))
	b = K.gen()
	C = [c.lift()(1/b) for c in self.C]
	if verb:
		print "C = %s"%C
	
	if verb: print K
	
	#détermine les places qu'il faut considérer
	parch = []
	for p in K.places(): #places archimédiennes
		if p(b).abs() > 1:
			parch += [p]
	pi = K.defining_polynomial()
	from sage.rings.arith import gcd
	pi = pi/gcd(pi.list()) #rend le polynôme à coefficients entiers et de contenu 1
	if verb: print "pi=%s"%pi
	lp = (Integer(pi.list()[0])).prime_divisors() #liste des nombres premiers concernés
	if verb: print "lp=%s"%lp
	pultra = [] #liste des places ultramétriques considérées
	for p in lp:
		#détermine toutes les places au dessus de p dans le corps de nombres K
		k = Qp(p)
		Kp = k['a']
		a = Kp.gen()
		for f in pi(a).factor():
			kp = f[0].root_field('e')
			if kp == k:
				c = f[0].roots(kp)[0][0]
			else:
				c = kp.gen()
			if verb: print "c=%s (abs=%s)"%(c, (c.norm().abs())**(1/f[0].degree()))
			if (c.norm().abs())**(1/f[0].degree()) > 1:
				pultra += [(c, f[0].degree())]
	
	if verb: print "places: "; print parch; print pultra
	
	self.parch = parch
	
	if (len(pultra) > 0):
		raise ValueError("Not implemented for b algebraic non-integer.")
	
	#calcule les bornes max pour chaque valeur absolue
	if Cd is None:
		Cd = Set([c-c2 for c in C for c2 in C])
	else:
		Cd = [K(c).lift()(1/b) for c in Cd]
	if verb: print "Cd = %s"%Cd
	m = dict([])
	for p in parch:
		m[p] = max([p(c).abs() for c in Cd])/abs(1-p(p.domain().gen()).abs())
	for p, d in pultra:
		m[p] = max([absp(c, p, d) for c in Cd])
		
	#convert the data in C
	n = K.degree()
	na = len(parch)
	ncmax = 100
	cdef InfoBetaAdic i
	if verb: print "alloc..."
	i = allocInfoBetaAdic(n, na, ncmax, verb)
	cdef int j
	#initialise bn
	if verb: print "init bn..."
	getElement(b**n, i.bn, n)
	#initialise les places
	if verb: print "init places..."
	for k in range(na):
		for j in range(n):
			i.p[k].c[j] = complex(parch[k](b**j))
	#initialise les chiffres et bornes
	if verb: print "init chiffres..."
	initCdInfoBetaAdic(self, &i, Cd=Cd, verb=verb)
	return i

cdef initCdInfoBetaAdic (self, InfoBetaAdic *i, Cd, verb=False):
	#recalcule les bornes max pour chaque valeur absolue
#	if Cd is None:
#		Cd = Set([c-c2 for c in self.C for c2 in self.C])
	Cd = list(Cd)
	if verb: print "Cd = %s"%Cd
	m = dict([])
	for p in self.parch:
		m[p] = max([p(c).abs() for c in Cd])/abs(1-p(p.domain().gen()).abs())
#	for p, d in self.pultra:
#		m[p] = max([absp(c, p, d) for c in Cd])
	#conversion en C
	i.nc = len(Cd)
	if i.nc > i.ncmax:
		raise ValueError("Trop de chiffres : %d > %d max."%(i.nc, i.ncmax))
	for j,c in enumerate(Cd):
		getElement(c, i.c[j], i.n)
	for j,p in enumerate(self.parch):
		i.cM[j] = m[p]**2

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
		Complexe* t #liste des translations
		int n		#nombre de translations
		Automate a
	cdef cppclass BetaAdic2:
		Complexe b
		Complexe* t #liste des translations
		int n		#nombre de translations
		Automate* a
		int na
	ctypedef Color* ColorList
		
	Surface NewSurface (int sx, int sy)
	void FreeSurface (Surface s)
	ColorList NewColorList (int n)
	void FreeColorList (ColorList l)
	Color randColor (int a)
	Automate NewAutomate (int n, int na)
	void FreeAutomate(Automate a)
	void FreeAutomates(Automate* a, int n)
	BetaAdic NewBetaAdic (int n)
	void FreeBetaAdic (BetaAdic b)
	BetaAdic2 NewBetaAdic2 (int n, int na)
	void FreeBetaAdic2 (BetaAdic2 b)
	void Draw (BetaAdic b, Surface s, int n, int ajust, Color col, int verb)
	void Draw2 (BetaAdic b, Surface s, int n, int ajust, Color col, int verb)
	void DrawList (BetaAdic2 b, Surface s, int n, int ajust, ColorList lc, double alpha, int verb)
	void print_word (BetaAdic b, int n, int etat)

cdef Complexe complex (c):
	cdef Complexe r
	r.x = c.real()
	r.y = c.imag()
	return r

cdef Color getColor (c):
	if len(c) < 4:
		raise ValueError("Colors must be defined by 4 float numbers between 0 and 1.")
	cdef Color r
	r.r = c[0]*255
	r.g = c[1]*255
	r.b = c[2]*255
	r.a = c[3]*255
	return r

cdef surface_to_img (Surface s):
	import numpy as np
	from PIL import Image	
	arr = np.zeros([s.sy, s.sx], dtype = [('r', 'uint8'), ('g', 'uint8'), ('b', 'uint8'), ('a', 'uint8')])
	cdef int x, y
	cdef Color c
	for x in range(s.sx):
		for y in range(s.sy):
			c = s.pix[x][s.sy -y-1]
			arr[y,x][0] = c.r
			arr[y,x][1] = c.g
			arr[y,x][2] = c.b
			arr[y,x][3] = c.a
	img = Image.fromarray(arr, 'RGBA')
	img.save("/Users/mercat/Desktop/output.png")
	img.save("output.png")

cdef Automate getAutomate (a, d, list C, iss=None, verb=False):
	cdef int i
	if verb:
		print "getAutomate %s..."%a
	cdef FastAutomaton fa
	if isinstance(a, FastAutomaton):
		fa = a
		fa.permut_op(C, verb=verb)
		#fa = fa.permut(C, verb=verb)
		return fa.a[0]
	#assume in the following that a is a Automaton
	lv = a.vertices()
	if hasattr(a, 'F'):
		F = a.F
	else:
		F = lv
	#alloue l'automate
	cdef Automate r = NewAutomate(a.num_verts(), len(a.Alphabet()))
	#réindice les sommets
	dv = {}
	for u,i in zip(lv, range(len(lv))):
		dv[u] = i
		if u in F:
			r.e[i].final = 1
	if verb:
		print len(lv)
	#copie l'automate en C
	le = a.edges()
	if verb:
		print "len(le)=%s"%len(le)
	for u,v,l in le:
		if d.has_key(l):
			#if dv.has_key(u) and dv.has_key(v):
			r.e[dv[u]].f[d[l]] = dv[v]
			#else:
			#   print "Erreur : pas de clef %s ou %s !"%(u,v)
		else:
			print "Erreur : pas de clef %s !"%l
	if verb:
		print "I..."
	if iss is not None:
		r.i = iss
	else:
		if hasattr(a, 'I') and len(a.I) > 0:
			r.i = dv[list(a.I)[0]]
		else:
			r.i = -1
			#raise ValueError("The initial state must be defined !")
	if verb:
		print "...getAutomate"
	return r

cdef BetaAdic getBetaAdic (self, prec=53, ss=None, tss=None, iss=None, transpose=True, add_letters=True, verb=False):
	from sage.rings.complex_field import ComplexField
	CC = ComplexField(prec)
	cdef BetaAdic b
	if ss is None:
		if hasattr(self, 'ss'):
			ss = self.ss
		else:
			ss = self.default_ss()
	else:
		if transpose and tss is None:
			if verb:
				print "Calcul de la transposée...\n"
			tss = ss.transpose().determinize()
			self.tss = tss
			if verb:
				print tss
	if transpose:
		if tss is None:
			if hasattr(self, 'tss'):
				tss = self.tss
			else:
				if verb:
					print "Calcul de la transposée...\n"
				tss = ss.transpose().determinize()
				self.tss = tss
				if verb:
					print tss
		a = tss
	else:
		a = ss
		
	if add_letters:
		C = set(self.C)
		C.update(a.Alphabet())
	else:
		C = a.Alphabet()
	C = list(C)
	
	b = NewBetaAdic(len(C))
	b.b = complex(CC(self.b))
	d = {}
	for i,c in zip(range(b.n), C):
		b.t[i] = complex(CC(c))
		d[c] = i
	#automaton
	#if isinstance(a, FastAutomaton):
	#	a = a.permut(C, verb=verb)
	b.a = getAutomate(a, d, C=C, iss=iss, verb=verb)
	return b

cdef BetaAdic2 getBetaAdic2 (self, la=None, ss=None, tss=None, prec=53, add_letters=True, verb=False):
	if verb:
		print "getBetaAdic %s"%self
	from sage.rings.complex_field import ComplexField
	from sage.combinat.words.automata import Automaton
	CC = ComplexField(prec)
	cdef BetaAdic2 b
	if la is None:
		la = self.get_la(ss=ss, tss=tss, verb=verb)
	  
	if add_letters:
		C = set(self.C)
		for a in la:
			C.update(a.Alphabet())
	C = list(C)
	
	b = NewBetaAdic2(len(C), len(la))
	b.b = complex(CC(self.b))
	d = {}
	for i,c in zip(range(b.n), C):
		b.t[i] = complex(CC(c))
		d[c] = i
	#automata
	for i in range(len(la)):
		#if isinstance(la[i], FastAutomaton):
		#	la[i] = la[i].permut(C, verb=verb)
		b.a[i] = getAutomate(la[i], d, C=C, iss=None, verb=verb);
	return b

def PrintWord (m, n):
	b = getBetaAdic(m, prec=53, ss=None, tss=None, iss=None, transpose=False, add_letters=True, verb=False)
	print_word(b, n, b.a.i)

class BetaAdicMonoid(Monoid_class):
	r"""
	``b``-adic monoid with numerals set ``C``.
	It is the beta-adic monoid generated by the set of affine transformations ``{x -> b*x + c | c in C}``.
	
	EXAMPLES::
	
		sage: m1 = BetaAdicMonoid(3, {0,1,3})
		sage: m2 = BetaAdicMonoid((1+sqrt(5))/2, {0,1})
		sage: b = (x^3-x-1).roots(ring=QQbar)[0][0]
		sage: m3 = BetaAdicMonoid(b, {0,1})
	"""
	def __init__ (self, b, C):
		r"""
		Construction of the b-adic monoid generated by the set of affine transformations ``{x -> b*x + c | c in C}``.
		
		EXAMPLES::
			sage: m1 = BetaAdicMonoid(3, {0,1,3})
			sage: m2 = BetaAdicMonoid((1+sqrt(5))/2, {0,1})
			sage: b = (x^3-x-1).roots(ring=QQbar)[0][0]
			sage: m3 = BetaAdicMonoid(b, {0,1})
		"""
		#print "init BAM with (%s,%s)"%(b,C)
		if b in QQbar:
			#			print b
			pi = QQbar(b).minpoly()
			K = NumberField(pi, 'b', embedding = QQbar(b))
		else:
			K = b.parent()
			try:
				K.places()
			except:
				print "b=%s must be a algebraic number, ring %s not accepted."%(b,K)
			
#		print K
		self.b = K.gen() #beta (element of an NumberField)
#		print "b="; print self.b
		self.C = Set([K(c) for c in C]) #set of numerals
#		print "C="; print self.C
	
	def gen (self, i):
		r"""
		Return the element of C of index i.
		"""
#		g(x) = self.b*x+self.C[i]
		return self.C[i]
	
	def ngens (self):
		r"""
		Return the number of elements of C.
		"""
		return len(self.C)
	
	def _repr_ (self):
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
				str=" with subshift of %s states"%self.ss.num_verts()
		
		from sage.rings.rational_field import QQ
		if self.b in QQ:
			return "Monoid of %s-adic expansion with numerals set %s"%(self.b,self.C) + str
		else:
			K = self.b.parent()
			if K.base_field() == QQ:
				return "Monoid of b-adic expansion with b root of %s and numerals set %s"%(self.b.minpoly(),self.C) + str
			else:
				if K.characteristic() != 0:
					return "Monoid of b-adic expansion with b root of %s and numerals set %s, in characteristic %s"%(self.b.minpoly(), self.C, K.characteristic()) + str
				else:
					return "Monoid of b-adic expansion with b root of %s and numerals set %s"%(K.modulus(),self.C) + str
		
	def default_ss (self, C=None):
		r"""
		Returns the full subshift (given by an Automaton) corresponding to the beta-adic monoid.

		EXAMPLES::

			sage: m=BetaAdicMonoid((1+sqrt(5))/2, {0,1})
			sage: m.default_ss()
			Finite automaton with 1 states
		"""
		if C is None:
			C = self.C
		ss = Automaton()
		ss.allow_multiple_edges(True)
		ss.allow_loops(True)
		ss.add_vertex(0)
		for c in C:
			ss.add_edge(0, 0, c)
		ss.I = [0]
		ss.F = [0]
		ss.A = C
		return ss
		
	#liste des automates donnant le coloriage de l'ensemble limite
	def get_la (self, ss=None, tss=None, verb=False):
		if hasattr(self, 'la'):
			return self.la
		if tss is None:
			if hasattr(self, 'tss'):
				tss = self.tss
				if not isinstance(tss, FastAutomaton):
					tss = FastAutomaton(tss)
		if ss is None:
			if hasattr(self, 'ss'):
				ss = self.ss
				if not isinstance(ss, FastAutomaton):
					ss = FastAutomaton(ss)
		if ss is None:
			if tss is None:
				tss = FastAutomaton(None).full(list(self.C))
				#raise ValueError("la, ss, or tss must be defined !")
			if verb:
				print "Compute the transposition of tss=%s..."%tss
			ss = tss.transpose().determinise() #ze()
			if verb:
				print ss
				print "simplify..."
			ss = ss.emonde_inf().emonde() #emonde0_simplify()
			if verb:
				print ss
		if tss is None:
			if ss is None:
				ss = FastAutomaton(None).full(list(self.C)) #self.default_ss()
				#raise ValueError("la, ss, or tss must be defined !")
			if verb:
				print "Compute the transposition of ss=%s..."%ss
			tss = ss.transpose().determinise() #ze()
			if verb:
				print tss
				print "simplify..."
			tss = tss.emonde_inf().emonde() #0_simplify()
			if verb:
				print tss

		#compute la
		a = {}
		for v in ss.states(): #ss.vertices():
			a[v] = ss.copy() #Automaton(ss) ##################################################################################
			a[v].set_final_states([v])
			if verb:
				print "Compute the transposition..."
			a[v] = a[v].transpose().determinise() #ze()
			if verb:
				print a[v]
				print "simplify..."
			a[v] = a[v].emonde_inf().emonde()
			if verb:
				print a[v]
		return [tss]+a.values()
	
	def points_exact (self, n=None, ss=None, iss=None):
		r"""
		Returns a set of exacts values (in the number field of beta) corresponding to the drawing of the limit set of the beta-adic monoid.

		INPUT:
		
		- ``n`` - integer (default: ``None``)
		  The number of iterations used to plot the fractal.
		  Default values: between ``5`` and ``16`` depending on the number of generators.
		
		- ``ss`` - Automaton (default: ``None``)
		  The subshift to associate to the beta-adic monoid for this drawing.
		
		- ``iss`` - set of initial states of the automaton ss (default: ``None``)
		
		OUTPUT:

			list of exact values

		EXAMPLES:

		#. The dragon fractal::

			sage: m=BetaAdicMonoid(1/(1+i), {0,1})
			sage: P = m.points_exact()
			sage: len(P)
			65536
		"""
		#global co
		
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
				print "Donner l'état de départ iss de l'automate ss !"
				return
		
		if n is None:
			if ng == 2:
				n = 16
			elif ng == 3:
				n = 9
			else:
				n = 5
		if n == 0:
			#donne un point au hasard dans l'ensemble limite
			#co = co+1
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
					#orbit_points = orbit_points.union([b*p+c for p in self.points_exact(n=n-1, ss=ss, iss=v)])
			#orbit_points0 = self.points_exact(n-1)
			#orbit_points = Set([])
			#for c in C:
			#	v = self.succ(i, c)
			#	if v is not None:
			#	   orbit_points = orbit_points.union(Set([b*p+c for p in orbit_points0]))
		#print "no=%s"%orbit_points.cardinality()
		return orbit_points
	
	def points (self, n=None, place=None, ss=None, iss=None, prec=53):
		r"""
		Returns a set of values (real or complex) corresponding to the drawing of the limit set of the beta-adic monoid.

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
		
		OUTPUT:

			list of real or complex numbers

		EXAMPLES:

		#. The dragon fractal::

			sage: m=BetaAdicMonoid(1/(1+I), {0,1})
			sage: P = m.points()
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
			n = int(5.2/-log(abs(self.b.N(prec=prec))))
		
		from sage.rings.complex_field import ComplexField
		CC = ComplexField(prec)
		if place is None:
			if abs(b) < 1:
				#garde la place courante
				#place = lambda x: CC(x.n())
				return [CC(c).N(prec) for c in self.points_exact(n=n, ss=ss, iss=iss)]
			else:
				#choisis une place
				places = K.places()
				place = places[0]
				for p in places:
					if abs(p(b)) < 1:
						place = p
						#break
		
	   #from sage.rings.qqbar import QQbar
		#from sage.rings.qqbar import QQbar, AA
		#if QQbar(self.b) not in AA:
		#	#print "not in AA !"
		#	return [(place(c).conjugate().N().real(), place(c).conjugate().N().imag()) for c in self.points_exact(n=n, ss=ss, iss=iss)]
		#else:
		#	#print "in AA !"
		#	return [place(c).conjugate().N() for c in self.points_exact(n=n, ss=ss, iss=iss)]
		return [place(c).N(prec) for c in self.points_exact(n=n, ss=ss, iss=iss)]
	
#		  if n == 0:
#			 #donne un point au hasard dans l'ensemble limite
#			 return [0]
#		 else:
#			 orbit_points0 = self.points(n-1)
#			 orbit_points = Set([])
#			 for c in C:
#				 orbit_points = orbit_points.union(Set([place(b)*p+place(c) for p in orbit_points0]))
#		 return orbit_points
	
	def plot2 (self, n=None, tss=None, ss=None, iss=None, sx=800, sy=600, ajust=True, prec=53, color=(0,0,0,255), method=0, add_letters=True, verb=False):
		r"""
		Draw the limit set of the beta-adic monoid (with or without subshift).

		INPUT:
		
		- ``n`` - integer (default: ``None``)
		  The number of iterations used to plot the fractal.
		  Default values: between ``5`` and ``16`` depending on the number of generators.
		
		- ``place`` - place of the number field of beta (default: ``None``)
		  The place used to evaluate elements of the number field.
		
		- ``ss`` - Automaton (default: ``None``)
		  The subshift to associate to the beta-adic monoid for this drawing.
		
		- ``iss`` - set of initial states of the automaton ss (default: ``None``)
		
		- ``sx, sy`` - dimensions of the resulting image (default : ``800, 600``)
		
		- ``ajust`` - adapt the drawing to fill all the image, with ratio 1 (default: ``True``)
		
		- ``prec`` - precision of returned values (default: ``53``)
		
		- ``color`` - list of three integer between 0 and 255 (default: ``(0,0,255,255)``)
		  Color of the drawing.
		
		- ``verb`` - bool (default: ``False``)
		  Print informations for debugging.
		
		OUTPUT:

			A Graphics object.

		EXAMPLES:

		#. The dragon fractal::

			sage: m=BetaAdicMonoid(1/(1+I), {0,1})
			sage: m.plot2()	 # long time

		#. The Rauzy fractal of the Tribonacci substitution::

			sage: s = WordMorphism('1->12,2->13,3->1')
			sage: m = s.rauzy_fractal_beta_adic_monoid()
			sage: m.plot2()	 # long time
		
		#. A non-Pisot Rauzy fractal::
			
			sage: s = WordMorphism({1:[3,2], 2:[3,3], 3:[4], 4:[1]})
			sage: m = s.rauzy_fractal_beta_adic_monoid()
			sage: m.b = 1/m.b
			sage: m.plot2(tss=m.ss)	 # long time
		
		#. The dragon fractal and its boundary::

			sage: m = BetaAdicMonoid(1/(1+I), {0,1})
			sage: ssi = m.intersection_words(w1=[0], w2=[1])	 # long time
			sage: m.plot2(tss=ssi)							   # long time
			sage: m.plot2()									  # long time
			
		#. The "Hokkaido" fractal and its boundary::
		  
			sage: s = WordMorphism('a->ab,b->c,c->d,d->e,e->a')
			sage: m = s.rauzy_fractal_beta_adic_monoid()
			sage: ssi = m.intersection_words(w1=[0], w2=[1])				  # long time
			sage: m.plot2(la=[la[0], ssi]+la[1:], colormap='gist_rainbow')	# long time
		
		#. A limit set that look like a tiling::
			
			sage: P=x^4 + x^3 - x + 1
			sage: b = P.roots(ring=QQbar)[2][0]
			sage: m = BetaAdicMonoid(b, {0,1})
			sage: m.plot2(19)								   # long time
		
		"""
		sig_on()
		cdef Surface s = NewSurface (sx, sy)
		cdef BetaAdic b
		b = getBetaAdic(self, prec=prec, tss=tss, ss=ss, iss=iss, add_letters=add_letters, transpose=True, verb=verb)
		#if verb:
		#	printAutomaton(b.a)
		#dessin
		cdef Color col
		col.r = color[0]
		col.g = color[1]
		col.b = color[2]
		col.a = color[3]
		if n is None:
			n = -1
		if method == 0:
			Draw(b, s, n, ajust, col, verb)
		elif method == 1:
			print "Not implemented !"
			return
			#lv = s.rauzy_fractal_projection_exact().values()
			#for i,v in zip(range(len(lv)),lv):
			#		b.t[i] = complex(CC(v))
			#Draw2(b, s, n, ajust, col, verb)
		#enregistrement du résultat
		surface_to_img(s)
		if verb:
			print "Free..."
		FreeSurface(s)
		if not isinstance(tss, FastAutomaton):
			FreeAutomate(b.a)
		FreeBetaAdic(b)
		sig_off()
		
	def plot3 (self, n=None, la=None, ss=None, tss=None, sx=800, sy=600, ajust=True, prec=53, colormap = 'hsv', backcolor=None, opacity = 1., add_letters=True, verb=False):
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
		
		- ``iss`` - set of initial states of the automaton ss (default: ``None``)
		
		- ``sx, sy`` - dimensions of the resulting image (default : ``800, 600``)
		
		- ``ajust`` - adapt the drawing to fill all the image, with ratio 1 (default: ``True``)
		
		- ``prec`` - precision of returned values (default: ``53``)
		
		- ``colormap`` - list of colors (default: ``hsv``)
		  Colors of the drawing.
		
		- ``opacity``- float (default: ``1.``)
		  Transparency of the drawing.
		
		- ``verb`` - bool (default: ``False``)
		  Print informations for debugging.
		
		OUTPUT:

			A Graphics object.

		EXAMPLES:

		#. The dragon fractal::

			sage: m=BetaAdicMonoid(1/(1+I), {0,1})
			sage: m.plot2()	 # long time

		#. The Rauzy fractal of the Tribonacci substitution::

			sage: s = WordMorphism('1->12,2->13,3->1')
			sage: m = s.rauzy_fractal_beta_adic_monoid()
			sage: m.plot3()	 # long time
		
		#. A non-Pisot Rauzy fractal::
			
			sage: s = WordMorphism({1:[3,2], 2:[3,3], 3:[4], 4:[1]})
			sage: m = s.rauzy_fractal_beta_adic_monoid()
			sage: m.b = 1/m.b
			sage: m.plot2(tss=m.ss)	 # long time
		
		#. The dragon fractal and its boundary::

			sage: m = BetaAdicMonoid(1/(1+I), {0,1})
			sage: ssi = m.intersection_words(w1=[0], w2=[1])							   # long time
			sage: m.plot2(tss=ssi) #plot the boundary									  # long time
			sage: m.plot3(la=[m.default_ss(), ssi], colormap=[(.5,.5,.5,.5), (0,0,0,1.)])  # long time
			
		#. The "Hokkaido" fractal and its boundary::
		  
			sage: s = WordMorphism('a->ab,b->c,c->d,d->e,e->a')
			sage: m = s.rauzy_fractal_beta_adic_monoid()
			sage: ssi = m.intersection_words(w1=[0], w2=[1])				# long time
			sage: m.plot3(la=[la[0], ssi]+la[1:], colormap='gist_rainbow')  # long time
		
		#. A limit set that look like a tiling::
			
			sage: P=x^4 + x^3 - x + 1
			sage: b = P.roots(ring=QQbar)[2][0]
			sage: m = BetaAdicMonoid(b, {0,1})
			sage: m.plot2(19)								   # long time
		
		"""
		cdef Surface s = NewSurface (sx, sy)
		cdef BetaAdic2 b
		sig_on()
		b = getBetaAdic2(self, la=la, ss=ss, tss=tss, prec=prec, add_letters=add_letters, verb=verb)
		sig_off()
		#dessin
		if n is None:
			n = -1
		
		# Manage colors
		if backcolor is None:
			backcolor = (.5, .5, .5, .5)
		cdef ColorList cl
		sig_on()
		cl = NewColorList(b.na)
		if isinstance(colormap, list):
			if b.na < len(colormap):
				raise ValueError("The list of color must contain at least %d elements."%b.na)
			for i in range(b.na):
				if i < len(colormap):
					cl[i] = getColor(colormap[i])
				else:
					cl[i] = randColor(255)
		elif isinstance(colormap, str):
			from matplotlib import cm
			if not colormap in cm.datad.keys():
				raise RuntimeError("Color map %s not known (type 'from matplotlib import cm' and look at cm for valid names)" % colormap)
			colormap = cm.__dict__[colormap]
			cl[0] = getColor(backcolor)
			for i in range(b.na-1):
				cl[i+1] = getColor(colormap(float(i)/float(b.na-1)))
		else:
			raise TypeError("Type of option colormap (=%s) must be list of colors or str" % colormap)
		sig_off()
		
		sig_on()
		DrawList(b, s, n, ajust, cl, opacity, verb)
		sig_off()
		#enregistrement du résultat
		sig_on()
		surface_to_img(s)
		sig_off()
		if verb:
			print "Free..."
		sig_on()
		FreeSurface(s)
		if la is None:
			FreeAutomates(b.a, b.na)
		else:
			for i,a in enumerate(la):
				if not isinstance(a, FastAutomaton):
					FreeAutomate(b.a[i])
		FreeBetaAdic2(b)
		FreeColorList(cl)
		sig_off()
		
	def plot (self, n=None, place=None, ss=None, iss=None, prec=53, point_size=None, color='blue', verb=False):
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

		EXAMPLES:

		#. The dragon fractal::

			sage: m=BetaAdicMonoid(1/(1+I), {0,1})
			sage: m.plot()	 # long time

		#. The Rauzy fractal of the Tribonacci substitution::

			sage: s = WordMorphism('1->12,2->13,3->1')
			sage: m = s.rauzy_fractal_beta_adic_monoid()
			sage: m.plot()	 # long time
		
		#. A non-Pisot Rauzy fractal::
			
			sage: s = WordMorphism({1:[3,2], 2:[3,3], 3:[4], 4:[1]})
			sage: m = s.rauzy_fractal_beta_adic_monoid()
			sage: m.b = 1/m.b
			sage: m.ss = m.ss.transpose().determinize().minimize()
			sage: m.plot()	 # long time
		
		#. The dragon fractal and its boundary::

			sage: m = BetaAdicMonoid(1/(1+I), {0,1})
			sage: p1 = m.plot()								  # long time
			sage: ssi = m.intersection_words(w1=[0], w2=[1])	 # long time
			sage: p2 = m.plot(ss = ssi, n=18)					# long time
			sage: p1+p2										  # long time
			
		#. The "Hokkaido" fractal and its boundary::
		  
			sage: s = WordMorphism('a->ab,b->c,c->d,d->e,e->a')
			sage: m = s.rauzy_fractal_beta_adic_monoid()
			sage: p1 = m.plot()									 # long time
			sage: ssi = m.intersection_words(w1=[0], w2=[1])		# long time
			sage: p2 = m.plot(ss=ssi, n=40)						 # long time
			sage: p1+p2											 # long time
		
		#. A limit set that look like a tiling::
			
			sage: P=x^4 + x^3 - x + 1
			sage: b = P.roots(ring=QQbar)[2][0]
			sage: m = BetaAdicMonoid(b, {0,1})
			sage: m.plot(18)									# long time
		
		"""
		
		global co
		
		co = 0
		orbit_points = self.points(n=n, place=place, ss=ss, iss=iss, prec=prec)
		if verb: print "co=%s"%co
		
		# Plot points size
		if point_size is None:
			point_size = 1
			
		# Make graphics
		from sage.plot.plot import Graphics
		G = Graphics()
		
		#dim = self.b.minpoly().degree()
		
		from sage.rings.qqbar import QQbar, AA
		if QQbar(self.b) not in AA: #2D plots
			from sage.all import points
			G = points(orbit_points, size=point_size, color=color)
		else: # 1D plots
			from sage.all import plot
			G += plot(orbit_points, thickness=point_size, color=color)
#			if plotbasis:
#				from matplotlib import cm
#				from sage.plot.arrow import arrow
#				canonicalbasis_proj = self.rauzy_fractal_projection(eig=eig, prec=prec)
#				for i,a in enumerate(alphabet):
#					x = canonicalbasis_proj[a]
#					G += arrow((-1.1,0), (-1.1,x[0]), color=cm.__dict__["gist_gray"](0.75*float(i)/float(size_alphabet))[:3])
#		else:
#			print "dimension too large !"
		G.set_aspect_ratio(1)

		return G
	
	def _relations_automaton_rec (self, current_state, di, parch, pultra, m, Cd, ext, verb=False, niter=-1):
		r"""
		Used by relations_automaton()
		"""
		
		if niter == 0:
			return di
		
		global count
		if verb and count%10000 == 0: print count
		if count == 0:
			return di
		count -= 1
		
		b = self.b
		if not di.has_key(current_state):
			di[current_state] = dict([])
		for c in Cd: #parcours les transitions partant de current_state
			#e = b*current_state + c #calcule l'état obtenu en suivant la transition c
			e = (current_state + c)/b #calcule l'état obtenu en suivant la transition c
			#if verb: print "b=%s, e=%s, cur=%s, c=%s, di=%s"%(b, e, current_state, c, di)
			if not di.has_key(e): #détermine si l'état est déjà dans le dictionnaire
				ok = True
				#calcule les valeurs abolues pour déterminer si l'état n'est pas trop grand
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
					continue #cesse de considérer cette transition
				for p, d in pultra:
					if absp(e, p, d) > m[p]:
						#if verb: print "abs(%s)=%s trop grand !"%(e, absp(e, p, d))
						ok = False
						break
				if ok:
					#on ajoute l'état et la transition à l'automate
					di[current_state][e] = c
					di = self._relations_automaton_rec (current_state=e, di=di, parch=parch, pultra=pultra, m=m, Cd=Cd, ext=ext, verb=verb, niter=niter-1)
			else:
				#ajoute la transition
				di[current_state][e] = c
		return di
	
	def relations_automaton (self, ext=False, ss=None, noss=False, Cd=None, verb=False, step=100, limit=None, niter=None):
		r"""
		Compute the relations automaton of the beta-adic monoid (with or without subshift).
		See http://www.latp.univ-mrs.fr/~paul.mercat/Publis/Semi-groupes%20fortement%20automatiques.pdf for a definition of such automaton (without subshift).
		
		INPUT:
		
		- ``ext`` - bool (default: ``False``)
		  If True, compute the extended relations automaton (which permit to describe infinite words in the monoid).
		
		- ``ss`` - Automaton (default: ``None``)
		  The subshift to associate to the beta-adic monoid for this operation.
		
		- ``noss`` - bool (default: ``False``)
		  
		
		- ``verb`` - bool (default: ``False``)
		  If True, print informations for debugging.
				
		- ``step`` - int (default: ``100``)
		  Stop to an intermediate state of the computing to verify that all is right.
		 
		- ``limit``- int (default: None)
		  Stop the computing after a number of states limit.
		
		OUTPUT:

		A Automaton.

		EXAMPLES::

			sage: m = BetaAdicMonoid(3, {0,1,3})
			sage: m.relations_automaton()
			Finite automaton with 3 states
			
			sage: b = (x^3-x-1).roots(ring=QQbar)[0][0]
			sage: m = BetaAdicMonoid(b, {0,1})
			sage: m.relations_automaton()
			Finite automaton with 179 states
		
		REFERENCES:

		..  [Me13] Mercat P.
			Bull. SMF 141, fascicule 3, 2013.
		
		"""
		
		if not noss:
			a = self.relations_automaton(ext=ext, ss=None, noss=True, Cd=Cd, verb=verb, step=step, limit=limit)
			if not step:
				return a
			step = step-1
			if ss is None:
				if hasattr(self, 'ss'):
					ss = self.ss
				else:
					return a #pas de sous-shift
			if not step:
				return ss
			step = step-1
			d=dict()
			for u in self.C:
				for v in self.C:
					if not d.has_key(u-v):
						d[u-v] = []
					d[u - v] += [(u,v)]
			if not step:
				return d
			step = step-1
			ss = ss.emonde0_simplify()
			P = ss.product(A=ss)
			#P = P.emonde0_simplify()
			if not step:
				return P
			step = step-1
			a.relabel2(d)
			if not step:
				return a
			step = step-1
			a = a.intersection(A=P)
			if not step:
				return a
			step = step-1
			a = a.emonde0_simplify()
			if not step:
				return a
			step = step-1
			if not ext:
				a.emondeF()
				if not step:
					return a
				step = step-1
			#a = a.determinize(A=a.A, noempty=True)
			#if not step:
			#	return a
			#step = step-1
			#return a
			return a.minimize()
		
		K = self.C[0].parent()
		b = self.b
		
		if verb: print K
		
		#détermine les places qu'il faut considérer
		parch = []
		for p in K.places(): #places archimédiennes
			if p(b).abs() < 1:
				parch += [p]
		pi = K.defining_polynomial()
		from sage.rings.arith import gcd
		pi = pi/gcd(pi.list()) #rend le polynôme à coefficients entiers et de contenu 1
#		den = pi.constant_coefficient().denominator()
#		lp = (pi.list()[pi.degree()].numerator()*den).prime_divisors() #liste des nombres premiers concernés
		lp = (Integer(pi.list()[0])).prime_divisors() #liste des nombres premiers concernés
		pultra = [] #liste des places ultramétriques considérées
		for p in lp:
			#détermine toutes les places au dessus de p dans le corps de nombres K
			k = Qp(p)
			Kp = k['a']
			a = Kp.gen()
			for f in pi(a).factor():
				kp = f[0].root_field('e')
#				c = kp.gen()
				if kp == k:
					c = f[0].roots(kp)[0][0]
				else:
					c = kp.gen()
				if verb: print "c=%s (abs=%s)"%(c, (c.norm().abs())**(1/f[0].degree()))
				if (c.norm().abs())**(1/f[0].degree()) < 1: #absp(c, c, f[0].degree()) > 1:
					pultra += [(c, f[0].degree())]
		
		if verb: print "places: "; print parch; print pultra
		
		#calcule les bornes max pour chaque valeur absolue
		if Cd is None:
			Cd = Set([c-c2 for c in self.C for c2 in self.C])
		if verb: print "Cd = %s"%Cd
		m = dict([])
		for p in parch:
			m[p] = max([p(c).abs() for c in Cd])/abs(1-p(p.domain().gen()).abs())
		for p, d in pultra:
			m[p] = max([absp(c, p, d) for c in Cd])
		
		if verb: print "m = %s"%m
		
		if verb: print K.zero().parent()
		
		global count
		#print limit
		if limit is None:
			count = -1
		else:
			count = limit
		if niter is None:
			niter = -1
		#print count
		if verb: print "Parcours..."
		di = self._relations_automaton_rec (current_state=K.zero(), di=dict([]), parch=parch, pultra=pultra, m=m, Cd=Cd, ext=ext, verb=verb, niter=niter)
		
		if count == 0:
			print "Nombre max d'états atteint."
		else:
			if verb:
				if limit is None:
					print "%s états parcourus."%(-1-count)
				else:
					print "%s états parcourus."%(limit-count)
		
		#a = Automaton([K.zero()], [K.zero()], di)
		
		#if verb: print "di = %s"%di
		
		res = Automaton(di, loops=True) #, multiedges=True)
		
		if verb: print "Avant emondation : %s"%res
		
		res.I = [K.zero()]
		res.A = Cd #Set([c-c2 for c in self.C for c2 in self.C])
		if verb: print "Emondation..."
		if not ext:
			res.F = [K.zero()]
			res.emonde()
		else:
			#res = res.emonde0_simplify() #pour retirer les états puits
			res.emonde0()
			res.F = res.vertices()
		return res
	
	def relations_automaton2 (self, verb=False, step=100, limit=None, niter=None):
		r"""
		
		Do the same as relations_automaton, but avoid recursivity in order to avoid the crash of sage.
		
		INPUT:
				
		- ``verb`` - bool (default: ``False``)
		  If True, print informations for debugging.
				
		- ``step`` - int (default: ``100``)
		  Stop to an intermediate state of the computing to verify that all is right.
		 
		- ``limit``- int (default: None)
		  Stop the computing after a number of states limit.
		
		OUTPUT:

		A Automaton.
		"""
		
		K = self.C[0].parent()
		b = self.b
		
		if verb: print K
		
		#détermine les places qu'il faut considérer
		parch = []
		for p in K.places(): #places archimédiennes
			if p(b).abs() < 1:
				parch += [p]
		pi = K.defining_polynomial()
		from sage.rings.arith import gcd
		pi = pi/gcd(pi.list()) #rend le polynôme à coefficients entiers et de contenu 1
		lp = (Integer(pi.list()[0])).prime_divisors() #liste des nombres premiers concernés
		pultra = [] #liste des places ultramétriques considérées
		for p in lp:
			#détermine toutes les places au dessus de p dans le corps de nombres K
			k = Qp(p)
			Kp = k['a']
			a = Kp.gen()
			for f in pi(a).factor():
				kp = f[0].root_field('e')
				if kp == k:
					c = f[0].roots(kp)[0][0]
				else:
					c = kp.gen()
				if verb: print "c=%s (abs=%s)"%(c, (c.norm().abs())**(1/f[0].degree()))
				if (c.norm().abs())**(1/f[0].degree()) < 1: #absp(c, c, f[0].degree()) > 1:
					pultra += [(c, f[0].degree())]
		
		if verb: print "places: "; print parch; print pultra
		
		#calcule les bornes max pour chaque valeur absolue
		Cd = Set([c-c2 for c in self.C for c2 in self.C])
		if verb: print "Cd = %s"%Cd
		m = dict([])
		for p in parch:
			m[p] = max([p(c).abs() for c in Cd])/abs(1-p(p.domain().gen()).abs())
		for p, d in pultra:
			m[p] = max([absp(c, p, d) for c in Cd])
		
		if verb: print "m = %s"%m
		
		if verb: print K.zero().parent()
		
		if limit is None:
			count = -1
		else:
			count = limit
		if niter is None:
			niter = -1

		if verb: print "Parcours..."
		
		di=dict([])
		S = [K.zero()] #set of states to look at
		iter = 0
		while len(S) != 0:
			if iter == niter:
				break
			for s in S:
				S.remove(s)
				if not di.has_key(s):
					di[s] = dict([])
					if count%10000 == 0:
						print count
					count-=1
					if count == 0:
						iter = niter-1 #to break the main loop
						break
				for c in Cd: #parcours les transitions partant de current_state
					e = (s + c)/b #calcule l'état obtenu en suivant la transition c
					#if verb: print "b=%s, e=%s, cur=%s, c=%s, di=%s"%(b, e, current_state, c, di)
					if not di.has_key(e): #détermine si l'état est déjà dans le dictionnaire
						ok = True
						#calcule les valeurs abolues pour déterminer si l'état n'est pas trop grand
						for p in parch:
							if p(e).abs() >= m[p]:
								ok = False
								break
						if not ok:
							continue #cesse de considérer cette transition
						for p, d in pultra:
							if absp(e, p, d) > m[p]:
								#if verb: print "abs(%s)=%s trop grand !"%(e, absp(e, p, d))
								ok = False
								break
						if ok:
							#on ajoute l'état et la transition à l'automate
							di[s][e] = c
							S.append(e)
						else:
							#ajoute la transition
							di[s][e] = c
							iter+=1
							
				if count == 0:
					print "Nombre max d'états atteint."
					return
				else:
					if verb:
						if limit is None:
							print "%s états parcourus."%(-1-count)
						else:
							print "%s états parcourus."%(limit-count)
							
		res = Automaton(di, loops=True) #, multiedges=True)
		
		if verb: print "Avant emondation : %s"%res
		
		res.I = [K.zero()]
		res.A = Set([c-c2 for c in self.C for c2 in self.C])
		res.F = [K.zero()]
		if verb: print "Emondation..."
		res.emonde()
		return res
	
	def relations_automaton3 (self, isvide=False, Cd=None, ext=False, verb=False):
		r"""
			Compute the relation automaton of the beta-adic monoid.
			For beta algebraic integer only.
			If isvide is True, it only checks if the automaton is trivial or not.
			Cd is the set of differences A-B where A and B are the alphabets to compare.
			ext : automate des relations à l'infini ou pas.
		"""
		if Cd is None:
			Cd = Set([c-c2 for c in self.C for c2 in self.C])
		Cd = list(Cd)
		sig_on()
		cdef InfoBetaAdic ib
		ib = initInfoBetaAdic(self, Cd=Cd, verb=verb)
		cdef Automate a
		a = RelationsAutomaton(ib, isvide, ext, verb)
		r = FastAutomaton(None)
		r.a[0] = a
		r.A = Cd
		freeInfoBetaAdic(ib)
		sig_off()
		if isvide:
			return a.na != 0
		if ext:
			r2 = r.emonde_inf()
			r2.set_final_states(r2.states())
		else:
			r2 = r.emonde()
		return r2
	
	def critical_exponent_aprox (self, niter=10, verb=False):
		b = self.b
		K = b.parent()
		C = self.C
		S = set([K.zero()])
		for i in range(niter):
			S2 = set([])
			for s in S:
				for c in C:
					S2.add((s+c)/b)
			#intervertit S et S2
			S3 = S2
			S2 = S
			S = S3
			if verb: print len(S)
		from sage.functions.log import log
		print "%s"%(log(len(S)).n()/(niter*abs(log(abs(b.n())))))
	
	def complexity (self, verb = False):
		r"""
		Return a estimation of an upper bound of the number of states of the relations automaton.
		
		INPUT:
		
		 - ``verb`` - Boolean (default: False) Display informations for debug.
		
		OUTPUT:

		A positive real number.

		EXAMPLES::

			sage: m = BetaAdicMonoid(3, {0,1,3})
			sage: m.complexity()
			3.0...
		"""
		K = self.C[0].parent()
		b = self.b
		
		if verb: print K
		
		#détermine les places qu'il faut considérer
		parch = K.places()
		r = len(parch)
		pi = K.defining_polynomial()
		from sage.rings.arith import gcd
		pi = pi/gcd(pi.list()) #rend le polynôme à coefficients entiers et de contenu 1
#		den = pi.constant_coefficient().denominator()
#		lp = (pi.list()[pi.degree()].numerator()*den).prime_divisors() #liste des nombres premiers concernés
		lp = (Integer(pi.list()[0])).prime_divisors() #liste des nombres premiers concernés
		pultra = [] #liste des places ultramétriques considérées
		for p in lp:
			#détermine toutes les places au dessus de p dans le corps de nombres K
			k = Qp(p)
			Kp = k['a']
			a = Kp.gen()
			for f in pi(a).factor():
				kp = f[0].root_field('e')
#				c = kp.gen()
				if kp == k:
					c = f[0].roots(kp)[0][0]
				else:
					c = kp.gen()
				if verb: print "c=%s (abs=%s)"%(c, (c.norm().abs())**(1/f[0].degree()))
				if c.norm().abs() != 1: #absp(c, c, f[0].degree()) > 1:
					pultra += [(c, f[0].degree())]
		
		if verb: print "places: "; print parch; print pultra
		
		#calcule les bornes max pour chaque valeur absolue
		Cd = Set([c-c2 for c in self.C for c2 in self.C])
		if verb: print "Cd = %s"%Cd
		vol = 1.
		#from sage.rings.real_mpfr import RR
		for p in parch:
			if (p(b)).imag() == 0:
				vol *= 2*max([p(c).abs() for c in Cd])/abs(1-p(b).abs())
				if verb: print "place réelle %s"%p
			else:
				vol *= 3.1415926535*(max([p(c).abs() for c in Cd])/abs(1-p(b).abs()))**2
				if verb: print "place complexe %s"%p
			#vol *= max([p(c).abs() for c in Cd])/abs(1-p(p.domain().gen()).abs())
			#vol *= max(1, max([p(c).abs() for c in Cd])/abs(1-p(p.domain().gen()).abs()))
		for p, d in pultra:
			vol *= max([(c.polynomial())(p).norm().abs() for c in Cd])
			#vol *= max([absp(c, p, d) for c in Cd])
			#vol *= max(1, max([absp(c, p, d) for c in Cd]))
		#from sage.functions.other import sqrt
		#return vol/(K.regulator()*(sqrt(r+1).n()))
		return vol
	
	#def infinite_relations_automaton (self, verb=False):
	#	a = self.relations_automaton (ext=True, verb=verb)
	#	#retire la composante fortement connexe de l'etat initial
	#	K = self.b.parent()
	#	for s in a.strongly_connected_components_subgraphs():
	#		if K.zero() in s:
	#			a.delete_vertices(a.strongly_connected_components_subgraphs()[0].vertices())
	#			return a
	
	def intersection (self, ss, ss2=None, Iss=None, Iss2=None, ext=True, verb = False): #calcule le sous-shift correspondant à l'intersection des deux monoïdes avec sous-shifts
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

		EXAMPLES:
			
			#. Compute the boundary of the dragon fractal (see intersection_words for a easier way) ::

				sage: m = BetaAdicMonoid(1/(1+I), {0,1})
				sage: ss=m.default_ss()
				sage: iss=ss.I[0]
				sage: ss0 = ss.prefix(w=[0], i=iss)
				sage: ss1 = ss.prefix(w=[1], i=iss)
				sage: ssi = m.intersection(ss=ss0, ss2=ss1)
				sage: ssd = ssi.determinize(A=m.C, noempty=True)
				sage: ssd = ssd.emonde0_simplify()
				sage: m.plot(ss = ssd, n=19)	 # long time
		"""
		
		m = None
		
		if ss2 is None:
			if hasattr(self, 'ss'):
				ss2 = self.ss
			else:
				raise ValueError("Only one sub-shift given, I need two !")
			if type(ss) == type(BetaAdicMonoid(2,{0,1})):
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
			if verb: print "Iss = %s, Iss2 = %s"%(Iss, Iss2)
		
		a = ss.product(ss2)
		if verb: print "Product = %s"%a
		
		ar = self.relations_automaton(ext=ext, noss=True)
		if verb: print "Arel = %s"%ar
		
		#maps actual edges to the list of corresponding couple
		m = dict([])
		for c in self.C:
			for c2 in self.C:
				if m.has_key(c-c2):
					m[c-c2] += [(c,c2)]
				else:
					m[c-c2] = [(c,c2)]
		if verb: print "m = %s"%m
		
		L = a.Alphabet() #a.edge_labels()
		LA = ar.Alphabet() #ar.edge_labels()
		d = dict([])
		for u, v in L:
			for ka in LA:
				for u2, v2 in m[ka]:
					if u == u2 and v == v2:
						d[((u,v), ka)] = u
						break
					else:
						d[((u,v), ka)] = None
		if verb: print "d = %s"%d
		p = a.product(A=ar, d=d)
		#I = [((i,i2),self.b.parent().zero()) for i in Iss for i2 in Iss2]
		#if verb: print "I = %s"%I
		#p.emondeI(I=I)
		#if verb: print "%s"%p
		#p.emonde0(I=I)
		p.I = [((i,i2), self.b.parent().zero()) for i,i2 in zip(Iss, Iss2)]
		p = p.emonde0_simplify()
		if m is not None:
			ssd = p.determinize2(A=self.C, noempty=True)
			ssd = ssd.emonde0_simplify()
			return ssd
		return p
	
	def intersection2 (self, ss1, ss2, verb = False): #calcule le sous-shift correspondant à l'intersection des deux monoïdes avec sous-shifts, utilise des FastAutomaton
		r"""
		Compute the intersection of two beta-adic monoid with subshifts given by FastAutomaton
		
		INPUT:

		- ``ss``- FastAutomaton (default: ``None``)
		  The first subshift to associate to the beta-adic monoid for this operation.
		
		- ``ss2``- FastAutomaton (default: ``None``)
		  The second subshift to associate to the beta-adic monoid for this operation.
		
		- ``verb``- bool (default: ``False``)
		  If True, print informations for debugging.
		
		OUTPUT:

		A FastAutomaton.

		EXAMPLES:
			
			#. Compute the boundary of the dragon fractal (see intersection_words for a easier way) ::

				sage: m = BetaAdicMonoid(1/(1+I), {0,1})
				sage: import sage.combinat.words.cautomata
				sage: from sage.combinat.words.cautomata import FastAutomaton
				sage: ss0 = FastAutomaton([(0,1,0)]+[(1,1,l) for l in m.C], i=0, final_states=[1])
				sage: ss1 = FastAutomaton([(0,1,1)]+[(1,1,l) for l in m.C], i=0, final_states=[1])
				sage: ssi = m.intersection2(ss0, ss1)
				sage: m.plot2(tss = ssi)	 # long time
		"""
		a = self.relations_automaton3()
		a = a.emonde_inf()
		a.set_final_states(a.states())
		ssp = ss1.product(ss2)
		ssp = ssp.emonde()
		d = {}
		for (la1,la2) in ssp.Alphabet():
			for lb in a.Alphabet():
				if lb == la1-la2:
					d[((la1,la2), lb)] = (la1,la2)
		ssi = ssp.product(a, d)
		ssi = ssi.emonde_inf()
		ssi = ssi.emonde()
		d = {}
		for (l1,l2) in ssi.Alphabet():
			d[(l1,l2)] = l1
		ssi = ssi.determinise_proj(d)
		ssi = ssi.emonde_inf()
		ssi = ssi.emonde()
		return ssi.minimise()
		
	def intersection_words (self, w1, w2, ss=None, iss=None):
		r"""
		Compute the intersection of two beta-adic monoid with subshifts corresponding to two prefix
		
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

		EXAMPLES:
			
			#. Compute the boundary of the dragon fractal::

				sage: m = BetaAdicMonoid(1/(1+I), {0,1})
				sage: m.intersection_words(w1=[0], w2=[1])
				Finite automaton with 21 states
			
			#. Draw the intersection of two sub-sets of a limit set::
				
				sage: m = BetaAdicMonoid(1/(1+I), {0,1})
				sage: ssi = m.intersection_words(w1=[0], w2=[1])
				sage: m.plot(n=10, ss=ssi)						# long time
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
		ssd = ssd.emonde0_simplify()
		return ssd
	
	
	def reduced_words_automaton2 (self, step=100, verb=False, transpose=False): #, FastAutomaton aut=None):
		r"""
		Compute the reduced words automaton of the beta-adic monoid (without considering the automaton of authorized words).
		See http://www.latp.univ-mrs.fr/~paul.mercat/Publis/Semi-groupes%20fortement%20automatiques.pdf for a definition of such automaton.
		Use fast C-functions but works only for algebraic integer.
		(Consider using reduced_words_automaton() if you're not in this case.)
		
		INPUT: 
		
		- ``verb`` - bool (default: ``False``)
		  If True, print informations for debugging.
		
#		- ``aut`` - FastAutomaton (default: ``None``, full language)
#		  Automaton describing the language in which we live.
		
		OUTPUT:

		FastAutomaton.
		"""
		
		#compute the relations automaton
		Cd = list(set([c-c2 for c in self.C for c2 in self.C]))
		Cdp = [k for k in range(len(Cd)) if Cd[k] in [self.C[j]-self.C[i] for i in range(len(self.C)) for j in range(i)]] #indices des chiffres strictements négatifs dans Cd
		arel = self.relations_automaton3(Cd=Cd, ext=False)
		arel = arel.emonde()
		if transpose:
			arel = arel.transpose_det()
		if verb:
			print "arel = %s"%arel
		if step == 1:
			return arel
		
		#add a new state
		cdef int ne, ei
		ei = arel.initial_state()
		ne = arel.n_states() #new added state
		arel.add_state(True)
		arel.set_final_state(ei, final=False) #it is the new final state
		if step == 2:
			return arel
		
		#add edges from the new state (copy edges from the initial state)
		for j in range(len(Cd)):
			arel.set_succ(ne, j, arel.succ(ei, j))
		if step == 3:
			return arel
		
		#suppress some edges from the initial state
		for j in Cdp:
			arel.set_succ(ei, j, -1)
		if step == 4:
			return arel
			
		#change edges that point to the initial state : make them point to the new state
		for e in arel.states():
			if e != ei:
				for j in range(len(Cd)):
					if arel.succ(e, j) == ei:
						arel.set_succ(e, j, ne)
		if step == 5:
			return arel
		
		#project, determinise and take the complementary
		d = {}
		for a in self.C:
			for b in self.C:
				if not d.has_key(a-b):
					d[a-b] = []
				d[a-b].append((a,b))
		if verb:
			print d
		arel = arel.duplicate(d) #replace differences with couples
		d = {}
		for j in self.C:
			for i in self.C:
				d[(i,j)] = i
		if verb:
			print d
			print arel
		arel = arel.determinise_proj(d, noempty=False, nof=True) #, verb=True) #project on the first value of the couple, determinise and take the complementary
		if verb:
			print arel
		arel = arel.emonde()
		arel = arel.emonde_inf()
		if step == 10:
			return arel
		return arel.minimise()
	
	
	
	
	def reduced_words_automaton (self, ss=None, Iss=None, ext=False, verb=False, step=None, arel=None):
		r"""
		Compute the reduced words automaton of the beta-adic monoid (with or without subshift).
		See http://www.latp.univ-mrs.fr/~paul.mercat/Publis/Semi-groupes%20fortement%20automatiques.pdf for a definition of such automaton (without subshift).
		
		WARNING: It seems there is a bug : result may be incorrect if ss is not None.
		
		INPUT:
		
		- ``ss``- Automaton (default: ``None``)
		  The first subshift to associate to the beta-adic monoid for this operation.
		
		- ``Iss``- set of states of ss (default: ``None``)
				
		- ``ext`` - bool (default: ``True``)
		  If True, compute the extended relations automaton (which permit to describe infinite words in the monoid).  
		
		- ``verb`` - bool (default: ``False``)
		  If True, print informations for debugging.
		
		- ``step`` - int (default: ``None``)
		  Stop to a intermediate state of the computing to make verifications.
		
		- ``arel`` - Automaton (default: ``None``)
		  Automaton of relations.
		
		OUTPUT:

		A Automaton.

		EXAMPLES:
			
			#. 3-adic expansion with numerals set {0,1,3}::
			
				sage: m = BetaAdicMonoid(3, {0,1,3})
				sage: m.reduced_words_automaton()
				Finite automaton with 2 states
			
			#. phi-adic expansion with numerals set {0,1}::
			
				sage: m = BetaAdicMonoid((1+sqrt(5))/2, {0,1})
				sage: m.reduced_words_automaton()
				Finite automaton with 3 states
			
			#. beta-adic expansion with numerals set {0,1} where beta is the plastic number::
				sage: b = (x^3-x-1).roots(ring=QQbar)[0][0]
				sage: m = BetaAdicMonoid(b, {0,1})
				sage: m.reduced_words_automaton()		# long time
				Finite automaton with 5321 states
		"""
		
		if ss is None:
			if hasattr(self, 'ss'):
				ss = self.ss
				if hasattr(self.ss, 'I'):
					Iss = self.ss.I
		
		if step is None:
			step = 1000
		
		K = self.C[0].parent()
		
		if verb: print "Calcul de l'automate des relations..."
		
		if arel is None:
			a = self.relations_automaton(noss=True)
		else:
			a = arel
		
		if verb: print " -> %s"%a
		
		if step == 1:
			return ("automate des relations", a)
		
		#add a state copy of K.0 (it will be the new initial state)
		a.add_vertex('O')
		
#		#add transitions to K.0 to 'O'
#		for f, d, l in a.incoming_edges(K.zero(), labels=True):
#			if f == K.zero():
#				a.add_edge('O', 'O', l)
#			else:
#				a.add_edge(f, 'O', l)
		
		#subset of positives labels
		Cdp = []
		for i in range(self.C.cardinality()):
			for j in range(i):
				Cdp += [self.C[i] - self.C[j]]
		
		#redirect positives transitions from K.0
		for f, d, l in a.outgoing_edges(K.zero(), labels=True):
			if l in Cdp:
#				a.delete_edge(K.zero(), d, l)
				#add the edge
				a.add_edge('O', d, l)
				
		a.add_edge('O', 'O', a.edge_label(K.zero(), K.zero()))
		
		#if verb: print a.incoming_edges(K.zero(), labels=True)
		
		#remove outgoing edges from K.0 (except from K.0 to K.0)
		for f, d, l in a.outgoing_edges(K.zero(), labels=True):
			if f != d:
				a.delete_edge(f, d, l)
		
		if step == 2:
			a.I = ['O']
			a.F = Set([K.zero()])
			return ("automate des relations ordonnées", a)
		
		a.emondeI(I=['O'])
		
		if step == 3:
			return ("automate des relations ordonnées, émondé", a)
		
		if ss is not None: #not full sub-shift
			if Iss is None:
				Iss = [ss.vertices()[0]]
			#maps actual edges to the list of corresponding couple
			m = dict([])
			for c in self.C:
				for c2 in self.C:
					if m.has_key(c-c2):
						m[c-c2] += [(c,c2)]
					else:
						m[c-c2] = [(c,c2)]
			#if verb: print "m=%s"%m
			
			#calculate the 'product to the right' of a with ss			
			d = dict([])
			La = a.edge_labels()
			Lss = ss.edge_labels()
			for ka in La: 
				for kss in Lss:
					d[(ka,kss)] = None
					for k in m[ka]:
						if k[1] == kss:
							d[(ka,kss)] = k[0]
							break
							
			#if verb: print "d=%s"%d
			if verb: print "avant produit : a=%s (%s etats)"%(a, a.num_verts())
			a = a.product(A=ss, d=d)
			if verb: print "après produit : a=%s"%a
			if step == 4:
				return ("automate des mots généraux non réduits", a)
			
			I = [('O',i) for i in Iss]
			nof = Set([(K.zero(),i) for i in ss.vertices()])
			
			#if verb: print "I=%s, F=%s"%(I, nof)
			
			if ext:
				#a.emondeI(I=I)
				#a.emonde0(I=I) #pour retirer les états puits
				a = a.emonde0_simplify(I=I)
			else:
				a = a.emonde0_simplify(I=I)
				a.emonde(I=I, F=nof)
			#a.emondeI(I=I)
			#a.emondeF(F=nof)
			#if step == 4:
			#	return ("automate des mots généraux non réduits, émondé", a)
			#a.emondeF(F=nof)	
			
			if verb: print "après émondation : a=%s"%a
			if step == 5:
				return ("automate des mots généraux non réduits, émondé", a)
			
			#return a
		else:
			#maps actual edges to element of self.C (the left part when writted c-c2)
			m = dict([])
			for c in self.C:
				for c2 in self.C:
					if m.has_key(c-c2):
						m[c-c2] += [c]
					else:
						m[c-c2] = [c]
			#if verb: print "m=%s"%m
			
			a.allow_multiple_edges(True)
			#replace each label by its mapping
			for f, d, l in a.edges():
				a.delete_edge(f, d, l)
				for l2 in m[l]:
					a.add_edge(f, d, l2)
			
			I = ['O']
			nof=Set([K.zero()])
		
		a.I = I
		a.F = nof
		a.C = self.C
		
		if verb: print "avant determinisation : a=%s"%a
		if step == 6:
			return ("automate des mots généraux non réduits, émondé", a)
		
		#rend l'automate plus simple
		a = a.emonde0_simplify()
		
		if verb: print "simplification : a=%s"%a
		
		if verb: print "Determinization..."
		#determinize
		ad = a.determinize2(nof=a.F)
		#ad = a.determinize(nof=a.F, verb=False)
		#ad = a.determinize(I, self.C, nof, verb=verb)
		
		if verb: print " -> %s"%ad
		if step == 7:
			return ("automate des mots généraux réduits", ad)
		
		if ss is not None: #not full sub-shift
			#calculate the intersection with ss
			ad = ad.emonde0_simplify()
			ad = ad.intersection(ss)
			if verb: print "apres intersection : a=%s"%ad
		
		if step == 8:
			return ("automate des mots réduits", ad)
		
		#F2=[e for e in a.vertices() nof in e[0]]
		#if verb: print "I2=%s"%I2 #, F2=%s"%(I2,F2)
		ad.A = self.C
		#ad.emondeI(I=I2) #, F=F2)
		ad = ad.emonde0_simplify()
		ad.F = ad.vertices()
		
		if verb: print "apres émondation : a=%s"%ad
		if step == 9:
			return ("automate des mots réduits, émondé", ad)
		
		return ad
	
	def critical_exponent_free (self, ss=None, prec = None, verb=False):
		r"""
		Compute the critical exponent of the beta-adic monoid (with or without subshift), assuming it is free.
		See http://www.latp.univ-mrs.fr/~paul.mercat/Publis/Semi-groupes%20fortement%20automatiques.pdf for a definition (without subshift). 
		
		INPUT:
		
		- ``ss``- Automaton (default: ``None``)
		  The first subshift to associate to the beta-adic monoid for this operation.
		
		- ``prec``- precision (default: ``None``)
				
		- ``verb``- bool (default: ``False``)
		  If True, print informations for debugging.
		
		OUTPUT:

		A real number.

		EXAMPLES:
			
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
			
				sage: m = BetaAdicMonoid(1/(1+I), {0,1})
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
				sage: m.ss = m.ss.transpose().determinize().minimize()
				sage: m.critical_exponent_free()
				log(y)/log(|1/2*b^2 - 1/2*b + 1/2|) where y is the max root of x^3 - x^2 + x - 2
				1.5485260383...
		"""
		
		if ss is None:
			if hasattr(self, 'ss'):
				ss = self.ss
		if ss is None:
			y = self.C.cardinality()
			print "log(%s)/log(|%s|)"%(y, self.b)
		else:
			if verb: print ""
			M = ss.adjacency_matrix()
			if verb: print "Valeurs propres..."
			e = M.eigenvalues()
			if verb: print "max..."
			y = max(e, key=abs)
			if verb: print ""
			print "log(y)/log(|%s|) where y is the max root of %s"%(self.b, QQbar(y).minpoly())
			y = y.N(prec)
		from sage.functions.log import log
		b = self.b.N(prec)
		if verb: print "y=%s, b=%s"%(y, b)
		return abs(log(y)/log(abs(b))).N()
		
	def critical_exponent (self, ss=None, prec = None, verb=False):
		r"""
		Compute the critical exponent of the beta-adic monoid (with or without subshift).
		See http://www.latp.univ-mrs.fr/~paul.mercat/Publis/Semi-groupes%20fortement%20automatiques.pdf for a definition (without subshift). 
		
		INPUT:
		
		- ``ss``- Automaton (default: ``None``)
		  The first subshift to associate to the beta-adic monoid for this operation.
		
		- ``prec``- precision (default: ``None``)
				
		- ``verb``- bool (default: ``False``)
		  If True, print informations for debugging.
		
		OUTPUT:

		A real number.

		EXAMPLES:
			
			#. Hausdorff dimension of limit set of 3-adic expansion with numerals set {0,1,3}::
			
				sage: m = BetaAdicMonoid(3, {0,1,3})
				sage: m.critical_exponent()
				log(y)/log(|3|) where y is the max root of x^2 - 3*x + 1
				0.8760357589...
			
			#. Hausdorff dimension of limit set of phi-adic expansion with numerals set {0,1}::
			
				sage: m = BetaAdicMonoid((1+sqrt(5))/2, {0,1})
				sage: m.critical_exponent()
				log(y)/log(|b|) where y is the max root of x^2 - x - 1
				1.0000000000...
			
			#. Critical exponent that is not the Hausdorff dimension of the limit set::
			
				sage: P = x^7 - 2*x^6 + x^3 - 2*x^2 + 2*x - 1
				sage: b = P.roots(ring=QQbar)[3][0]
				sage: m = BetaAdicMonoid(b, {0,1})
				sage: m.critical_exponent()					# long time
				log(y)/log(|b|) where y is the max root of x^11 - 2*x^10 - 4*x^2 + 8*x + 2
				3.3994454205...
			
			#. See more examples with doc critical_exponent_free()
			
		"""
		#	 #. Critical exponent that is not the Hausdorff dimension of the limit set::
		#
		#	 sage:

		if ss is None:
			if hasattr(self, 'ss'):
				ss = self.ss
		if verb:
			print "Calcul de l'automate des mots réduits...\n"
		a = self.reduced_words_automaton(ss=ss, verb=verb)
		return self.critical_exponent_free (prec=prec, ss=a, verb=verb)
	
	#test if 0 is an inner point of the limit set
	def ZeroInner (self, verb=False):
	
		if not hasattr(self, 'ss'):
			self.ss = self.default_ss()
		
		if verb: print "relations automaton..."
		
		ar = self.relations_automaton(ext=True)
		ar.complementary()
		
		if verb: print "complementaire : %s"%ar
		#return ar
		
		a = self.default_ss().product(self.ss)
		
		if verb: print "a = %s"%a
		#return a
		
		#a = ar.intersection(a)
		
		if verb: print "product..."
		
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
		
		if verb: print "product = %s"%a2
		#return a2
		
		#test if there is a cycle in the graph
		if a2.is_directed_acyclic():
			print "0 is not an inner point."
		else:
			ss = self.ss
			self.ss = None
			print "Zero is an inner point iff the %s has non-empty interior."%self
			self.ss = ss
		
	#complete the langage of a
	def complete (self, FastAutomaton a, C=None, ext=False, arel=None, verb=False):
		r"""
		Return an automaton that recognize the language of all words over C that represent elements recognized by a.
		If ext is True, this also include words equal at infinity.
	
		INPUT:
	
		- ``a`` - A FastAutomaton.
	
		- ``C`` - list of digits (default : ``None``).
		"""
		if C is None:
			C = self.C
		ap = FastAutomaton(None).full(list(C)).product(a)
		if ext:
			ap = ap.prefix_closure()
		if verb:
			if ap.n_states() < 100:
				ap.plot()
			print "ap=%s"%ap
		d = dict()
		Cd = [c - c2 for c in C for c2 in a.A]
		for c in Cd:
			d[c] = []
		for c in C:
			for c2 in a.A:
				d[c - c2].append((c,c2))
		if arel is None:
			arel = self.relations_automaton3(Cd=Cd, ext=ext).duplicate(d)
		else:
			arel = arel.duplicate(d)
		if verb:
			if arel.n_states() < 100:
				arel.plot()
			print "arel=%s"%arel
		ai = ap.intersection(arel)
		if ext:
			ai = ai.emonde_inf()
		ai = ai.emonde()
		ai = ai.minimise()
		if verb:
			if ai.n_states() < 100:
				ai.plot()
			print "ai=%s"%ai
		d = dict()
		for c in C:
			for c2 in a.A:
				d[(c,c2)] = c
		ac = ai.determinise_proj(d)
		if ext:
			ac = ac.emonde_inf()
		ac = ac.emonde()
		ac = ac.minimise()
		return ac
	
	#donne l'automate décrivant l'adhérence de l'ensemble limite avec un nouvel alphabet C
	def adherence (self, tss=None, C=None, C2=None, ext=False, verb=False, step=None):
		if tss is None:
			if hasattr(self, 'tss'):
				tss = self.tss
			else:
				tss = self.default_ss()
		if C is None:
			C = list(set(self.C))
		if C2 is None:
			C2 = list(set(tss.Alphabet()))
		if verb:
			print "Calcul de l'automate des relations..."
		Cd = list(set([c1-c2 for c1 in C2 for c2 in C]))
		if verb:
			print "Cd=%s"%Cd
		a = self.relations_automaton3(Cd=Cd, ext=ext)
		if verb:
			print " -> %s"%a
		if step == 1:
			return a
		if ext:
			a = a.emonde_inf()
		a = a.emonde()
		if verb:
			print " Après émondation : %s"%a
		if step == 2:
			return a
		d={}
		for c1 in C2:
			for c2 in C:
				if not d.has_key(c1-c2):
					d[c1-c2] = []
				d[c1-c2].append((c1,c2))
		if verb:
			print d
		a2 = a.duplicate(d, verb=verb)
		if verb:
			print a2.Alphabet()
			print a2
		if step == 3:
			return a2
		ap = tss.product(FastAutomaton(self.default_ss(C=C)), verb=verb)
		if ext:
			ap = ap.prefix_closure()
		if step == 4:
			return ap
		a2 = ap.intersection(a2)
		if step == 5:
			return a2
		if ext:
			a2 = a2.emonde_inf()
		a2 = a2.emonde()
		if step == 6:
			return a2
		a2 = a2.minimise()
		if step == 7:
			return a2
		if verb:
			print "déterminise..."
		d={}
		for c1,c2 in a2.Alphabet():
			d[(c1,c2)] = c2
		a2 = a2.determinise_proj(d, verb=verb)
		if step == 8:
			return a2
		if verb:
			print " -> %s"%a2
		if ext:
			a2 = a2.emonde_inf()
		a2 = a2.emonde()
		if step == 9:
			return a2
		if verb:
			print "Après simplification : %s"%a2
		return a2.minimise()
		
	#donne l'automate décrivant le translaté de +t, avec les chiffres C
	def move (self, t, FastAutomaton tss=None, list C=None, step=None):
		if tss is None:
			if hasattr(self, 'tss'):
				if isinstance(self.tss, FastAutomaton):
					tss = self.tss
				else:
					tss = FastAutomaton(self.tss)
			else:
				tss = FastAutomaton(self.default_ss())
		if C is None:
			C = list(set(self.C))
		
		A = tss.Alphabet()
		k = self.b.parent()
		nA = list(set([k(a+t2) for a in A for t2 in [0,t]]))
		a = tss.bigger_alphabet(nA)
		
		#add a new state
		cdef int ne, ei
		ei = a.initial_state()
		ne = a.n_states() #new added state
		a.add_state(a.is_final(ei))
		a.set_initial_state(ne) #it is the new initial state
		if step == 2:
			return a
		
		#add edges from the new state (copy edges from the initial state and move them)
		cdef s
		for j in range(len(A)):
			a.set_succ(ne, nA.index(A[j] + t), tss.succ(ei, j))
		if step == 3:
			return a

		#compute the adherence of the new automaton
		return self.adherence(tss=a, C=C, C2=nA)
	
	#calcule l'intersection des ensembles limites
	def intersection2 (self, FastAutomaton a, FastAutomaton b, ext=True):
		a2 = self.complete(a, ext=ext)
		b2 = self.complete(b, ext=ext)
		return a2.intersection(b2).emonde()
	
	#determine if the limit sets intersect
	def intersect (self, FastAutomaton a, FastAutomaton b, ext=True, verb=False):
		a2 = self.complete(a, ext=ext)
		if verb: print "a2=%s"%a2
		b2 = self.complete(b, ext=ext)
		if verb: print "b2=%s"%b2
		return not a2.intersection(b2).is_empty(ext)
	
	#Dit si toutes les tuiles du pavages autosimilaire sont connexes ou pas
	def is_all_connected (self, FastAutomaton a=None, ext=True, verb=False):
		if a is None:
			if hasattr(self, 'tss'):
				if isinstance(self.tss, FastAutomaton):
					a = m.tss
				else:
					a = FastAutomaton(self.tss)
			else:
				a = FastAutomaton(None).full(list(self.C))
		
		from sage.graphs.graph import Graph
		n = a.n_states()
		na = len(a.Alphabet())
		d = dict([]) #dictionnaire des automates complétés
		if verb: print "Automate des relations..."
		arel = self.relations_automaton3(ext=ext)
		for i in range(n):
			if verb: print "morceau %s"%i
			g = Graph({j:{} for j in range(na) if a.succ(i, j) != -1})
			#compute the neighboorhood graph of the piece i
			for j in g.vertices():
				for k in g.vertices():
					if k >= j:
						continue
					if verb: print " intersection %s et %s..."%(j,k)
					la = []
					la.append(a.piece(j, e=i).minimise())
					la.append(a.piece(k, e=i).minimise())
					for l in range(2):
						a2 = d.get(la[l]) #récupère le complété dans le dictionnaire s'il y est
						if a2 is None:
							if verb: print "  complete %s..."%la[l]
							a2 = self.complete(la[l], ext=ext)
							d[la[l]] = a2
							la[l] = a2
						else:
							la[l] = a2
							if verb: print "  déjà calculé !"
					if verb: print "  intersect..."
					#if self.intersect(a.piece(j, e=i), a.piece(k, e=i), ext=ext, verb=verb):
					if la[0].intersect(la[1], ext=ext):
						if verb: print "  oui !"
						g.add_edge(j,k)
			if not g.is_connected():
				return False
		return True
	
	#Dit si l'ensemble limite est connexe ou pas
	def is_connected (self, FastAutomaton a=None):
		if a is None:
			if hasattr(self, 'tss'):
				a = FastAutomaton(self.tss)
			else:
				a = FastAutomaton(None).full(self.C)
		
		n = a.n_states()
		na = len(a.Alphabet())
		rules = [[[l] for l in a.Alphabet()] for i in range(n)]
		gvois = [Graph(na) for i in range(n)] #graphe des morceaux voisins
		gnvois = [Graph(na) for i in range(n)] #graphe des morceaux non voisins
		
		#liste des morceaux dont il faut tester la connexité du graphe de voisinage
		m = [a.initial_state()]
		
		while len(m) > 0:
			i = m.pop()
			gvois[i].connected_component_containing_vertex(i)
			
		raise ValueError("Not implemented !")

	#Dit si l'ensemble limite est simplement connexe ou pas
	def is_simply_connected (self, FastAutomaton a=None):
		if a is None:
			if hasattr(self, 'tss'):
				a = FastAutomaton(self.tss)
			else:
				a = FastAutomaton(None).full(self.C)
		
		#TODO!
		
		raise ValueError("Not implemented !")




