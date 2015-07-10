# coding=utf8
"""
Finite state machines using C

AUTHORS:

- Paul Mercat
"""

#*****************************************************************************
#	   Copyright (C) 2014 Paul Mercat <mercatp@icloud.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#	This code is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#	General Public License for more details.
#
#  The full text of the GPL is available at:
#
#				  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.stdlib cimport malloc, free

cimport sage.combinat.words.cautomata
include "sage/ext/interrupt.pxi"

#ctypedef Automate Automaton

cdef extern from "automataC.h":
	ctypedef Automate Automaton
	ctypedef NAutomate NAutomaton
	cdef cppclass Dict:
		int* e
		int n
	cdef cppclass InvertDict:
		Dict* d
		int n
	
	Automaton NewAutomaton (int n, int na)
	void FreeAutomaton (Automaton *a)
	int hashAutomaton (Automaton a)
	void FreeNAutomaton (NAutomaton *a)
	Automaton CopyAutomaton (Automaton a, int nalloc, int naalloc)
	Automaton PieceAutomaton (Automaton a, int *w, int n, int e)
	void init (Automaton *a)
	void printAutomaton (Automaton a)
	void plotTikZ (Automaton a, const char **labels, const char *graph_name, double sx, double sy)
	void NplotTikZ (NAutomaton a, const char **labels, const char *graph_name, double sx, double sy)
	Automaton Product(Automaton a1, Automaton a2, Dict d)
	Automaton Determinise (Automaton a, Dict d, bool noempty, bool onlyfinals, bool nof, bool verb)
	Automaton DeterminiseN (NAutomaton a, bool puits)
	Automaton emonde_inf (Automaton a, bool verb)
	Automaton emonde (Automaton a, bool verb)
	Automaton emondeI (Automaton a, bool verb)
	bool equalsAutomaton (Automaton a1, Automaton a2)
	Dict NewDict (int n)
	void FreeDict (Dict d)
	void printDict (Dict d)
	InvertDict NewInvertDict (int n)
	void FreeInvertDict (InvertDict id)
	void printInvertDict (InvertDict id)
	Automaton Duplicate (Automaton a, InvertDict id, int na2, bool verb)
	Automaton TransposeDet (Automaton a)
	NAutomaton Transpose (Automaton a)
	int StronglyConnectedComponents (Automaton a, int *res)
	Automaton SubAutomaton (Automaton a, Dict d, bool verb)
	Automaton Permut (Automaton a, int *l, int na, bool verb)
	void PermutOP (Automaton a, int *l, int na, bool verb)
	Automaton Minimise (Automaton a, bool verb)
	void DeleteVertexOP (Automaton* a, int e)
	Automaton DeleteVertex (Automaton a, int e)
	bool equalsLangages (Automaton *a1, Automaton *a2, Dict a1toa2, bool minimized)
	bool emptyLangage (Automaton a)
	void AddEtat (Automaton *a, bool final)
	bool IsCompleteAutomaton (Automaton a)
	bool CompleteAutomaton (Automaton *a)
	Automaton BiggerAlphabet (Automaton a, Dict d, int nna) #copy the automaton with a new bigger alphabet
	void Test ()

#dictionnaire numérotant l'alphabet projeté
cdef imagDict (dict d, list A, list A2=[]):
	d1 = {}
	i = 0
	for a in A:
		if d.has_key(a):
			if not d1.has_key(d[a]):
				d1[d[a]] = i
				A2.append(d[a])
				i+=1
	return d1

#dictionnaire numérotant le nouvel alphabet
cdef imagDict2 (dict d, list A, list A2=[]):
	#print "d=%s, A=%s"%(d,A)
	d1 = {}
	i = 0
	for a in A:
		if d.has_key(a):
			for v in d[a]:
				if not d1.has_key(v):
					d1[v] = i
					A2.append(v)
					i+=1
	return d1

cdef Dict getDict (dict d, list A, dict d1=None):
	A = list(A)
	cdef Dict r
	r = NewDict(len(A))
	cdef int i
	if d1 is None:
		d1 = imagDict(d, A)
	#print d1
	for i in range(r.n):
		if d.has_key(A[i]):
			r.e[i] = d1[d[A[i]]]
		else:
			r.e[i] = -1
	return r
	
cdef Dict list_to_Dict (list l):
	cdef Dict d = NewDict(len(l))
	cdef int i
	for i in range(len(l)):
		d.e[i] = l[i]
	return d

cdef InvertDict getDict2 (dict d, list A, dict d1=None):
	A = list(A)
	cdef InvertDict r
	r = NewInvertDict(len(A))
	cdef int i
	if d1 is None:
		d1 = imagDict2(d, A)
	#print d1
	for i in range(r.n):
		if d.has_key(A[i]):
			r.d[i] = NewDict(len(d[A[i]]))
			for j in range(r.d[i].n):
				r.d[i].e[j] = d1[d[A[i]][j]]
		else:
			r.d[i].n = 0;
	return r

#dictionnaire numérotant l'alphabet projeté
cdef imagProductDict (dict d, list A1, list A2, list Av=[]):
	dv = {}
	i = 0
	for a1 in A1:
		for a2 in A2:
			if d.has_key((a1,a2)):
				if not dv.has_key(d[(a1,a2)]):
					dv[d[(a1,a2)]] = i
					Av.append(d[(a1,a2)])
					i+=1
	return dv

cdef Dict getProductDict (dict d, list A1, list A2, dict dv=None, verb=True):
	cdef Dict r
	d1 = {}
	d2 = {}
	cdef int i, n1, n2
	n1 = len(A1)
	for i in range(n1):
		d1[A1[i]] = i
	if verb:
		print d1
	n2 = len(A2)
	for i in range(n2):
		d2[A2[i]] = i
	if verb:
		print d2
	if dv is None:
		dv = imagProductDict(d, A1, A2)
	r = NewDict(n1*n2)
	Keys = d.keys()
	if verb:
		print "Keys=%s"%Keys
	for (a1,a2) in Keys:
		if d1.has_key(a1) and d2.has_key(a2):
			r.e[d1[a1]+d2[a2]*n1] = dv[d[(a1,a2)]]
	return r

def TestAutomaton (a):
	cdef Automaton r
	#d = {}
	#da = {}
	r = getAutomaton(a) #, d, da)
	printAutomaton(r)
	#print d, da
	print a.vertices(), list(a.Alphabet())

def TestProduct (a1, a2, di):
	cdef Automaton a,b,c
	a = getAutomaton(a1)
	b = getAutomaton(a2)
	printAutomaton(a)
	print a1.vertices(), a1.Alphabet()
	printAutomaton(b)
	print a2.vertices(), a2.Alphabet()
	cdef Dict d
	d = getProductDict(di, list(a1.Alphabet()), list(a2.Alphabet()))
	print "dictionnaire du produit :"
	printDict(d)
	c = Product(a, b, d)
	print "résultat :"
	printAutomaton(c)

#def TestDeterminise (a, d, noempty=True, verb=True):
#	cdef Dict di = getDict(d, a.Alphabet())
#	cdef Automaton au = getAutomaton(a)
#	if verb:
#		printDict(di)
#	if verb:
#		printAutomaton(au)
#	cdef Automaton r = Determinise(au, di, noempty, verb)
#	printAutomaton(r)

#def TestDeterminiseEmonde (a, d, noempty=True, verb=True):
#	cdef Dict di = getDict(d, a.Alphabet())
#	cdef Automaton au = getAutomaton(a)
#	if verb:
#		printDict(di)
#	if verb:
#		printAutomaton(au)
#	cdef Automaton r = Determinise(au, di, noempty, verb)
#	print "Avant émondation :"
#	printAutomaton(r)
#	cdef Automaton r2 = emonde_inf(r)
#	print "Après émondation :"
#	printAutomaton(r2)
#	if equalsAutomaton(r, r2):
#		print "equals !"
#	else:
#		print "differents !"

def TestEmonde (a, noempty=True, verb=True):
	cdef Automaton au = getAutomaton(a)
	if verb:
		print "Avant émondation :"
		printAutomaton(au)
	cdef Automaton r = emonde_inf(au, verb)
	if verb:
		print "Après émondation :"
		printAutomaton(r)
	if equalsAutomaton(r, au):
		print "equals !"
	else:
		print "differents !"
	return AutomatonGet(r)
	

cdef Automaton getAutomaton (a, initial=None, F=None, A=None):
	sig_on()
	d = {}
	da = {}
	if F is None:
		if not hasattr(a, 'F'):
			F = a.vertices()
		else:
			F = a.F
	cdef Automaton r
	if A is None:
		A = list(a.Alphabet())
	V = list(a.vertices())
	cdef int n = len(V)
	cdef int na = len(A)
	r = NewAutomaton(n, na)
	init(&r)
	for i in range(na):
		da[A[i]] = i
	for i in range(n):
		r.e[i].final = 0
		d[V[i]] = i
	for v in F:
		if not d.has_key(v):
			FreeAutomaton(&r)
			r = NewAutomaton(0,0)
			print "Error : Incorrect set of final states."
			return r;
		r.e[d[v]].final = 1
	if initial is None:
		if not hasattr(a, 'I'):
			I = []
			#raise ValueError("I must be defined !")
		else:
			I = list(a.I)
		if len(I) > 1:
			print "L'automate doit être déterministe ! (I=%s)"%a.I
		if len(I) >= 1:
			r.i = d[I[0]]
		else:
			r.i = -1
	else:
		r.i = initial
	for e,f,l in a.edges():
		r.e[d[e]].f[da[l]] = d[f]
	sig_off()
	return r

cdef AutomatonGet (Automaton a, A=None):
	from sage.combinat.words.automata import Automaton
	r = Automaton(multiedges=True, loops=True)
	cdef int i,j
	r.F = []
	if A is None:
		A = [i for i in range(a.na)]
	for i in range(a.n):
		for j in range(a.na):
			if a.e[i].f[j] != -1:
				r.add_edge((i, a.e[i].f[j], A[j]))
		if a.e[i].final:
			r.F.append(i)
	r.I = [a.i]
	return r

#cdef initFA (Automaton *a):
#	*a = NewAutomaton(1,1)

cdef Bool (int x):
	if x:
		return True
	return False

cdef class NFastAutomaton:

#   cdef NAutomaton* a
#	cdef list A
	
	def __cinit__ (self):
		#print "cinit"
		self.a = <NAutomaton *>malloc(sizeof(NAutomaton))
		#initialise
		self.a[0].e = NULL
		self.a[0].n = 0
		self.a[0].na = 0
		self.A = []
	
	def __init__(self, a, i=None, F=None, A=None):
		#print "init"
		if a is None:
			return
		else:
			raise ValueError("Cannot construct directly a NFastAutomaton for the moment.")
	
	def __dealloc__ (self):
		#print "free"
		FreeNAutomaton(self.a)
		free(self.a)
	
	def __repr__ (self):
		return "NFastAutomaton with %d states and an alphabet of %d letters"%(self.a.n, self.a.na)
	
	def determinise (self, puits=False):
		cdef Automaton a
		sig_on()
		r = FastAutomaton(None)
		a = DeterminiseN(self.a[0], puits)
		r.a[0] = a
		r.A = self.A
		sig_off()
		return r
	
	def plot (self, int sx=10, int sy=8):
		sig_on()
		cdef char** ll
		ll = <char **>malloc(sizeof(char*)*self.a.na)
		cdef int i
		strA = []
		for i in range(self.a.na):
			strA.append(str(self.A[i]))
			ll[i] = strA[i]
		NplotTikZ(self.a[0], ll, "Automaton", sx, sy)
		free(ll);
		sig_off()

cdef class FastAutomaton:
	
#	cdef Automaton* a
#	cdef list A

	def __cinit__ (self):
		#print "cinit"
		self.a = <Automaton *>malloc(sizeof(Automaton))
		#initialise
		self.a[0].e = NULL
		self.a[0].n = 0
		self.a[0].na = 0
		self.a[0].i = -1
		self.A = []
	
	def __init__(self, a, i=None, final_states=None, A=None):
		#print "init"
		if a is None:
			return
		from sage.graphs.digraph import DiGraph
		if isinstance(a, list):
			a = DiGraph(a, multiedges=True, loops=True)
		if isinstance(a, DiGraph):
			if A is None:
				if hasattr(a, 'A'):
					A = list(a.A)
				else:
					A = list(set(a.edge_labels()))
			self.A = A
			sig_on()
			self.a[0] = getAutomaton(a, initial=i, F=final_states, A=self.A)
			sig_off()
		else:
			raise ValueError("Cannot convert the input to FastAutomaton.")
	
	def __dealloc__ (self):
		#print "free"
		sig_on()
		FreeAutomaton(self.a)
		free(self.a)
		sig_off()
	
	def __repr__ (self):
		return "FastAutomaton with %d states and an alphabet of %d letters"%(self.a.n, self.a.na)
	
	def __hash__ (self):
		h = 3;
		for a in self.A:
			h += hash(a)
			h = (h*19) % 1000000007
		h += hashAutomaton(self.a[0])
		#print "hash=%s"%h
		return h
	
	#######
	#def __richcmp__(left, right, int op):
	#######
	
	def __cmp__ (self, FastAutomaton other):
		#if type(other) != FastAutomaton:
		#	return 1
		cdef int r
		sig_on()
		r = equalsAutomaton(self.a[0], other.a[0]) and self.A == other.A
		sig_off()
		#print "cmp %s"%r
		return (r == 0)

	def Automaton (self):
		return AutomatonGet(self.a[0], self.A)
	
	#give a FastAutomaton recognizing the full language over A.
	def full (self, list A):
		cdef Automaton a
		sig_on()
		r = FastAutomaton(None)
		a = NewAutomaton(1, len(A))
		for i in range(len(A)):
			a.e[0].f[i] = 0
		a.e[0].final = True
		a.i = 0
		r.a[0] = a
		r.A = A
		sig_off()
		return r

	
	def plot (self, int sx=10, int sy=8):
		sig_on()
		cdef char** ll
		ll = <char **>malloc(sizeof(char*)*self.a.na)
		cdef int i
		strA = []
		for i in range(self.a.na):
			strA.append(str(self.A[i]))
			ll[i] = strA[i]
		plotTikZ(self.a[0], ll, "Automaton", sx, sy)
		free(ll);
		sig_off()
		#self.Automaton().plot2()
	
	def Alphabet (self):
		return self.A
	
	def setAlphabet (self, list A):
		self.A = A
	
	def initial_state (self):
		return self.a.i
	
	def set_initial_state (self, int i):
		self.a.i = i
	
	def final_states (self):
		l = []
		for i in range(self.a.n):
			if self.a.e[i].final:
				l.append(i)
		return l
	
	def states (self):
		return range(self.a.n)
	
	def set_final_states (self, lf):
		cdef int f
		for f in range(self.a.n):
			self.a.e[f].final = 0
		for f in lf:
			if f < 0 or f >= self.a.n:
				raise ValueError("%d is not a state !"%f)
			self.a.e[f].final = 1
	
	def is_final (self, int e):
		if e >= 0 and e < self.a.n:
			return Bool(self.a.e[e].final)
		else:
			return False
	
	def set_final_state (self, int e, final=True):
		self.a.e[e].final = final
	
	def succ (self, int i, int j):
		if i<0 or i>=self.a.n or j<0 or j>=self.a.na:
			return -1
		return self.a.e[i].f[j]
	
	#suit le chemin étiqueté par l
	def path (self, list l, i=None):
		if i is None:
			i = self.a.i
		for j in l:
			i = self.succ(i, j)
		return i
	
	def set_succ (self, int i, int j, int k):
		if i<0 or i>=self.a.n or j<0 or j>=self.a.na:
			raise ValueError("set_succ(%s, %s) : index out of bounds !"%(i,j))
		self.a.e[i].f[j] = k
	
	def emonde_inf (self, verb=False):
		sig_on()
		cdef Automaton a
		r = FastAutomaton(None)
		a = emonde_inf(self.a[0], verb)
		r.a[0] = a
		r.A = self.A
		sig_off()
		return r
	
	def emonde_i (self, verb=False):
		sig_on()
		cdef Automaton a
		r = FastAutomaton(None)
		a = emondeI(self.a[0], verb)
		r.a[0] = a
		r.A = self.A
		sig_off()
		return r
	
	def emonde (self, verb=False):
		sig_on()
		cdef Automaton a
		r = FastAutomaton(None)
		a = emonde(self.a[0], verb)
		r.a[0] = a
		r.A = self.A
		sig_off()
		return r
	
#	def equals (self, FastAutomaton b):
#		return Bool(equalsAutomaton(self.a[0], b.a[0]))
	
	#assume that the dictionnary d is injective !!!
	def product (self, FastAutomaton b, dict d=None, verb=False):
		if d is None:
			d = {}
			for la in self.A:
				for lb in b.A:
					d[(la,lb)] = (la,lb)
		sig_on()
		cdef Automaton a
		cdef Dict dC
		r = FastAutomaton(None)
		Av = []
		dv = imagProductDict(d, self.A, b.A, Av=Av)
		if verb:
			print "Av=%s"%Av
			print "dv=%s"%dv
		dC = getProductDict(d, self.A, b.A, dv=dv, verb=verb)
		if verb:
			print "dC="
			printDict(dC)
		a = Product(self.a[0], b.a[0], dC)
		FreeDict(dC)
		r.a[0] = a
		r.A = Av
		sig_off()
		return r
	
	def intersection (self, FastAutomaton a, verb=False):
		d = {}
		for l in self.A:
			if l in a.A:
				d[(l,l)] = l
		if verb:
			print "d=%s"%d
		return self.product(a, d, verb=verb)
	
	#determine if the automaton is complete (i.e. with his hole state)
	def is_complete (self):
		sig_on()
		res = IsCompleteAutomaton(self.a[0])
		sig_off()
		return Bool(res)
	
	#give a complete automaton (i.e. with his hole state)
	def complete (self):
		sig_on()
		res = CompleteAutomaton(self.a)
		sig_off()
		return Bool(res)
	
	#give the smallest language stable by prefix containing the language of self
	#i.e. every states begin finals
	def prefix_closure (self):
		cdef int i
		cdef Automaton a
		r = FastAutomaton(None)
		sig_on()
		a = emonde(self.a[0], False)
		sig_off()
		r.a[0] = a
		r.A = self.A
		for i in range(a.n):
			a.e[i].final = True
		return r
	
	def union (self, FastAutomaton a, verb=False):
		#complete the automata
		sig_on()
		CompleteAutomaton(self.a)
		CompleteAutomaton(a.a)
		sig_off()
		
		#make the product
		d = {}
		for l in self.A:
			if l in a.A:
				d[(l,l)] = l
		
		cdef Automaton ap
		cdef Dict dC
		r = FastAutomaton(None)
		Av = []
		sig_on()
		dv = imagProductDict(d, self.A, a.A, Av=Av)
		sig_off()
		if verb:
			print "Av=%s"%Av
			print "dv=%s"%dv
		sig_on()
		dC = getProductDict(d, self.A, a.A, dv=dv, verb=verb)
		sig_off()
		if verb:
			print "dC="
			printDict(dC)
		sig_on()
		ap = Product(self.a[0], a.a[0], dC)
		FreeDict(dC)
		sig_off()
		
		#set final states
		cdef int i,j
		cdef n1 = self.a.n
		for i in range(n1):
			for j in range(a.a.n):
				ap.e[i+n1*j].final = self.a.e[i].final or a.a.e[j].final
		
		r.a[0] = ap
		r.A = Av
		#if verb:
		#	print r
		#r = r.emonde()
		#r = r.minimise()
		return r.emonde().minimise()
	
	#split the automaton with respect to a
	def split (self, FastAutomaton a, verb=False):
		#complete the automaton a
		sig_on()
		CompleteAutomaton(a.a)
		sig_off()
		
		#make the product
		d = {}
		for l in self.A:
			if l in a.A:
				d[(l,l)] = l
		
		cdef Automaton ap
		cdef Dict dC
		r = FastAutomaton(None)
		r2 = FastAutomaton(None)
		Av = []
		sig_on()
		dv = imagProductDict(d, self.A, a.A, Av=Av)
		sig_off()
		if verb:
			print "Av=%s"%Av
			print "dv=%s"%dv
		sig_on()
		dC = getProductDict(d, self.A, a.A, dv=dv, verb=verb)
		sig_off()
		if verb:
			print "dC="
			printDict(dC)
		sig_on()
		ap = Product(self.a[0], a.a[0], dC)
		FreeDict(dC)
		sig_off()
		
		#set final states for the intersection
		cdef int i,j
		cdef n1 = self.a.n
		for i in range(n1):
			for j in range(a.a.n):
				ap.e[i+n1*j].final = self.a.e[i].final and a.a.e[j].final
		
		#complementary of a in self
		cdef Automaton ap2
		ap2 = CopyAutomaton(ap, ap.n, ap.na)
		#set final states
		for i in range(n1):
			for j in range(a.a.n):
				ap2.e[i+n1*j].final = self.a.e[i].final and not a.a.e[j].final
		
		r.a[0] = ap
		r.A = Av
		r2.a[0] = ap2
		r2.A = Av
		return [r.emonde().minimise(), r2.emonde().minimise()]
	
	#return the automaton recognizing the langage shifted
	def shift (self, verb=False):
		#détermine la liste des successeurs de l'état initial
		cdef int i
		cdef int e
		l = set()
		for j in range(self.a.na):
			e = self.a.e[self.a.i].f[j]
			if e != -1:
				l.add(e)
		l = list(l)
		if verb:
			print "états à considérer : %s"%l
		#calcule l'union des automates
		a = self.copy()
		a.set_initial_state(l[0])
		a = a.emonde().minimise()
		for i in range(1, len(l)):
			b = self.copy()
			b.set_initial_state(l[i])
			b = b.emonde().minimise()
			a = a.union(b)
		return a
	
	#modify the automaton to recognize the langage shifted by a (letter given by its index)
	def shift1OP (self, int a, verb=False):
		if self.a.i != -1:
			self.a.i = self.a.e[self.a.i].f[a]
	
	def unshift1 (self, int a):
		r = FastAutomaton(None)
		sig_on()
		cdef Automaton aut
		aut = CopyAutomaton(self.a[0], self.a.n+1, self.a.na)
		cdef int i
		cdef int ne = self.a.n
		for i in range(aut.na):
			aut.e[ne].f[i] = -1
		aut.e[ne].f[a] = self.a.i
		aut.i = ne
		r.a[0] = aut
		r.A = self.A
		return r
		
	def determinise_proj (self, dict d, noempty=True, onlyfinals=False, nof=False, verb=False):
		sig_on()
		cdef Automaton a
		cdef Dict dC
		r = FastAutomaton(None)
		A2 = []
		d1 = imagDict(d, self.A, A2=A2)
		if verb:
			print "d1=%s, A2=%s"%(d1,A2)
		dC = getDict(d, self.A, d1=d1)
		a = Determinise (self.a[0], dC, noempty, onlyfinals, nof, verb)
		FreeDict(dC)
		#FreeAutomaton(self.a[0])
		r.a[0] = a
		r.A = A2
		sig_off()
		return r
	
	#change les lettres selon d, en dupliquant les arêtes si nécessaire
	#the result is assumed deterministic !!!
	def duplicate (self, dict d, verb=False):
		sig_on()
		cdef Automaton a
		cdef InvertDict dC
		r = FastAutomaton(None)
		A2 = []
		d1 = imagDict2(d, self.A, A2=A2)
		if verb:
			print "d1=%s, A2=%s"%(d1,A2)
		dC = getDict2(d, self.A, d1=d1)
		if verb:
			printInvertDict(dC)
		a = Duplicate (self.a[0], dC, len(A2), verb)
		if verb:
			print "end..."
		FreeInvertDict(dC)
		r.a[0] = a
		r.A = A2
		sig_off()
		return r
	
	#change les lettres
	#le dictionnaire est supposé bijectif de A dans le nouvel alphabet
	#opération sur place !
	def relabel (self, dict d):
		self.A = [d[c] for c in self.A] 
	   
	#permute les lettres
	#A = liste des lettres dans le nouvel ordre (il peut y en avoir moins)
	def permut (self, list A, verb=False):
		if verb:
			print "A=%s"%A
		sig_on()
		cdef Automaton a
		r = FastAutomaton(None)
		cdef int *l = <int*>malloc(sizeof(int)*len(A))
		cdef int i
		for i in range(self.a.na):
			l[i] = -1
		d = {}
		for i,c in enumerate(self.A):
			d[c] = i
		for i,c in enumerate(A):
			if d.has_key(c):
				l[i] = d[c] #l gives the old index from the new one
		if verb:
			str = "l=["
			for i in range(len(A)):
				str+=" %s"%l[i]
			str+=" ]"
			print str 
		a = Permut (self.a[0], l, len(A), verb)
		free(l)
		r.a[0] = a
		r.A = A
		sig_off()
		return r
	
	#permute les lettres SUR PLACE
	#A = liste des lettres dans le nouvel ordre (il peut y en avoir moins)
	def permut_op (self, list A, verb=False):
		if verb:
			print "A=%s"%A
		sig_on()
		cdef int *l = <int*>malloc(sizeof(int)*len(A))
		cdef int i
		for i in range(self.a.na):
			l[i] = -1
		d = {}
		for i,c in enumerate(self.A):
			d[c] = i
		for i,c in enumerate(A):
			if d.has_key(c):
				l[i] = d[c] #l gives the old index from the new one
		if verb:
			str = "l=["
			for i in range(len(A)):
				str+=" %s"%l[i]
			str+=" ]"
			print str 
		PermutOP (self.a[0], l, len(A), verb)
		free(l)
		self.A = A
		sig_off()
	
	def transpose_det (self):
		sig_on()
		r = FastAutomaton(None)
		r.a[0] = TransposeDet(self.a[0])
		r.A = self.A
		sig_off()
		return r
	
	def transpose (self):
		sig_on()
		r = NFastAutomaton(None)
		r.a[0] = Transpose(self.a[0])
		r.A = self.A
		sig_off()
		return r
	
	def strongly_connected_components (self):
		sig_on()
		cdef int* l = <int*>malloc(sizeof(int)*self.a.n)
		cdef int ncc = StronglyConnectedComponents(self.a[0], l)
		#inverse la liste
		l2 = {}
		cdef int i
		for i in range(self.a.n):
			if not l2.has_key(l[i]):
				l2[l[i]] = []
			l2[l[i]].append(i)
		free(l)
		sig_off()
		return l2.values()
	
	def sub_automaton (self, l, verb=False):
		sig_on()
		r = FastAutomaton(None)
		r.a[0] = SubAutomaton(self.a[0], list_to_Dict(l), verb)
		r.A = self.A
		sig_off()
		return r
	
	def minimise (self, verb=False):
		sig_on()
		r = FastAutomaton(None)
		r.a[0] = Minimise(self.a[0], verb)
		r.A = self.A
		sig_off()
		return r
	
	def adjacency_matrix (self, sparse=None):
		if sparse is None:
			if self.a.n <= 128:
				sparse=False
			else:
				sparse=True
		
		d = {}
		cdef int i,j,f
		for i in range(self.a.n):
			for j in range(self.a.na):
				f = self.a.e[i].f[j]
				if f != -1:
					if d.has_key((i,f)):
						d[(i,f)] += 1
					else:
						d[(i,f)] = 1
		from sage.matrix.constructor import matrix
		from sage.rings.integer_ring import IntegerRing
		return matrix(IntegerRing(), self.a.n, self.a.n, d, sparse=sparse)
	
	def delete_vertex (self, int i):
		sig_on()
		r = FastAutomaton(None)
		r.a[0] = DeleteVertex(self.a[0], i)
		r.A = self.A
		sig_off()
		return r
	
	def delete_vertex_op (self, int i):
		sig_on()
		DeleteVertexOP(self.a, i)
		sig_off()
	
	def spectral_radius (self, only_non_trivial=False, verb=False):
		sig_on()
		a = self.minimise()
		if verb:
			print "Automate minimal : %s"%a
		l = a.strongly_connected_components()
		if verb:
			print "%s composantes fortement connexes."%len(l)
		r = 0 #valeur propre maximale trouvée
		for c in l:
			if not only_non_trivial or len(c) > 1:
				if verb:
					print "composante ayant %s états..."%len(c)
				b = a.sub_automaton(c)
				m = b.adjacency_matrix()
				cp = m.charpoly()
				fs = cp.factor()
				if verb:
					print fs
				for f in fs:
					if verb:
						print f
					from sage.functions.other import real_part
					from sage.rings.qqbar import AlgebraicRealField
					r = max([ro[0] for ro in f[0].roots(ring=AlgebraicRealField())]+[r])
		sig_off()
		return r
		
	def test (self):
		Test()
	
	def copy (self):
		sig_on()
		r = FastAutomaton(None)
		r.a[0] = CopyAutomaton(self.a[0], self.a.n, self.a.na)
		r.A = self.A
		sig_off()
		return r
	
	def has_empty_langage (self):
		sig_on()
		res = emptyLangage(self.a[0])
		sig_off()
		return Bool(res)
	
	def equals_langages (self, FastAutomaton a2, minimized=False):
		sig_on()
		cdef Dict d = NewDict(self.a.na)
		cdef int i,j
		for i in range(self.a.na):
			for j in range(a2.a.na):
				if self.A[i] == a2.A[j]:
					d.e[i] = j
					break
		res = equalsLangages(self.a, a2.a, d, minimized)
		sig_off()
		return Bool(res)

	def add_state (self, bool final):
		sig_on()
		AddEtat(self.a, final)
		sig_off()

	def n_states (self):
		return self.a.n

	def bigger_alphabet (self, list nA):
		cdef Dict d
		d = NewDict(self.a.na)
		for i in range(self.a.na):
			d.e[i] = nA.index(self.A[i])
		r = FastAutomaton(None)
		sig_on()
		r.a[0] = BiggerAlphabet(self.a[0], d, len(nA))
		sig_off()
		r.A = nA
		return r
	
	def complementaryOP (self):
		self.complete()
		cdef i
		for i in range(self.a.n):
			self.a.e[i].final = not self.a.e[i].final
			
	def included (self, FastAutomaton a, bool verb=False, step=None):
		d = {}
		for l in self.A:
			if l in a.A:
				d[(l,l)] = l
		if verb:
			print "d=%s"%d
		a.complete()
		cdef FastAutomaton p = self.product(a, d, verb=verb)
		
		#set final states
		cdef int i,j
		cdef n1 = self.a.n
		for i in range(n1):
			for j in range(a.a.n):
				p.a.e[i+n1*j].final = self.a.e[i].final and not a.a.e[j].final
		
		if step == 1:
			return p;
		
		return p.has_empty_langage()
	
	#donne un automate reconnaissant w(w^(-1)L) où L est le langage de a partant de e
	def piece (self, w, e=None):
		cdef int* l = <int*>malloc(sizeof(int)*self.a.n)
		cdef int i
		if type(w) != list:
			w = [int(w)]
		for i in range(len(w)):
			l[i] = w[i]
		if e is None:
			e = self.a.i
		r = FastAutomaton(None)
		sig_on()
		r.a[0] = PieceAutomaton(self.a[0], l, len(w), e)
		sig_off()
		free(l)
		r.A = self.A
		return r

	#tell if the language of the automaton is empty
	#(this function is not very efficient)
	def is_empty (self, ext=True):
		if ext:
			return self.emonde().emonde_inf().n_states() == 0
		else:
			return self.emonde().n_states() == 0
	
	#determine if the languages intersect
	def intersect (self, FastAutomaton b, ext=True):
		return not self.intersection(b).is_empty(ext)

