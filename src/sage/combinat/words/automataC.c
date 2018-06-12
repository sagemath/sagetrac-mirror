#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Automaton.h"
#include "automataC.h"

typedef Automate Automaton;

//test if the dot command is installed in the system
bool DotExists ()
{
    bool res = true;
    system("command -v dot > testdot");
    FILE *f = fopen("testdot","r");
    fseek(f, 0, SEEK_END);
    if (ftell(f) == 0)
        res = false;
    fclose(f);
    return res;
}

Dict NewDict (int n)
{
	Dict r;
	r.n = n;
	if (n == 0)
		return r;
	r.e = (int *)malloc(sizeof(int)*n);
	if (!r.e)
	{
		printf("Out of memory !\n");
		exit(15);
	}
	int i;
	for (i=0;i<n;i++)
	{
		r.e[i] = -1;
	}
	return r;
}

void FreeDict (Dict *d)
{
	if (d->n == 0)
		return;
	free(d->e);
	d->n = 0;
	d->e = NULL;
}

void printDict (Dict d)
{
	int i;
	printf("[ ");
	for (i=0;i<d.n;i++)
	{
		printf("%d ", d.e[i]);
	}
	printf("]\n");
}

//Add an element to the dictionnary (even if already in)
void dictAdd (Dict *d, int e)
{
	d->n++;
	if (d->n == 1)
		d->e = (int *)malloc(sizeof(int));
	else
		d->e = (int *)realloc(d->e, sizeof(int)*d->n);
	if (!d->e)
	{
		printf("Out of memory !");
		exit(1);
	}
	d->e[d->n-1] = e;
}

#define DISP_MEMORY	false

Automaton NewAutomaton (int n, int na)
{
#if DISP_MEMORY
	printf("New %d...\n", n);
#endif
	Automaton a;
	a.n = n;
	a.na = na;
	a.i = -1;
	if (n == 0)
	{
		a.e = NULL;
		return a;
	}
	a.e = (Etat *)malloc(sizeof(Etat)*n);

#if DISP_MEMORY
	printf("new aut %ld\n", a.e);
#endif
	
	//a.nalloc = n;
	if (!a.e)
	{
		printf("Out of memory !");
		exit(6);
	}
	int i,j;
	for (i=0;i<n;i++)
	{
		a.e[i].f = (int *)malloc(sizeof(int)*na);
		if (!a.e[i].f)
		{
			printf("Out of memory !");
			exit(7);
		}
		for (j=0;j<na;j++)
		{
			a.e[i].f[j] = -1;
		}
	}
	return a;
}

void FreeAutomaton(Automaton *a)
{
#if DISP_MEMORY
	printf("Free %d ...\n", a->n);
#endif
	//if (a.nalloc == 0 || a.n==0)
	if (a->n == 0 || a->e == NULL)
		return;
	int i;
	for (i=0;i<a->n;i++)
	{
		free(a->e[i].f);
	}
	//printf("free %ld\n", a->e);
	free(a->e);
	//a.nalloc = 0;
	a->n = 0;
}

void FreeAutomates (Automate* a, int n)
{
	int j;
	for (j=0;j<n;j++)
	{
		FreeAutomaton(&a[j]);
	}
}

int hashAutomaton (Automaton a)
{
	int h = 3;
	h += a.n + 1009*a.na + 1000003*a.i;
	int i, j;
	for (i=0;i<a.n;i++)
	{
		for (j=0;j<a.na;j++)
		{
			h += a.e[i].f[j];
			h = (h*1009) % 1000000007;
		}
		h += (a.e[i].final == 0);
		h = (h*1009) % 1000000007;
	}
	return h;
}

NAutomaton NewNAutomaton(int n, int na)
{
	NAutomaton a;
	a.n = n;
	a.na = na;
	if (n == 0)
	{
		a.e = NULL;
		return a;
	}
	a.e = (NEtat *)malloc(sizeof(NEtat)*n);
	if (!a.e)
	{
		printf("Out of memory !");
		exit(6);
	}
	int i;
	for (i=0;i<n;i++)
	{
		a.e[i].a = NULL;
		a.e[i].n = 0;
	}
	return a;
}

void ReallocNAutomaton (NAutomaton *a, int n)
{
	if (a->n)
		a->e = (NEtat*)realloc(a->e, sizeof(NEtat)*n);
	else
		a->e = (NEtat*)malloc(sizeof(NEtat)*n);
	if (a->n < n)
	{
		int i;
		for (i=a->n;i<n;i++)
		{
			a->e[i].a = NULL;
			a->e[i].n = 0;
			a->e[i].initial = false;
			a->e[i].final = false;
		}
	}
	a->n = n;
}

void FreeNAutomaton (NAutomaton *a)
{
	if (a->n == 0)
		return;
	int i;
	for (i=0;i<a->n;i++)
	{
		if (a->e[i].n > 0)
		{
			free(a->e[i].a);
		}
	}
	free(a->e);
	a->n = 0;
}

//add an edge on the NFastAutomaton
void AddEdgeN (NAutomaton *a, int e, int f, int l)
{
	a->e[e].n++;
	if (a->e[e].n > 1)
		a->e[e].a = (Arete *)realloc(a->e[e].a, sizeof(Arete)*a->e[e].n);
	else
		a->e[e].a = (Arete *)malloc(sizeof(Arete));
	a->e[e].a[a->e[e].n - 1].e = f;
	a->e[e].a[a->e[e].n - 1].l = l;
}

//Add a path between states e and f of NFastAutomaton a
void AddPathN (NAutomaton *a, int e, int f, int *l, int len, bool verb)
{
	int n = a->n;
	if (len > 1)
	{
		//Add new states
		if (verb)
			printf("realloc to size %d\n", a->n+len-1);
		ReallocNAutomaton(a, a->n+len-1);
		int i;
		if (verb)
			printf("add edge %d --%d--> %d\n", e, l[0], n+i);
		AddEdgeN(a, 0, n, l[0]);
		for (i=1;i<len-1;i++)
		{
			if (verb)
				printf("add edge %d --%d--> %d\n", n+i-1, l[i], n+i);
			AddEdgeN(a, n+i-1, n+i, l[i]);
		}
		if (verb)
			printf("add edge %d --%d--> %d\n", n+len-2, l[len-1], f);
		AddEdgeN(a, n+len-2, f, l[len-1]);
	}else
	{
		if (len == 1)
		{
			AddEdgeN(a, e, f, l[0]);
		}
	}
}

void ReallocAutomaton (Automaton *a, int n, bool init)
{
	if (a->n > n)
	{
		//free states to suppress
		int i;
		for (i=n;i<a->n;i++)
		{
			free(a->e[i].f);
		}
	}
	if (a->n)
		a->e = (Etat*)realloc(a->e, sizeof(Etat)*n);
	else
		a->e = (Etat*)malloc(sizeof(Etat)*n);
	if (a->n < n)
	{
		int i;
		for (i=a->n;i<n;i++)
		{
			a->e[i].f = (int *)malloc(sizeof(int)*a->na);
			if (init)
			{
				int j;
				for (j=0;j<a->na;j++)
				{
					a->e[i].f[j] = -1;
				}
			}
		}
	}
	a->n = n;
}

Automaton CopyAutomaton(Automaton a, int nalloc, int naalloc)
{
	//a.n, a.na
	Automaton r = NewAutomaton(nalloc, naalloc);
	int i,j;
	for (i=0;i<a.n;i++)
	{
		r.e[i].final = a.e[i].final;
		for (j=0;j<a.na;j++)
		{
			r.e[i].f[j] = a.e[i].f[j];
		}
	}
	r.i = a.i;
	return r;
}

//give an automaton recognizing w(w^(-1)L) where L is the language of a starting from state e
Automaton PieceAutomaton (Automaton a, int *w, int n, int e)
{
	int i, j, f;
	Automaton r = NewAutomaton(a.n+n, a.na);
	//put the states for the word w
	for (i=0;i<n;i++)
	{
		for (j=0;j<a.na;j++)
		{
			r.e[i].f[j] = -1;
		}
		r.e[i].f[w[i]] = i+1;
		r.e[i].final = false;
		f = a.e[e].f[w[i]];
		e = f;
	}
	if (n > 0)
		r.e[n-1].f[w[n-1]] = e+n;
	//test if empty
	if (n > 0 && f == -1)
	{
		FreeAutomaton(&r);
		return NewAutomaton(0, a.na);
	}
	//put the states of the automaton a
	for (i=0;i<a.n;i++)
	{
		for (j=0;j<a.na;j++)
		{
			if (a.e[i].f[j] == -1)
				r.e[i+n].f[j] = -1;
			else
				r.e[i+n].f[j] = a.e[i].f[j]+n;
		}
		r.e[i+n].final = a.e[i].final;
	}
	if (n > 0)
		r.i = 0;
	else
		r.i = a.i;
	return r;
}

void init (Automaton *a)
{
	int i,j;
	a->i = -1;
	for (i=0;i<a->n;i++)
	{
		a->e[i].final = false;
		for (j=0;j<a->na;j++)
		{
			a->e[i].f[j] = -1;
		}
	}
}

void printAutomaton (Automaton a)
{
	printf("Automaton with %d states, %d letters.\n", a.n, a.na);
	int i, j;
	for (i=0;i<a.n;i++)
	{
		for (j=0;j<a.na;j++)
		{
			if (a.e[i].f[j] != -1)
			{
				printf("%d --%d--> %d\n", i, j, a.e[i].f[j]);
			}
		}
	}
	printf("initial State %d.\n", a.i);
}

void plotDot (const char *file, Automaton a, const char **labels, const char *graph_name, double sx, double sy, const char **vlabels, bool html, bool verb, bool run_dot)
{
	char tamp[1024];
	FILE *f = fopen(file, "w");
	if (!f)
	{
		printf("Unable to open file %s !\n", file);
		return;
	}
	
	if (verb)
		printf("start...\n");
	fprintf(f, "digraph %s\n{\n"\
	"	node[fontsize=20]"\
	"	edge[fontsize=20, arrowhead = open]"\
	"	rankdir = LR;\n"\
	"	size = \"%lf, %lf\";\n"\
	"	center = 1;\n"\
	"	nodesep = \"0.2\"\n", graph_name, sx, sy);
//	"	ranksep = \"0.4 equally\";\n", graph_name);
//	"	rotate = -90\n"\
//	"	orientation=landscape\n"\
//	"orientation = Landscape\n");
	if (verb)
		printf("write...\n");
	
	fprintf(f, "	\n");
	int i,j;
	for (i=0;i<a.n;i++)
	{
		fprintf(f, "	%d [shape=", i);
		if (a.e[i].final)
			fprintf(f, "doublecircle");
		else
			fprintf(f, "circle");
		fprintf(f, ", style=");
		if (i == a.i)
			fprintf(f, "bold");
		else
			fprintf(f, "solid");
		if (vlabels)
		{
			if (verb)
			{
				printf("ptz : i=%d : %s\n", i, vlabels[i]);
			}
			fprintf(f, ", label=%s", vlabels[i]);
		}
		fprintf(f, ", fontsize=20, margin=0]\n");
	}
	fprintf(f, "	\n");
	for (i=0;i<a.n;i++)
	{
		for (j=0;j<a.na;j++)
		{
			if (a.e[i].f[j] != -1)
			{
			    const char *ptr = labels[j];
			    bool dotted = false;
			    if (labels[j][0] == '.' && labels[j][1] == '.' && labels[j][2] == '.')
			    {
			        ptr += 3;
			        dotted = true;
				}
				if (html)
    				fprintf(f, "	%d -> %d [label=<%s>", i, a.e[i].f[j], ptr);
    			else
    			    fprintf(f, "	%d -> %d [label=\"%s\"", i, a.e[i].f[j], ptr);
    			if (dotted)
    			{
    			    fprintf(f, ", style=\"dotted\"");
    			}
    			fprintf(f, "]\n");
			}
		}
	}
	fprintf(f, "}\n");
	
	fclose(f);
	if (run_dot)
	{
        if (verb)
            printf("draw...\n");
        sprintf(tamp, "dot %s -Gname -Tpng > %s.png", file, file);
        system(tamp);
	}
}

void NplotDot (const char *file, NAutomaton a, const char **labels, const char *graph_name, double sx, double sy, bool run_dot)
{
	bool verb = false;
	char tamp[1024];
	//FILE *f = fopen(temp_dot_file_name, "w");
	FILE *f = fopen(file, "w");
	if (!f)
	{
		printf("Unable to open file a.dot !\n");
		return;
	}
	
	if (verb)
		printf("start...\n");
	fprintf(f, "digraph %s\n{\n"\
	"	node[fontsize=20]"\
	"	edge[fontsize=20, arrowhead = open]"\
	"	rankdir = LR;\n"\
	"	size = \"%lf, %lf\";\n"\
	"	center = 1;\n"\
	"	nodesep = \"0.2\"\n", graph_name, sx, sy);
//	"	ranksep = \"0.4 equally\";\n", graph_name);
//	"	rotate = -90\n"\
//	"	orientation=landscape\n"\
//	"orientation = Landscape\n");
	if (verb)
		printf("write...\n");
	
	fprintf(f, "	\n");
	int i,j;
	for (i=0;i<a.n;i++)
	{
		fprintf(f, "	%d [shape=", i);
		if (a.e[i].final)
			fprintf(f, "doublecircle");
		else
			fprintf(f, "circle");
		fprintf(f, ", style=");
		if (a.e[i].initial)
			fprintf(f, "bold");
		else
			fprintf(f, "solid");
		fprintf(f, ", fontsize=20, margin=0]\n");
	}
	fprintf(f, "	\n");
	for (i=0;i<a.n;i++)
	{
		for (j=0;j<a.e[i].n;j++)
		{
			if (a.e[i].a[j].e != -1)
			{
				if (a.e[i].a[j].l != -1)
					fprintf(f, "	%d -> %d [label=\"%s\"]\n", i, a.e[i].a[j].e, labels[a.e[i].a[j].l]);
				else
					fprintf(f, "	%d -> %d [label=\"Ɛ\"]\n", i, a.e[i].a[j].e);
			}
		}
	}
	fprintf(f, "}\n");
	
	fclose(f);
	if (run_dot)
	{
        if (verb)
            printf("draw...\n");
        sprintf(tamp, "dot %s -Gname -Tpng > %s.png", file, file);
        system(tamp);
	}
	//sprintf(tamp, "scp %s %s &> /dev/null", file, temp_dot_file_name); //copy the file in the choosen place
	//system(tamp);
}

//determine if the automaton is complete (i.e. with his hole state)
bool IsCompleteAutomaton (Automaton a)
{
	int i,j;
	for (i=0;i<a.n;i++)
	{
		for (j=0;j<a.na;j++)
		{
			if (a.e[i].f[j] == -1)
				return false;
		}
	}
	return true;
}

//complete the automaton (i.e. add a hole state if necessary)
//return true iff a state was added
bool CompleteAutomaton (Automaton *a)
{
	int ne = a->n; //new state
	int i,j;
	bool add_etat = false;
	for (i=0;i<ne;i++)
	{
		for (j=0;j<a->na;j++)
		{
			if (a->e[i].f[j] == -1)
			{
				a->e[i].f[j] = ne;
				add_etat = true;
			}
		}
	}
	if (a->i == -1)
	{
		a->i = ne;
		add_etat = true;
	}
	if (!add_etat)
		return false;
	AddEtat(a, false); //add the hole state
	for (j=0;j<a->na;j++)
	{
		a->e[ne].f[j] = ne;
	}
	if (a->i == -1)
		a->i = ne;
	return true;
}

//determine if the two automata are the same (different if permuted states)
bool equalsAutomaton(Automaton a1, Automaton a2)
{
	if (a1.n != a2.n || a1.na != a2.na || a1.i != a2.i)
		return false;
	int i, j;
	for (i=0;i<a1.n;i++)
	{
		for (j=0;j<a1.na;j++)
		{
			if (a1.e[i].f[j] != a2.e[i].f[j])
				return false;
		}
		if (a1.e[i].final != a2.e[i].final)
			return false;
	}
	return true;
}

//used by equalsLanguages
//determine if the languages of states e1 of a1 and state e2 of a2 are the same
bool equalsLanguages_rec (Automaton a1, Automaton a2, Dict a1toa2, Dict a2toa1, int e1, int e2, bool verb)
{
	if ((a1.e[e1].final & 1) != (a2.e[e2].final & 1))
		return false; //one of the states is final but not the other one
	if (a1.e[e1].final & 2 && a2.e[e2].final & 2)
		return true; //state already seen
	//indicate that the state a has been seen
	a1.e[e1].final |= 2;
	a2.e[e2].final |= 2;
	//browse the son of e1 in a1
	int i;
	for (i=0;i<a1.na;i++)
	{
		if (a1.e[e1].f[i] != -1)
		{//this edge exists in a1
			if (a1toa2.e[i] != -1)
			{
				if (a2.e[e2].f[a1toa2.e[i]] == -1)
				{//this edge doesn't correspond to an edge in a2
					if (verb)
						printf("%d -%d-> exists in a1 but %d -%d-> doesn't exists in a2.", e1, i, e2, a1toa2.e[i]);
					return false;
				}
				if (!equalsLanguages_rec(a1, a2, a1toa2, a2toa1, a1.e[e1].f[i], a2.e[e2].f[a1toa2.e[i]], verb))
				{
					return false;
				}
			}else
			{
				if (verb)
					printf("%d -%d-> exists in a1 but %d -%d-> doesn't exists in a2.", e1, i, e2, a1toa2.e[i]);
				return false;
			}
		}else
		{
			if (a1toa2.e[i] != -1)
			{
				if (a2.e[e2].f[a1toa2.e[i]] != -1)
				{				
					if (verb)
						printf("%d -%d-> doesn't exists in a1 but %d -%d-> exists in a2.", e1, i, e2, a1toa2.e[i]);
					return false;
				}
			}
		}
	}
	return true;
}

//determine if the languages of the two automata are the same
//the dictionnary gives the letters of a2 depending of the ones of a1 (-1 if the letter of a1 doesn't correspond to any letter for a2).
//This dictionnary is assumed to be invertible.
//if minimized is true, the automaton a1 and a2 are assumed to be minimal.
bool equalsLanguages(Automaton *a1, Automaton *a2, Dict a1toa2, bool minimized, bool pruned, bool verb)
{
	int i;
	if (!pruned)
	{
		if (verb)
			printf("Emonde...\n");
		//prun the automata
		Automaton a3 = prune(*a1, false);
		FreeAutomaton(a1);
		*a1 = a3;
		a3 = prune(*a2, false);
		FreeAutomaton(a2);
		*a2 = a3;
	}
	if (!minimized)
	{
		if (verb)
			printf("Minimise...\n");
		//minimise les automates
		Automaton a3 = Minimise(*a1, false);
		FreeAutomaton(a1);
		*a1 = a3;
		a3 = Minimise(*a2, false);
		FreeAutomaton(a2);
		*a2 = a3;
	}
	if (verb)
	{
		printf("Automata : ");
		printAutomaton(*a1);
		printAutomaton(*a2);
	}
	//inverse the dictionnary
	Dict a2toa1 = NewDict(a2->na);
	for (i=0;i<a1toa2.n;i++)
	{
		a2toa1.e[a1toa2.e[i]] = i;
	}
	if (verb)
		printDict(a2toa1);
	//
	bool res;
	if (a1->i == -1 || a2->i == -1)
		res = (a1->i == a2->i);
	else
		res = equalsLanguages_rec(*a1, *a2, a1toa2, a2toa1, a1->i, a2->i, verb);
	//put back final states
	for (i=0;i<a1->n;i++)
	{
		a1->e[i].final &= 1;
	}
	for (i=0;i<a2->n;i++)
	{
		a2->e[i].final &= 1;
	}
	return res;
}

//used by emptyLanguage
//determine if the language of the state e is empty
bool emptyLanguage_rec (Automaton a, int e)
{
	if (a.e[e].final)
		return false;
	//indicate that the state has been seen
	a.e[e].final |= 2;
	//browse the sons
	int i;
	for (i=0;i<a.na;i++)
	{
		if (a.e[e].f[i] != -1)
		{
			if (a.e[a.e[e].f[i]].final & 2)
				continue; //this son was already seen
			if (!emptyLanguage_rec(a, a.e[e].f[i]))
				return false;
		}
	}
	return true;
}

//determine if the language of the automaton is empty
bool emptyLanguage (Automaton a)
{
	if (a.i == -1)
		return true;
	bool res = emptyLanguage_rec(a, a.i);
	//put back final states
	int i;
	for (i=0;i<a.n;i++)
	{
		a.e[i].final &= 1;
	}
	return res;
}

bool findWord_rec (Automaton a, int e, int n, Dict *w, bool verb)
{
	if (a.e[e].final)
	{
		if (verb)
			printf("Allocated one word of size %d.\n", n);
		*w = NewDict(n);
		return true;
	}
	//indicate that the state a has been seen
	a.e[e].final |= 2;
	//browse sons
	int i;
	for (i=0;i<a.na;i++)
	{
		if (a.e[e].f[i] != -1)
		{
			if (a.e[a.e[e].f[i]].final & 2)
				continue; //this son was already seen
			if (findWord_rec(a, a.e[e].f[i], n+1, w, verb))
			{
				if (verb)
				{
					printf("w[%d] = %d\n", n, i);
				}
				w->e[n] = i;
				return true;
			}
		}
	}
	return false;
}

//return a word in the language of a
bool findWord (Automaton a, Dict *w, bool verb)
{
	if (a.i == -1)
		return false;
	bool res = findWord_rec(a, a.i, 0, w, verb);
	//put back final states
	int i;
	for (i=0;i<a.n;i++)
	{
		a.e[i].final &= 1;
	}
	return res;
}

//return a shortest word of the language
bool shortestWord (Automaton a, Dict *w, int init, int end, bool verb)
{
	if (a.i == -1)
		return false;
	//Dijkstra's algorithm
	int *d = malloc(sizeof(int)*a.n); //table of distances from initial state
	int *prec = malloc(sizeof(int)*a.n); //table of predecessors
	int *vu = malloc(sizeof(int)*a.n);
	int i, j, imin, f;
	//initialization
	for (i=0;i<a.n;i++)
	{
		d[i] = a.n;
		prec[i] = -1;
		vu[i] = 0;
	}
	d[init] = 0;
	//Dijkstra (not very optimized)
	for (i=0;i<a.n;i++)
	{
		//find the nearest state not seen yet (could be optimized)
		imin = 0;
		while(vu[imin])
			imin++;
		for (j=0;j<a.n;j++)
		{
			if (d[j] < d[imin] && !vu[j])
			{
				imin = j;
			}
		}
		vu[imin] = 1;
		//browse sons
		for (j=0;j<a.na;j++)
		{
			f = a.e[imin].f[j];
			if (f != -1)
			{
				if (d[f] > d[imin]+1)
				{
					d[f] = d[imin]+1;
					prec[f] = imin;
				}
			}
		}
	}
	free(vu);
	//Determine the nearest final state (or end if end != -1)
	f = end;
	if (f == -1)
	{
		for (i=0;i<a.n;i++)
		{
			if (a.e[i].final)
			{
				if (f == -1)
					f = i;
				else if (d[i] < d[f])
				{
					f = i;
				}
			}
		}
	}
	free(d);
	if (f == -1)
	{
		free(prec);
		return false;
	}
	//find the shortest path to this state
	i = 0;
	imin = f;
	while(prec[f] != -1)
	{
		i++;
		f = prec[f];
	}
	//alloc the word
	*w = NewDict(i);
	f = imin;
	while(prec[f] != -1)
	{
		i--;
		//find a letter to go from one state to the other
		for (j=0;j<a.na;j++)
		{
			if (a.e[prec[f]].f[j] == f)
			{
				w->e[i] = j;
				break;
			}
		}
		f = prec[f];
	}
	free(prec);
	return true;
}

//return the shortest words up to each state
bool shortestWords (Automaton a, Dict *w, int init, bool verb)
{
	if (a.i == -1)
		return false;
	//Dijkstra's algorithm
	int *d = malloc(sizeof(int)*a.n); //table of distances from the initial state
	int *prec = malloc(sizeof(int)*a.n); //table of predecessors
	int *vu = malloc(sizeof(int)*a.n);
	int i, j, k, imin, f;
	//initialization
	for (i=0;i<a.n;i++)
	{
		d[i] = a.n;
		prec[i] = -1;
		vu[i] = 0;
	}
	d[init] = 0;
	if (verb)
		printf("Dijkstra...\n");
	//Dijkstra (not very well optimized)
	for (i=0;i<a.n;i++)
	{
		//find the nearest state not already seeen (could be optimized)
		imin = 0;
		while(vu[imin])
			imin++;
		for (j=0;j<a.n;j++)
		{
			if (d[j] < d[imin] && !vu[j])
			{
				imin = j;
			}
		}
		vu[imin] = 1;
		//browse sons
		for (j=0;j<a.na;j++)
		{
			f = a.e[imin].f[j];
			if (f != -1)
			{
				if (d[f] > d[imin]+1)
				{
					d[f] = d[imin]+1;
					prec[f] = imin;
				}
			}
		}
	}
	free(vu);
	if (verb)
	{
		for (i=0;i<a.n;i++)
		{
			printf("prec[%d] = %d, d[%d] = %d\n", i, prec[i], i, d[i]);
		}
	}
	if (verb)
		printf("Fill the list...\n");
	//browse the states
	for (k=0;k<a.n;k++)
	{
		//find the shortest path to this state
		//find the length of the path
		i = 0;
		f = k;
		while(prec[f] != -1)
		{
			i++;
			f = prec[f];
		}
		if (verb)
			printf("alloc %d (size %d)...\n", k, i);
		//alloc the word
		w[k] = NewDict(i);
		f = k;
		while(prec[f] != -1)
		{
			i--;
			//find a letter to go from one state to the next one
			for (j=0;j<a.na;j++)
			{
				if (a.e[prec[f]].f[j] == f)
				{
					w[k].e[i] = j;
					break;
				}
			}
			f = prec[f];
		}
	}
	if (verb)
		printf("end...\n");
	free(prec);
	return true;
}

//check that the word w is recognized by the automaton
bool rec_word (Automaton a, Dict d)
{
	int i;
	int e = a.i;
	for (i=0;i<d.n;i++)
	{
		e = a.e[e].f[d.e[i]];
		if (e == -1)
		{
			return false;
		}
	}
	return a.e[e].final;
}

int contract (int i1, int i2, int n1)
{
	return i1+n1*i2;
}

int geti1 (int c, int n1)
{
	return c%n1;
}

int geti2 (int c, int n1)
{
	return c/n1;
}

void Product_rec(Automaton r, int i1, int i2, Automaton a1, Automaton a2, Dict d)
{
	int i,j;
	int e1, e2;
	int a;
	Etat *current = &r.e[contract(i1, i2, a1.n)];
	int *next;
	current->final = true; //indicate that the state has been visited
	for (i=0;i<a1.na;i++)
	{
		e1 = a1.e[i1].f[i];
		if (e1 < 0)
			continue;
		for (j=0;j<a2.na;j++)
		{
			e2 = a2.e[i2].f[j];
			a = d.e[contract(i,j,a1.na)];
			next = &current->f[a];
			if (a != -1)
			{
				if (e2 < 0)
					*next = -1;
				else
				{
					*next = contract(e1, e2, a1.n);
					if (!r.e[*next].final)
						Product_rec(r, e1, e2, a1, a2, d);
				}
			}
		}
	}
}

Automaton Product(Automaton a1, Automaton a2, Dict d, bool verb)
{
	if (verb)
	{
		printAutomaton(a1);
		printAutomaton(a2);
	}
	//count the number of letters of the final alphabet
	int i, na=0;
	for (i=0;i<d.n;i++)
	{
		if (d.e[i] >= na)
			na = d.e[i]+1;
	}
	Automaton r = NewAutomaton(a1.n*a2.n, na);
	init(&r);
	if (a1.i == -1 || a2.i == -1)
	{
		r.i = -1;
		return r;
	}
	r.i = contract(a1.i, a2.i, a1.n);
	Product_rec(r, a1.i, a2.i, a1, a2, d);
	//put the final states
	for (i=0;i<r.n;i++)
	{
		r.e[i].final = a1.e[geti1(i, a1.n)].final && a2.e[geti2(i, a1.n)].final;
	}
	return r;
}

bool Intersect_rec(int i1, int i2, Automaton a1, Automaton a2, bool *vu, bool verb)
{
	if (verb)
		printf("Intersect_rec %d %d...\n", i1, i2);
	int i,j;
	int e1, e2;
	if (a1.e[i1].final && a2.e[i2].final)
		return true;
	int a = contract(i1, i2, a1.n);
	if (vu[a])
		return false;
	vu[a] = true; //indicate that the state has been visited
	for (i=0;i<a1.na;i++)
	{
		e1 = a1.e[i1].f[i];
		if (e1 < 0)
			continue;
		e2 = a2.e[i2].f[i];
		if (e2 < 0)
			continue;
		if (Intersect_rec(e1, e2, a1, a2, vu, verb))
			return true;
	}
	return false;
}

//determine if the intersection is empty
bool Intersect (Automaton a1, Automaton a2, bool verb)
{
	if (verb)
	{
		printAutomaton(a1);
		printAutomaton(a2);
	}
	if (a1.i == -1 || a2.i == -1)
	{
		if (verb)
			printf("One of the automata doesn't have an initial state !\n");
		return false;
	}
	bool *vu = (bool *)malloc(sizeof(bool)*a1.n*a2.n);
	int i;
	for (i=0;i<a1.n*a2.n;i++)
		vu[i] = false;
	bool res = Intersect_rec(a1.i, a2.i, a1, a2, vu, verb);
	free(vu);
	return res;
}

bool Included_rec(int i1, int i2, Automaton a1, Automaton a2, bool *vu)
{
	int i,j;
	int e1, e2;
	if (a1.e[i1].final && !a2.e[i2].final)
		return false;
	int a = contract(i1, i2, a1.n);
	if (vu[a])
		return true;
	vu[a] = true; //indicate that the state has been visited
	for (i=0;i<a1.na;i++)
	{
		e1 = a1.e[i1].f[i];
		if (e1 < 0)
			continue;
		if (i >= a2.na)
			return false;
		e2 = a2.e[i2].f[i];
		if (e2 < 0)
			return false;
		if (!Included_rec(e1, e2, a1, a2, vu))
			return false;
	}
	return true;
}

//determine if the language of a1 is included in the language of a2
bool Included(Automaton a1, Automaton a2, bool pruned, bool verb)
{
	if (!pruned)
	{
		a1 = prune(a1, verb);
	}
	if (verb)
	{
		printAutomaton(a1);
		printAutomaton(a2);
	}
	if (a1.i == -1)
		return true;
	if (a2.i == -1)
		return emptyLanguage(a1);
	bool *vu = (bool *)malloc(sizeof(bool)*a1.n*a2.n);
	int i;
	for (i=0;i<a1.n*a2.n;i++)
		vu[i] = false;
	bool res = Included_rec(a1.i, a2.i, a1, a2, vu);
	free(vu);
	return res;
}

void AddEtat (Automaton *a, bool final)
{
	a->n++;
	if (a->n == 1)
		a->e = (Etat *)malloc(sizeof(Etat));
	else
		a->e = (Etat *)realloc(a->e, sizeof(Etat)*a->n);
	if (!a->e)
	{
		printf("Out of memory !");
		exit(2);
	}
	a->e[a->n-1].f = (int *)malloc(sizeof(int)*a->na);
	if (!a->e[a->n-1].f)
	{
		printf("Out of memory !");
		exit(3);
	}
	int i;
	for (i=0;i<a->na;i++)
	{
		a->e[a->n-1].f[i] = -1;
	}
	a->e[a->n-1].final = final;
}

Etats NewEtats (int n)
{
	Etats e;
	e.n = n;
	e.e = (int *)malloc(sizeof(int)*n);
	if (!e.e)
	{
		printf("Out of memory !");
		exit(7);
	}
return e;
}

void FreeEtats (Etats e)
{
	free(e.e);
}

void initEtats (Etats e)
{
	int i;
	for (i=0;i<e.n;i++)
	{
		e.e[i] = 0;
	}
}

void printEtats (Etats e)
{
	int i;
	printf("[ ");
	for (i=0;i<e.n;i++)
	{
		printf("%d ", e.e[i]);
		fflush(stdout);
	}
	printf("]\n");
}

bool equals (Etats e1, Etats e2)
{
	if (e1.n != e2.n)
		return false;
	int i;
	for (i=0;i<e1.n;i++)
	{
		if (e1.e[i] != e2.e[i])
			return false;
	}
	return true;
}

Etats copyEtats(Etats e)
{
	Etats r = NewEtats(e.n);
	int i;
	for (i=0;i<e.n;i++)
	{
		r.e[i] = e.e[i];
	}
	return r;
}

void printListEtats (ListEtats l)
{
	int i;
	for (i=0;i<l.n;i++)
	{
		printf("%d : ", i);
		printEtats(l.e[i]);
	}
}

//add an element if not already in the list
bool AddEl (ListEtats *l, Etats e, int* res)
{
	int i;
	for (i=0;i<l->n;i++)
	{
		if (equals(l->e[i], e))
		{
			if (res)
				*res = i;
			return false;
		}
	}
	//add the element
	l->n++;
	if (l->n == 1)
		l->e = (Etats*)malloc(sizeof(Etats));
	else
		l->e = (Etats*)realloc(l->e, sizeof(Etats)*l->n);
	if (!l->e)
	{
		printf("Out of memory !");
		exit(4);
	}
	l->e[l->n-1] = copyEtats(e);
	if (res)
		*res = l->n-1;
	return true;
}

//add an element even if already in the list
void AddEl2 (ListEtats *l, Etats e)
{
	//add the element
	l->n++;
	if (l->n == 1)
		l->e = (Etats*)malloc(sizeof(Etats));
	else
		l->e = (Etats*)realloc(l->e, sizeof(Etats)*l->n);
	if (!l->e)
	{
		printf("Out of memory !");
		exit(5);
	}
	l->e[l->n-1] = copyEtats(e);
}

///////////////////////////////////////////////////////////////////
Etats2 NewEtats2 (int n)
{
	Etats2 e;
	e.n = n;
	e.e = (uint64 *)malloc(sizeof(uint64)*((n+63)/64));
	if (!e.e)
	{
		printf("Out of memory !");
		exit(7);
	}
return e;
}

void FreeEtats2 (Etats2 e)
{
	free(e.e);
}

void initEtats2 (Etats2 e)
{
	int i;
	int n = (e.n+63)/64;
	for (i=0;i<n;i++)
	{
		e.e[i] = 0;
	}
}

void printEtats2 (Etats2 e)
{
	int i, j;
	int n = (e.n+63)/64;
	printf("[ ");
	for (i=0;i<n;i++)
	{
		for (j=0;j<64;j++)
		{
			if (e.e[i] & (uint64)1<<j)
			{
				printf("%d ", 64*i+j);
				fflush(stdout);
			}
		}
	}
	printf("]\n");
}

bool isNullEtats2 (Etats2 e)
{
	int i;
	int n = (e.n+63)/64;
	for (i=0;i<n;i++)
	{
		if (e.e[i] != 0)
			return false;
	}
	return true;
}

bool equalsEtats2 (Etats2 e1, Etats2 e2)
{
	if (e1.n != e2.n)
		return false;
	int i;
	int n = (e1.n+63)/64;
	for (i=0;i<n;i++)
	{
		if (e1.e[i] != e2.e[i])
			return false;
	}
	return true;
}

bool hasEtats2(Etats2 e, uint64 i)
{
	return (e.e[i/64] & ((uint64)1<<(i%64))) != 0;
}

Etats2 copyEtats2(Etats2 e)
{
	Etats2 r = NewEtats2(e.n);
	int i;
	int n = (e.n+63)/64;
	for (i=0;i<n;i++)
	{
		r.e[i] = e.e[i];
	}
	return r;
}

void addEtat (Etats2 *e, uint64 i)
{
	e->e[i/64] |= ((uint64)1<<(i%64));
}

ListEtats2 NewListEtats2 (int n, int na)
{
	ListEtats2 r;
	r.n = n;
	r.na = na;
	r.e = (Etats2*)malloc(sizeof(Etats2)*na);
	return r;	
}

void ReallocListEtats2(ListEtats2 *l, int n, bool marge)
{
	if (!marge)
	{
		if (l->na)
			l->e = (Etats2*)realloc(l->e, sizeof(Etats2)*n);
		else
			l->e = (Etats2*)malloc(sizeof(Etats2)*n);
		l->na = n;
	}else
	{
		if (n > l->na)
		{
			if (l->na)
				l->e = (Etats2*)realloc(l->e, sizeof(Etats2)*n*2);
			else
				l->e = (Etats2*)malloc(sizeof(Etats2)*n*2);
			l->na = n*2;
		}
	}
	l->n = n;
}

void FreeListEtats2 (ListEtats2* l)
{
	int i;
	for (i=0;i<l->n;i++)
	{
		FreeEtats2(l->e[i]);
	}
	free(l->e);
	l->n = 0;
	l->na = 0;
}

void printListEtats2 (ListEtats2 l)
{
	int i;
	for (i=0;i<l.n;i++)
	{
		printf("%d : ", i);
		printEtats2(l.e[i]);
	}
	printf("(%d allocated states2 )\n", l.na);
}

///////////////////////////////////////////////////////////////////

Dict* lhash;
const int nhash = 10000019;
bool allocated = false;
bool filled = false;

void AllocHash ()
{
	if (!allocated)
	{
		lhash = (Dict*)malloc(sizeof(Dict)*nhash);
		if (!lhash)
		{
				printf("Out of memory !");
			exit(8);
		}
		allocated = true;
	}
	if (!filled)
	{
		int i;
		for (i=0;i<nhash;i++)
		{
			lhash[i].e = NULL;
			lhash[i].n = 0;
		}
		filled = true;
	}
}

void FreeHash ()
{
	free(lhash);
	allocated = false;
	filled = false;
}

int hashEtats (Etats e)
{
	int i;
	int h = 1;
	for (i=0;i<e.n;i++)
	{
		h *= 2;
		h += e.e[i];
		h %= nhash;
	}
	return h;
}

int hash2 (Etats2 e)
{
	int i;
	uint64 h = 1;
	int n = (e.n+63)/64;
	for (i=0;i<n;i++)
	{
		h *= 2;
		h += e.e[i];
		h %= nhash;
	}
	return h;
}

//add the element if not already in the hash table
bool addEtats2 (const ListEtats2* l, Etats2 e, int *k)
{
	int h = hash2(e);
	int i, v;
	for (i=0;i<lhash[h].n;i++)
	{
		////////check
		if (lhash[h].e[i] >= l->n)
		{
			printf("***************\nError : element of the hash table too big !!!\n****************\n");
		}
		/////////////
		v = lhash[h].e[i];
		if (k)
			*k = v;
		if (equalsEtats2(l->e[v], e))
		{
			//printf("equals !\n");
			return false;
		}
	}
	//add the element
	if (k)
		*k = l->n;
	dictAdd(&lhash[h], l->n);
	return true;
}

//add the element if not already in the hash table
bool addH (const ListEtats *l, Etats e, int* nf)
{
	int h = hashEtats(e);
	
	int i, v;
	for (i=0;i<lhash[h].n;i++)
	{
		////////verif
		if (lhash[h].e[i] >= l->n)
		{
			printf("***************\nError : element of the hash table too big !!!\n****************\n");
		}
		/////////////
		v = lhash[h].e[i];
		if (nf)
			*nf = v;
		if (v < 0)
			return false;
		if (equals(l->e[v], e))
		{
			//printf("equals !\n");
			return false;
		}
	}
	//add the element
	if (nf)
		*nf = l->n;
	dictAdd(&lhash[h], l->n);
	return true;
}

InvertDict NewInvertDict(int n)
{
	InvertDict r;
	r.n = n;
	if (n == 0)
		return r;
	r.d = (Dict*)malloc(sizeof(Dict)*n);
	if (!r.d)
	{
		printf("Out of memory !");
		exit(9);
	}
	return r;
}

InvertDict invertDict(Dict d)
{
	//count the number of different values (assumed to be consecutive)
	int i;
	int nv = 0;
	for (i=0;i<d.n;i++)
	{
		if (d.e[i] >= nv)
			nv = d.e[i]+1;
	}
	//alloc the inverse
	InvertDict r;
	r.n = nv;
	r.d = (Dict *)malloc(sizeof(Dict)*nv);
	if (!r.d)
	{
		printf("Out of memory !");
		exit(10);
	}
	//initialize
	for (i=0;i<nv;i++)
	{
		r.d[i].n = 0;
		r.d[i].e = NULL;
	}
	//fill the dictionnary
	for (i=0;i<d.n;i++)
	{
		if (d.e[i] != -1)
			dictAdd(&r.d[d.e[i]], i);
	}
	return r;
}

void FreeInvertDict (InvertDict id)
{
	if (id.n == 0)
		return;
	int i;
	for (i=0;i<id.n;i++)
	{
		FreeDict(&id.d[i]);
	}
	free(id.d);
}

void printInvertDict (InvertDict id)
{
	int i;
	for (i=0;i<id.n;i++)
	{
		printf("%d : ", i);
		printDict(id.d[i]);
	}
}

////////////////////////////////// to improve, with a hash table !!!!
void putEtat (Etats *f, int ef)
{
	int i;
	for (i=0;i<f->n;i++)
	{
		if (f->e[i] == ef)
			return;
	}
	f->e[f->n] = ef;
	f->n++;
}

//function used by Determinize()
//Etats : list of states of a
void Determinize_rec (Automaton a, InvertDict id, Automaton *r, ListEtats* l, bool onlyfinals, bool nof, int niter)
{
	int current = l->n-1;
	Etats c = l->e[current];
	//Browse sons
	Etats f = NewEtats(a.n);
	int nf;
	Etat e;
	int i,j,k;
	int ef;
	bool final;
	
	for (i=0;i<id.n;i++) //browse letters of the new alphabet
	{
		//fill f (list of states of a corresponding to the state of r where we come)
		f.n = 0;
		final = false;
		for (j=0;j<c.n;j++) //browse the states of the list
		{
			e = a.e[c.e[j]]; //corresponding state
			for (k=0;k<id.d[i].n;k++) //browse the original letters, corresponding to the new chosen letter
			{
				ef = e.f[id.d[i].e[k]];
				if (ef != -1)
				{
					//check that the state of a is not already in the list
					putEtat(&f, ef);
					if (a.e[ef].final)
						final = true;
				}
			}
		}
		if (onlyfinals && !final)
			continue;
		if (nof && final)
			continue;
		//test if the state has been seen, and otherwise mark it as seen
		if (addH(l, f, &nf)) //ajoute l'état à la table de hachage si nouveau
		{	
			//add the state to the list
			AddEl2(l, f);
			//add the state to r
			if (nof)
				final = true;
			AddEtat(r, final);
			//induction
			Determinize_rec(a, id, r, l, onlyfinals, nof, niter+1);
		}
		if (nf != -1)
		{
			//add the edge
			r->e[current].f[i] = nf;
		}
	}
	FreeEtats(f);
}

//Determinize the automaton obtained by changing the alphabet, using the Power Set Construction.
Automaton Determinize(Automaton a, Dict d, bool noempty, bool onlyfinals, bool nof, bool verb)
{
	int i;
	
	Automaton r;
	if (a.i == -1)
	{
		//compute the size of the alphabet
		int nv = 0;
		for (i=0;i<d.n;i++)
		{
			if (d.e[i] >= nv)
				nv = d.e[i]+1;
		}
		//
		if (verb)
			printf("No initial state !\n");
		if (nof)
		{
			r = NewAutomaton(1, nv);
			r.i = 0;
			r.e[0].final = true;
			for (i=0;i<nv;i++)
			{
				r.e[0].f[i] = 0;
			}
		}else
			r = NewAutomaton(0, nv);
		return r;
	}
	
/*
	//increase the stack size
	const rlim_t kStackSize = 32 * 1024 * 1024;
	struct rlimit rl;
	int result;
	result = getrlimit(RLIMIT_STACK, &rl);
	if (result == 0)
	{
		if (rl.rlim_cur < kStackSize)
		{
			if (verb)
				printf("limite : %d -> %d\n", rl.rlim_cur, kStackSize);
			rl.rlim_cur = kStackSize;
			result = setrlimit(RLIMIT_STACK, &rl);
			if (result != 0)
			{
				fprintf(stderr, "setrlimit returned result = %d\n", result);
			}
		}
	}
	//
*/
	
	if (verb)
	{
		if (onlyfinals)
			printf("onlyfinals\n");
		if (nof)
			printf("nof\n");
		if (noempty)
			printf("noempty\n");
		printf("Dictionary : ");
		printDict(d);
	}
	
	//compute the inverse of the dictionnary
	InvertDict id = invertDict(d);
	if (verb && id.n == d.n)
	{
		printf("The dictionary is inversible : trivial determination !\n");
	}
	
	if (verb)
	{
		printf("Inverse Dictionary :\n");
		printInvertDict(id);
	}
	
	//alloc the hash table
	AllocHash();
	//put the empty set in the hash table if we don't want this state in r
	if (noempty)
	{
		Etats e = NewEtats(0); //empty set
		int h = hashEtats(e);
		if (verb)
			printf("hash empty : %d\n", h);
		dictAdd(&lhash[h], -1);
	}
	
	//initialize the result automaton with just the initial state
	if (verb)
		printf("Init r...\n");
	r.n = 1;
	r.na = id.n;
	r.i = 0;
	r.e = (Etat *)malloc(sizeof(Etat));
	if (!r.e)
	{
		printf("Out of memory !");
		exit(11);
	}
	if (nof)
		r.e[0].final = true;
	else
		r.e[0].final = a.e[a.i].final;
	r.e[0].f = (int *)malloc(sizeof(int)*r.na);
	if (!r.e[0].f)
	{
		printf("Out of memory !");
		exit(12);
	}
	for (i=0;i<r.na;i++)
	{
		r.e[0].f[i] = -1;
	}
	if (verb)
		printAutomaton(r);
	
	//initialize the list of states of the new automaton and the hash table
	//(this list is used to number with consecutive numbers the states of the new automaton that are list of states of a)
	if (verb)
		printf("Init l...\n");
	ListEtats l;
	l.n = 0;
	l.e = NULL;
	Etats e = NewEtats(1);
	e.e[0] = a.i; //the initial state of the new automaton is the list of initial states (here with cardinality 1)
	addH(&l, e, NULL);
	AddEl2(&l, e);
	
	//initialize the hash table
	bool b = addH(&l, l.e[0], NULL); //add the initial state to the hash table
	if (verb)
		printListEtats(l);
	
	if (verb)
		printf("Induction...\n");
	
	Determinize_rec(a, id, &r, &l, onlyfinals, nof, 0);
	
	if (verb)
		printf("Free...\n");
	
	//put back to zero the element of the hash table
	for (i=0;i<l.n;i++)
	{
		FreeDict(&lhash[hashEtats(l.e[i])]);
	}
	
	for (i=0;i<l.n;i++)
	{
		FreeEtats(l.e[i]);
	}
	free(l.e);
	
	FreeInvertDict(id);
	//FreeHash(); //do not free the memory, in order to use it faster after
	
	return r;
}

NAutomaton Concat (Automaton a, Automaton b, bool verb)
{
	int i, j, f;
	// !!! should take care of the fact that alphabet could be differents !!!
	// We assume for the moment that a and b share the same alphabet
	NAutomaton r = NewNAutomaton(a.n+b.n, a.na);
	//copie a
	for (i=0;i<a.n;i++)
	{
		//count the edges
		int nv = 0, rnv;
		for (j=0;j<a.na;j++)
		{
			f = a.e[i].f[j];
			if (f != -1)
				nv++;
		}
		if (a.e[i].final)
			nv++; //add an edge to the other automaton
		rnv = nv;
		//alloc new edges
		r.e[i].a = (Arete *)malloc(sizeof(Arete)*nv);
		r.e[i].n = nv;
		//copy edges
		nv = 0;
		for (j=0;j<a.na;j++)
		{
			f = a.e[i].f[j];
			if (f != -1)
			{
				r.e[i].a[nv].l = j;
				r.e[i].a[nv].e = f;
				nv++;
			}
		}
		if (a.e[i].final)
		{
			//add an espilon-transition to the initial state of the other automaton
			r.e[i].a[nv].l = -1;
			r.e[i].a[nv].e = a.n+b.i;
			nv++;
		}
		if (rnv != nv)
		{
			printf("Error : nv has not the right value !!!\n");
			exit(1);
		}
		r.e[i].final = false;
		r.e[i].initial = (a.i == i);
	}
	//copy b
	for (i=0;i<b.n;i++)
	{
		//count the edges
		int nv = 0;
		for (j=0;j<b.na;j++)
		{
			if (b.e[i].f[j] != -1)
				nv++;
		}
		//alloc the new edges
		r.e[a.n+i].a = (Arete *)malloc(sizeof(Arete)*nv);
		r.e[a.n+i].n = nv;
		//copy the edges
		nv = 0;
		for (j=0;j<b.na;j++)
		{
			if (b.e[i].f[j] != -1)
			{
				r.e[a.n+i].a[nv].l = j;
				r.e[a.n+i].a[nv].e = a.n + b.e[i].f[j];
				nv++;
			}
		}
		r.e[a.n+i].final = b.e[i].final;
		r.e[a.n+i].initial = false;
	}
	return r;
}

//convert a deterministic automaton to a non-deterministic one
NAutomaton CopyN(Automaton a, bool verb)
{
	int i,j;
	NAutomaton r = NewNAutomaton(a.n, a.na);
	for (i=0;i<a.n;i++)
	{
		//count the edges
		int nv = 0;
		for (j=0;j<a.na;j++)
		{
			if (a.e[i].f[j] != -1)
				nv++;
		}
		//alloc the new edges
		r.e[i].a = (Arete *)malloc(sizeof(Arete)*nv);
		r.e[i].n = nv;
		//copy the edges
		nv = 0;
		for (j=0;j<a.na;j++)
		{
			if (a.e[i].f[j] != -1)
			{
				r.e[i].a[nv].l = j;
				r.e[i].a[nv].e = a.e[i].f[j];
				nv++;
			}
		}
		r.e[i].final = a.e[i].final;
		r.e[i].initial = (a.i == i);
	}
	return r;
}

//change the alphabet of the automaton
NAutomaton Proj (Automaton a, Dict d, bool verb)
{
	//compute the size of the new alphabet
	int nv = 0;
	int i,j;
	for (i=0;i<d.n;i++)
	{
		if (d.e[i] >= nv)
			nv = d.e[i]+1;
	}
	//
	NAutomaton r = NewNAutomaton(a.n, nv);
	for (i=0;i<a.n;i++)
	{
		//count the edges
		nv = 0;
		for (j=0;j<a.na;j++)
		{
			if (a.e[i].f[j] != -1 && d.e[j] != -1)
				nv++;
		}
		//alloc new edges
		r.e[i].a = (Arete *)malloc(sizeof(Arete)*nv);
		r.e[i].n = nv;
		//copy edges
		nv = 0;
		for (j=0;j<a.na;j++)
		{
			if (a.e[i].f[j] != -1 && d.e[j] != -1)
			{
				r.e[i].a[nv].l = d.e[j];
				r.e[i].a[nv].e = a.e[i].f[j];
				nv++;
			}
		}
		r.e[i].final = a.e[i].final;
		r.e[i].initial = (a.i == i);
	}
	return r;
}

//Add to e all states reached by an epsilon-transition from state i
void EpsilonParcours (NAutomaton a, uint i, Etats2 e)
{
	if (!hasEtats2(e, i))
	{
		addEtat(&e, i);
		uint j;
		for (j=0;j<a.e[i].n;j++)
		{
			if (a.e[i].a[j].l == -1) //it is an epsilon-transition
			{
				EpsilonParcours(a, a.e[i].a[j].e, e);
			}
		}
	}
}

//compute the epsilon-closure of e (result in ec)
void EpsilonCloture (NAutomaton a, Etats2 e, Etats2 ec)
{
	int j;
	initEtats2(ec);
	for (j=0;j<a.n;j++)
	{
		if (hasEtats2(e, j))
		{
			EpsilonParcours(a, j, ec);
		}
	}
}

//Determinize a non-deterministic automaton
Automaton DeterminizeN (NAutomaton a, bool puits, int verb)
{
	if (verb)
		printf("allocation...\n");
	
	Automaton r = NewAutomaton(a.n, a.na);
	
	if (a.n == 0)
		return r;
	
	if (verb >= 20)
		printf("allocates the hash table...\n");
	
	AllocHash();
	
	if (verb >= 20)
		printf("allocates states...\n");
	
	int i,j,k,u;
	ListEtats2 l = NewListEtats2(1, 1024);
	Etats2* e = (Etats2*)malloc(sizeof(Etats2)*a.na);
	for (i=0;i<a.na;i++)
	{
		e[i] = NewEtats2(a.n);
	}
	
	Etats2 ec = NewEtats2(a.n);
	
	if (verb >= 20)
		printf("allocates the first state...\n");
	
	//initialize the first state
	r.i = 0;
	l.e[0] = NewEtats2(a.n);
	initEtats2(l.e[0]);
	for (i=0;i<a.n;i++)
	{
		if (a.e[i].initial)
			addEtat(&l.e[0], i);
	}
	l.n = 0;
	addEtats2(&l, l.e[0], NULL); //add to the hash table
	l.n++;
	
	if (verb)
		printf("way...\n");
	
	for (i=0;i<l.n;i++)
	{ //browse states of the new automaton
		
		if (verb >= 20)
		{
			printf("state %d : ", i);
			printEtats2(l.e[i]);
		}
		
		r.e[i].final = false;
		for (j=0;j<a.na;j++)
		{
			initEtats2(e[j]);
		}
		//compute the closure by epsilon-transition of l.e[i]
		EpsilonCloture(a, l.e[i], ec);
		//browse all edges leaving from states of ec
		for (j=0;j<a.n;j++)
		{ //browse states of a being in the current state
			if (hasEtats2(ec, j))
			{
				r.e[i].final = r.e[i].final || a.e[j].final;
				for (u=0;u<a.e[j].n;u++)
				{ //browse edges leaving from state j of a
					if (a.e[j].a[u].l >= 0)
						addEtat(&e[a.e[j].a[u].l], a.e[j].a[u].e);
				}
			}
		}
		if (verb > 20)
			printf(" -> reached states :\n");
		for (j=0;j<a.na;j++)
		{ //browse reached states
			if (verb > 20)
			{
				printf("	");
				printEtats2(e[j]);
			}
			if (puits || !isNullEtats2(e[j]))
			{
				EpsilonCloture(a, e[j], ec);
				//determine if the state is new or not
				if (addEtats2(&l, ec, &k))
				{
					ReallocListEtats2(&l, l.n+1, true);
					l.e[k] = copyEtats2(ec);
					//add a state to the automaton
					if (l.n > r.n)
					{ //alloc the memory if needed
						ReallocAutomaton(&r, 2*l.n, false);
					}
				}
				//add a edge labeled by j toward k in the new automaton
				r.e[i].f[j] = k;
			}else
			{
				r.e[i].f[j] = -1;
			}
		}
	}
	
	if (verb)
		printf("free...\n");
	
	//libère la mémoire
	ReallocAutomaton(&r, l.n, false);
	
	//put to zero elements of the hash table
	for (i=0;i<l.n;i++)
	{
		FreeDict(&lhash[hash2(l.e[i])]);
	}
	
	FreeListEtats2(&l);
	for (i=0;i<a.na;i++)
	{
		FreeEtats2(e[i]);
	}
	free(e);
	//FreeHash(); //do not free memory in order to use it faster later
	return r;
}

//change the alphabet by duplicating edges if needed
//the result is assumed to be deterministic !!!!
Automaton Duplicate (Automaton a, InvertDict id, int na2, bool verb)
{
	if (verb)
	{
		printf("NewAutomaton(%d, %d)\n", a.n, na2);
	}
	Automaton r = NewAutomaton(a.n, na2);

	if (verb)
	{
		printf("NewAutomaton(%d, %d) done.\n", a.n, na2);
	}
	
	int i,j,k;
	r.i = a.i;
	for (i=0;i<r.n;i++)
	{
		r.e[i].final = a.e[i].final;
		for (j=0;j<r.na;j++)
		{
			r.e[i].f[j] = -1;
		}
		for (j=0;j<a.na;j++)
		{
			for (k=0;k<id.d[j].n;k++)
			{
				r.e[i].f[id.d[j].e[k]] = a.e[i].f[j];
			}
		}
	}
	return r;
}

/*
//stabilization of the language by prefix (i.e. 
void PrefixStabilize (Automaton *a)
{
		
}
*/

void ZeroComplete_rec(Automaton *a, int etat, bool *vu, int l0, bool verb)
{
	if (verb)
		printf("state %d ..\n", etat);
	vu[etat] = true;
	int i, e;
	for (i=0;i<a->na;i++)
	{
		e = a->e[etat].f[i];
		if (e != -1 && e < a->n)
		{
			if (!vu[e])
				ZeroComplete_rec(a, e, vu, l0, verb);
			if (i == l0 && a->e[e].final)
				a->e[etat].final = true;
		}
	}
}

void ZeroComplete(Automaton *a, int l0, bool verb)
{
	if (verb)
		printf("l0 = %d\n", l0);
	if (a->i == -1)
		return;
	bool *vu = (bool *)malloc(sizeof(bool)*a->n); //list of seen states
	if (!vu)
	{
		printf("Out of memory !\n");
		exit(25);
	}
	int i;
	for (i=0;i<a->n;i++)
		vu[i] = false;
	ZeroComplete_rec(a, a->i, vu, l0, verb);
	free(vu);
}

Automaton ZeroComplete2 (Automaton *a, int l0, bool etat_puits, bool verb)
{
	NAutomaton r = NewNAutomaton(a->n+1, a->na);
		
	int i,j,k;
	for (i=0;i<a->n;i++)
	{
		if (i == a->i)
			r.e[i].initial = true;
		else
			r.e[i].initial = false;
		r.e[i].final = a->e[i].final;
		r.e[i].n = 0;
		//compte les arêtes
		for (j=0;j<a->na;j++)
		{
			if (a->e[i].f[j] != -1)
				r.e[i].n++;
		}
		if (a->e[i].final)
			r.e[i].n++; //edge 0 added
		//alloc
		r.e[i].a = (Arete *)malloc(sizeof(Arete)*r.e[i].n);
		//fill
		k = 0;
		for (j=0;j<a->na;j++)
		{
			if (a->e[i].f[j] != -1)
			{
				r.e[i].a[k].e = a->e[i].f[j];
				r.e[i].a[k].l = j; 
				k++;
			}
		}
		if (a->e[i].final)
		{
			r.e[i].a[k].e = a->n; //add the edge to the state that recognize 0
			r.e[i].a[k].l = l0;
		}
	}
	r.e[a->n].n = 1;
	r.e[a->n].a = (Arete *)malloc(sizeof(Arete));
	r.e[a->n].a[0].e = a->n;
	r.e[a->n].a[0].l = l0;
	r.e[a->n].initial = false;
	r.e[a->n].final = true;
	
	return DeterminizeN(r, etat_puits, 0);
}

Automaton ZeroInv(Automaton *a, int l0)
{
	NAutomaton r = NewNAutomaton(a->n+1, a->na);
		
	int i,j,k;
	for (i=0;i<a->n;i++)
	{
		if (i == a->i)
			r.e[i].initial = true;
		else
			r.e[i].initial = false;
		r.e[i].final = a->e[i].final;
		r.e[i].n = 0;
		//count the edges
		for (j=0;j<a->na;j++)
		{
			if (a->e[i].f[j] != -1)
				r.e[i].n++;
		}
		//alloc
		r.e[i].a = (Arete *)malloc(sizeof(Arete)*r.e[i].n);
		//fill
		k = 0;
		for (j=0;j<a->na;j++)
		{
			if (a->e[i].f[j] != -1)
			{
				r.e[i].a[k].e = a->e[i].f[j];
				r.e[i].a[k].l = j; 
				k++;
			}
		}
	}
	r.e[a->n].n = 2;
	r.e[a->n].a = (Arete *)malloc(sizeof(Arete)*2);
	r.e[a->n].a[0].e = a->n;
	r.e[a->n].a[0].l = l0;
	r.e[a->n].a[1].e = a->i;
	r.e[a->n].a[1].l = l0;
	r.e[a->n].initial = true;
	r.e[a->n].final = true;
	
	return DeterminizeN(r, false, 0);
}

int countEtats = 0;
bool prune_inf_rec(Automaton a, int etat)
{
	int i, f;
	bool cycle = false;
	a.e[etat].final = 1; //mark the state as currently seen
	for (i=0;i<a.na;i++)
	{
		f = a.e[etat].f[i];
		if (f == -1)
			continue;
		if (a.e[f].final == 1)
			cycle = true; //the state is part of a cycle
		if (a.e[f].final == 0)
		{
			if (prune_inf_rec(a, f))
				cycle = true; //the state permits to reach a cycle
		}
	}
	if (!cycle)
		a.e[etat].final = 2; //indicate that we won't keep the state (but it has been seen)
	else
		countEtats++; //count the states that we keep
	return cycle;
}

//remove every state from which there is no infinite path
Automaton prune_inf(Automaton a, bool verb)
{
	int *l = (int *)malloc(sizeof(int)*a.n);
	if (!l)
	{
		printf("Out of memory !\n");
		exit(25);
	}
	
	//start by computing the list of states of a to keep
	int i;
	int *finaux = (int *)malloc(sizeof(int)*a.n);
	if (!finaux)
	{
		printf("Out of memory !\n");
		exit(13);
	}
	for (i=0;i<a.n;i++)
	{
		finaux[i] = a.e[i].final;
		a.e[i].final = 0; //states not seen
	}
	if (verb)
		printf("induction...\n");
	countEtats = 0;
	if (a.i != -1)
		prune_inf_rec (a, a.i);
	
	if (verb)
	{
		printf("States counter = %d\n", countEtats);
		printf("count...\n");
	}
	
	//count the states to keep
	int cpt = 0;
	for (i=0;i<a.n;i++)
	{
		if (a.e[i].final & 1)
		{ //new state to add
			l[i] = cpt;
			cpt++;
		}else
			l[i] = -1;
	}
	
	if (verb)
		printf("cpt = %d\n", cpt);
	
	//create the new automaton
	int j, f;
	Automaton r = NewAutomaton(cpt, a.na);
	for (i=0;i<a.n;i++)
	{
		if (l[i] == -1)
			continue;
		for (j=0;j<a.na;j++)
		{
			f = a.e[i].f[j];
			r.e[l[i]].f[j] = -1;
			if (f == -1)
				continue;
			if (l[f] != -1)
			{
				r.e[l[i]].f[j] = l[f];
			}
		}
	}
	
	if (verb)
		printf("final states...\n");
	
	//put back final states as before
	for (i=0;i<a.n;i++)
	{
		a.e[i].final = finaux[i];
		if (l[i] != -1)
			r.e[l[i]].final = finaux[i];
	}
	
	//initial state
	if (a.i != -1)
		r.i = l[a.i];
	else
		r.i = -1;
	
	free(finaux);
	free(l);
	return r;
}

//Compute the transposition, assuming it is deterministic
Automaton TransposeDet(Automaton a)
{
	Automaton r = NewAutomaton(a.n, a.na);
	int i,j;
	for (i=0;i<a.n;i++)
	{
		if (a.e[i].final)
			r.i = i;
		if (i == a.i)
			r.e[i].final = true;
		else
			r.e[i].final = false;
		for (j=0;j<a.na;j++)
		{
			r.e[i].f[j] = -1;
		}
	}
	int f;
	for (i=0;i<a.n;i++)
	{
		for (j=0;j<a.na;j++)
		{
			f = a.e[i].f[j];
			if (f != -1)
			{
				//printf("%d --%d--> %d\n", i, j, f);
				r.e[f].f[j] = i;
			}
		}
	}
	return r;
}

//Compute the transposition
NAutomaton Transpose(Automaton a)
{
	NAutomaton r = NewNAutomaton(a.n, a.na);
		
	int i,j;
	for (i=0;i<a.n;i++)
	{
		if (i == a.i)
			r.e[i].final = true;
		else
			r.e[i].final = false;
		r.e[i].initial = a.e[i].final;
		r.e[i].a = NULL;
		r.e[i].n = 0;
	}
		
	int f;
	for (i=0;i<a.n;i++)
	{
		for (j=0;j<a.na;j++)
		{
			f = a.e[i].f[j];
			if (f != -1)
			{
				//add an edge from f to i labeled by j
				r.e[f].n++;
				if (r.e[f].n == 1)
					r.e[f].a = (Arete *)malloc(sizeof(Arete));
				else
					r.e[f].a = (Arete *)realloc(r.e[f].a, sizeof(Arete)*r.e[f].n);
				r.e[f].a[r.e[f].n-1].l = j;
				r.e[f].a[r.e[f].n-1].e = i;
			}
		}
	}
	return r;
}

int min (int a, int b)
{
	if (a < b)
		return a;
	return b;
}

int count2;
void StronglyConnectedComponents_rec(Automaton a, int etat, int *pile, int *m, int *res)
{
	int j,f,c;
	pile[countEtats] = etat;
	m[etat] = countEtats;
	c = countEtats;
	a.e[etat].final |= 2; //mark the state as seen
	countEtats++;
	for (j=0;j<a.na;j++)
	{
		f = a.e[etat].f[j];
		if (f == -1)
			continue;
		if (!(a.e[f].final & 2))
		{
			StronglyConnectedComponents_rec(a, f, pile, m, res);
			m[etat] = min(m[etat], m[f]);
		}else if (res[f] == -1)
		{
			m[etat] = min(m[etat], m[f]);
		}
	}
	if (m[etat] == c)
	{ //we have a strongly connected component
	    //pops the component
		do
		{
			countEtats--;
			res[pile[countEtats]] = count2;
		}while(pile[countEtats] != etat);
		count2++;
	}
}

//Tarjan's algorithm
int StronglyConnectedComponents (Automaton a, int *res)
{
	int *m = (int *)malloc(sizeof(int)*a.n);
	int *pile = (int *)malloc(sizeof(int)*a.n);
	int i;
	for (i=0;i<a.n;i++)
	{
		res[i] = -1;
	}
	countEtats = 0; //count the states added to the stack
	count2 = 0; //count the strongly connected components
	for (i=0;i<a.n;i++)
	{
		if (res[i] == -1)
			StronglyConnectedComponents_rec(a, i, pile, m, res);
	}
	//put back final states
	for (i=0;i<a.n;i++)
	{
		a.e[i].final &= 1;
	}
	free(pile);
	free(m);
	return count2;
}

//determine reachable and co-reachable states
void prune_rec(Automaton a, int *l, InvertDict id, int etat)
{
	//printf("prune_rec %d...\n", etat);
	int i, j, f;
	a.e[etat].final |= 2; //mark that this state is currently seen
	for (i=0;i<a.na;i++)
	{
		f = a.e[etat].f[i];
		if (f == -1)
			continue;
		if (!(a.e[f].final & 2))
		{ //the state has not been seen
			prune_rec(a, l, id, f);
		}
		if ((a.e[f].final & 4) && !(a.e[etat].final & 4))
		{ //we get to a co-final state that is not yet marked as co-final
		    //propagate the information to the strongly connected component
			for (j=0;j<id.d[l[etat]].n;j++)
			{
				a.e[id.d[l[etat]].e[j]].final |= 4;
			}
		}
	}
}

//remove every non-reachable and non-co-reachable states
//
// not really efficient : to improve !!!!!!!!!!!!!!!!!
//
Automaton prune(Automaton a, bool verb)
{
	int i,j,f;
	
	//determine reachable and co-reachable states
	int *l = (int *)malloc(sizeof(int)*a.n);
	if (!l)
	{
		printf("Out of memory !\n");
		exit(15);
	}
	int ncc = StronglyConnectedComponents(a, l);
	if (verb)
	{
		printf("%d components : [", ncc);
		for (i=0;i<a.n;i++)
		{
			printf(" %d", l[i]);
		}
		printf(" ]\n");
	}
	InvertDict id = NewInvertDict(ncc);
	for (i=0;i<ncc;i++)
	{
		id.d[i] = NewDict(0);
	}
	for (i=0;i<a.n;i++)
	{
		if (a.e[i].final)
			a.e[i].final = 1;
		dictAdd(&id.d[l[i]], i);
	}
	if (verb)
		printInvertDict(id);
	//propagate final states to strongly connected components
	for (i=0;i<a.n;i++)
	{
		if (a.e[i].final & 1)
		{ //the state i is final and co-reachable but not already marked as co-reachable
			for (j=0;j<id.d[l[i]].n;j++)
			{
				a.e[id.d[l[i]].e[j]].final |= 4;
				if (verb)
					printf("%d co-acc\n", id.d[l[i]].e[j]);
			}
		}
	}
	//
	if (verb)
		printf("rec...\n");
	if (a.i != -1)
		prune_rec(a, l, id, a.i);
	
	//count the states to keep
	int cpt = 0;
	for (i=0;i<a.n;i++)
	{
		if ((a.e[i].final & 2) && (a.e[i].final & 4))
		{ //new state to add
			l[i] = cpt;
			cpt++;
		}else
			l[i] = -1;
	}
	
	if (verb)
	{
		printf("l : [");
		for (i=0;i<a.n;i++)
		{
			printf(" %d(%d)", l[i], a.e[i].final);
		}
		printf(" ]\n");
		
		printf("create the new automaton %d %d...\n", cpt, a.na);
	}
	//create the new automaton
	Automaton r = NewAutomaton(cpt, a.na);
	for (i=0;i<a.n;i++)
	{
		if (l[i] == -1)
		{
			if (verb)
				printf("pass %d\n", i);
			continue;
		}
		for (j=0;j<a.na;j++)
		{
			f = a.e[i].f[j];
			r.e[l[i]].f[j] = -1;
			if (f == -1)
				continue;
			if (l[f] != -1)
			{
				r.e[l[i]].f[j] = l[f];
			}
		}
	}
	
	//put back the final states of a
	if (verb)
	{
		printf("deleted States : [");
		fflush(stdout);
	}
	for (i=0;i<a.n;i++)
	{
		if (verb)
		{
			if (l[i] == -1)
			{
				printf(" %d(", i);
				if (!(a.e[i].final & 2))
					printf(" non-acc");
				if (!(a.e[i].final & 4))
					printf(" non-co-acc");
				printf(" )");
			}
		}
		a.e[i].final &= 1;
		if (l[i] != -1)
			r.e[l[i]].final = a.e[i].final;
	}
	if (verb)
		printf(" ]\n");
	
	//initial state
	if (a.i != -1)
		r.i = l[a.i];
	else
		r.i = -1;
	
	FreeInvertDict(id);
	free(l);
	return r;
}

bool AccCoAccRec (Automaton *a, int *coa, int e)
{
	if (coa[e] == 2)
	{
		int i, f;
		bool coacc = a->e[e].final;
		coa[e] = 0; //indicate that the state has been seen
		for (i=0;i<a->na;i++)
		{
			f = a->e[e].f[i];
			if (f != -1)
				if (AccCoAccRec(a, coa, f))
					coacc = true;
		}
		if (coacc)
			coa[e] = 1;
		return coacc;
	}else
		return (coa[e] == 1);
}

//determine the reachable and co-reachable states
// 0 : non co-reachable but reachable
// 1 : co-reachable and reachable
// 2 : non-reachable (but we don't know the co-reachability)
void AccCoAcc(Automaton *a, int *coa)
{
	int i;
	for (i=0;i<a->n;i++)
	{
		coa[i] = 2; //indicate states not seen yet
	}
	if (a->i != -1)
		AccCoAccRec(a, coa, a->i);
}

bool CoAccRec (Automaton *a, int *coa, int e)
{
	if (coa[e] == 2)
	{
		int i, f;
		bool coacc = a->e[e].final;
		coa[e] = 0; //indicate that the state has been seen
		for (i=0;i<a->na;i++)
		{
			f = a->e[e].f[i];
			if (f != -1)
				if (CoAccRec(a, coa, f))
					coacc = true;
		}
		if (coacc)
			coa[e] = 1;
		return coacc;
	}else
		return (coa[e] == 1);
}

//determine the co-reachable states
// 0 : non co-reachable
// 1 : reachable
void CoAcc (Automaton *a, int *coa)
{
	int i;
	for (i=0;i<a->n;i++)
	{
		coa[i] = 2; //indicate states non seen yet
	}
	if (a->i != -1)
		CoAccRec(a, coa, a->i);
	for (i=0;i<a->n;i++)
	{
		if (coa[i] == 2)
			CoAccRec(a, coa, i);
	}
}

//determine reachable states
void pruneI_rec(Automaton a, int etat)
{
	int i, f;
	a.e[etat].final |= 2; //mark that the state is currently seen
	for (i=0;i<a.na;i++)
	{
		f = a.e[etat].f[i];
		if (f == -1)
			continue;
		if (!(a.e[f].final & 2))
		{ //the state has not been seen yet
			pruneI_rec(a, f);
		}
	}
}

//remove every non-reachable states
Automaton pruneI(Automaton a, bool verb)
{

	int i,j,f;
	//determine reachable states
	if (a.i != -1)
		pruneI_rec(a, a.i);

	int *l = (int *)malloc(sizeof(int)*a.n);
	if (!l)
	{
		printf("Out of memory !\n");
		exit(15);
	}

	//count the states to keep
	int cpt = 0;
	for (i=0;i<a.n;i++)
	{
		if (a.e[i].final & 2)
		{ //new state to add
			l[i] = cpt;
			cpt++;
		}else
			l[i] = -1;
	}
	
	//create the new automaton
	Automaton r = NewAutomaton(cpt, a.na);

	for (i=0;i<a.n;i++)
	{
		if (l[i] == -1)
			continue;
		for (j=0;j<a.na;j++)
		{
			f = a.e[i].f[j];
			r.e[l[i]].f[j] = -1;
			if (f == -1)
				continue;
			if (l[f] != -1)
			{
				r.e[l[i]].f[j] = l[f];
			}
		}
	}

	//put back the final states of a
	if (verb)
	{
		printf("deleted States : [");
    	fflush(stdout);
	}
	for (i=0;i<a.n;i++)
	{
		if (verb)
		{
			if (l[i] == -1)
			{
				printf(" %d(", i);
				if (!(a.e[i].final & 2))
					printf(" not-acc");
				if (!(a.e[i].final & 4))
					printf(" not-co-acc");
				printf(" )");
			}
		}
		a.e[i].final &= 1;
		r.e[l[i]].final = a.e[i].final;
	}
	if (verb)
		printf(" ]\n");

	//initial state
	if (a.i != -1)
		r.i = l[a.i];
	else
		r.i = -1;

	free(l);

	return r;
}

Automaton SubAutomaton(Automaton a, Dict d, bool verb)
{
	if (verb)
	{
		printf("dict = ");
		printDict(d);
	}
	Automaton r = NewAutomaton(d.n, a.na);
	int *l = (int *)malloc(sizeof(int)*a.n);
	int i,j;
	for (i=0;i<a.n;i++)
	{
		l[i] = -1;
	}
	for (i=0;i<d.n;i++)
	{
		l[d.e[i]] = i;
	}
	
	if (verb)
	{
		printf("l = [");
		for (i=0;i<a.n;i++)
		{
			printf(" %d", l[i]);
		}
		printf(" ]\n");
	}
	
	for (i=0;i<r.n;i++)
	{
		for (j=0;j<r.na;j++)
		{
			r.e[i].f[j] = -1;
		}
	}
	for (i=0;i<a.n;i++)
	{
		if (l[i] != -1)
		{
			r.e[l[i]].final = a.e[i].final;
			for (j=0;j<a.na;j++)
			{
			
				if (a.e[i].f[j] != -1)
					r.e[l[i]].f[j] = l[a.e[i].f[j]];
				else
					r.e[l[i]].f[j] = -1;
			}
		}
	}
	r.i = -1;
	if (a.i != -1)
	{
		r.i = l[a.i];
	}
	free(l);
	return r;
}

/////////////////////////
// 
//  Tout le code qui suit est à tester !!!!!!!!!!!!!!!!!!!
//
/////////////////////////

//permute les labels des arêtes
//l donne les anciens indices à partir des nouveaux
Automaton Permut (Automaton a, int *l, int na, bool verb)
{
	if (verb)
	{
		int i;
		printf("l = [ ");
		for (i=0;i<na;i++)
		{
			printf("%d ", l[i]);
		}
		printf("]\n");
	}
	Automaton r = NewAutomaton(a.n, na);
	int i,j;
	for (i=0;i<a.n;i++)
	{
		for (j=0;j<na;j++)
		{		
			if (l[j] != -1)
				r.e[i].f[j] = a.e[i].f[l[j]];
			else
				r.e[i].f[j] = -1;
		}
		r.e[i].final = a.e[i].final;
	}
	r.i = a.i;
	return r;
}

//permute les labels des arêtes SUR PLACE
//l donne les anciens indices à partir des nouveaux
void PermutOP (Automaton a, int *l, int na, bool verb)
{
	if (verb)
	{
		int i;
		printf("l = [ ");
		for (i=0;i<na;i++)
		{
			printf("%d ", l[i]);
		}
		printf("]\n");
	}
	int *lf = (int*)malloc(sizeof(int)*a.na);
	int i,j;
	for (i=0;i<a.n;i++)
	{
		//sauvegarde les arêtes
		for (j=0;j<a.na;j++)
		{
			lf[j] = a.e[i].f[j];
			a.e[i].f[j] = -1;
		}
		//met les nouvelles
		for (j=0;j<na;j++)
		{		
			if (l[j] != -1)
				a.e[i].f[j] = lf[l[j]];
		}
	}
	free(lf);
}

/////////////////////// the following is an implementation of Hopcroft's algorithm minimization

typedef int Couple[2];

int *partition; //state --> indice
int *partitioni; //indice --> state
int *class; //classe de chaque state
Couple *class_indices; //intervalle d'indices de la classe
int nclass = 0; //nb de classes
Dict **transitioni; //inverse des transitions de l'automate : state, lettre --> liste d'states
int *L; //liste des classes par rapport auxquelles il faut raffiner
int nL; //nb d'éléments de L
int *pt_visited_class; //premier indice non rencontré dans la classe
int *visited_class; //liste des classes visitées dernièrement (utilisé dans split)
int *etats; //states à parcourir

int global_n = 0;
void print_partition ()
{
	int i;
	printf("partition = [");
	for (i=0;i<global_n+1;i++)
	{
		printf(" %d", partition[i]);
	}
	printf(" ]\n");
	printf("partitioni = [");
	for (i=0;i<global_n+1;i++)
	{
		printf(" %d", partitioni[i]);
	}
	printf(" ]\n");
}

void print_classes ()
{
	//affiche la liste des classes
	int l,h,i,j;
	for (i=0;i<nclass;i++)
	{
		printf("class %d : ", i);
		l = class_indices[i][0];
		h = class_indices[i][1];
		for (j=l;j<h;j++)
		{
			printf("%d ", partitioni[j]);
		}
		printf("\n");
	}
}

//échange les states i et j
void swap (int i, int j)
{
	if (i == j)
		return;
	int k = partition[i];
	partition[i] = partition[j];
	partition[j] = k;
	partitioni[k] = j;
	partitioni[partition[i]] = i;
}

void split (int C, int a, bool verb)
{
	//compute the préimage of C
	int i,j,l,h, e, p, lp, ep, cp;
	int nrc = 0; //nombre de classes rencontrées
	l = class_indices[C][0];
	h = class_indices[C][1];
	//copie la liste des sommets à parcourir (au cas où celle-ci soit modifiée pendant le parcours)
	for (i=l;i<h;i++)
	{ //parcours la classe C
		etats[i] = partitioni[i]; //state d'indice i
	}
	for (i=l;i<h;i++)
	{ //parcours la classe C
		e = etats[i]; //state d'indice i
		for (j=0;j<transitioni[e][a].n;j++)
		{ //parcours l'image inverse de l'state e par la lettre a
			p = transitioni[e][a].e[j]; //parent
			cp = class[p]; //classe de p
			if (!pt_visited_class[cp])
			{ //la classe de p n'a pas encore été vue dans cet appel de split
				if (verb)
					printf("new visited class : %d (%d parent of %d)\n", cp, p, e);
				visited_class[nrc] = cp;
				pt_visited_class[cp] = class_indices[cp][0]; //lowest indice of the class of p
				nrc++;
			}else
			{
				if (verb)
					printf("re-visited class : %d (%d parent of %d)\n", cp, p, e);
			}
			ep = pt_visited_class[cp]; //indice de l'élément à permuter avec p
			if (ep > partition[p])
			{
				if (verb)
					printf("vertex %d already seen\n", p);
				continue; //on a déjà vu l'state p
			}
			ep = partitioni[ep]; //élément à permuter avec p
			swap(ep, p);
			pt_visited_class[cp]++;
			//if (verb)
			//	print_partition();
		}
	}
	
	if (verb)
	{
		//print_partition();
		print_classes();
		printf("%d class encountered\n", nrc);
	}
	
	//create new classes
	for (i=0;i<nrc;i++)
	{
		cp = visited_class[i];
		
		/////only for verification : to be avoided
		if (pt_visited_class[cp] > class_indices[cp][1])
		{
			printf("***********\nError !!!\n***********\n");
		}
		////
		
		l = class_indices[cp][0];
		h = class_indices[cp][1];
		j = pt_visited_class[cp];
		
		if (verb)
			printf("class %d : l = %d %d %d = h\n", cp, l, j, h);
		
		if (j < h)
		{ //on doit ajouter une nouvelle classe
			//choisi la plus petite classe
			if (h - j > j - l)
			{ //on choisi la partie gauche
				class_indices[cp][0] = j; //l'ancienne classe devient la partie droite
				class_indices[nclass][0] = l;
				class_indices[nclass][1] = j;
			}else
			{ //on choisit la partie droite
				class_indices[cp][1] = j; //l'ancienne classe devient la partie gauche
				class_indices[nclass][0] = j;
				class_indices[nclass][1] = h;
			}
			//met à jour les classes des sommets
			for (j=class_indices[nclass][0];j<class_indices[nclass][1];j++)
			{
				class[partitioni[j]] = nclass;
			}
			L[nL] = nclass; //ajoute la nouvelle classe à L
			nL++;
			nclass++;
		}
		
		pt_visited_class[cp] = 0; //remet à 0
	}
}

//minimisation par l'algo d'Hopcroft
//voir "Around Hopcroft’s Algorithm" de Manuel BACLET and Claire PAGETTI
Automaton Minimise(Automaton a, bool verb)
{
	if (verb)
		global_n = a.n;
	//allocations
	transitioni = (Dict **)malloc(sizeof(Dict *)*(a.n+1)); //inverse des partitions
	partition = (int *)malloc(sizeof(int)*(a.n+1));
	partitioni = (int *)malloc(sizeof(int)*(a.n+1));
	class = (int *)malloc(sizeof(int)*(a.n+1)); //classe d'un state
	nclass = 0;
	class_indices = (Couple *)malloc(sizeof(Couple)*(a.n+1));
	visited_class = (int *)malloc(sizeof(int)*(a.n+1));
	pt_visited_class = (int *)malloc(sizeof(int)*(a.n+1));
	etats = (int *)malloc(sizeof(int)*(a.n+1));
	L = (int *)malloc(sizeof(int)*(a.n+1));	 //liste des classes à partir desquelles raffiner
	nL = 0;
	int i,j,f;
	//initialise
	for (i=0;i<a.n+1;i++)
	{
		partition[i] = i;
		partitioni[i] = i;
		pt_visited_class[i] = 0;
	}
	
	//if (verb)
	//	print_partition();
	
	//initialise l'inverse des transitions
	for (i=0;i<a.n+1;i++)
	{
		transitioni[i] = (Dict *)malloc(sizeof(Dict)*a.na);
		for (j=0;j<a.na;j++)
		{
			transitioni[i][j].e = NULL;
			transitioni[i][j].n = 0;
		}
	}
	for (i=0;i<a.n;i++)
	{
		for (j=0;j<a.na;j++)
		{
			f = a.e[i].f[j];
			if (f != -1)
			{
				dictAdd(&transitioni[f][j], i);
			}else
			{
				dictAdd(&transitioni[a.n][j], i); //state puits
			}
		}
	}
	for (j=0;j<a.na;j++)
	{
		dictAdd(&transitioni[a.n][j], a.n); //transitions de l'state puits
	}
	
	/**/
	if (verb)
	{
		for (i=0;i<a.n+1;i++)
		{
			for (j=0;j<a.na;j++)
			{
				printf("transition i[%d][%d] = [", i, j);
				for (f=0;f<transitioni[i][j].n;f++)
				{
					printf(" %d", transitioni[i][j].e[f]);
				}
				printf(" ]\n");
			}
		}
	}
	/**/
	
	//commence par séparer states finaux et non-finaux
	f = 0; //count du nombre d'states finaux
	for (i=0;i<a.n;i++)
	{
		if (a.e[i].final)
		{
			class[i] = 0;
			//printf("swap %d %d\n", partitioni[f], i);
			swap(partitioni[f], i);
			f++;
			/*
			if (verb)
			{
				printf("%d final\n", i);
				print_partition();
			}
			
			*/
		}else
			class[i] = 1;
	}
	class[a.n] = 1;
	//classe 0 : final states, classe 1 : the rest
	class_indices[0][0] = 0;
	class_indices[0][1] = f;
	class_indices[1][0] = f;
	class_indices[1][1] = a.n+1;
	visited_class[0] = 0;
	visited_class[1] = 0;
	nclass = 2;
	
	if (verb)
		print_partition();
	
	if (verb)
	{
		printf("Initial partition :\n");
		print_classes();
	}
	
	//choisi la classe la plus petite
	if (f <= (a.n+1)/2)
		L[0] = 0;
	else
		L[0] = 1;
	nL = 1;
	
	//algo
	int C; //current class
	while (nL)
	{
		//retire la première classe de la liste L
		nL--;
		C = L[nL];
		//partionne selon cette classe
		for (j=0;j<a.na;j++)
		{
			if (verb)
				printf("split %d %d...\n", C, j);
			split(C, j, verb);
		}
	}
	
	if (verb)
	{
		printf("Final partition :\n");
		print_classes();
	}
	
	//créé le nouvel automate
	int e;
	Automate r = NewAutomaton(nclass, a.na);
	for (i=0;i<nclass;i++)
	{
		e = partitioni[class_indices[i][0]]; //un state de la classe
		if (e >= a.n)
		{ //state puits
			for (j=0;j<a.na;j++)
				r.e[i].f[j] = -1;
			r.e[i].final = false;
			continue;
		}
		for (j=0;j<a.na;j++)
		{
			if (a.e[e].f[j] != -1)
			{
				f = class[a.e[e].f[j]];
				r.e[i].f[j] = f;
			}else
				r.e[i].f[j] = -1;
		}
		r.e[i].final = a.e[e].final;
	}
	
	if (verb)
	{
		printf("a.i = %d", a.i);
		if (a.i != -1)
		{
			printf(" class %d", class[a.i]);
		}
		printf("\n");
	}
	
	if (a.i != -1)
		r.i = class[a.i];
	else
		r.i = -1;
	
	//retire l'state puits si pas présent dans l'automate initial
	i = class[a.n];
	if (class_indices[i][1] == class_indices[i][0]+1)
	{ //il faut retirer l'state puits
		if (verb)
			printf("removes the hole state  %d...\n", i);
		DeleteVertexOP(&r, i);
	}
	
	//libère la mémoire
	free(transitioni);
	free(partition);
	free(partitioni);
	free(class);
	free(class_indices);
	free(visited_class);
	free(pt_visited_class);
	free(etats);
	free(L);
	
	return r;
}

//////////////////////////////////////////////////////////////////////////////////////////////

/*
int sign (int a)
{
	if (a > 0)
		return 1;
	if (a < 0)
		return -1;
	return 0;
}

int delta (int a)
{
	if (a)
		return 1;
	return 0;
}
*/

void DeleteVertexOP(Automaton *a, int e)
{
	if (e < 0 || e >= a->n)
		printf("The states %d is not in the automaton !\n", e);
	int i,j,f;
	a->n--;
	if (!a->n)
	{
		free(a->e);
		a->e = NULL;
		return;
	}
	for (i=0;i<a->n;i++)
	{
		for (j=0;j<a->na;j++)
		{
			f = a->e[i+(i>=e)].f[j];
			if (f != e)
				a->e[i].f[j] = f-(f>=e);
			else
				a->e[i].f[j] = -1;
		}
		a->e[i].final = a->e[i+(i>=e)].final;
	}
	if (e == a->i)
		a->i = -1;
	else
		a->i = a->i - (a->i>=e);
}

Automaton DeleteVertex (Automaton a, int e)
{
	if (e < 0 || e >= a.n)
		printf("The state %d is not in the automaton !\n", e);
	Automaton r = NewAutomaton(a.n-1, a.na);
	int i,j,f;
	for (i=0;i<a.n-1;i++)
	{
		for (j=0;j<a.na;j++)
		{
			f = a.e[i+(i>=e)].f[j];
			if (f != e)
				r.e[i].f[j] = f-(f>=e);
			else
				r.e[i].f[j] = -1;
		}
		r.e[i].final = a.e[i+(i>=e)].final;
	}
	r.i = a.i - (a.i>=e);
	return r;
}

Automaton BiggerAlphabet (Automaton a, Dict d, int nna)
{
	if (d.n != a.na)
	{
		printf("BA Error : the dictionnary must be of the same size as the previous alphabet.\n");
		return NewAutomaton(0,0);
	}
	Automaton r = NewAutomaton(a.n, nna);
	init(&r);
	int i,j;
	for (i=0;i<a.n;i++)
	{
		for (j=0;j<a.na;j++)
		{
			r.e[i].f[d.e[j]] = a.e[i].f[j];
		}
		r.e[i].final = a.e[i].final;
	}
	r.i = a.i;
	return r;
}

void Test ()
{
	printf("sizeof(Automaton)=%ld\n", sizeof(Automaton));
}
