#include <stdlib.h>
#include "complex.h"
#include "Automaton.h"
#include "relations.h"

Automate NewAutomaton (int n, int na)
{
	Automate a;
	a.n = n;
	a.na = na;
	if (n == 0)
	{
		a.i = -1;
		a.e = NULL;
		return a;
	}
	a.e = (Etat *)malloc(sizeof(Etat)*n);
	//a.nalloc = n;
	if (!a.e)
	{
		printf("Out of memory !");
		exit(6);
	}
	int i;
	for (i=0;i<n;i++)
	{
		a.e[i].f = (int *)malloc(sizeof(int)*na);
		if (!a.e[i].f)
		{
			printf("Out of memory !");
			exit(7);
		}
	}
	return a;
}

#define nhash 1000003
struct Elements
{
	Element *e;
//	int *ne; //nombre de fois qu'a été vu chaque élément
	int *ind; //indice de chaque élément
	int n;
};
typedef struct Elements Elements;

InfoBetaAdic iba; //variable globale contenant les infos sur le développement en base bêta nécessaires au calcul de l'automate des relations

Element NewElement (int n)
{
	Element e;
	e.c = (int*)malloc(sizeof(int)*n);
	return e;
}

void FreeElement (Element e)
{
	free(e.c);
}

Elements *thash; //table de hachage
	
Element *pile; //pile utilisée pendant le parcours
int npile;

InfoBetaAdic allocInfoBetaAdic (int n, int na, int ncmax, bool verb)
{
	if (verb)
	{
		printf("alloc n=%d na=%d ncmax=%d...\n", n, na, ncmax);
	}
	iba.n = n;
	iba.ncmax = ncmax;
	iba.nc = 0;
	iba.bn = NewElement(n);
	iba.c = (Element *)malloc(sizeof(Element)*ncmax);
	int i;
	for (i=0;i<ncmax;i++)
	{
		iba.c[i] = NewElement(n);
	}
	iba.na = na;
	iba.p = (PlaceArch *)malloc(sizeof(PlaceArch)*na);
	for (i=0;i<na;i++)
	{
		iba.p[i].c = (Complexe *)malloc(sizeof(Complexe)*n);
	}
	iba.cM = (double *)malloc(sizeof(double)*na);
	
	thash = (Elements*)malloc(sizeof(Elements)*nhash);
	
	npile = 1000000;
	pile = (Element *)malloc(sizeof(Element)*npile);
	return iba;
}

void initHash ()
{
	int i;
	for (i=0;i<nhash;i++)
	{
		thash[i].e = NULL;
		thash[i].n = 0;
	}
}

void freeInfoBetaAdic (InfoBetaAdic iba)
{
	int i;
	for (i=0;i<iba.ncmax;i++)
	{
		FreeElement(iba.c[i]);
	}
	FreeElement(iba.bn);
	free(iba.c);
	free(iba.p);
}

//calcul le fils de l'élément e pour le chiffre i
void succ (Element e, int i, Element *r)
{
	int j;
	//multiplication par b
	for (j=0;j<iba.n-1;j++)
	{
		r->c[j+1] = e.c[j];
	}
	r->c[0] = 0;
	//ajout de e.c[n-1]*b^n + iba.c[i]
	for (j=0;j<iba.n;j++)
	{
		r->c[j] += e.c[iba.n-1]*iba.bn.c[j] + iba.c[i].c[j];
	}
}

//évalue l'élément dans la place p
Complexe eval (Element e, int p)
{
	Complexe c = zero();
	int i;
	for (i=0;i<iba.n;i++)
	{
		addOP(&c, mul_i(iba.p[p].c[i], e.c[i]));
	}
	return c;
}

//détermine si un élément doit être gardé
//en l'évaluant dans les différentes places et en comparant aux valeurs max autorisées
bool keep (Element e)
{
	int i;
	for (i=0;i<iba.na;i++)
	{
		if (cnorm(eval(e, i)) >= iba.cM[i])
			return false;
	}
}

//haché d'un élément
int hash (Element e)
{
	int i, j;
	int h = 1;
	for (i=0;i<iba.n;i++)
	{
		j = e.c[i];
		do
		{
			h = (h*256 + j%256)%nhash;
			j /= 256;
		}while(j);
	}
	return h;
}

Element zeroElement ()
{
	Element e = NewElement(iba.n);
	int i;
	for (i=0;i<iba.n;i++)
	{
		e.c[i] = 0;
	}
	return e;
}

bool equalsElements (Element e, Element f)
{
	int i;
	for (i=0;i<iba.n;i++)
	{
		if (e.c[i] != f.c[i])
			return false;
	}
	return true;
}

void copy (Element f, Element d)
{
	int i;
	for (i=0;i<iba.n;i++)
	{
		d.c[i] = f.c[i];
	}
}

int compteur = 0; //nombre d'états de l'automate
bool vide = 1;

bool inHash (Element e)
{
	int h = hash(e);
	int i;
	for (i=0;i<thash[h].n;i++)
	{
		if (equalsElements(e, thash[h].e[i]))
		{
			//if (!thash[h].ne[i])
			//	compteur++;
			vide = false; //l'automate n'est pas trivial
			//thash[h].ne[i]++; //note que l'état a été revu
			return true;
		}
	}
	//ajoute l'élément
	thash[h].e = (Element *)realloc(thash[h].e, thash[h].n+1);
	//thash[h].ne = (int *)realloc(thash[h].ne, thash[h].n+1);
	thash[h].ind = (int *)realloc(thash[h].ind, thash[h].n+1);
	thash[h].e[thash[h].n] = NewElement(iba.n);
	copy(e, thash[h].e[thash[h].n]);
	//thash[h].ne = 0;
	thash[h].ind = compteur;
	thash[h].n++;
	return false;
}

//trouve l'indice d'un élément
int indElement (Element e)
{
	int h = hash(e);
	int i;
	for (i=0;i<thash[h].n;i++)
	{
		if (equalsElements(e, thash[h].e[i]))
		{
			return thash[h].ind[i];
		}
	}
	return -1;
}

void printElement (Element e)
{
	int i;
	for (i=0;i<iba.n;i++)
	{
		printf("%d ", e.c[i]);
	}
}

//calcule l'automate des relations
Automate RelationsAutomaton (InfoBetaAdic iba2, bool isvide, bool verb)
{
	iba = iba2;
	
	////afficher les données : chiffres, places, bornes, etc... pour vérif
	
	if (verb)
		printf("init hash...\n");
	//table de hachage servant à repérer les éléments déjà rencontrés
	initHash();
	
	if (verb)
		printf("parcours...\n");
	int n = 1; //nombre d'éléments sur la pile
	compteur = 1; //nombre d'états de l'automate
	int i,j;
	pile[0] = zeroElement();
	Element e = NewElement(iba.n);
	Element s = NewElement(iba.n); //fils
	vide = true;
	while (n)
	{
		//parcours le dernier élément mis sur la pile
		n--; //dépile
		copy(pile[n], e);
		FreeElement(pile[n]);
		if (verb)
		{
			printf("état ");
			printElement(e);
			printf("vu\n");
		}
		for (i=0;i<iba.na;i++)
		{
			succ(e, i, &s);
			if (verb)
			{
				printf("succ %d : ", i);
				printElement(e);
				printf("\n");
			}
			if (keep(s))
			{ //l'élément est dans l'automate
				//teste si l'élément a déjà été vu et l'ajoute si non
				if (inHash(s))
				{ //l'élement est nouveau et a été ajouté à la table de hachage
					if (!vide)
					{ //l'automate n'est pas trivial
						return NewAutomaton(1,0);
					}
					//empile
					pile[n] = NewElement(iba.n);
					copy(e, pile[n]);
					n++;
					compteur++;
					if (n > npile)
					{
						printf("Erreur : dépassement de la pile !!!\n");
					}
				}
			}
		}
	}
	if (verb)
		printf("..fini !\n");
	FreeElement(e);
	
	if (verb)
		printf("%d états rencontrés.\n", compteur);
	
	//créé l'automate
	Automate r = NewAutomaton(compteur, iba.nc);
	for (i=0;i<r.n;i++)
	{
		for (j=0;j<r.na;j++)
		{
			r.e[i].f[j] = -1;
		}
	}
	int k, ind;
	for (i=0;i<nhash;i++)
	{
		for (j=0;j<thash[i].n;j++)
		{
			e = thash[i].e[j];
			r.e[thash[i].ind[j]].final = false;
			for (k=0;k<iba.na;k++)
			{
				succ(e, k, &s);
				ind = indElement(s);
				if (ind != -1)
				{ //ajoute la transition
					r.e[thash[i].ind[j]].f[k] = ind;
				}
			}
			FreeElement(e);
		}
	}
	if (verb)
		printf("free...\n");
	FreeElement(s);
	//états initiaux et finaux : zéro
	e = zeroElement();
	ind = indElement(e);
	FreeElement(e);
	r.i = ind;
	r.e[ind].final = true;
	return r;
}


