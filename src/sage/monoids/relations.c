#include <stdlib.h>
#include "complex.h"
#include "Automaton.h"
#include "relations.h"
//#include "../combinat/words/automataC.h"

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
	iba.b1 = NewElement(n);
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
	free(pile);
	free(thash);
	free(iba.cM);
	for (i=0;i<iba.ncmax;i++)
	{
		FreeElement(iba.c[i]);
	}
	FreeElement(iba.bn);
	FreeElement(iba.b1);
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

//calcul le fils de l'élément e pour le chiffre i (autre sens)
void succ2 (Element e, int i, Element *r)
{
	//printf("succ2 : ");
	//printElement(e);
	//printf("\n");
	int j;
	//retranche iba.c[i]
	for (j=0;j<iba.n;j++)
	{
		r->c[j] = e.c[j] - iba.c[i].c[j];
	}
	//printElement(*r);
	//printf("\n");
	//division par b
	int r0 = r->c[0]; //retient le coefficient de degré 0
	for (j=0;j<iba.n-1;j++)
	{
		r->c[j] = r->c[j+1];
	}
	r->c[iba.n-1] = 0;
	//printElement(*r);
	//printf("\n");
	//ajout de r0*(1/b)
	for (j=iba.n-1;j>=0;j--)
	{
		r->c[j] += r0*iba.b1.c[j];
	}
	//printElement(*r);
	//printf("\n");
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
		if (cnorm(eval(e, i)) - .0000001 > iba.cM[i])
			return false;
	}
	return true;
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

bool isNull (Element e)
{
	int i;
	for (i=0;i<iba.n;i++)
	{
		if (e.c[i] != 0)
			return false;
	}
	return true;
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

void copy (const Element f, Element d)
{
	int i;
	for (i=0;i<iba.n;i++)
	{
		d.c[i] = f.c[i];
	}
}

void printElement (Element e)
{
	int i;
	for (i=0;i<iba.n;i++)
	{
		printf("%d ", e.c[i]);
	}
}

int compteur = 0; //nombre d'états de l'automate
bool vide = 1;

bool inHash (Element e)
{
	int h = hash(e);
	/*
	printf("hash ");
	printElement(e);
	printf(": %d\n", h);
	*/
	int i;
	for (i=0;i<thash[h].n;i++)
	{
		/*
		printf("equals ");
		printElement(e);
		printElement(thash[h].e[i]);
		printf("...\n");
		*/
		if (equalsElements(e, thash[h].e[i]))
		{
			//printf(" -> equals !\n");
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
	if (!thash[h].e || !thash[h].ind)
	{
		printf("Out of memory !!!\n");
		exit(-2);
	}
	thash[h].e[thash[h].n] = NewElement(iba.n);
	copy(e, thash[h].e[thash[h].n]);
	//thash[h].ne = 0;
	thash[h].ind[thash[h].n] = compteur;
	thash[h].n++;
	compteur++;
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

//calcule l'automate des relations
Automate RelationsAutomaton (InfoBetaAdic iba2, bool isvide, bool ext, bool verb)
{
	int i,j;
	
	iba = iba2;
	
	////affiche les données : chiffres, places, bornes pour vérif
	if (verb)
	{
		for (i=0;i<iba.nc;i++)
		{
			printf("chiffre %d : ", i);
			printElement(iba.c[i]);
			printf("\n");
		}
		for (i=0;i<iba.na;i++)
		{
			printf("place %d : ", i);
			for (j=0;j<iba.n;j++)
			{
				printf("(%lf, %lf) ", iba.p[i].c[j].x, iba.p[i].c[j].y);
			}
			printf("borne %lf\n", iba.cM[i]);
		}
	}

	if (verb)
		printf("init hash...\n");
	//table de hachage servant à repérer les éléments déjà rencontrés
	initHash();
	
	if (verb)
		printf("parcours...\n");
	int n = 1; //nombre d'éléments sur la pile
	compteur = 0; //nombre d'états de l'automate
	//état initial 
	pile[0] = zeroElement();
	inHash(pile[0]); //ajoute l'élément à la table de hachage
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
		for (i=0;i<iba.nc;i++)
		{
			succ(e, i, &s);
			if (verb)
			{
				printf("succ %d/%d : ", i, iba.nc);
				printElement(s);
				printf("\n");
			}
			if (keep(s))
			{ //l'élément est dans l'automate
				//teste si l'élément a déjà été vu et l'ajoute si non
				if (!inHash(s))
				{ //l'élement est nouveau et a été ajouté à la table de hachage
					/*
					if (isvide && !vide)
					{ //l'automate n'est pas trivial
						return NewAutomaton(1,0);
					}
					*/
					//empile
					pile[n] = NewElement(iba.n);
					copy(s, pile[n]);
					n++;
					if (n > npile)
					{
						printf("Erreur : dépassement de la pile !!!\n");
					}
				}else
				{//on retombe sur un état déjà vu
					if (isvide)
					{
						//if (!isNull(e) && 
						if (ext || isNull(s))
						{ //l'automate émondé inf ou émondé n'est pas vide
							return NewAutomaton(1,0);
						}
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
		if (verb)
		{
			if (thash[i].n > 0)
				printf("hash %d : %d éléments.\n", i, thash[i].n);
		}
		for (j=0;j<thash[i].n;j++)
		{
			e = thash[i].e[j];
			if (verb)
			{
				printf("Element ");
				printElement(e);
				printf("indice %d\n", thash[i].ind[j]);			
			}
			r.e[thash[i].ind[j]].final = false;
			for (k=0;k<iba.nc;k++)
			{
				succ(e, k, &s);
				ind = indElement(s);
				/*
				printf("indice de ");
				printElement(s);
				printf(": %d\n", ind);
				*/
				if (ind != -1)
				{ //ajoute la transition
					r.e[thash[i].ind[j]].f[k] = ind;
				}
			}
		}
	}
	//états initiaux et finaux : zéro
	e = zeroElement();
	ind = indElement(e);
	r.i = ind;
	r.e[ind].final = true;
	if (verb)
		printf("free...\n");
	//libère les éléments de la table de hachage
	for (i=0;i<nhash;i++)
	{
		for (j=0;j<thash[i].n;j++)
		{
			FreeElement(thash[i].e[j]);
		}
	}
	FreeElement(s);
	FreeElement(e);
	return r;
}

//calcule l'automate des relations avec translation
Automate RelationsAutomatonT (InfoBetaAdic iba2, Element t, bool isvide, bool ext, bool verb)
{
	int i,j;
	
	iba = iba2;
	
	////affiche les données : chiffres, places, bornes pour vérif
	if (verb)
	{
		printf("isvide = %d\n", isvide);
		printf("translation : ");
		printElement(t);
		printf("\n");
		for (i=0;i<iba.nc;i++)
		{
			printf("chiffre %d : ", i);
			printElement(iba.c[i]);
			printf("\n");
		}
		for (i=0;i<iba.na;i++)
		{
			printf("place %d : ", i);
			for (j=0;j<iba.n;j++)
			{
				printf("(%lf, %lf) ", iba.p[i].c[j].x, iba.p[i].c[j].y);
			}
			printf("borne %lf\n", iba.cM[i]);
		}
	}

/*	
	//teste si la translation n'est pas trop grande
	if (!keep(t))
	{
		return NewAutomaton(0,0);
	}
*/

	if (verb)
		printf("init hash...\n");
	//table de hachage servant à repérer les éléments déjà rencontrés
	initHash();
	
	if (verb)
		printf("parcours...\n");
	int n = 1; //nombre d'éléments sur la pile
	compteur = 0; //nombre d'états de l'automate
	//état initial 
	pile[0] = NewElement(iba.n);
	copy(t, pile[0]);
	inHash(pile[0]); //ajoute l'élément à la table de hachage
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
		for (i=0;i<iba.nc;i++)
		{
			succ2(e, i, &s);
			if (verb)
			{
				printf("succ %d/%d : ", i, iba.nc);
				printElement(s);
				printf("\n");
			}
			if (keep(s))
			{ //l'élément est dans l'automate
				//teste si l'élément a déjà été vu et l'ajoute si non
				if (!inHash(s))
				{ //l'élement est nouveau et a été ajouté à la table de hachage
					/*
					if (isvide && !vide)
					{ //l'automate n'est pas trivial
						return NewAutomaton(1,0);
					}
					*/
					//empile
					pile[n] = NewElement(iba.n);
					copy(s, pile[n]);
					n++;
					if (n > npile)
					{
						printf("Erreur : dépassement de la pile !!!\n");
					}
				}else
				{//on retombe sur un état déjà vu
					if (isvide && ext)
					{
						//l'automate émondé inf n'est pas vide
						return NewAutomaton(1,0);
					}
				}
				if (isvide && isNull(s))
				{ //l'automate émondé n'est pas vide
					return NewAutomaton(1,0);
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
		if (verb)
		{
			if (thash[i].n > 0)
				printf("hash %d : %d éléments.\n", i, thash[i].n);
		}
		for (j=0;j<thash[i].n;j++)
		{
			e = thash[i].e[j];
			if (verb)
			{
				printf("Element ");
				printElement(e);
				printf("indice %d\n", thash[i].ind[j]);			
			}
			r.e[thash[i].ind[j]].final = false;
			for (k=0;k<iba.nc;k++)
			{
				succ2(e, k, &s);
				ind = indElement(s);
				/*
				printf("indice de ");
				printElement(s);
				printf(": %d\n", ind);
				*/
				if (ind != -1)
				{ //ajoute la transition
					r.e[thash[i].ind[j]].f[k] = ind;
				}
			}
		}
	}
	//états initiaux et finaux
	r.i = indElement(t);
	e = zeroElement();
	ind = indElement(e);
	r.e[ind].final = true;
	if (verb)
		printf("free...\n");
	//libère les éléments de la table de hachage
	for (i=0;i<nhash;i++)
	{
		for (j=0;j<thash[i].n;j++)
		{
			FreeElement(thash[i].e[j]);
		}
	}
	FreeElement(s);
	FreeElement(e);
	return r;
}


