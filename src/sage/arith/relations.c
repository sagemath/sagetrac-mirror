#include <stdio.h>
#include <stdlib.h>
#include "complex.h"
#include "Automaton.h"
#include "relations.h"
#include "automataC.h"
//#include "../combinat/words/automataC.h"

#define MEM_DEBUG	yes

//#define nhash 1000003

InfoBetaAdic iba; //variable globale contenant les infos sur le développement en base bêta nécessaires au calcul de l'automate des relations

Element NewElement (int n)
{
	//printf("NewElement(%d)\n", n);
	Element e;
	e.c = (coeff *)malloc(sizeof(coeff)*n);
	return e;
}

void FreeElement (Element e)
{
	free(e.c);
}

InfoBetaAdic allocInfoBetaAdic(int n, int na, int ncmax, int nhash, bool verb)
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
	
	iba.nhash = nhash;
	iba.thash = (Elements*)malloc(sizeof(Elements)*nhash);
	
	iba.npile = 256; //1000000;
	iba.pile = (Element *)malloc(sizeof(Element)*iba.npile);
	return iba;
}

void initHash (InfoBetaAdic *iba)
{
	int i;
	for (i=0;i<iba->nhash;i++)
	{
		//thash[i].e = NULL;
		iba->thash[i].n = 0;
	}
}

void freeInfoBetaAdic (InfoBetaAdic *iba)
{
	int i,j;
	free(iba->pile);
	for (i=0;i<iba->nhash;i++)
	{
		for (j=0;j<iba->thash[i].n;j++)
		{
			FreeElement(iba->thash[i].e[j]);
		}
		if (iba->thash[i].n)
		{
			free(iba->thash[i].e);
			free(iba->thash[i].ind);
		}
	}
	free(iba->thash);
	free(iba->cM);
	for (i=0;i<iba->ncmax;i++)
	{
		FreeElement(iba->c[i]);
	}
	FreeElement(iba->bn);
	FreeElement(iba->b1);
	free(iba->c);
	free(iba->p);
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
			h = (h*256 + j%256)%iba.nhash;
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
		printf("%ld ", e.c[i]);
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
	for (i=0;i<iba.thash[h].n;i++)
	{
		/*
		printf("equals ");
		printElement(e);
		printElement(thash[h].e[i]);
		printf("...\n");
		*/
		if (equalsElements(e, iba.thash[h].e[i]))
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
	if (iba.thash[h].n)
	{
		iba.thash[h].e = (Element *)realloc(iba.thash[h].e, sizeof(Element)*(iba.thash[h].n+1));
		//thash[h].ne = (int *)realloc(thash[h].ne, thash[h].n+1);
		iba.thash[h].ind = (int *)realloc(iba.thash[h].ind, sizeof(int)*(iba.thash[h].n+1));
	}else
	{
		iba.thash[h].e = (Element *)malloc(sizeof(Element));
		//thash[h].ne = (int *)realloc(thash[h].ne, thash[h].n+1);
		iba.thash[h].ind = (int *)malloc(sizeof(int));
	}
	if (!iba.thash[h].e || !iba.thash[h].ind)
	{
		printf("Out of memory !!!\n");
		exit(-2);
	}
	iba.thash[h].e[iba.thash[h].n] = NewElement(iba.n);
	copy(e, iba.thash[h].e[iba.thash[h].n]);
	//thash[h].ne = 0;
	iba.thash[h].ind[iba.thash[h].n] = compteur;
	iba.thash[h].n++;
	compteur++;
	return false;
}

//trouve l'indice d'un élément
int indElement (Element e)
{
	int h = hash(e);
	int i;
	for (i=0;i<iba.thash[h].n;i++)
	{
		if (equalsElements(e, iba.thash[h].e[i]))
		{
			return iba.thash[h].ind[i];
		}
	}
	return -1;
}

//calcule l'automate des relations avec translation
Automaton RelationsAutomatonT(InfoBetaAdic *iba2, Element t, bool isvide, bool ext, bool verb)
{
	int i,j;
	
	iba = *iba2;
	
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
	initHash(&iba);
	
	if (verb)
		printf("parcours...\n");
	int n = 1; //nombre d'éléments sur la pile
	compteur = 0; //nombre d'états de l'automate
	//état initial 
	iba.pile[0] = NewElement(iba.n);
	copy(t, iba.pile[0]);
	inHash(iba.pile[0]); //ajoute l'élément à la table de hachage
	Element e = NewElement(iba.n);
	Element s = NewElement(iba.n); //fils
	vide = true;
	while (n)
	{
		//parcours le dernier élément mis sur la pile
		n--; //dépile
		copy(iba.pile[n], e);
		FreeElement(iba.pile[n]);
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
					iba.pile[n] = NewElement(iba.n);
					copy(s, iba.pile[n]);
					n++;
					if (n >= iba.npile)
					{
						iba.npile *= 2;
						iba.pile = (Element *)realloc(iba.pile, sizeof(Element)*iba.npile);
						if (!iba.pile)
						{
						    printf("Error: failed to reallocate the stack (size %d).\n", iba.npile);
						    exit(EXIT_FAILURE);
						}
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
	
	if (verb)
		printf("NewAutomaton(%d, %d)\n", compteur, iba.nc);
	//créé l'automate
	Automaton r = NewAutomaton(compteur, iba.nc);
	int k, ind;
	for (i=0;i<iba.nhash;i++)
	{
		if (verb)
		{
			if (iba.thash[i].n > 0)
				printf("hash %d : %d éléments.\n", i, iba.thash[i].n);
		}
		for (j=0;j<iba.thash[i].n;j++)
		{
			e = iba.thash[i].e[j];
			if (verb)
			{
				printf("Element ");
				printElement(e);
				printf("indice %d\n", iba.thash[i].ind[j]);			
			}
			r.e[iba.thash[i].ind[j]].final = false;
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
					r.e[iba.thash[i].ind[j]].f[k] = ind;
				}
			}
		}
	}
	//états initiaux et finaux
	r.i = indElement(t);
	e = zeroElement();
	ind = indElement(e);
	if (ind != -1)
		r.e[ind].final = true;
	else
	{
		//Automate émondé trivial
	}
	if (verb)
		printf("free...\n");
	//libère les éléments de la table de hachage
	for (i=0;i<iba.nhash;i++)
	{
		for (j=0;j<iba.thash[i].n;j++)
		{
			FreeElement(iba.thash[i].e[j]);
		}
		if (iba.thash[i].n)
		{
			free(iba.thash[i].e);
			free(iba.thash[i].ind);
		}
		iba.thash[i].n = 0;
	}
	//if (verb)
	//	printf("..free..\n");
	FreeElement(s);
	FreeElement(e);
	//if (verb)
	//	printf("done.\n");
	return r;
}


