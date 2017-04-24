
typedef char bool;
#define true 1
#define false 0

struct Etat
{
	int *f; //liste des na fils
	bool final;
};
typedef struct Etat Etat;

struct Automate
{
	Etat *e; //états
	int n; //nombre d'états
	int na; //nombre de lettres
	int i; //état initial
	//
//	int nalloc; //usage interne, mémoire allouée
};
typedef struct Automate Automate;
typedef struct Automate Automaton;

//Automaton NewAutomaton (int n, int na);
//void FreeAutomaton (Automaton *a);

///////////////////////////////////////////////////////
//Automates non déterministe
///////////////////////////////////////////////////////

struct Arete
{
	int l; //label (-1 : epsilon-transition)
	int e; //état d'arrivée
};
typedef struct Arete Arete;

struct NEtat
{
	Arete *a;
	int n;
	bool final;
	bool initial;
};
typedef struct NEtat NEtat;

struct NAutomate
{
	NEtat *e; //états
	int n; //nombre d'états
	int na; //nombre de lettres
};
typedef struct NAutomate NAutomate;
