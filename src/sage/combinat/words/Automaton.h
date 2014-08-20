
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
