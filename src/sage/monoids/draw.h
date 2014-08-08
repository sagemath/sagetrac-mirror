
//typedef int *Etat;
struct Etat
{
	int* f; //liste des na fils indexés par les lettres
	int final; //l'état est final ?
};
typedef struct Etat Etat;

struct Automate
{
	Etat* e; //liste des états (l'état initial est le premier)
	int n; //nombre d'états
	int na; //nombre de lettres différentes
	int i; //état initial
};
typedef struct Automate Automate;

typedef unsigned char uint8;

struct Color
{
	uint8	r,g,b,a;
};
typedef struct Color Color;

struct Surface
{
	Color** pix;
	int sx, sy; //dimensions de l'image
};
typedef struct Surface Surface;

struct BetaAdic
{
	Complexe b; //
	Complexe *t; //liste des translations
	int n; //nombre de translations
	Automate a;
};
typedef struct BetaAdic BetaAdic;

struct BetaAdic2
{
	Complexe b; //
	Complexe *t; //liste des translations
	int n; //nombre de translations
	Automate *a; //liste des couleurs
	int na; //nombre d'automates
};
typedef struct BetaAdic2 BetaAdic2;
