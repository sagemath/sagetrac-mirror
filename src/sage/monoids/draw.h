
typedef unsigned char uint8;

struct Color
{
	uint8	r,g,b,a;
};
typedef struct Color Color;

typedef struct Color* ColorList;

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

ColorList NewColorList (int n);
void FreeColorList (ColorList l);
Surface NewSurface (int sx, int sy);
void FreeSurface (Surface s);
void Fill (Surface s, Color c);
Automate NewAutomate (int n, int na);
void FreeAutomate (Automate a);
void FreeAutomates (Automate* a, int n);
BetaAdic NewBetaAdic (int n);
void FreeBetaAdic (BetaAdic b);
BetaAdic2 NewBetaAdic2 (int n, int na);
void FreeBetaAdic2 (BetaAdic2 b);
inline Color moy (Color a, Color b, double ratio);
void set_pix (Surface s, Complexe p);
void print_word (BetaAdic b, int n, int etat);
Color randColor (int a);
void Draw (BetaAdic b, Surface s, int n, int ajust, Color col, int verb);
void Draw2 (BetaAdic b, Surface s, int n, int ajust, Color col, int verb);
void printAutomate (Automate a);
void DrawList (BetaAdic2 b, Surface s, int n, int ajust, ColorList cl, double alpha, int verb);
