
struct Element
{
	long int *c; //liste des n coeffs
};
typedef struct Element Element;

struct PlaceArch
{
	Complexe *c; //1, b, b^2, ... pour cette place
};
typedef struct PlaceArch PlaceArch;

Element NewElement (int n);
void FreeElement (Element e);

//structure contenant les infos nécessaires pour calculer l'automate des relations
struct InfoBetaAdic
{
	int n; //degré
	Element bn; //expression de b^n comme un polynome en b de degré < n
	Element b1; //expression de 1/b comme un polynôme en b de degré < n
	Element *c; //liste des chiffres utilisés pour le calcul de l'automate des relations
	int nc; //nombre de chiffre
	int ncmax; //nombre de chiffres alloués
	PlaceArch *p; //liste des na places
	double *cM; //carré des valeurs absolues max
	int na; //nombre de places
};
typedef struct InfoBetaAdic InfoBetaAdic;

InfoBetaAdic allocInfoBetaAdic (int n, int na, int ncmax, bool verb);
void freeInfoBetaAdic (InfoBetaAdic iba);

//calcule l'automate des relations
Automate RelationsAutomaton (InfoBetaAdic iba2, bool isvide, bool ext, bool verb);

//calcule l'automate des relations avec translation
Automate RelationsAutomatonT (InfoBetaAdic iba2, Element t, bool isvide, bool ext, bool verb);

