#include <stdint.h>

typedef uint64_t uint64;
typedef unsigned int uint;

typedef Automate Automaton;
typedef NAutomate NAutomaton;

struct Dict
{
	int *e;
	int n;
};
typedef struct Dict Dict;

bool DotExists ();

int hashAutomaton (Automaton a);
bool findWord (Automaton a, Dict *w, bool verb); //return a word of the language of a
bool rec_word (Automaton a, Dict d); //check that the word w is recognized by the automaton a
bool shortestWord (Automaton a, Dict *w, int i, int f, bool verb); //return a shortest word of the language of a
bool shortestWords (Automaton a, Dict *w, int i, bool verb); //return shortest words toward each state
Dict NewDict (int n);
void FreeDict (Dict *d);
void printDict (Dict d);
void dictAdd (Dict *d, int e); //add the element to the disctionnary (even if already present)
Automaton NewAutomaton (int n, int na);
void ReallocNAutomaton (NAutomaton *a, int n);
void FreeAutomaton (Automaton *a);
void FreeAutomates (Automate* a, int n);
void FreeNAutomaton (NAutomaton *a);
void AddEdgeN (NAutomaton *a, int e, int f, int l);
void AddPathN (NAutomaton *a, int e, int f, int *l, int len, bool verb);
Automaton CopyAutomaton (Automaton a, int nalloc, int naalloc);
Automaton PieceAutomaton (Automaton a, int *w, int n, int e); //gives an automaton recognizing w(w^(-1)L) where L is the language of a starting from e
void init (Automaton *a);
void printAutomaton (Automaton a);
void plotDot (const char *file, Automaton a, const char **labels, const char *graph_name, double sx, double sy, const char **vlabels, bool html, bool verb, bool run_dot);
bool equalsAutomaton (Automaton a1, Automaton a2); //determine if automata are the same (differents if permuted states)
int contract (int i1, int i2, int n1);
int geti1 (int c, int n1);
int geti2 (int c, int n1);
Automaton Product (Automaton a1, Automaton a2, Dict d, bool verb);
void AddEtat (Automaton *a, bool final);

struct Etats
{
	int *e;
	int n;	
};
typedef struct Etats Etats;

Etats NewEtats (int n);
void FreeEtats (Etats e);
void initEtats (Etats e);
void printEtats (Etats e);
bool equals (Etats e1, Etats e2);
Etats copyEtats (Etats e);

struct ListEtats
{
	Etats *e;
	int n;
};
typedef struct ListEtats ListEtats;

void printListEtats (ListEtats l);
bool AddEl (ListEtats *l, Etats e, int* res); //add an element if not already in the list
void AddEl2 (ListEtats *l, Etats e); //add an element even if already in the list

////////////////
struct Etats2
{
	uint n;
	uint64 *e;
};
typedef struct Etats2 Etats2;

Etats2 NewEtats2 (int n);
void FreeEtats2 (Etats2 e);
void initEtats2 (Etats2 e);
void printEtats2 (Etats2 e);
bool isNullEtats2 (Etats2 e);
bool equalsEtats2 (Etats2 e1, Etats2 e2);
bool hasEtats2 (Etats2 e, uint64 i);
Etats2 copyEtats2 (Etats2 e);
void addEtat (Etats2 *e, uint64 i);

struct ListEtats2
{
	Etats2 *e;
	int n; //number of states
	int na; //memory allocated
};
typedef struct ListEtats2 ListEtats2;

ListEtats2 NewListEtats2(int n, int na);
void ReallocListEtats2(ListEtats2* l, int n, bool marge);
void FreeListEtats2 (ListEtats2* l);
void printListEtats2 (ListEtats2 l);

//inverse of a dictionnary
struct InvertDict
{
	Dict *d;
	int n;
};
typedef struct InvertDict InvertDict;

InvertDict NewInvertDict (int n);
InvertDict invertDict (Dict d);
void FreeInvertDict (InvertDict id);
void printInvertDict (InvertDict id);
void putEtat (Etats *f, int ef); /////////////to improve !!!!
void Determinize_rec (Automaton a, InvertDict id, Automaton* r, ListEtats* l, bool onlyfinals, bool nof, int niter);
Automaton Determinize (Automaton a, Dict d, bool noempty, bool onlyfinals, bool nof, bool verb);
NAutomaton Concat (Automaton a, Automaton b, bool verb);
NAutomaton CopyN (Automaton a, bool verb);
NAutomaton Proj (Automaton a, Dict d, bool verb);

Automaton DeterminizeN (NAutomaton a, bool puits, int verb);

//change the alphabet, duplicating edges if necessary
//the result is assumed deterministic !!!!
Automaton Duplicate (Automaton a, InvertDict id, int na2, bool verb);

//add all the words that are in the language if we remove some ending zeroes
//zero is the letter of the alphabet of index l0
//i.e. the result has the language L(l0*)^(-1), if L is the language of a
void ZeroComplete (Automaton *a, int l0, bool verb);

//add all the words that can be completed to a word of the language by adding some ending zeroes
//zero is the letter of index l0
//i.e. the result has the language L(l0*), if L is tha language of a
Automaton ZeroComplete2 (Automaton *a, int l0, bool etat_puits, bool verb);

//Compute an automaton recognizing the language (l0*)L, where L is the language of a
Automaton ZeroInv (Automaton *a, int l0);

//retire tous les états à partir desquels il n'y a pas de chemin infini
Automaton prune_inf (Automaton a, bool verb);

//Compute the transposition, assuming it is deterministic
Automaton TransposeDet (Automaton a);

//Compute the transposition
NAutomaton Transpose (Automaton a);

//Tarjan algorithm
int StronglyConnectedComponents (Automaton a, int *res);

//determine accessible and co-accessible states
void AccCoAcc (Automaton *a, int *coa);

//determine co-accessible states
void CoAcc (Automaton *a, int *coa);

//retire tous les états non accessible ou non co-accessible
Automaton prune (Automaton a, bool verb);

//retire tous les états non accessible
Automaton pruneI (Automaton a, bool verb);

Automaton SubAutomaton (Automaton a, Dict d, bool verb);

//permute les labels des arêtes
//l donne les anciens indices à partir des nouveaux
Automaton Permut (Automaton a, int *l, int na, bool verb);
//idem mais SUR PLACE
void PermutOP (Automaton a, int *l, int na, bool verb);

//minimisation par l'algo d'Hopcroft
//voir "Around Hopcroft’s Algorithm" de Manuel BACLET and Claire PAGETTI
Automaton Minimise (Automaton a, bool verb);

void DeleteVertexOP (Automaton *a, int e);
Automaton DeleteVertex (Automaton a, int e);

//détermine si les langages des automates sont les mêmes
//le dictionnaires donne les lettres de a2 en fonction de celles de a1 (-1 si la lettre de a1 ne correspond à aucune lettre de a2). Ce dictionnaire est supposé inversible.
//if minimized is true, the automaton a1 and a2 are assumed to be minimal.
bool equalsLanguages (Automaton *a1, Automaton *a2, Dict a1toa2, bool minimized, bool pruned, bool verb);

//détermine si le langage de l'automate a1 est inclus dans celui de a2
//le dictionnaires donne les lettres de a2 en fonction de celles de a1 (-1 si la lettre de a1 ne correspond à aucune lettre de a2). Ce dictionnaire est supposé inversible.
//if pruned is true, the automaton a1 and a2 are assumed to be pruned.


//détermine si les langages des automates ont une intersection non vide
//le dictionnaires donne les lettres de a2 en fonction de celles de a1 (-1 si la lettre de a1 ne correspond à aucune lettre de a2). Ce dictionnaire est supposé inversible.
//if pruned is true, the automaton a1 and a2 are assumed to be pruned otherwise it prunes.
//bool intersectLanguage (Automaton *a1, Automaton *a2, Dict a1toa2, bool pruned, bool verb);

//détermine si l'intersection est vide ou non
bool Intersect (Automaton a1, Automaton a2, bool verb);

//détermine si l'on a inclusion des langages
bool Included (Automaton a1, Automaton a2, bool pruned, bool verb);

//détermine si le langage de l'automate est vide
bool emptyLanguage (Automaton a);

//determine if the automaton is complete (i.e. with his hole state)
bool IsCompleteAutomaton (Automaton a);

//complete the automaton (i.e. add a hole state if necessary)
bool CompleteAutomaton (Automaton *a);

//copy the automaton with a new bigger alphabet
Automaton BiggerAlphabet (Automaton a, Dict d, int nna);


