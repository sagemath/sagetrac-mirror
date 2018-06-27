
typedef int bool;
#define true 1
#define false 0

struct State
{
	int *f; //list of na sons
	bool final;
};
typedef struct State State;

struct Automate
{
	State *e; //states
	int n; //number of states
	int na; //number of letters
	int i; //initial state
	//
//	int nalloc; //internal usage, allocated memory
};
typedef struct Automate Automate;
typedef struct Automate Automaton;

//Automaton NewAutomaton (int n, int na);
//void FreeAutomaton (Automaton *a);

///////////////////////////////////////////////////////
//Non Deterministaic Automata
///////////////////////////////////////////////////////

struct Transition
{
	int l; //label (-1 : epsilon-transition)
	int e; //arrival state
};
typedef struct Transition Transition;

struct NState
{
	Transition *a;
	int n;
	bool final;
	bool initial;
};
typedef struct NState NState;

struct NAutomate
{
	NState *e; //states
	int n; //number of states
	int na; //number of letters
};
typedef struct NAutomate NAutomate;
