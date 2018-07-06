
typedef int bool;
#define true 1
#define false 0

struct State
{
	int *f; //list of na sons
	bool final;
};
typedef struct State State;

struct Automaton
{
	State *e; //states
	int n; //number of states
	int na; //number of letters
	int i; //initial state
	//
//	int nalloc; //internal usage, allocated memory
};
typedef struct Automaton Automaton;
typedef struct Automaton Automaton;

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

struct NAutomaton
{
	NState *e; //states
	int n; //number of states
	int na; //number of letters
};
typedef struct NAutomaton NAutomaton;
