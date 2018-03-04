
struct Complexe
{
	double x,y;
};
typedef struct Complexe Complexe;

Complexe prod (Complexe a, Complexe b);
Complexe mul_i (Complexe a, int i);
Complexe zero ();
Complexe un ();
Complexe powC (Complexe a, int n);
Complexe sub (Complexe a, Complexe b);
Complexe add (Complexe a, Complexe b);
void addOP (Complexe *a, Complexe b);
double carre (double x);
double cnorm (Complexe c);
Complexe inv (Complexe c);
