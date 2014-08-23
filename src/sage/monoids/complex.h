
struct Complexe
{
	double x,y;
};
typedef struct Complexe Complexe;

Complexe prod (Complexe a, Complexe b);
Complexe mul_i (Complexe a, int i);
Complexe zero ();
Complexe un ();
Complexe add (Complexe a, Complexe b);
void addOP (Complexe *a, Complexe b);
inline double carre (double x);
double cnorm (Complexe c);
Complexe inv (Complexe c);
