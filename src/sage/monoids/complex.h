
struct Complexe
{
	double x,y;
};
typedef struct Complexe Complexe;

Complexe prod (Complexe a, Complexe b)
{
	Complexe r;
	r.x = a.x*b.x - a.y*b.y;
	r.y = a.x*b.y + a.y*b.x;
	return r;
}

Complexe zero ()
{
	Complexe r;
	r.x = 0;
	r.y = 0;
	return r;
}

Complexe un ()
{
	Complexe r;
	r.x = 1;
	r.y = 0;
	return r;
}

Complexe add (Complexe a, Complexe b)
{
	Complexe r;
	r.x = a.x + b.x;
	r.y = a.y + b.y;
	return r;
}

inline double carre (double x)
{
	return x*x;
}

double cnorm (Complexe c)
{
	return carre(c.x) + carre(c.y);
}

Complexe inv (Complexe c)
{
	Complexe r;
	double cn = cnorm(c);
	r.x = c.x/cn;
	r.y = -c.y/cn;
	return r;
}
