#include <stdlib.h>
#include "complex.h"
#include "draw.h"

double mx = -2, my = -2, Mx = 2, My = 2; //zone de dessin

double mx2 = 1000000, my2 = 1000000, Mx2 = -1000000, My2 = -1000000; //extremum observés

Color color0; //couleur du fond
Color color; //couleur de dessin
Color* colors; //liste de couleurs de dessin

Surface NewSurface (int sx, int sy)
{
	Surface s;
	s.sx = sx;
	s.sy = sy;
	s.pix = (Color **)malloc(sizeof(Color *)*sx);
	int x;
	for (x=0;x<sx;x++)
	{
		s.pix[x] = (Color *)malloc(sizeof(Color)*sy);
	}
	return s;
}

void FreeSurface (Surface s)
{
	int x;
	for (x=0;x<s.sx;x++)
	{
		free(s.pix[x]);
	}
	free(s.pix);
}

void Fill (Surface s, Color c)
{
	int x,y;
	for (y=0;y<s.sy;y++)
	{
		for (x=0;x<s.sx;x++)
		{
			s.pix[x][y] = c;
		}
	}
}

Automate NewAutomate (int n, int na)
{
	Automate a;
	a.n = n;
	a.na = na;
	a.e = (Etat *)malloc(sizeof(Etat)*n);
	int i, j;
	for (i=0;i<n;i++)
	{
		a.e[i].f = (int *)malloc(sizeof(int)*na);
		for (j=0;j<na;j++)
		{
			a.e[i].f[j] = -1;
			a.e[i].final = 0;
		}
	}
	return a;
}

void FreeAutomate (Automate a)
{
	int i;
	for (i=0;i<a.n;i++)
	{
		free(a.e[i].f);
	}
	free(a.e);
}

BetaAdic NewBetaAdic (int n)
{
	BetaAdic b;
	b.n = n;
	b.t = (Complexe *)malloc(sizeof(Complexe)*n);
	return b;
}

void FreeBetaAdic (BetaAdic b)
{
	free(b.t);
}

inline Color moy (Color a, Color b, double ratio)
{
	Color r;
	//ratio = 1.-ratio;
	r.r = a.r*(1.-ratio) + b.r*ratio;
	r.g = a.g*(1.-ratio) + b.g*ratio;
	r.b = a.b*(1.-ratio) + b.b*ratio;
	r.a = a.a*(1.-ratio) + b.a*ratio;
	return r;
}

void set_pix (Surface s, Complexe p)
{
	if (p.x < mx || p.x >= Mx || p.y < my || p.y >= My)
		return;
	double fx = (p.x - mx)*s.sx/(Mx-mx);
	double fy = (p.y - my)*s.sy/(My-my);
	unsigned int x = fx;
	unsigned int y = fy;
	if (x < s.sx && y < s.sy)
	{
		if (x+1 < s.sx && y+1 < s.sy)
		{
			s.pix[x][y] = moy(s.pix[x][y], color, (1.-fx+x)*(1.-fy+y));
			s.pix[x+1][y] = moy(s.pix[x+1][y], color, (fx-x)*(1.-fy+y));
			s.pix[x][y+1] = moy(s.pix[x][y+1], color, (1.-fx+x)*(fy-y));
			s.pix[x+1][y+1] = moy(s.pix[x+1][y+1], color, (fx-x)*(fy-y));
		}else
			s.pix[x][y] = color;
	}
}

void print_word (BetaAdic b, int n, int etat)
{
	if (etat < 0)
		return;
	if (n == 0)
		printf("(%d)", etat);
	else
	{
		int i;
		for (i=0;i<b.a.na;i++)
		{
			print_word(b, n-1, b.a.e[etat].f[i]);
		}
	}
}

//int niter;

Complexe pos;

void Draw_rec2 (BetaAdic b, Surface s, int n, int etat)
{
	if (etat < 0)
		return;
	if (n == 0)
	{
		if (etat >= 0 && etat < b.a.n)
		{
			if (!b.a.e[etat].final)
			{
				//printf("%d pas final !", etat);
				return;
			}
		}else
		{
			printf("état %d !\n", etat);
			return;
		}
		if (pos.x < mx2)
			mx2 = pos.x;
		if (pos.x > Mx2)
			Mx2 = pos.x;
		if (pos.y < my2)
			my2 = pos.y;
		if (pos.y > My2)
			My2 = pos.y;
		pos = add(pos, b.t[etat]);
		set_pix (s, pos);
	}else
	{
		int i;
		for (i=0;i<b.a.na;i++)
		{
			Draw_rec2 (b, s, n-1, b.a.e[etat].f[i]);
		}
	}
}

void Draw_rec (BetaAdic b, Surface s, int n, Complexe p, Complexe bn, int etat)
{
	if (n == 0)
	{
		if (etat >= 0 && etat < b.a.n)
		{
			if (!b.a.e[etat].final)
			{
				//printf("%d pas final !", etat);
				return;
			}
		}else
		{
			printf("état %d !\n", etat);
			return;
		}
		if (p.x < mx2)
			mx2 = p.x;
		if (p.x > Mx2)
			Mx2 = p.x;
		if (p.y < my2)
			my2 = p.y;
		if (p.y > My2)
			My2 = p.y;
		set_pix (s, p);
	}else
	{
		int i;
		for (i=0;i<b.a.na;i++)
		{
			//if (n == niter-1)
			//{
		    //	color = colors[i];
	        //}
			if (b.a.e[etat].f[i] != -1)
				Draw_rec (b, s, n-1, add(p, prod(bn, b.t[i])), prod(bn, b.b), b.a.e[etat].f[i]);
		}
	}
}

int max (int a, int b)
{
	if (a < b)
		return b;
	return a;
}

Color randCol (int a)
{
	Color c;
	c.r = rand();
	c.g = rand();
	c.b = rand();
	c.a = a;
	return c;
}

void Draw (BetaAdic b, Surface s, int n, int ajust, Color col, int verb)
{
	int auto_n = (n < 0);
	//set global variables
	mx2 = 1000000, my2 = 1000000, Mx2 = -1000000, My2 = -1000000; //extremum observés
	color0.r = color0.g = color0.b = 255;
	color0.a = 0;
	color = col;
	colors = (Color *)malloc(sizeof(Color)*b.a.n);
	int i;
	for (i=0;i<b.a.n;i++)
	{
		colors[i] = randCol(255);
	}
	if (verb)
	{
		printf("%d translations, %d lettres, %d états.\n", b.n, b.a.na, b.a.n);
		printf("couleur de fond : %d %d %d %d\n", color0.r, color0.g, color0.b, color0.a);
		printf("couleur de dessin : %d %d %d %d\n", color.r, color.g, color.b, color.a);
		printf("état initial : %d\n", b.a.i);
		printf("états finaux : ");
		for (i=0;i<b.a.n;i++)
		{
			if (b.a.e[i].final)
				printf("%d ", i);
		}
		printf("\n");
	}
	//ajust the window of the drawing
	if (ajust)
	{
		//ajuste le cadre
		if (auto_n)
		{
			//printf("max = %d\n", max(s.sx, s.sy));
			//printf("maj = %lf\n", (1.-sqrt(cnorm(b.b))));
			n = -2.*log(max(s.sx, s.sy)*(1.-sqrt(cnorm(b.b))))/log(cnorm(b.b));
			if (verb)
				printf("n = %d\n", n);
		}
		//niter = n;
		Draw_rec (b, s, n, zero(), un(), b.a.i);
		mx = mx2 - (Mx2 - mx2)/100;
		my = my2 - (My2 - my2)/100;
		Mx = Mx2 + (Mx2-mx2)/s.sx + (Mx2 - mx2)/100;;
		My = My2 + (My2-my2)/s.sy + (My2 - my2)/100;
		//preserve le ratio
		double delta = (Mx - mx)*s.sy - (My - my)*s.sx;
		if (delta > 0)
		{
			My = My + delta/(2*s.sx);
			my = my - delta/(2*s.sx);
		}else
		{
			Mx = Mx - delta/(2*s.sy);
			mx = mx + delta/(2*s.sy);
		}
	}
	if (verb)
	{
		printf("Zone de dessin : (%lf, %lf) (%lf, %lf)\n", mx, my, Mx, My);
	}
	Fill(s, color0);
	if (auto_n)
	{
		n = .5 - 2.*log(max(3.*s.sx/(Mx-mx), 3.*s.sy/(My-my)))/log(cnorm(b.b));
		if (verb)
			printf("n = %d\n", n);
	}
	//niter = n;
	Draw_rec (b, s, n, zero(), un(), b.a.i);
}

void Draw2 (BetaAdic b, Surface s, int n, int ajust, Color col, int verb)
{
	int auto_n = (n < 0);
	//set global variables
	mx2 = 1000000, my2 = 1000000, Mx2 = -1000000, My2 = -1000000; //extremum observés
	color0.r = color0.g = color0.b = 255;
	color0.a = 0;
	color = col;
	colors = (Color *)malloc(sizeof(Color)*b.a.n);
	int i;
	for (i=0;i<b.a.n;i++)
	{
		colors[i] = randCol(255);
	}
	if (verb)
	{
		printf("%d translations, %d lettres, %d états.\n", b.n, b.a.na, b.a.n);
		printf("couleur de fond : %d %d %d %d\n", color0.r, color0.g, color0.b, color0.a);
		printf("couleur de dessin : %d %d %d %d\n", color.r, color.g, color.b, color.a);
		printf("état initial : %d\n", b.a.i);
		printf("états finaux : ");
		for (i=0;i<b.a.n;i++)
		{
			if (b.a.e[i].final)
				printf("%d ", i);
		}
		printf("\n");
		printf("translations : ");
		for (i=0;i<b.n;i++)
		{
			printf("(%lf, %lf) ", b.t[i].x, b.t[i].y);
		}
		printf("\n");
	}
	//ajust the window of the drawing
	if (ajust)
	{
		//ajuste le cadre
		if (auto_n)
		{
			//printf("max = %d\n", max(s.sx, s.sy));
			//printf("maj = %lf\n", (1.-sqrt(cnorm(b.b))));
			n = -2.*log(max(s.sx, s.sy)*(1.-sqrt(cnorm(b.b))))/log(cnorm(b.b));
			if (verb)
				printf("n = %d\n", n);
		}
		pos = zero();
		Draw_rec2 (b, s, n, b.a.i);
		mx = mx2 - (Mx2 - mx2)/100;
		my = my2 - (My2 - my2)/100;
		Mx = Mx2 + (Mx2-mx2)/s.sx + (Mx2 - mx2)/100;;
		My = My2 + (My2-my2)/s.sy + (My2 - my2)/100;
		//preserve le ratio
		double delta = (Mx - mx)*s.sy - (My - my)*s.sx;
		if (delta > 0)
		{
			My = My + delta/(2*s.sx);
			my = my - delta/(2*s.sx);
		}else
		{
			Mx = Mx - delta/(2*s.sy);
			mx = mx + delta/(2*s.sy);
		}
	}
	if (verb)
	{
		printf("Zone de dessin : (%lf, %lf) (%lf, %lf)\n", mx, my, Mx, My);
	}
	Fill(s, color0);
	if (auto_n)
	{
		n = .5 - 2.*log(max(3.*s.sx/(Mx-mx), 3.*s.sy/(My-my)))/log(cnorm(b.b));
		if (verb)
			printf("n = %d\n", n);
	}
	pos = zero();
	Draw_rec2 (b, s, n, b.a.i);
}
