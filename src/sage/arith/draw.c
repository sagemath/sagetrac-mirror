#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include "complex.h"
#include "Automaton.h"
#include "automataC.h"
#include "numpy/ndarraytypes.h"
#include "draw.h"

//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

double mx = -2, my = -2, Mx = 2, My = 2; //zone de dessin

double mx2 = 1000000, my2 = 1000000, Mx2 = -1000000, My2 = -1000000; //extremum observés

Color color0; //couleur du fond
Color color; //couleur de dessin
Color* colors; //liste de couleurs de dessin

//rend une SDL_Surface contenant l'image
void* OpenImage (const char *file_name)
{
	return IMG_Load(file_name);
}

bool InImage (void* img, int x, int y)
{
	SDL_Surface *s = (SDL_Surface *)img;
	if (x < 0 || y < 0 || x >= s->w || y >= s->h)
		return false;
	Uint8 r,g,b,a;
	SDL_GetRGBA(*((Uint32 *)s->pixels + x + (s->pitch/4)*y), s->format, &r, &g, &b, &a);
	return a >= 128;
}

int ImageWidth (void *img)
{
	SDL_Surface *s = (SDL_Surface *)img;
	return s->w;
}

int ImageHeight (void *img)
{
	SDL_Surface *s = (SDL_Surface *)img;
	return s->h;
}

void CloseImage (void* img)
{
	SDL_Surface *s = (SDL_Surface *)img;
	SDL_FreeSurface(s);
}

//////////////////////////////TEST(
#define WIDTH	800
#define HEIGHT	600

void DrawRond (int x, int y, SDL_Surface *s)
{
	Uint32 *pix = s->pixels;
	int i,j;
	int size = 4;
	for (i=x-size;i<=x+size;i++)
	{
		for (j=y-size;j<=y+size;j++)
		{
			if ((i-x)*(i-x)+(j-y)*(j-y) < size*size)
			{
				if (i >= 0 && j >= 0 && i < WIDTH && j < HEIGHT)
					*(pix+i+(s->pitch/4)*j) = SDL_MapRGB(s->format, 255, 0, 0);
			}
		}
	}
}

void TestSDL()
{
	SDL_Window* win;
	SDL_Surface * s;
    int i, j;
    if (SDL_Init(SDL_INIT_VIDEO) == -1)
    {
        printf("Error during usage of SDL: %s\n", SDL_GetError());
        return;
    }
    
    win = SDL_CreateWindow("Test SDL 2.0", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 800, 600, SDL_WINDOW_SHOWN);
    if( win == NULL )
	{
		printf( "Window could not be created! SDL_Error: %s\n", SDL_GetError() );
		exit(1);
	}else
	{
		//Get window surface
		s = SDL_GetWindowSurface( win );
	}

    printf("Video Mode: %dx%d %d bits/pixel\n", s->w, s->h, s->format->BitsPerPixel);
           
    SDL_FillRect(s, NULL, SDL_MapRGB(s->format, 0x00, 0xff, 0xff));
    SDL_UpdateWindowSurface(win);
    
    Uint32 rmask, gmask, bmask, amask;

    /* SDL interprets each pixel as a 32-bit number, so our masks must depend
       on the endianness (byte order) of the machine */
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
    rmask = 0xff000000;
    gmask = 0x00ff0000;
    bmask = 0x0000ff00;
    amask = 0x000000ff;
#else
    rmask = 0x000000ff;
    gmask = 0x0000ff00;
    bmask = 0x00ff0000;
    amask = 0xff000000;
#endif

    //SDL_Surface *s = SDL_GetVideoSurface();SDL_CreateRGBSurface(0, WIDTH, HEIGHT, 32, rmask, gmask, bmask, amask);
    /*
    if(s == NULL)
    {
        fprintf(stderr, "CreateRGBSurface failed: %s\n", SDL_GetError());
        exit(1);
    }*/
    
    Uint32 *pix = s->pixels;
    int x, y;
    for (y=0;y<HEIGHT;y++)
    {
    	for (x=0;x<WIDTH;x++)
    	{
    		*pix = SDL_MapRGB(s->format, x, y, x+y);
    		pix++;
    	}
    	pix += (s->pitch/4-WIDTH);
    }
    
    SDL_UpdateWindowSurface(win);
    
    int quit = 0;
    pix = s->pixels;
	SDL_Event event;
	for(;;)
	{
		SDL_WaitEvent(&event); // Récupération des actions de l'utilisateur
		switch(event.type)
		{
			case SDL_QUIT: // Clic sur la croix
				quit=1;
				break;
			case SDL_KEYUP: // Relâchement d'une touche
				if ( event.key.keysym.sym == SDLK_f ) // Touche f
				{
				    
				}
				break;
			case SDL_MOUSEMOTION:
				if (event.motion.state & SDL_BUTTON_LMASK)
				{
					x = event.motion.x;
					y = event.motion.y;
					DrawRond(x,y,s);
					SDL_UpdateWindowSurface(win);
				}
				break;
		}
		if (quit)
			break;
	}            
	
	//SDL_FreeSurface(s);
    SDL_Quit();
}
//////////////////////////////)TEST

void *GetSDL_SurfaceFromNumpy (PyArrayObject *o)
{
    //PyArrayObject *o = (PyArrayObject *)np;
    if (o->nd != 2)
    {
        printf("Error: numpy array must be two-dimensional (here %d-dimensional).", o->nd);
        return NULL;
    }
    if (o->strides[0] != 4)
    {
        printf("Error: pixels must be stored with 4 bytes (RGBA format). Here %ld bytes/pixel.", o->strides[0]);
    }
    
    Uint8 *data = (Uint8 *)o->data;
    
    Uint32 rmask, gmask, bmask, amask;
    /* SDL interprets each pixel as a 32-bit number, so our masks must depend
       on the endianness (byte order) of the machine */
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
    rmask = 0xff000000;
    gmask = 0x00ff0000;
    bmask = 0x0000ff00;
    amask = 0x000000ff;
#else
    rmask = 0x000000ff;
    gmask = 0x0000ff00;
    bmask = 0x00ff0000;
    amask = 0xff000000;
#endif
    int sx = o->dimensions[1];
    int sy = o->dimensions[0];
	SDL_Surface *r = SDL_CreateRGBSurface(0, sx, sy, 32, rmask, gmask, bmask, amask);
	int x,y;
	Uint32 *ptr = r->pixels;
	for (y=0;y<sy;y++)
	{
		for (x=0;x<sx;x++)
		{
			// *ptr = SDL_MapRGBA(r->format, x, y, x-y, 255);
			*ptr = SDL_MapRGBA(r->format, *data, *(data+1), *(data+2), *(data+3));
			data+=4;
			ptr++;
		}
		ptr += (r->pitch/4) - sx;
	}
	return (void *)r;
}

//dessine la surface dans la SDL_Surface
SDL_Surface *GetSurface (Surface s)
{
	Uint32 rmask, gmask, bmask, amask;

    /* SDL interprets each pixel as a 32-bit number, so our masks must depend
       on the endianness (byte order) of the machine */
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
    rmask = 0xff000000;
    gmask = 0x00ff0000;
    bmask = 0x0000ff00;
    amask = 0x000000ff;
#else
    rmask = 0x000000ff;
    gmask = 0x0000ff00;
    bmask = 0x00ff0000;
    amask = 0xff000000;
#endif
	SDL_Surface *r = SDL_CreateRGBSurface(0, s.sx, s.sy, 32, rmask, gmask, bmask, amask);
	int x,y;
	Uint32 *ptr = r->pixels;
	for (y=0;y<s.sy;y++)
	{
		for (x=0;x<s.sx;x++)
		{
			// *ptr = SDL_MapRGBA(r->format, x, y, x-y, 255);
			*ptr = SDL_MapRGBA(r->format, s.pix[x][y].r, s.pix[x][y].g, s.pix[x][y].b, s.pix[x][y].a);
			ptr++;
		}
		ptr += (r->pitch/4) - s.sx;
	}
	return r;
}

Complexe getComplexe (int x, int y, int sx, int sy)
{
	Complexe r;
	r.x = x*(Mx-mx)/sx + mx;
	r.y = y*(My-my)/sy + my;
	return r;
}

void ComplexeToPoint (Complexe c, int *x, int *y, int sx, int sy)
{
	*x = (c.x - mx)*sx/(Mx - mx);
	*y = (c.y - my)*sy/(My - my);
}

//dessine la transfomée inverse de l'image s
void drawTransf (SDL_Surface *s, SDL_Surface *screen, Complexe m, Complexe t, Color col)
{
	int x, y, i, j, k;
	Uint32 *f = s->pixels;
	Uint32 *d = screen->pixels;
	Complexe c;
	Uint8 r, g, b, a;
	Uint8 r2, g2, b2, a2;
	//calcule un encadrement
	int xmin = screen->w-1, ymin = screen->h-1, xmax = 0, ymax = 0;
	Complexe im = inv(m);
	int cx[4];
	int cy[4];
	cx[0] = cy[0] = 0;
	cx[1] = 0;
	cy[1] = s->h-1;
	cx[2] = s->w-1;
	cy[2] = 0;
	cx[3] = s->w-1;
	cy[3] = s->h-1;
	for (k=0;k<4;k++)
	{
		ComplexeToPoint(sub(prod(getComplexe(cx[k], cy[k], s->w, s->h), m), t), &i, &j, screen->w, screen->h);
		if (i < xmin)
			xmin = i;
		if (i > xmax)
			xmax = i;
		if (j < ymin)
			ymin = j;
		if (j > ymax)
			ymax = j;
	}
	if (xmin < 0)
		xmin = 0;
	if (xmax >= screen->w)
		xmax = screen->w-1;
	if (ymin < 0)
		ymin = 0;
	if (ymax >= screen->h)
		ymax = screen->h-1;
	//
	for (y=ymin;y<=ymax;y++)
	{
		for (x=xmin;x<=xmax;x++)
		{
			c = getComplexe(x, y, screen->w, screen->h);
			ComplexeToPoint(prod(m, add(c, t)), &i, &j, s->w, s->h);
			if (i >= 0 && j >= 0 && i < s->w && j < s->h)
			{
				d = ((Uint32 *)screen->pixels) + x + (screen->pitch/4)*y;
				SDL_GetRGBA(*(f+i+(s->pitch/4)*j), s->format, &r, &g, &b, &a);
				SDL_GetRGBA(*d, screen->format, &r2, &g2, &b2, &a2);
				r = ((Uint32)col.r*a+(Uint32)r2*(255-a))/256;
				g = ((Uint32)col.g*a+(Uint32)g2*(255-a))/256;
				b = ((Uint32)col.b*a+(Uint32)b2*(255-a))/256;
				*d = SDL_MapRGBA(screen->format, r, g, b, a2);
			}
			//d++;
		}
		//d += (screen->pitch/4 - screen->w);
	}
}

int lt[256]; //liste des indices des translations du morceau courant

Complexe barycentre;

//cherche la translation donnant le morceau le plus proche du point
//morceau de taille b^n
Complexe FindTr (int n, Complexe c, BetaAdic b, SDL_Surface *s, bool *ok, bool verb)
{
	Complexe r = zero();
	int i, j;
	double nn, nmax;
	int imax;
	Complexe ib = inv(b.b);
	Complexe m = un();
	Complexe bb = prod(barycentre, b.b);
	int x,y;
	Uint8 r0,g0,b0,a;
	if (ok)
		*ok = false;
	for (j=0;j<n;j++)
	{
		nmax = -1;
		for (i=0;i<b.n;i++)
		{
			//calcule la distance entre c et b.t[i]
			nn = cnorm(sub(c, add(bb, b.t[i])));
			if (nmax == -1 || nn < nmax)
			{
				//teste si l'on est bien inclus dans le morceau
				ComplexeToPoint(prod(sub(c, b.t[i]), ib), &x, &y, s->w, s->h);
				if (x < 0 || y < 0 || x >=s->w || y >= s->h)
					continue; //le point n'est pas dedans
				SDL_GetRGBA(*((Uint32 *)s->pixels+x+(s->pitch/4)*y), s->format, &r0, &g0, &b0, &a);
				if (a >= 50)
				{
					if (ok)
						*ok = true;
					imax = i;
					nmax = nn;
				}
			}
		}
		if (nmax == -1)
		{
			if (ok)
				*ok = false;
			return r;
		}
		lt[j] = imax;
		r = add(r, prod(b.t[imax], inv(m)));
		c = prod(sub(c, b.t[imax]), ib);
		//if (verb)
		//	printf("%d ", imax);
		m = prod(m, ib);
	}
	//if (verb)
	//	printf("\n");
	r.x = -r.x;
	r.y = -r.y;
	return r;
}

bool addA(Automaton *a, int n)
{
	int e = a->i;
	int i;
	bool r = false;
	for (i=0;i<n-1;i++)
	{
		if (a->e[e].f[lt[i]] == -1)
		{
			r = true;
			AddState(a, false);
			a->e[e].f[lt[i]] = a->n-1;
		}
		e = a->e[e].f[lt[i]];
	}
	if (a->e[e].f[lt[i]] != 1)
	{
		r = true;
		a->e[e].f[lt[i]] = 1;
	}
	return r;
}

Automaton UserDraw (BetaAdic b, int sx, int sy, int n, int ajust, Color col, int only_pos, int verb)
{
	Automaton r = NewAutomaton(2, b.n);
	r.i = 0;
	int i;
	if (!only_pos)
	{
        for (i=0;i<b.n;i++)
        {
            r.e[1].f[i] = 1; //état reconnaissant tout
        }
	}
	r.e[0].final = false;
	r.e[1].final = true;

	if (SDL_Init(SDL_INIT_VIDEO) == -1)
    {
        printf("Erreur lors de l'initialisation de SDL: %s\n", SDL_GetError());
        return r;
    }

	Surface s0 = NewSurface(sx, sy);
	if (verb)
	{
		printf("Dessin de la fractale...\n");
		printf("n=%d, ajust=%d, verb=%d\n", n, ajust, verb);
	}
	Draw(b, s0, n, ajust, col, 3., verb);
	if (verb)
		printf("Conversion en SDL...\n");
	SDL_Surface *s = GetSurface(s0);	//utilisé pour dessiner les transformées
	
	    Uint32 rmask, gmask, bmask, amask;

    // SDL interprets each pixel as a 32-bit number, so our masks must depend
    //  on the endianness (byte order) of the machine
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
    rmask = 0xff000000;
    gmask = 0x00ff0000;
    bmask = 0x0000ff00;
    amask = 0x000000ff;
#else
    rmask = 0x000000ff;
    gmask = 0x0000ff00;
    bmask = 0x00ff0000;
    amask = 0xff000000;
#endif
    
	if (verb)
	{
		printf("s de taille %dx%d\n", s->w, s->h);
		printf("Ouverture de la fenêtre...\n");
    }
    
    SDL_Window* win;
	SDL_Surface *screen;

    // Création de la fenêtre
    win = SDL_CreateWindow("Fenetre user_draw", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, sx, sy, SDL_WINDOW_SHOWN);
    if( win == NULL )
	{
		printf( "Window could not be created! SDL_Error: %s\n", SDL_GetError() );
		exit(1);
	}else
	{
		//Get window surface
		screen = SDL_GetWindowSurface( win );
	}
    if (verb)
    	printf("Mode vidéo: %dx%d %d bits/pixel\n", screen->w, screen->h, screen->format->BitsPerPixel);
    
    //surface dans laquelle on dessine le résultat
    SDL_Surface *sf = SDL_CreateRGBSurface(0, screen->w, screen->h, 32, rmask, gmask, bmask, amask);
    if(sf == NULL)
    {
        fprintf(stderr, "CreateRGBSurface failed: %s\n", SDL_GetError());
        exit(1);
    }
    //SDL_SetAlpha(sf, 0, 255);
         
    //SDL_FillRect(screen, NULL, SDL_MapRGB(screen->format, 0x00, 0x00, 0x00));
    SDL_FillRect(sf, NULL, SDL_MapRGB(sf->format, 0x00, 0x00, 0x00));
    col.r = 80;
    col.g = 80;
    col.b = 80;
    drawTransf(s, sf, un(), zero(), col);
    if (SDL_BlitSurface(sf, NULL, screen, NULL) < 0)
    {
        printf("Error : %s\n", SDL_GetError());
    }
    int np = 4;
    Complexe f = powC(b.b, -np);
    Complexe t;
    Complexe rt;
    t.x = -b.t[0].x;
    t.y = -b.t[0].y;
    Color colf;
    col.r = 100;
    col.g = 200;
    col.b = 50;
    colf.r = 150;
    colf.g = 100;
    colf.b = 200;
    drawTransf(s, screen, f, t, col);
    SDL_UpdateWindowSurface(win);
    
    int quit = 0;
    int x, y;
	SDL_Event event;
	bool ok;
	bool clic = false;
	for(;;)
	{
		SDL_WaitEvent(&event); // Récupération des actions de l'utilisateur
		switch(event.type)
		{
			case SDL_QUIT: // Clic sur la croix
				quit=1;
				break;
			case SDL_KEYUP: // Relâchement d'une touche
				if ( event.key.keysym.sym == SDLK_p ) // Touche p
				{
					if (np < 255)
					{
				    	np++;
				    	f = powC(b.b, -np);
				    	if (verb)
					    	printf("np = %d\n", np);
				    }
				}
				if ( event.key.keysym.sym == SDLK_m ) // Touche m
				{
				    np--;
				    f = powC(b.b, -np);
				    if (verb)
					    printf("np = %d\n", np);
				}
				if ( event.key.keysym.sym == SDLK_ESCAPE )
				{
					quit = 1;
				}
				break;
			case SDL_MOUSEBUTTONUP:
				clic = false;
				break;
			case SDL_MOUSEBUTTONDOWN:
			case SDL_MOUSEMOTION:
				x = event.motion.x;
				y = event.motion.y;
				Complexe c = getComplexe(x, y, screen->w, screen->h);
				//cherche le morceau le plus proche du point c
				rt = t;
				t = FindTr(np, c, b, s, &ok, verb);
				if (!ok)
					break;
				//
				if (event.type == SDL_MOUSEBUTTONDOWN || (clic && event.motion.state & SDL_BUTTON_LMASK))
				{ //clic
					clic = true;
					if (addA(&r, np)) //ajoute le morceau à l'automate
					{ //si morceau ajouté
						drawTransf(s, sf, f, t, colf);
						if (SDL_BlitSurface(sf, NULL, screen, NULL) < 0)
                        {
                            printf("Error : %s\n", SDL_GetError());
                        }
						ComplexeToPoint(zero(), &x, &y, screen->w, screen->h);
						DrawRond(x, y, screen);
						SDL_UpdateWindowSurface(win);
					}
				}else
				{
					if (rt.x != t.x || rt.y != t.y)
					{
						if (SDL_BlitSurface(sf, NULL, screen, NULL) < 0)
                        {
                            printf("Error : %s\n", SDL_GetError());
                        }
						drawTransf(s, screen, f, t, col);
						ComplexeToPoint(zero(), &x, &y, screen->w, screen->h);
						DrawRond(x, y, screen);
						SDL_UpdateWindowSurface(win);
					}
				}
				break;
		}
		if (quit)
			break;
	}            
	
	SDL_FreeSurface(sf);
	SDL_FreeSurface(s);
    SDL_Quit();
    return r;
}

void invert(const SDL_PixelFormat* format, Uint32 *p)
{
	Uint8 r,g,b;
	SDL_GetRGB(*p, format, &r, &g, &b);
	*p = SDL_MapRGB(format, 255-r, 255-g, 255-b);
}

int mini (int a, int b)
{
	if (a < b)
		return a;
	return b;
}

int max (int a, int b)
{
	if (a > b)
		return a;
	return b;
}

void invertRect(SDL_Surface *s, int x1, int y1, int x2, int y2)
{
	if (x1 < 0 || x2 < 0 || y1 < 0 || y2 < 0)
		return;
	if (x1 >= s->w || x2 >= s->w || y1 >= s->h || y2 >= s->h)
		return;
	int x,y;
	Uint32 *ptr = s->pixels;
	for (x=mini(x1,x2);x<max(x1,x2);x++)
	{
		invert(s->format, ptr+x+y1*s->pitch/4);
		invert(s->format, ptr+x+y2*s->pitch/4);
	}
	for (y=mini(y1,y2);y<max(y1,y2);y++)
	{
		invert(s->format, ptr+x1+y*s->pitch/4);
		invert(s->format, ptr+x2+y*s->pitch/4);
	}
}

int sign(int a)
{
	if (a > 0)
		return 1;
	if (a == 0)
		return 0;
	return -1;
}

double mousex = 0, mousey = 0;

double *Rmaj = NULL; //liste de majorants mesurés pour chaque état

//ouvre une fenêtre avec la fractale sur laquelle on peut zoomer
//coeff = coefficient de calcul du nombre d'itérations
void DrawZoom (BetaAdic b, int sx, int sy, int n, int ajust, Color col, double coeff, int verb)
{
	if (SDL_Init(SDL_INIT_VIDEO) == -1)
    {
        printf("Erreur lors de l'initialisation de SDL: %s\n", SDL_GetError());
        return;
    }
    
    Surface s0 = NewSurface(sx, sy);
    
    Rmaj = NULL;
    if (ajust)
    {
    	//réinitialise les extremum observés
    	mx2 = 1000000; my2 = 1000000; Mx2 = -1000000; My2 = -1000000;
    	mx = -1000000; my = -1000000; Mx = 1000000; My = 1000000;
    	/*
    	if (b.a.n < 10)
    	{
    		Rmaj = (double *)malloc(sizeof(double)*b.a.n);
    		if (verb)
    			printf("Calcule des bornes pour chaque état...\n");
    		//calcule les majorants
    		int i;
    		for (i=0;i<b.a.n;i++)
    		{
    			printf("Calul du majorant de %d...\n", i);
    			mx2 = 1000000; my2 = 1000000; Mx2 = -1000000; My2 = -1000000; //réinitialise les extremum observés
    			Draw(b, s0, n, ajust, col, verb);
    			Rmaj[i] = sqrt((Mx2-mx2)*(Mx2-mx2) + (My2-my2)*(My2-my2))/2;
    			printf("	-> %lf\n", Rmaj[i]);
    		}
    	}
    	*/
    }
    
	if (verb)
	{
		printf("Dessin de la fractale...\n");
		printf("n=%d, ajust=%d, verb=%d\n", n, ajust, verb);
	}
	Draw(b, s0, n, ajust, col, coeff, verb);
	if (verb)
		printf("Conversion en SDL...\n");
	SDL_Surface *s = GetSurface(s0);	//utilisé pour dessiner les transformées
    
	if (verb)
	{
		printf("s de taille %dx%d\n", s->w, s->h);
		printf("Ouverture de la fenêtre...\n");
    }
    
    SDL_Window* win;
	SDL_Surface *screen;

    // Création de la fenêtre
    win = SDL_CreateWindow("Fenetre Draw zoom", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, sx, sy, SDL_WINDOW_SHOWN);
    if( win == NULL )
	{
		printf( "Window could not be created! SDL_Error: %s\n", SDL_GetError() );
		exit(1);
	}else
	{
		//Get window surface
		screen = SDL_GetWindowSurface( win );
	}
    if (verb)
    	printf("Mode vidéo: %dx%d %d bits/pixel\n", screen->w, screen->h, screen->format->BitsPerPixel);
         
    //SDL_FillRect(screen, NULL, SDL_MapRGB(screen->format, 0x00, 0x00, 0x00));
    SDL_FillRect(screen, NULL, SDL_MapRGB(screen->format, 0xFF, 0xFF, 0xFF));
    SDL_BlitSurface(s, NULL, screen, NULL);
    SDL_UpdateWindowSurface(win);
       
    int quit = 0;
    int x, y, rx, ry, x1, y1, x2, y2;
	SDL_Event event;
	bool clic = false;
	bool redraw = false;
	for(;;)
	{
		SDL_WaitEvent(&event); // Récupération des actions de l'utilisateur
		switch(event.type)
		{
			case SDL_QUIT: // Clic sur la croix
				quit=1;
				break;
			case SDL_KEYUP: // Relâchement d'une touche
				if (event.key.keysym.sym == SDLK_ESCAPE )
				{
					quit = 1;
				}
				if (event.key.keysym.sym == SDLK_RIGHT)
				{
					double px = (Mx - mx)/10;
					Mx += px;
					mx += px;
					redraw = true;
				}
				if (event.key.keysym.sym == SDLK_LEFT)
				{
					double px = (Mx - mx)/10;
					Mx -= px;
					mx -= px;
					redraw = true;
				}
				if (event.key.keysym.sym == SDLK_UP)
				{
					double py = (My - my)/10;
					My -= py;
					my -= py;
					redraw = true;
				}
				if (event.key.keysym.sym == SDLK_DOWN)
				{
					double py = (My - my)/10;
					My += py;
					my += py;
					redraw = true;
				}
				if (event.key.keysym.sym == SDLK_m)
				{
					double m = 1.4;
					double m2;
					m2 = (Mx+mx - (Mx-mx)*m)/2;
					Mx = (Mx+mx + (Mx-mx)*m)/2;
					mx = m2;
					m2 = (My+my - (My-my)*m)/2;
					My = (My+my + (My-my)*m)/2;
					my = m2;
					redraw = true;
				}
				if (event.key.keysym.sym == SDLK_p)
				{
					double m = 1/1.4;
					double m2;
					m2 = (Mx+mx - (Mx-mx)*m)/2;
					Mx = (Mx+mx + (Mx-mx)*m)/2;
					mx = m2;
					m2 = (My+my - (My-my)*m)/2;
					My = (My+my + (My-my)*m)/2;
					my = m2;
					redraw = true;
				}
				if (event.key.keysym.sym == SDLK_r)
				{
					//réinitialise les extremum observés
			    	mx2 = 1000000; my2 = 1000000; Mx2 = -1000000; My2 = -1000000;
			    	mx = -1000000; my = -1000000; Mx = 1000000; My = 1000000;
    				Draw(b, s0, n, true, col, coeff, verb);
					redraw = true;
				}
				break;
			case SDL_MOUSEBUTTONUP:
			case SDL_MOUSEBUTTONDOWN:
			case SDL_MOUSEMOTION:
				x = event.motion.x;
				y = event.motion.y;
				//Complexe c = getComplexe(x, y, screen->w, screen->h);
				if (clic)
				{
					if (abs(x-x1)*sy > abs(y-y1)*sx)
					{
						x = x1 + sign(x-x1)*abs(y-y1)*sx/sy;
					}else
					{
						y = y1 + sign(y-y1)*abs(x-x1)*sy/sx;
					}
				}else if (event.type == SDL_MOUSEBUTTONDOWN)
				{
					clic = true;
					x1 = x;
					y1 = y;
					rx = -1;
				}
				if (event.type == SDL_MOUSEBUTTONUP)
				{
					x2 = rx;
					y2 = ry;
					// change la fenetre de dessin et redessine
					Mx2 = mx + (double)max(x1,x2)*(Mx-mx)/sx;
					My2 = my + (double)max(y1,y2)*(My-my)/sy;
					mx = mx + (double)mini(x1,x2)*(Mx-mx)/sx;
					my = my + (double)mini(y1,y2)*(My-my)/sy;
					Mx = Mx2;
					My = My2;
					//invertRect(screen, x1, y1, rx, ry); //éfface le précédent dessin
					redraw = true;
					clic = false;
				}
				if (clic && (x != rx || y != ry))
				{
					if (rx != -1)
						invertRect(screen, x1, y1, rx, ry); //éfface le précédent rectangle
					invertRect(screen, x1, y1, x, y); //dessine le nouveau rectangle
					SDL_UpdateWindowSurface(win);
				}
				rx = x;
				ry = y;
				break;
		}
		if (redraw)
		{
			redraw = false;
			if (verb)
			{
				printf("Dessin de la fractale...\n");
				printf("zone (%lf ... %lf)x(%lf ... %lf)\n", mx, Mx, my, My);
				printf("n=%d, ajust=%d, verb=%d\n", n, ajust, verb);
			}
			Draw(b, s0, n, false, col, coeff, verb);
			s = GetSurface(s0);	//utilisé pour dessiner les transformées
		    SDL_FillRect(screen, NULL, SDL_MapRGB(screen->format, 0xFF, 0xFF, 0xFF));
		    SDL_BlitSurface(s, NULL, screen, NULL);
		    SDL_UpdateWindowSurface(win);
		}
		if (quit)
			break;
	}            
	
	if (Rmaj)
		free(Rmaj);
	SDL_FreeSurface(s);
    SDL_Quit();
}

ColorList NewColorList (int n)
{
	ColorList l = malloc(sizeof(Color)*n);
	return l;
}

void FreeColorList (ColorList l)
{
	free(l);
}

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

/*
Automate NewAutomate (int n, int na)
{
	Automate a;
	a.n = n;
	a.na = na;
	if (n == 0)
	{
		a.e = NULL;
		return a;
	}
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
	if (a.n)
		free(a.e);
}

void FreeAutomates (Automate* a, int n)
{
	int j;
	for (j=0;j<n;j++)
	{
		FreeAutomaton(&a[j]);
	}
}
*/

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

BetaAdic2 NewBetaAdic2 (int n, int na)
{
	BetaAdic2 b;
	b.n = n;
	b.na = na;
	b.t = (Complexe *)malloc(sizeof(Complexe)*n);
	b.a = (Automaton *)malloc(sizeof(Automaton)*na);
	return b;
}

void FreeBetaAdic2 (BetaAdic2 b)
{
	free(b.t);
	free(b.a);
}

Color moy (Color a, Color b, double ratio)
{
	Color r;
	//ratio = 1.-ratio;
	r.r = a.r*(1.-ratio) + b.r*ratio;
	r.g = a.g*(1.-ratio) + b.g*ratio;
	r.b = a.b*(1.-ratio) + b.b*ratio;
	r.a = a.a*(1.-ratio) + b.a*ratio;
	return r;
}

int ns;
Complexe cs;

bool set_pix (Surface s, Complexe p)
{
	if (p.x < mx || p.x >= Mx || p.y < my || p.y >= My)
		return false;
	cs.x += p.x;
	cs.y += p.y;
	ns++;
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
		return true;
	}
	return false;
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
double Maj = 1000; //majorant
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

void DrawList_rec (BetaAdic2 b, Surface s, int n, Complexe p, Complexe bn, int *etat)
{
	if (etat[0] == -1)
		return;
	int i;
	if (n == 0)
	{
		if (!b.a[0].e[etat[0]].final)
		{
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
		color = colors[0];
		for (i=b.na-1;i>0;i--)
		{
			if (etat[i] != -1)
			{
				if (b.a[i].e[etat[i]].final)
					color = colors[i];
			}
		}
		set_pix (s, p);
	}else
	{
	    if (n > 5)
		{
			//teste si l'on sort de la zone de dessin
			//sous-arbre inclus dans le disque de centre p et de rayon abs(bn)*M
			double Mn = Maj*sqrt(cnorm(bn));
			if (p.x + Mn > mx && p.x - Mn < Mx && p.y + Mn > my && p.y - Mn < My)
			{
				/*
				if (Rmaj != NULL)
				{
					///////////TODO !!!
				}
				*/
			}else
				return; //intersection des rectangles vide
		}
		int *etat2 = (int *)malloc(sizeof(int)*b.na);
		int j;
		for (i=0;i<b.n;i++)
		{
			for (j=0;j<b.na;j++)
			{
				if (etat[j] != -1)
					etat2[j] = b.a[j].e[etat[j]].f[i];
				else
					etat2[j] = -1;
			}
			DrawList_rec (b, s, n-1, add(p, prod(bn, b.t[i])), prod(bn, b.b), etat2);
		}
		free(etat2);
	}
}

//compute the (almost) longest word w such that w(w^{-1}L) and L give the same drawing in the zone
//where L is the language given by b
//NON UTILISE : NE FONCTIONNE PAS BIEN
void WordZone (BetaAdic b, int *word, int nmax)
{
	Complexe p = zero();
	Complexe bn = un();
	int e = b.a.i;
	int i, n, f;
	if (cnorm(b.b) < 1)
	{
		//calcul du majorant
		Maj = 0;
		double m;
		for (i=0;i<b.n;i++)
		{
			m = cnorm(b.t[i]);
			if (Maj < m)
				Maj = m;
		}
		Maj = sqrt(Maj)/(1. - sqrt(cnorm(b.b)));
	}
	//parcours
	bool ok;
	for (n=0;n<nmax;n++)
	{
		ok = false;
		//parcours les fils
		for (i=0;i<b.a.na;i++)
		{
			f = b.a.e[e].f[i];
			if (f != -1)
			{
				//teste si l'on est dans la zone de dessin ou pas
				double Mn = Maj*sqrt(cnorm(bn));
				if (p.x + Mn > mx && p.x - Mn < Mx && p.y + Mn > my && p.y - Mn < My)
				{
					if (ok)
					{ //il y a déjà un autre fils dans la zone de dessin
						word[n] = -1;
						return;
					}
					ok = true;
					word[n] = i;
				}
			}
		}
		if (!ok)
		{ //il n'y a rien dans la zone de dessin
			word[n] = -1;
			return;
		}
	}
}

int word[1024];

int *WordDrawn ()
{
	word[1023] = -1;
	return word;
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
		if (set_pix(s, p))
			word[0] = -2;
	}else
	{
		if (n > 5)
		{
			//teste si l'on sort de la zone de dessin
			//sous-arbre inclus dans le disque de centre p et de rayon abs(bn)*M
			double Mn = Maj*sqrt(cnorm(bn));
			if (p.x + Mn > mx && p.x - Mn < Mx && p.y + Mn > my && p.y - Mn < My)
			{
				/*
				if (Rmaj != NULL)
				{
					///////////TODO !!!
				}
				*/
			}else
				return; //intersection des rectangles vide
		}
		int i;
		for (i=0;i<b.a.na;i++)
		{
			//if (n == niter-1)
			//{
			//	color = colors[i];
			//}
			if (b.a.e[etat].f[i] != -1)
			{
				Draw_rec (b, s, n-1, add(p, prod(bn, b.t[i])), prod(bn, b.b), b.a.e[etat].f[i]);
				if (n < 1024)
				{
					if (word[n-1] == -2)
					{
						word[n-1] = i;
						word[n] = -2;
					}
				}
			}
		}
	}
}

Color randColor (int a)
{
	Color c;
	c.r = rand();
	c.g = rand();
	c.b = rand();
	c.a = a;
	return c;
}

double absd (double f)
{
	if (f < 0)
		return -f;
	return f;
}

void Draw (BetaAdic b, Surface s, int n, int ajust, Color col, double coeff, int verb)
{
	int auto_n = (n < 0);
	//set global variables
	mx2 = 1000000, my2 = 1000000, Mx2 = -1000000, My2 = -1000000; //extremum observés
	if (ajust)
	{
		mx = -1000000; my = -1000000; Mx = 1000000; My = 1000000;
	}
	color0.r = color0.g = color0.b = 255;
	color0.a = 0;
	color = col;
	int i, j;
	/*
	colors = (Color *)malloc(sizeof(Color)*b.a.n);
	for (i=0;i<b.a.n;i++)
	{
		colors[i] = randCol(255);
	}
	*/
	
	word[0] = -1; //initialise le mot
	
	if (cnorm(b.b) < 1)
	{
		//calcul du majorant
		Maj = 0;
		double m;
		for (i=0;i<b.n;i++)
		{
			m = cnorm(b.t[i]);
			if (Maj < m)
				Maj = m;
		}
		Maj = sqrt(Maj)/(1. - sqrt(cnorm(b.b)));
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
		for (i=0;i<b.a.n;i++)
		{
			for (j=0;j<b.a.na;j++)
			{
				if (b.a.e[i].f[j] != -1)
					printf("%d --%d--> %d, ", i, j, b.a.e[i].f[j]);
			}
		}
		printf("\n");
	}
	//ajust the window of the drawing
	if (ajust)
	{
		//ajuste le cadre
		if (auto_n)
		{
			if (verb)
			{
				printf("max = %d\n", max(s.sx, s.sy));
				printf("maj = %lf\n", (1.-sqrt(cnorm(b.b))));
				printf("abs(b) = %lf\n", sqrt(cnorm(b.b)));
			}
			if (cnorm(b.b) < 1)
				n = .5 + -3.*log(max(s.sx, s.sy)*absd(1.-sqrt(cnorm(b.b))))/log(cnorm(b.b));
			else
				n = .5 + 3.*log(max(s.sx, s.sy)*absd(1.-1./sqrt(cnorm(b.b))))/log(cnorm(b.b));
			
			if (n < 0)
				n = 0;
			if (n > 1000000)
				n = 100;
				
			if (verb)
				printf("n = %d\n", n);
		}
		//niter = n;
		ns = 0;
		cs.x = 0;
		cs.y = 0;
		Draw_rec(b, s, n, zero(), un(), b.a.i);
		barycentre.x = cs.x/ns;
		barycentre.y = cs.y/ns;
		mx = mx2 - (Mx2 - mx2)/100;
		my = my2 - (My2 - my2)/100;
		Mx = Mx2 + (Mx2-mx2)/s.sx + (Mx2 - mx2)/20;
		My = My2 + (My2-my2)/s.sy + (My2 - my2)/20;
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
		if (cnorm(b.b) < 1)
			n = coeff + 2.*absd(log(max(s.sx/(Mx-mx), s.sy/(My-my)))/log(cnorm(b.b)));
		else
		{
			double ratio = pow(cnorm(b.b), n/2);
			n = .75 + coeff*absd(log(ratio*max(3.*s.sx/(Mx-mx), 3.*s.sy/(My-my)))/log(cnorm(b.b)));
			//change d'échelle pour dessiner la fractale à la bonne échelle
			ratio = pow(cnorm(b.b), n/2)/ratio;
			Mx = Mx*ratio;
			mx = mx*ratio;
			My = My*ratio;
			my = my*ratio;
			if (verb)
			{
				printf("Zone de dessin : (%lf, %lf) (%lf, %lf)\n", mx, my, Mx, My);
			}
		}
		
		if (n < 0)
			n = 0;
		if (n > 1000000)
			n = 100;
		
		if (verb)
			printf("n = %d\n", n);
	}
	//niter = n;
	Draw_rec(b, s, n, zero(), un(), b.a.i);
}

void Draw2 (BetaAdic b, Surface s, int n, int ajust, Color col, int verb)
{
	int auto_n = (n < 0);
	//set global variables
	mx2 = 1000000, my2 = 1000000, Mx2 = -1000000, My2 = -1000000; //extremum observés
	color0.r = color0.g = color0.b = 255;
	color0.a = 0;
	color = col;
	int i, j;
	/*
	colors = (Color *)malloc(sizeof(Color)*b.a.n);
	for (i=0;i<b.a.n;i++)
	{
		colors[i] = randCol(255);
	}
	*/
	
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
		for (i=0;i<b.a.n;i++)
		{
			for (j=0;j<b.a.na;j++)
			{
				if (b.a.e[i].f[j] != -1)
					printf("%d --%d--> %d, ", i, j, b.a.e[i].f[j]);
			}
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
			if (cnorm(b.b) < 1)
				n = -2.*log(max(s.sx, s.sy)*(1.-sqrt(cnorm(b.b))))/log(cnorm(b.b));
			else
				n = 2.*log(max(s.sx, s.sy)*(1.-1./sqrt(cnorm(b.b))))/log(cnorm(b.b));
			if (n < 0)
				n = 0;
			if (n > 1000000)
				n = 100;
			if (verb)
				printf("n = %d\n", n);
		}
		pos = zero();
		Draw_rec2 (b, s, n, b.a.i);
		mx = mx2 - (Mx2 - mx2)/100;
		my = my2 - (My2 - my2)/100;
		Mx = Mx2 + (Mx2-mx2)/s.sx + (Mx2 - mx2)/20;
		My = My2 + (My2-my2)/s.sy + (My2 - my2)/20;
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
		if (cnorm(b.b) < 1)
			n = .75 + 2.*absd(log(max(3.*s.sx/(Mx-mx), 3.*s.sy/(My-my)))/log(cnorm(b.b)));
		else
			n = .75 + 2.*absd(log(max(.3*s.sx/(Mx-mx), .3*s.sy/(My-my)))/log(cnorm(b.b)));
		if (n < 0)
			n = 0;
		if (n > 1000000)
			n = 100;
		if (verb)
			printf("n = %d\n", n);
	}
	pos = zero();
	Draw_rec2 (b, s, n, b.a.i);
}

void printAutomate(Automaton a)
{
	printf("Automate ayant %d sommets, %d lettres.\n", a.n, a.na);
	int i,j;
	for (i=0;i<a.n;i++)
	{
		for (j=0;j<a.na;j++)
		{
			if (a.e[i].f[j] != -1)
				printf("%d -(%d)-> %d, ", i, j, a.e[i].f[j]);
		}
	}
	printf("état initial %d.\n", a.i);
}

/*
#used by Approx
int Approx_rec (Automaton *a, test, f, x, int n, int n2)
{
	if (n == 0)
	{
		if test(x):
			return f
		else:
			return -1
	}else:
		e1 = self.Approx_rec(a, test, f, x, n-1, n2)
		e2 = self.Approx_rec(a, test, f, x + self.b**(n2-n), n-1, n2)
		if e1 != -1 or e2 != -1:
			e3 = a.add_state(0)
			if e1 != -1:
				a.add_edge(e3, 0, e1)
			if e2 != -1:
				a.add_edge(e3, 1, e2)
			return e3
		return -1
}
	
#gives a automaton describing a approximation of a set defined by the characteritic function test
Automaton ApproxImage (BetaAdic b, SDLImage s, int n)
{
	Automaton r;
	f = a.add_state(1)
	e = Approx_rec(&r, test, f, 0, n, n)
	a.add_edge(f, 0, f)
	a.add_edge(f, 1, f)
	a.set_initial_state(e)
	return r;
}
*/

void DrawList (BetaAdic2 b, Surface s, int n, int ajust, ColorList cl, double alpha, int verb)
{
	int auto_n = (n < 0);
	//set global variables
	mx2 = 1000000, my2 = 1000000, Mx2 = -1000000, My2 = -1000000; //extremum observés
	color0.r = color0.g = color0.b = 255;
	color0.a = 0;
	colors = (Color *)malloc(sizeof(Color)*b.na);
	int i;
	for (i=0;i<b.na;i++)
	{
		colors[i] = cl[i]; //randCol(255);
		colors[i].a = cl[i].a*alpha;
		if (verb)
		{
			printf("couleur %d : %d %d %d %d\n", i, cl[i].r, cl[i].g, cl[i].b, cl[i].a);
		}
	}
	
	if (cnorm(b.b) < 1)
	{
		//calcul du majorant
		Maj = 0;
		double m;
		for (i=0;i<b.n;i++)
		{
			m = cnorm(b.t[i]);
			if (Maj < m)
				Maj = m;
		}
		Maj = sqrt(Maj)/(1. - sqrt(cnorm(b.b)));
	}
	
	if (verb)
	{
		printf("%d translations, %d automates.\n", b.n, b.na);
		printf("couleur de fond : %d %d %d %d\n", color0.r, color0.g, color0.b, color0.a);
		printf("translations : ");
		for (i=0;i<b.n;i++)
		{
			printf("(%lf, %lf) ", b.t[i].x, b.t[i].y);
		}
		printf("\n");
		printf("Automates :\n");
		for (i=0;i<b.na;i++)
		{
			printAutomate(b.a[i]);
		}
	}
	int *etat = (int *)malloc(sizeof(int)*b.na);
	//ajust the window of the drawing
	if (ajust)
	{
		//ajuste le cadre
		if (auto_n)
		{
			//printf("max = %d\n", max(s.sx, s.sy));
			//printf("maj = %lf\n", (1.-sqrt(cnorm(b.b))));
			if (cnorm(b.b) < 1)
				n = .5 + -3.*log(max(s.sx, s.sy)*(1.-sqrt(cnorm(b.b))))/log(cnorm(b.b));
			else
				n = .5 + 3.*log(max(s.sx, s.sy)*(1.-1./sqrt(cnorm(b.b))))/log(cnorm(b.b));
			if (n < 0)
				n = 0;
			if (n > 1000000)
				n = 100;
			if (verb)
				printf("n = %d\n", n);
		}
		pos = zero();
		
		for (i=0;i<b.na;i++)
		{
			etat[i] = b.a[i].i;
		}
		DrawList_rec (b, s, n, zero(), un(), etat);
		mx = mx2 - (Mx2 - mx2)/20;
		my = my2 - (My2 - my2)/20;
		Mx = Mx2 + (Mx2-mx2)/s.sx + (Mx2 - mx2)/20;
		My = My2 + (My2-my2)/s.sy + (My2 - my2)/20;
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
		if (cnorm(b.b) < 1)
			n = .5 + 2.*absd(log(max(3.*s.sx/(Mx-mx), 3.*s.sy/(My-my)))/log(cnorm(b.b)));
		else
			n = .5 + 2.*absd(log(max(.3*s.sx/(Mx-mx), .3*s.sy/(My-my)))/log(cnorm(b.b)));
		if (n < 0)
			n = 0;
		if (n > 1000000)
			n = 100;
		if (verb)
			printf("n = %d\n", n);
	}
	pos = zero();
	for (i=0;i<b.na;i++)
	{
		etat[i] = b.a[i].i;
	}
	DrawList_rec (b, s, n, zero(), un(), etat);
	if (verb)
		printf("Fin de DrawList...\n");
	free(etat);
}
