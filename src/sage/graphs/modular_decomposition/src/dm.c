/******************************************************

Copyright 2004, 2010 Fabien de Montgolfier
fm@liafa.jussieu.fr

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
**********************************************************/

/********************************************************

    DECOMPOSITION MODULAIRE DE GRAPHES NON-ORIENTES

Cet algorithme construit l'arbre de decomposition modulaire
d'un graphe donne sous forme d'une matrice d'adjacence.
Il s'effectue en temps O(m log n) temps et $O(m)$ espace.
Il est la concatenation de deux algorithmes distincts.
Le premier realise une permutation factorisante des sommets du graphe
(pour cette notion, cf these de Christian Capelle)
grace a une technique d'affinage de partitions (cf Habib, Paul & Viennot)
Le second construit l'arbre de decomposition modulaire a partir de cette
+permutation, cf
[Anne Bergeron, Cedric Chauve, Fabien de Montgolfier, Mathieu Raffinot,
 Computing Common Intervals of K Permutations, with Applications to Modular
 Decomposition of Graphs. SIAM J. Discrete Math. 22(3): 1022-1039 (2008)]

Montpellier, decembre 2000 et Paris, novembre 2010
********************************************************/

#include "dm_english.h"
#include <stdio.h>
#include <stdlib.h>

#define DEBUG 0                /* si 0 aucune sortie graphique!! */

/* dm.h definit les constantes FEUILLE, MODULE, etc...
ainsi que les structures noeud et fils. Les autres
structures n'ont pas a etre connues par les programmes
exterieurs et sont donc definies ici. */


/* un sommet du graphe
(utilise dans la premiere partie seulement, ainsi que ce qui suit)*/
typedef struct Sommet {
  int place;/* numero du sommet dans la NOUVELLE numerotation */
  int nom;  /* numero du sommet dans l'ANCIENNE numerotation */
  /* On a donc sigma(nom)=place */
  struct Sadj *adj;
  struct SClasse *classe;
} sommet;

/* liste d'adjacence d'un sommet, DOUBLEMENT chainee*/
typedef struct Sadj {
  struct Sommet *pointe;
  struct Sadj *suiv;
  struct Sadj *prec;
} sadj;

/* classe de la partition courante,
   organisees en liste chainnee selon l'ordre de la partition */
typedef struct SClasse {
    int debut;
    int fin;
    struct Sommet *firstpivot;
    int inpivot;                /*indice de la classe dans le tableau pivot */
    int inmodule;                /* (resp module); -1 si non present */
    int whereXa;                /* lie le couple X/Xa: vaut
                                   0 si X n'est actuellement lie a aucun Xa
                                   -1 si Xa est a gauche
                                   +1 si Xa est a droite */
    struct SClasse *suiv;        /* forment une liste chainee */
    struct SClasse *prec;        /*...doublement */
} sclasse;

/* plein de parametres statiques que algo1() donne a Raffine() */
typedef struct Info {
  sclasse **pivot;
  int *ipivot;
  sclasse **module;
  int *imodule;
  int *numclasse;
  int *n;
} info;

/* clef a deux entrees utilisee pour le tri lineaire
   represente l'arrete ij */
typedef struct Clef2{
  int i; //sommet pointeur
  int nom; // nom du sommet pointe
  int place; //place du sommet pointe
} clef2;

/*************************************************************
utilitaires
*************************************************************/
void *fabmalloc(size_t s)
/* malloc sans erreur */
{
  void *p;
  p=malloc(s);
  if(p==NULL)
    {
      perror("Erreur de malloc!\n");
      exit(1);
    }
  return p;
}

int min(int a, int b)
{
  return (a<b) ? a : b;
}

int max(int a, int b)
{
    return (a > b) ? a : b;
}
/**************************************************************
Premiere passe de l'algorithme: il s'agit de trouver une
permutation factorisante du graphe. Nous utilisons les
techniques de raffinement de partition. Tout cela est
explique dans l'article de Habib, Viennot & Paul, dont je ne
fais ici que transcrire le travail.
****************************************************************/

void printS(sommet ** S, int n)
{
  /* imprimme S selon S et selon les classes */
  int i;
  sclasse *s;

  for (s = S[0]->classe; s != NULL; s = s->suiv) {
    printf("[ ");
    for (i = s->debut; i <= s->fin; i++)
      printf("%i ", S[i]->nom);
    printf("] ");
  }
  printf("\n");
}

sclasse *nouvclasse(sclasse * un, sclasse * deux)
{
  /* cree une nouvelle classe et l'insere entre un et deux;
     on suppose que si un et deux sont pas NULL alors
     FORCEMENT un=deux->suiv */

  sclasse *nouv;
  nouv = (sclasse *) fabmalloc(sizeof(sclasse));
  nouv->whereXa = 0;
  nouv->inpivot = -1;
  nouv->inmodule = -1;
  nouv->firstpivot = NULL;
  nouv->prec = un;
  if (un != NULL)                /* accroche pas en tete de chaine */
    un->suiv = nouv;
  nouv->suiv = deux;
  if (deux != NULL)                /* pas en queue de chaine */
    deux->prec = nouv;

    /* debut et fin ne sont PAS initialises! */
  return nouv;
}

void permute(sommet ** S, int a, int b)
{
  /* transpose les sommets a et b dans S */
  /* ne touche pas aux classes! */
  sommet *tmp;
  S[a]->place = b;
  S[b]->place = a;
  tmp = S[a];
  S[a] = S[b];
  S[b] = tmp;
}

void Raffiner(sommet ** S, sommet * p, sommet * centre, info * I)
{
  /* melange raffiner, pivotset, insertright et addpivot */
  sadj *a;                        /* parcours l'adjacence du pivot */
  sommet *x;                        /* sommet quiva changer de classe */
  sclasse *X, *Xa;                /* x in X; Xa nouv classe de x */
  sclasse *Z;
  sclasse **pivot;
  sclasse **module;
  int *ipivot, *imodule, *numclasse, n;

  if (DEBUG)
    printf("Raffinage avec le pivot %i\n", p->nom);
  pivot = I->pivot;
  module = I->module;
  ipivot = I->ipivot;
  imodule = I->imodule;
  numclasse = I->numclasse;
  n = *(I->n);

  for (a = p->adj; a != NULL; a = a->suiv) {
    x = a->pointe;
    X = x->classe;
    if (X == p->classe)
      continue;                /* on raffine pas la classe du pivot! */

    if (X->whereXa == 0) {
      /* c'est la premiere fois qu'on trouve un x
         appartenant a X lors de cet appel a raffiner */

      if ((centre->place < x->place && x->place < p->place)
          || (p->place < x->place && x->place < centre->place)) {
        /* insere a gauche */
        Xa = nouvclasse(X->prec, X);
        (*numclasse)++;
        permute(S, x->place, X->debut);
        X->debut++;
        X->whereXa = -1;
        Xa->whereXa = 1;        /* besoin dans le second tour */
      }
      else {                /* insere a droite */

        Xa = nouvclasse(X, X->suiv);
        (*numclasse)++;
        permute(S, x->place, X->fin);
        X->fin--;
        X->whereXa = 1;
        Xa->whereXa = -1;
      }
      x->classe = Xa;
      Xa->debut = x->place;
      Xa->fin = x->place;
    }
    else {
      if (X->whereXa == -1) {
        Xa = X->prec;
        permute(S, x->place, X->debut);
        X->debut++;
        Xa->fin++;
      }
      else {
        Xa = X->suiv;
        permute(S, x->place, X->fin);
        X->fin--;
        Xa->debut--;
      }
      x->classe = Xa;
    }
  }

  for (a = p->adj; a != NULL; a = a->suiv)
    /* deuxieme couche! Maintenant on va faire les addpivot,
           et remettre les whereXa a 0
           Noter qu'on lit les Xa et plus les X */
    {
      x = a->pointe;
      Xa = x->classe;
      if (Xa->whereXa == 0)
        continue;                /* deja remis a zero! */
      if (Xa->whereXa == -1)
        X = Xa->prec;
      else
        X = Xa->suiv;

      if (X->debut > X->fin) {
        /*on a trop enleve! X est vide
          -> on va le supprimer mechamment */

        (*numclasse)--;
        if (Xa->whereXa == 1) {        /*deconnecte */
          Xa->suiv = X->suiv;
          if (Xa->suiv != NULL)
            Xa->suiv->prec = Xa;
        }
        else {
          Xa->prec = X->prec;
          if (Xa->prec != NULL)
            Xa->prec->suiv = Xa;
        }
        Xa->inpivot = X->inpivot;
        if (X->inpivot != -1)        /* ecrase X dans pivot */
          pivot[X->inpivot] = Xa;
        Xa->inmodule = X->inmodule;
        if (X->inmodule != -1)        /* ecrase X dans pivot */
          module[X->inmodule] = Xa;

        Xa->whereXa = 0;
        continue;
      }

      /* Maintenant on fait addpivot(X,Xa)
         noter que X et Xa sont non vides */

      if (X->inpivot == -1) {
        if ((X->inmodule != -1)
            && (X->fin - X->debut < Xa->fin - Xa->debut)) {
          /* remplace X par Xa dans module */
          module[X->inmodule] = Xa;
          Xa->inmodule = X->inmodule;
          X->inmodule = -1;
          if (DEBUG)
            printf("Dans module %i-%i ecrase %i-%i\n",
                    S[Xa->debut]->nom, S[Xa->fin]->nom,
                    S[X->debut]->nom, S[X->fin]->nom);
        }
        else {
          if (X->inmodule == -1) {
            if (X->fin - X->debut < Xa->fin - Xa->debut)
              Z = Xa;
            else
              Z = X;
            /* ajoute Z (=max(X,Xa)) a module */
            module[(*imodule)] = Z;
            Z->inmodule = (*imodule);
            (*imodule)++;
            if (DEBUG)
              printf("module empile:%i-%i\n",
                      S[Z->debut]->nom, S[Z->fin]->nom);
          }
        }
      }

      if (X->inpivot != -1)
        Z = Xa;
      else if (X->fin - X->debut < Xa->fin - Xa->debut)
        Z = X;
      else
        Z = Xa;
      /* on empile Z dans pivot */
      pivot[(*ipivot)] = Z;
      Z->inpivot = (*ipivot);
      (*ipivot)++;
      if (DEBUG)
        printf("pivot empile: %i-%i\n", S[Z->debut]->nom,
                S[Z->fin]->nom);
      X->whereXa = 0;
      Xa->whereXa = 0;
    }
  if (DEBUG) {
    printS(S, n);
    printf("\n");
  }
}

sommet **algo1(graphe G)
     /* Entree: un graphe G
        Sortie: une permutation factorisante de G,
        donnee sous la forme d'un tableau de structures Sommet ordonnees selon sigma.
        d'apres le travail de Habib/Paul/Viennot */
{
  int n; // nombre de sommets de G

  sclasse **pivot;                /*pile des pivots */
  int ipivot = 0;                /*indice sur la precedante */

    sclasse **module;                /*idem, modules */
    int imodule = 0;

    sclasse *singclasse;
    /*invariant: toute classe avant singclasse dans la chaine */
    /*a un seul element */
    int numclasse;                /* quand vaut n, on a fini!! */

    sclasse *C1;                /*premiere classe, tete de la chaine */
    sclasse *Y;                        /*classe qui raffine */
    sclasse *X;                        /*classe raffinee */
    sclasse *Xc;                /* morceau utilise de X */
    /* sclasse *Xa, *Xc; morceau inutilise de X */
    sommet *x;                        /* x in X */
    sommet *y;                        /* y in Y */
    sommet *centre;                /* le centre du raffinage actuel */

    sommet **S;                        /*la permutation factorisante ! */

    int i, j;                        /*divers indices */
    sommet *scourant;                /* pour l'init */
    sadj *nextadj;                /*sommet adjacent suivant */
    adj *nextadj2;              /* idem mais de type adj */
    info Inf;                        /* diverses info a passer a raffiner */

    /* debut des initialisations */
    n=G.n;
    /*initialisation des tableaux */
    module = (sclasse **) fabmalloc(n * sizeof(sclasse *));
    pivot = (sclasse **) fabmalloc(n * sizeof(sclasse *));
    S = (sommet **) fabmalloc(n * sizeof(sommet *));
    /* on va initialiser la permutation factorisante,
       ainsi que chaque structure sommet */
    C1 = nouvclasse(NULL, NULL);
    numclasse = 1;
    singclasse = C1;
    C1->debut = 0;
    C1->fin = n - 1;
    for (i = 0; i < n; i++) {
        /* initialisation des sommets */
        /* notre bebe est le sommet i dans M */
        scourant = (sommet *) fabmalloc(sizeof(struct Sommet));
        scourant->nom = i;
        scourant->place = i;        /* a ce point S=identite */
        scourant->adj = NULL; /* pas encore d'adjacence */
        scourant->classe = C1;
        S[i] = scourant;
    }
    for (i = 0; i < n; i++)
      {
        nextadj2 = G.G[i];
        while(nextadj2 != NULL)
          {
            j=nextadj2->s; //numero du sommet pointe
            if((j<0)||(j>=n))
              {
                perror("Graphe invalide (numero de sommet erronne)!\n");
                exit(1);
              }
            nextadj = (sadj *) fabmalloc(sizeof(struct Sadj));
            //un nouveau sadj
            nextadj->pointe = S[j];
            nextadj->suiv = S[i]->adj; //tete de liste
            if(nextadj->suiv!=NULL)
              nextadj->suiv->prec=nextadj;
            nextadj->prec=NULL;
            S[i]->adj = nextadj;        /*et le tour est joue */
            nextadj2=nextadj2->suiv;
          }
        }
    /* NB: module et pivot sont vides */
    Inf.pivot = pivot;
    Inf.ipivot = &ipivot;
    Inf.module = module;
    Inf.imodule = &imodule;
    Inf.numclasse = &numclasse;
    Inf.n = &n;
    /* init terminnee */

    while (1) {
        while (ipivot > 0 || imodule > 0) {
            while (ipivot > 0) {
                /*cette boucle raffine selon tous les sommets
                   de la premiere classe dans pivot */

                Y = pivot[ipivot - 1];
                ipivot--;
                Y->inpivot = -1;

                for (i = Y->debut; i <= Y->fin; i++)
                    Raffiner(S, S[i], centre, &Inf);

                /* une optimisation de la fin de l'algo */
                if (numclasse == n)
                    return (S);
            }
            /*maintenant pivot est vide, mais peut-etre pas module */
            if (imodule > 0) {
                /* relance par un sommet (pas au pif...) */
                /* de chaque module qui le represente */
                Y = module[imodule - 1];
                imodule--;
                Y->inmodule = -1;
                y = S[Y->debut];        /* le firstpivot sera toujours... */
                Y->firstpivot = y;        /* le premier!! */
                if (DEBUG)
                    printf("module-pivot %i-%i: sommet %i\n",
                            S[Y->debut]->nom, S[Y->fin]->nom,
                            y->nom);
                Raffiner(S, y, centre, &Inf);
            }
        }
        /* a ce point, pivot et module sont vides...
           pas de pb! On va faire initpartition HERE */
        if (DEBUG)
            printf("\nInit Partition\n");
        /**** ajoute ici pour debbugger, mais moche!! */
        singclasse = S[0]->classe;
        while ((singclasse != NULL) &&
               (singclasse->debut == singclasse->fin))
          {
            singclasse = singclasse->suiv;
          }
        /* singclasse est la premiere classe
           non singlette, sauf si: */
        if (singclasse == NULL)
            /* on a n classes singlettes? ben c'est gagne! */
          {
            return (S);
          }
        if (singclasse == NULL && numclasse < n) {
            perror("c'est pas normal! Ca termine trop vite!\n");
            exit(1);
        }

        X = singclasse;
        x = X->firstpivot;
        if (x == NULL)
            x = S[X->debut];
        else                        /* remet firstpivot a NULL!! */
            X->firstpivot = NULL;

        if (DEBUG)
            printf("Relance dans le module %i-%i avec le sommet %i\n",
                    S[X->debut]->nom, S[X->fin]->nom, x->nom);

        centre = x;                /*important! */
        /* astuce: on place {x} en tete de X
           ensuite, on raffine S selon x -> seule X est coupee
           il y a alors {x} X Xa
           -> on met {x} en queue de X et c'est bon!
           ainsi on a bien nonvoisins-x-voisons */
        Xc = nouvclasse(X->prec, X);
        numclasse++;
        x->classe = Xc;
        permute(S, x->place, X->debut);
        X->debut++;
        Xc->debut = x->place;
        Xc->fin = x->place;
        Raffiner(S, x, x, &Inf);
        /* X existe-il encore? */
        if (X->debut > X->fin)
            continue;
        /* echange de x et {x}. Init: -{x}-X- */
        Xc->suiv = X->suiv;
        if (X->suiv != NULL)
            X->suiv->prec = Xc;
        X->prec = Xc->prec;
        if (Xc->prec != NULL)
            Xc->prec->suiv = X;
        X->suiv = Xc;
        Xc->prec = X;
        permute(S, x->place, X->fin);
        Xc->debut = x->place;
        Xc->fin = x->place;
        X->debut--;
        X->fin--;
        //antibug?
        singclasse=X;
        /* now -X-{x}- */
        if (DEBUG)
            printS(S, n);
    }
  free(module);
  free(pivot);
}

/***************************************************************
Etape intermediaire: trier toutes les listes d'adjacence
selon S. ce sont les listes de type sadj qui sont concernees
***************************************************************/
int Calculm(graphe G)
/* compte le nombre d'arretes du graphe */
{
  int i,r; adj *a;
  r=0;
  for(i=0;i<G.n;i++)
    {
      a=G.G[i];
      while(a!=NULL)
        {
          a=a->suiv;
          r++;
        }
    }
  if(r%2!=0)
    {
      perror("Erreur: nombre impaire d'arrete, graphe non-oriente??\n");
      exit(1);
    }
  return r/2; // G symetrique!
}

void TrierTous(sommet **S, int n, int m)
/* trie chaque liste d'adjacence de S*/
{
  //n sommets, m arretes
  int i; // numero du sommet courant
  sadj *a,*atmp;// parcours sa liste d'adjacence
  clef2 *c; // enregistrement a trier
  int *tab1; clef2 **tab2; //tableaux du tri par seaux
  tab1=(int *)fabmalloc(n*sizeof(int));
  tab2=(clef2 **)fabmalloc(m * 2 * sizeof(clef2 *));
  for(i=0; i<n; i++)
    tab1[i]=0;

  // premiere passe: construit tab1:
  // tab[i] est la frequence du ieme (selon S) sommet
  for(i=0; i<n; i++)
    {
      a=S[i]->adj;
      while(a!=NULL)
        {
          tab1[i]++;
          a=a->suiv;
        }
    }
  //deuxieme passe: frequences cumulees a rebours
  // (car les listes d'adjacences se construisent a l'envers
  //tab1[n-1]--; // a cause des indices de tableau qui commence a zero
  //for(i=n-1;i>0;i--)
  //  tab1[i-1]+=tab1[i];

  //deuxieme passe: frequences cumulees
  for(i=1;i<n;i++)
    tab1[i]+=tab1[i-1];

  //troisieme passe: liste double
  for(i=0; i<n; i++)
    {
      a=S[i]->adj;
      while(a!=NULL)
        {
          /* cree un nouveau record */
          c=(clef2 *)fabmalloc(sizeof(struct Clef2));
          c->i=i;
          c->nom=a->pointe->nom;
          c->place=a->pointe->place;
          /* le place bien dans tab2 */
          tab1[c->place]--;
          tab2[tab1[c->place]]=c;
          /*et on continue */
          a=a->suiv;
        }
    }

  //quatrieme passe: detruit les vielles listes d'adjacence
  for(i=0; i<n; i++)
    {
      a=S[i]->adj;
      while(a!=NULL)
        {
          atmp=a->suiv;
          free(a);
          a=atmp;
        }
      S[i]->adj=NULL;
    }

  //derniere passe: reconstruit les listes d'adjacence
  for(i=0;i<2*m;i++)
    {
      c=tab2[i];
      a=(sadj *)fabmalloc(sizeof(struct Sadj));
      a->pointe=S[c->i];
      a->suiv=S[c->place]->adj; //insere en tete
      if(a->suiv!=NULL)
        a->suiv->prec=a;
      a->prec=NULL;
      S[c->place]->adj=a;
      //nettoie
      free(c);
   }
  free(tab1);
  free(tab2);
}


/***************************************************************
 Maintenant, la deuxieme partie de l'aglorithme
 On va, etant donne la matrice M construite a l'etape precedante,
 etablir l'arbre de decomposition modulaire.
 Tous les details sont dans mon memoire de DEA
****************************************************************/
noeud *nouvnoeud(int type, noeud * pere, int sommet, int fv, int lv)
{
    /* cree un nouveau noeud. Noter que l'on est oblige de passer n
       comme parametre car les bords et separateurs droits doivent
       etre initilises avec des valeurs >n */
    noeud *nn;
    static int compteur = 0;
    /*pour donner un ID unique aux noeuds. juste pour debug */

    nn = (noeud *) fabmalloc(sizeof(noeud));
    nn->type = type;
    nn->pere = pere;
    /* nn->fpere ne peut etre deja mis a jour... */
    nn->nom = sommet;
    nn->fv=fv;
    nn->lv=lv;
    nn->fils = NULL;
    nn->lastfils = NULL;
    nn->id = compteur;
    compteur++;
    return nn;
}

void ajoutfils(noeud * pere, noeud * nfils)
{
    fils *nf;
    /* noter que c'est un ajout en queue! */
    nf = (fils *) fabmalloc(sizeof(fils));
    nf->pointe = nfils;
    nf->suiv = NULL;
    if (pere->fils == NULL)
        pere->fils = nf;        /* on cree le premier fils */
    else
        pere->lastfils->suiv = nf;        /* on ajoute nf a la chaine */
    pere->lastfils = nf;
    nfils->pere = pere;                /* normalement: redondant,mais... */
    nfils->fpere = nf;
}

void printnoeud(noeud * N, int level)
{
    /* imprime recursivement l'arbre par parcours en profondeur */
    fils *ffils;
    noeud *nfils;
    int i;
    ffils = N->fils;

    for (i = 0; i < level - 1; i++)
        printf("  |");
    if (N->pere == NULL)
        printf(" ");
    else
        printf("  +-");
    switch (N->type) {
    case UNKN:
        printf("Noeud\n");
        break;
    case MODULE:
        printf("Module\n");
        break;
    case ARTEFACT:
        printf("Artefact\n");
        break;
    case SERIE:
        printf("Serie \n");
        break;
    case PARALLELE:
        printf("Parallele \n");
        break;
    case PREMIER:
        printf("Premier \n");
        break;
    }

    do {
        nfils = ffils->pointe;
        if (nfils->type == FEUILLE) {
            for (i = 0; i < level; i++)
                printf("  |");
            printf("  +--");
            printf("%i\n", nfils->nom);
        }
        else {
            printnoeud(nfils, level + 1);
        }
        ffils = ffils->suiv;
    }
    while (ffils != NULL);
}

void printarbre(noeud * N)
{
    printnoeud(N, 0);
}


void typenoeud(noeud * N, int *L2, int *R2, sommet **S)
{
  // type recursivement l'arbre
  if(N->type == FEUILLE)
    return;
  else if (N->type != UNKN)
    {
      printf("Erreur typage!\n");
      return;
    }
  // maintenant N est supposé avoir deux fils
  noeud *premier, *second;
  premier = N->fils->pointe;
  if(premier==NULL)
    {
      printf("Erreur typage!\n");
      return;
    }
  second = N->fils->suiv->pointe;
  if(second==NULL)
    {
      printf("Erreur typage!\n");
      return;
    }

  sadj *a1;
  int i = premier->fv;
  int j = second->lv;
  if(L2[j]<=i && j<= R2[i])
    // ah ah les deux premiers fils forment un module
    // on recherche si adjacents
    {
      N->type = PARALLELE; // sera modifie si adjacence
      a1=S[i]->adj;
      while(a1!=NULL)
        {
          if(a1->pointe->place == S[j]->place)
            // adjacence trouvee !
            {
              N->type = SERIE;
              break;
            }
          a1=a1->suiv;
        }
    }
  else
    N->type = PREMIER;

    fils *ffils;
    ffils = N->fils;
    do {
      typenoeud( ffils -> pointe, L2, R2, S);
      ffils = ffils->suiv;
    }
    while (ffils != NULL);
}


noeud *algo2(graphe G, sommet **S)
{
/* algorithme de decomposition modulaire, deuxieme passe
entree: le graphe G, et sa permutation factorisante S.
sortie: un pointeur sur un arbre de decomposition modulaire
*/
    /* debug: S n'est utilise que pour mettre vrainom a jour */
    int n; //nombre de sommets du graphe

    int *ps;                        /* ps[i]=premier separateur de (i,i+1) */
    int *ds;

    int i, j;                        /*indices de paires ou de sommets */

    sadj *a1, *a2;                /* parcours de liste d'adjacence */


    int *R2, *L2; // le generateur des intervalle-modules produit Capelle-like en phase 2
    int *R, *L; // le  generateur canonique des intervalle-modules
    int *SupR, *SupL; // le  support pour passer de (R2,L2) a (R,L)

    int *pile;  // pile utilisee pour calcule generateur
    int sommet; // sommet de ladite pile. Indique la derniere place PLEINE

    int *tab1,*tab2; // utilise pour le tri lineaire avant algo 6
    int *fortg, *fortd; // bornes de mes modules forts
    int nb_forts=0; // nombre de modules forts

    int *newfortg, *newfortd; // tableau de modules forts utilises temporairement pour bucket sort

    /*PROPHASE : initialisations */
    n=G.n;
    ps = (int *) fabmalloc(n * sizeof(int));
    ds = (int *) fabmalloc(n * sizeof(int));
    R2 = (int *)fabmalloc(n * sizeof(int));
    L2 = (int *)fabmalloc(n * sizeof(int));
    SupR = (int *)fabmalloc(n * sizeof(int));
    SupL = (int *)fabmalloc(n * sizeof(int));
    R = (int *)fabmalloc(n * sizeof(int));
    L = (int *)fabmalloc(n * sizeof(int));
    pile = (int *)fabmalloc(n * sizeof(int));
    tab1 = (int *)fabmalloc(4 * n * sizeof(int));
    tab2 = (int *)fabmalloc(2 * n * sizeof(int));
    fortg=fabmalloc(2 * n * sizeof(int));
    fortd=fabmalloc(2 * n * sizeof(int));
    newfortg= (int *)fabmalloc(2 * n * sizeof(int));
    newfortd=(int *) fabmalloc(2 * n * sizeof(int));

    /* PREMIERE PASSE
       on va parentheser la permutation factorisante.
       complexite: O(m) */
    // utilise ps ds

    for (i = 0; i < n - 1; i++) {
      /*recherche de ps(i,i+1) */
      a1=S[i]->adj;
      a2=S[i+1]->adj;
      while((a1!=NULL) && (a2!=NULL) && (a1->pointe->place<i) &&
            (a2->pointe->place<i) && (a1->pointe->place == a2->pointe->place))
        {
          a1=a1->suiv;
          a2=a2->suiv;
        }

      //arbre de decision complique pour trouver le premier separateur!
      if( ((a1==NULL) && (a2==NULL))
          ||((a1==NULL) &&(a2->pointe->place >= i))
          ||((a2==NULL) && (a1->pointe->place >= i))
          ||((a1!=NULL) && (a2!=NULL) && (a1->pointe->place >= i) && (a2->pointe->place >= i)))
        //pas de separateur
        ps[i]=-1;
      else
        {
          if((a1==NULL) || (a1->pointe->place >= i))
            ps[i]=a2->pointe->place;
          else if((a2==NULL) || (a2->pointe->place >= i))
            ps[i]=a1->pointe->place;
          else
            {
              if((a1->suiv!=NULL)&&(a1->suiv->pointe->place == a2->pointe->place))
                ps[i]=a1->pointe->place;
              else if ((a2->suiv!=NULL)&&(a2->suiv->pointe->place == a1->pointe->place))
                ps[i]=a2->pointe->place;
              else
                ps[i]=min(a1->pointe->place , a2->pointe->place);
            }
        }
      if (DEBUG)
        printf("ps(%i,%i)=%i ps(%i,%i)=%i\n", i , i+1, ps[i],
               S[i]->nom, S[i + 1]->nom, ps[i]==-1?-1:S[ps[i]]->nom );
        /*recherche de ds(i,i+1)
          plus penible encore!*/
      a1=S[i]->adj;
      if(a1!=NULL) // se place en queue de liste.
        while(a1->suiv!=NULL)
          a1=a1->suiv;
      a2=S[i+1]->adj;
      if(a2!=NULL)
        while(a2->suiv!=NULL)
          a2=a2->suiv;
      while((a1!=NULL) && (a2!=NULL) && (a1->pointe->place > i+1) &&
            (a2->pointe->place > i+1) && (a1->pointe->place == a2->pointe->place))
        {
          a1=a1->prec;
          a2=a2->prec;
        }
      if( ((a1==NULL) && (a2==NULL))
          ||((a1==NULL) && (a2->pointe->place <= i+1))
          ||((a2==NULL) && (a1->pointe->place <= i+1))
          ||((a1!=NULL) && (a2!=NULL) && (a1->pointe->place <= i+1) && (a2->pointe->place <= i+1)))
        //pas de separateur
        ds[i]=-1;
      else
        {
          if((a1==NULL) || (a1->pointe->place <= i+1))
            ds[i]=a2->pointe->place;
          else if((a2==NULL) || (a2->pointe->place <= i+1))
            ds[i]=a1->pointe->place;
          else
            {
              if((a1->prec!=NULL)&&(a1->prec->pointe->place == a2->pointe->place))
                ds[i]=a1->pointe->place;
              else if((a2->prec!=NULL)&&(a2->prec->pointe->place == a1->pointe->place))
                ds[i]=a2->pointe->place;
              else
                ds[i]=max(a1->pointe->place , a2->pointe->place);
            }


          //ds[i] = j;
        }
      if (DEBUG)
        printf("ds(%i,%i)=%i ds(%i,%i)=%i\n", i,i+1,ds[i],
               S[i]->nom, S[i + 1]->nom, ds[i]==-1?-1:S[ds[i]]->nom );
    }

    /* DEUXIEMME PASSE
       Calcule le generateur des interval-modules comme algo 9 du papier DAM
    */
    int v;
    // utilise ps ds L2 R2 pile

    // algo 9
    // Attention s[v] dans l'algo est  ds[v-1] !!
    sommet = -1;

    for(v = n-1; v>=0; v--)
      if(ds[v-1] != -1)
        {
          L2[v]=v;
          while( pile[sommet] < ds[v-1])
            L2[pile[sommet--]]=v;
        }
      else
        pile[++sommet]=v;

    while(sommet>=0)
      L2[pile[sommet--]]=0;

    sommet = -1;
    // algo 9 symétrique
    //re-attention là c'est ps[v]
    for(v = 0 ; v<=n-1; v++)
      if(ps[v] != -1)
        {
          R2[v]=v;
          while(  pile[sommet] > ps[v])
            R2[pile[sommet--]]=v;
        }
      else
        pile[++sommet]=v;

    while(sommet>=0)
      R2[pile[sommet--]]=n-1;

    if (DEBUG)
      {    printf("L2 = [");
        for(v = 0;v<n;v++)
          printf("%i ",L2[v]);
        printf("]\nR2 = [");
        for(v = 0;v<n;v++)
          printf("%i ",R2[v]);
        printf("]\n");
        for(i=0;i<n;i++)
          for(j=i+1;j<n;j++)
            if(L2[j]<=i && j<= R2[i])
              printf("[%i-%i] module\n",i,j);

      }



    /* TROISIEME PASSE
       generateur canonique. Algos 3 et 5 du papier de DAM
    */

    // Algo 3, R
    pile[0]=0;
    sommet =0;
    for(i=1;i<n;i++)
      {
        while(R2[pile[sommet]] < i)
          sommet --;
        SupR[i]=pile[sommet];
        pile[++sommet]=i;
      }

    // Algo 3, L

    pile[0]=n-1;
    sommet =0;
    for(i=n-2;i>=0;i--)
      {
        while(L2[pile[sommet]] > i)
          sommet --;
        SupL[i]=pile[sommet];
        pile[++sommet]=i;
      }

    if (DEBUG)
      {
        printf("SupL = [");
        for(v = 0;v<n;v++)
          printf("%i ",SupL[v]);
        printf("]\nSupR = [");
        for(v = 0;v<n;v++)
          printf("%i ",SupR[v]);
        printf("]\n");
      }

    // Algo 5, R
    int k;
    R[0]=n-1;
    for(k=1;k<n;k++)
      R[k]=k;
    for(k=n-1;k>=1;k--)
      if(L2[R[k]] <= SupR[k] && R[k] <= R2[SupR[k]])
        R[SupR[k]] = max(R[k], R[SupR[k]]);

    // Algo 5, L
    L[n-1]=0;
    for(k=0;k<n-1;k++)
      L[k]=k;
    for(k=0;k<n-1;k++)
      if(L2[SupL[k]] <= L[k] && SupL[k] <= R2[L[k]])
        L[SupL[k]] = min(L[k], L[SupL[k]]);

    if (DEBUG)
      {    printf("L = [");
        for(v = 0;v<n;v++)
          printf("%i ",L[v]);
        printf("]\nR = [");
        for(v = 0;v<n;v++)
          printf("%i ",R[v]);
        printf("]\n");
        for(i=0;i<n;i++)
          for(j=i+1;j<n;j++)
            if(L[j]<=i && j<= R[i])
              printf("[%i-%i] module\n",i,j);

      }

    /* QUATRIEME PASSE
       strong modules. Algo 6 du papier de DAM
    */
    //1 remplit le tableau des bornes
    // astuce : on met 2a pour a borne gauche et 2a+1 pour a borne droite. Puis tri lineaire

    for(i = 0; i<n; i++)
      {
        // intervalle (L[i]..i)
        tab1[4*i] = 2*L[i];
        tab1[4*i+1] = 2*i+1;
        // intervalle (i..R[i])
        tab1[4*i+2] = 2*i;
        tab1[4*i+3] = 2*R[i]+1;
      }
    // tableau des fréquences

    for( i=0;i<2*n;i++)
      tab2[i]=0;
    for( i=0;i<4*n;i++)
      tab2[tab1[i]]++;
    // tri
    j = 0;
    for( i=0;i<2*n;i++)
      while(tab2[i]-->0)
        {
          tab1[j++] = i;
        }

    // Algo 6 du papier de DAM
    // j'utilise maintenant tab2 comme stack de sommet : sommet
    sommet = -1;


    for( i=0;i<4*n;i++)
      {
        int ai,s; //meme noms de variables que papier
        ai=tab1[i];
        if((ai%2)==0)
          tab2[++sommet] = ai;
        else
          {
            ai/=2;
            s = tab2[sommet--]/2;

            if(nb_forts == 0 || (fortg[nb_forts - 1] != s) || (fortd[nb_forts - 1] != ai))
              {
                fortg[nb_forts] = s;
                fortd[nb_forts] = ai;
                nb_forts++;
              }
          }
      }
    if(DEBUG)
      {
        printf("Modules forts :\n");
        for( i=0;i<nb_forts;i++)
          printf("[%i %i] ", fortg[i], fortd[i]);
        printf("\n");
      }

    // TRI des modules forts par bornes droites decroissantes
    for( i=0;i<n;i++)
      tab1[i]=0;
    for( i=0;i<nb_forts;i++)
      tab1[fortd[i]]++;
    for( i=1;i<n;i++)
      tab1[i]+=tab1[i-1];

    // tab1 = tableau des frequences cumulées
    for( i=0;i<nb_forts;i++)
      {
        int pos = (nb_forts-1)-(tab1[fortd[i]]-1);
        newfortg[ pos ] = fortg[i];
        newfortd[ pos ] = fortd[i];
        tab1[fortd[i]]--;
      }
    if(DEBUG)
      {
        printf("Modules forts :\n");
        for( i=0;i<nb_forts;i++)
          printf("[%i %i] ", newfortg[i], newfortd[i]);
        printf("\n");
      }

   // now tri par bornes gauches croissantes
    for( i=0;i<n;i++)
      tab1[i]=0;
    for( i=0;i<nb_forts;i++)
      tab1[newfortg[i]]++;
    for( i=1;i<n;i++)
      tab1[i]+=tab1[i-1];

    // tab1 = tableau des frequences cumulées
    for( i = nb_forts-1 ; i>=0 ; i--)// parcours de droite a gauche car on retourne tout
      {
        int pos = tab1[newfortg[i]]-1;
        // printf("[%i %i] en pos %i\n",newfortg[i],newfortd[i],pos);
        fortg[ pos ] = newfortg[i];
        fortd[ pos ] = newfortd[i];
        tab1[newfortg[i]]--;
      }
    if(DEBUG)
      {
        printf("Modules forts :\n");
        for( i=0;i<nb_forts;i++)
          printf("[%i %i] ", fortg[i], fortd[i]);
        printf("\n");
      }

    // Finally algo 7 du papier : construit l'arbre
    noeud* racine = nouvnoeud( UNKN, NULL, -1, 0, n-1);

    noeud *F = racine;
    noeud *Nouv;
    for(k=1; k< nb_forts; )
      {
        if(F->fv <= fortg[k] && F->lv >= fortd[k] )
          {
            if(fortg[k]==fortd[k]) // on a une feuille
              Nouv = nouvnoeud(FEUILLE, F, S[fortg[k]]->nom, fortg[k], fortd[k]);
            else
              Nouv = nouvnoeud(UNKN, F, -1, fortg[k], fortd[k]);
            ajoutfils(F, Nouv);
            // printf("nouv [%i %i] de pere [%i %i]\n",
            // fortg[k], fortd[k], F->fv, F->lv);
            F = Nouv;
            k++;

          }
        else
          {
            // printf("remonte car [%i %i] pas fils de [%i %i]\n",
            // fortg[k], fortd[k], F->fv, F->lv);
            F = F->pere;
          }
      }
    if(DEBUG)
      {
        printf("arbre avant typage :\n");
        printarbre(racine);
      }

    /* CINQUIEME PASSE : typage */
    // Methode : on teste grace au generateur la modularite des deux premiers fils.
    // Ensuite si module est serie si adjacence parallele sinon
    typenoeud(racine, L2, R2, S);

    if(DEBUG)
      {
        printf("\n arbre final :\n");
        printarbre(racine);
      }
    /* nettoie un peu */
    free(ps);
    free(ds);
    free(R2);
    free(L2);
    free(SupR);
    free(SupL);
    free(R);
    free(L);
    free(pile);
    free(tab1);
    free(tab2);
    free(fortg);
    free(fortd);
    free(newfortg);
    free(newfortd);


    return racine;
}

void PrintG(graphe G)
/* affiche le graphe */
{
  int i,r; adj *a;
  r=0;
  for(i=0;i<G.n;i++)
    {
      printf("%i : ",i);
      a=G.G[i];
      while(a!=NULL)
        {
          printf("%i ", a->s);
          a=a->suiv;
        }
      printf("\n");
    }
}
void PrintGS(sommet **S, int n)
/* affiche le graphe trie selon S (celui utilise par algo2)*/
{
  int i; sadj *a;
  for(i=0;i<n;i++)
    {
      printf("%i : ",i);
      a=S[i]->adj;
      while(a!=NULL)
        {
          printf("%i ", a->pointe->place);
          a=a->suiv;
        }
      printf("\n");
    }
}

void PrintS2(sommet **S, int n)
     /* affiche la permutation factorisante */
{
  int i;
  printf("Attention dans debug maintenant usage des NOUVEAUX noms\n");
  printf(  "Place (nouveaux noms, nouvelle numerotation) ");
  for(i=0;i<n;i++)
    printf("%3i ",S[i]->place);
   printf("\nNom (ancienne num) :                         ");
  for(i=0;i<n;i++)
    printf("%3i ",S[i]->nom);
  printf("\n");
}

void nettoie(sommet **S, int n)
{
  // Detruit la copie locale du graphe
  // normalement, apres, tous les mallocs sont
  // annules sauf ceux de l'arbre (le resultat)
  int i;
  sommet * scourant;
  sadj *a,*atmp;
  for (i = 0; i < n; i++) {
        /* initialisation des sommets */
        /* notre bebe est le sommet i dans M */
    scourant = S[i];
    a=scourant->adj;
      while(a!=NULL)
        {
          atmp=a->suiv;
          free(a);
          a=atmp;
        }
      free(scourant->classe); // efface toutes les classes car sont singletons
      free(scourant);
  }
  free(S);
}

/* la fonction principale; qui fait pas grand'chose....*/
noeud *decomposition_modulaire(graphe G)
{

    sommet **S;                        /* la permutation factorisante */
    noeud *Racine;                /* le futur arbre de decomposition */

    setbuf(stdout,NULL);

    S = algo1(G);               /* premiere partie: calcul
                                   de la permutation factorisante */

    TrierTous(S,G.n,Calculm(G));/* Trie les listes d'adjacence selon S
                                 */
    if(DEBUG)
      {
        PrintGS(S,G.n);
        PrintS2(S,G.n);
      }
    Racine = algo2(G, S);       /* deuxieme partie: calcul de l'arbre */

    nettoie(S, G.n);

    return Racine;
}
