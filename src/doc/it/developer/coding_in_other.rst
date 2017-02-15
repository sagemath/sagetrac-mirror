.. _chapter-other:

===================================
Usare librerie esterne e interfacce
===================================

Quando si scrive codice per Sage, si usa Python per la struttura di base 
e l'interfaccia. Per velocit\`a, efficienza, o convenienza, puoi implementare 
parti del codice usando uno qualunque dei seguenti linguaggi: :ref:`Cython
<chapter-cython>`, C/C++, Fortran 95, GAP, Common Lisp, Singular, e 
PARI/GP. Puoi anche usare tute le librerie C/C++ incluse in Sage [SageComponents]_. 
E se il tuo codice funziona ma dipende su dei pacchetti opzionali di Sage, puoi 
usare Octave, o anche Magma, Mathematica, o Maple.

In questo capitolo discutiamo le interfacce fra Sage e :ref:`PARI
<section-pari-library>`, :ref:`section-gap` e 
:ref:`section-singular`.


.. _section-pari-library:

La libreria di interfaccia C di PARI
====================================

Ecco una guida passo-passo per aggiungere uove funzioni PARI a Sage. 
Usiamo come esempio la forma di Frobenius di una matrice. Alcuni algoritmi 
pesanti per le matrici di interi sono implementati usando la libreria PARI. 
Per calcolare la forma di Frobenius in PARI, si usa la funzione ``matfrobenius``.

Ci sono 2 modi di interagire con la libreria PARI da Sage. L'interfaccia gp usa 
l'interprete gp. L'interfaccia PARI usa delle chiamate dirette lee funzioni C di 
PARI---questo \`e il modo preferito perch\`e \`e il pi\`u veloce. Quindi questa 
sezione si focalizza sull'utilizzo di PARI.

Aggiungeremo un nuovo metodo alla classe ``gen``. Questa \`e la rappresentazione 
astratta di tuuti gli oggetti della libreria PARI. Ci\`o significa che una volta 
che abbiamo aggiunto un metodo a questa classe, ogni oggetto PARI, che sia un numero, 
un polinomio o una matrice, avr\`a il nostro nuovo metodo. Quindi puoi fare 
``pari(1).matfrobenius()``, ma poich\`e PARI vuole applicare ``matfrobenius`` alle 
matrici, non ai numeri, riceverai un ``PariError`` in questo caso.

La classe ``gen`` \`e definita in
:file:`SAGE_ROOT/src/sage/libs/pari/gen.pyx`, e qui \`e dove aggiungiamo il metodo 
``matfrobenius``::

    def matfrobenius(self, flag=0):
        r"""
        M.matfrobenius(flag=0): Return the Frobenius form of the square
        matrix M. If flag is 1, return only the elementary divisors (a list
        of polynomials). If flag is 2, return a two-components vector [F,B]
        where F is the Frobenius form and B is the basis change so that
        `M=B^{-1} F B`.

        EXAMPLES::

            sage: a = pari('[1,2;3,4]')
            sage: a.matfrobenius()
            [0, 2; 1, 5]
            sage: a.matfrobenius(flag=1)
            [x^2 - 5*x - 2]
            sage: a.matfrobenius(2)
            [[0, 2; 1, 5], [1, -1/3; 0, 1/3]]
        """
        sig_on()
        return self.new_gen(matfrobenius(self.g, flag, 0))

Nota l'uso di :ref:`sig_on() statement <section_sig_on>`.

La chiamata ``matfrobenius`` \`e solo una chiamata alla funzione ``matfrobenius`` 
della libreria C di PARI C con i parametri appropriati.

La chiamata ``self.new_gen(GEN x)`` costruisce un nuovo oggetto Sage ``gen`` 
da un dato ``GEN`` di PARI, dove il ``GEN`` \`e memorizzato come l'attrributo 
``.g``.  A parte questo, ``self.new_gen()`` chiama una macro ``sig_off()`` di 
chiusura ed inoltre pulisce lo stack di PARI, quindi \`e molto conveniente usarla 
in un'istruzione ``return`` come illustrato sopra. Quindi dopo ``self.new_gen()``, 
tuute le ``GEN`` di PARI che non sono convertite in ``gen`` di Sage sono perdute. 
C'\`e anche ``self.new_gen_noclear(GENx)`` che fa la stessa cosa di ``self.new_gen(GEN x)`` 
eccetto che *non* chiama ``sig_off()`` e non pulisce lo stack di PARI.

L'informazione su quale funzione chiamare e come chiamarla possono essere 
recuperate dal manuale utenete di PARI (nota: Sage include la versione in 
sviluppo di PARI, quindi cerca quella versione del manuale dell'utente). 
Cercando la stringa ``matfrobenius`` puoi trovare:

    La sintassi della libreria \`e ``GEN matfrobenius(GEN M, long flag, long v
    = -1)``, dove ``v`` \`e un numero variabile.

In caso ce fossi familiare con gp, nota che la funzione C di PARI pu\`o avere un 
nome differente dalla funzione gp corrispondente (ad esempio, vedi ``mathnf``), 
quindi verifica sempre il manuale.

Possiamo anche aggiungere un metodo ``frobenius(flag)`` alla classe ``matrix_integer``
dove chiamiamo il metodo ``matfrobenius()`` sull'oggetto di PARI associato alla 
matrice dopo aver fatto qualche verifica di correttezza. Poi convertiamo l'output da 
oggetti PARI ad oggetti Sage::

    def frobenius(self, flag=0, var='x'):
        """
        Return the Frobenius form (rational canonical form) of this
        matrix.

        INPUT:

        -  ``flag`` -- 0 (default), 1 or 2 as follows:

            -  ``0`` -- (default) return the Frobenius form of this
               matrix.

            -  ``1`` -- return only the elementary divisor
               polynomials, as polynomials in var.

            -  ``2`` -- return a two-components vector [F,B] where F
               is the Frobenius form and B is the basis change so that
               `M=B^{-1}FB`.

        -  ``var`` -- a string (default: 'x')

        ALGORITHM: uses PARI's matfrobenius()

        EXAMPLES::

            sage: A = MatrixSpace(ZZ, 3)(range(9))
            sage: A.frobenius(0)
            [ 0  0  0]
            [ 1  0 18]
            [ 0  1 12]
            sage: A.frobenius(1)
            [x^3 - 12*x^2 - 18*x]
            sage: A.frobenius(1, var='y')
            [y^3 - 12*y^2 - 18*y]
        """
        if not self.is_square():
            raise ArithmeticError("frobenius matrix of non-square matrix not defined.")

        v = self._pari_().matfrobenius(flag)
        if flag==0:
            return self.matrix_space()(v.python())
        elif flag==1:
            r = PolynomialRing(self.base_ring(), names=var)
            retr = []
            for f in v:
                retr.append(eval(str(f).replace("^","**"), {'x':r.gen()}, r.gens_dict()))
            return retr
        elif flag==2:
            F = matrix_space.MatrixSpace(QQ, self.nrows())(v[0].python())
            B = matrix_space.MatrixSpace(QQ, self.nrows())(v[1].python())
            return F, B



.. _section-gap:

GAP
===

Incapsulare una funzione GAP in Sage \`e questione di scrivere un programma 
in Python che usa l'interfaccia pexpect per passare vari comandi a GAP 
e ricevere l'input all'indietro verso Sage. Questo a volte \`e facile, altre 
meno.

Ad esempio, supposiamo che vogliamo fare un wrapper per il calcolo della matrice 
di Cartan di una semplice algebra di Lie. La matrice di Cartan di `G_2` \`e 
disponibile in GAP usando i comandi::

    gap> L:= SimpleLieAlgebra( "G", 2, Rationals );
    <Lie algebra of dimension 14 over Rationals>
    gap> R:= RootSystem( L );
    <root system of rank 2>
    gap> CartanMatrix( R );

In Sage si pu\`o accede a questi comandi digitando::

    sage: L = gap.SimpleLieAlgebra('"G"', 2, 'Rationals'); L
    Algebra( Rationals, [ v.1, v.2, v.3, v.4, v.5, v.6, v.7, v.8, v.9, v.10,
      v.11, v.12, v.13, v.14 ] )
    sage: R = L.RootSystem(); R
    <root system of rank 2>
    sage: R.CartanMatrix()
    [ [ 2, -1 ], [ -3, 2 ] ]

Nota la ``'"G"'`` che \`e valuatata in GAP come stringa ``"G"``.

Lo scopo di questa sezione \`e usare questo esempio per mostrare come si 
possa scrivere un programma Python/Sage il cui input \`e, diciamo, ``('G',2)`` 
ed il cui output sia la matrice suddetta (ma come tipo matrice di Sage---vedi il 
codice nella directory :file:`SAGE_ROOT/src/sage/matrix/` e le correspondenti 
parti nel manuale di riferimento di Sage).

Innanzitutto l'input dev'essere convertito in stringhe che consistono di 
comandi GAP consentiti. Poi occorre fare il parse dell'output di GAP, che \`e 
anch'esso una stringa, e va convertito se possibile al corrispondente oggetto di 
Sage/Python.

.. skip

::

    def cartan_matrix(type, rank):
        """
        Return the Cartan matrix of given Chevalley type and rank.

        INPUT:
            type -- a Chevalley letter name, as a string, for
                    a family type of simple Lie algebras
            rank -- an integer (legal for that type).

        EXAMPLES:
            sage: cartan_matrix("A",5)
            [ 2 -1  0  0  0]
            [-1  2 -1  0  0]
            [ 0 -1  2 -1  0]
            [ 0  0 -1  2 -1]
            [ 0  0  0 -1  2]
            sage: cartan_matrix("G",2)
            [ 2 -1]
            [-3  2]
        """

        L = gap.SimpleLieAlgebra('"%s"'%type, rank, 'Rationals')
        R = L.RootSystem()
        sM = R.CartanMatrix()
        ans = eval(str(sM))
        MS  = MatrixSpace(QQ, rank)
        return MS(ans)

L'output ``ans`` \`e una lista di Python. Le ultime 2 linee convertono quella 
lista ad un'istanza della classe ``Matrix`` di Sage.

In alternativa si pu\`o rimpiazzare la prima linea della suddetta funzione 
con questa::

        L = gap.new('SimpleLieAlgebra("%s", %s, Rationals);'%(type, rank))

Definire "facile" e "difficile" \`e soggettivo, ma qui c'\`e una definizione.
Fare un wrapper di una funzione GAP \`e "facile" se c'\`e gi\`a una classe 
corrispondente in Python o Sage per il tipo di dato di output della funzione 
GAP che stai cercando di incapsulare. Ad esempio, incapsulare una qualunque 
funzione GUAVA (il pacchetto di codici di correzione di GAP) \`e "facile" poich\`e 
i codici di correzione sono spazi vettoriali su campi finiti e le funzioni GUAVA 
restituiscono uno dei seguenti tipi di dato:

- vettori su campi finiti,

- polinomi su campi finiti,

- matricio su campi finiti,

- gruppi di permutazioni o loro elementi,

- interi.


Sage dispone gi\`a di classi per ciascunodei suddetti.

Un esempio "difficile" \`e lasciato come esercizio! Ecco alcune idee.

- Scrivi un wrapper per la funzione ``FreeLieAlgebra`` di GAP (o, pi\`u in 
  generale, tutte le funzioni di GAP per le algebre di Lie finitamente rappresentate). 
  Questo richieder\`a il creare nuovi oggetti Python.

- Scrivi un wrapper per la funizone ``FreeGroup`` di GAP (o, pi\`u in 
  generale, tutte le funzioni di GAP per gruppi finitamente rappresentati). Questo 
  richieder\`a lo scrivere alcuni nuovi oggetti di Python.

- Scrivi un wrapper per le tabelle di caratteri di GAP. Sebbene questo possa essere 
  fatto senza creare nuovi oggetti di Python, per utilizzare al massimo tali tabelle, 
  \`e probabilmente meglio avere dei nuovi oggetti Python per farlo.


.. _section_libgap:

LibGAP
======

Lo svantaggio di usare altri programmi attraverso delle interfacce \`e che 
c'\`e una certa latenza inevitabile (dell'ordine di 10ms) dovuta all'inviare 
l'input e ricevere il risultato. Se deve chiamare delle funzioni in un piccolo 
ciclo questo pu\`o essere inaccettabilmente lento. Richiamare da una libreria 
condivisa ha molta meno latenza ed inoltre evita di avere da convertire tutto 
in una stringa nel passaggio. \`E per questo che Sage include una versione a 
libreria condivisa del kernel di GAP, disponibile come `libgap` in Sage. 
L'analogo del primo esempio usando libgap in :ref:`section-gap` \`e::

    sage: SimpleLieAlgebra = libgap.function_factory('SimpleLieAlgebra')
    sage: L = SimpleLieAlgebra('G', 2, QQ)
    sage: R = L.RootSystem();  R
    <root system of rank 2>
    sage: R.CartanMatrix()    # output is a GAP matrix
    [ [ 2, -1 ], [ -3, 2 ] ]
    sage: matrix(R.CartanMatrix())   # convert to Sage matrix
    [ 2 -1]
    [-3  2]


.. _section-singular:

Singular
========

Usare funzioni Singular da Sage non \`e molto differente concettualmente 
dall'usare funzioni GAP da Sage. Come con GAP, questo varia da facile a 
difficile, in base a quanto della struttura dati dell'output della funzione 
di Singular \`e gi\`a presente in Sage.

Innanzitutto, un po' di terminologia. Per noi, una *curva* `X` su un campo finito 
`F` \`e un'equazione della forma `f(x,y) = 0`, dove `f \in F[x,y]` \`e un 
polinomio. Pu\`o essere o meno singolare. Un *luogo di grado* `d` \`e un'orbita 
di Galois di `d` punti in `X(E)`, dove `E/F` \`e di grado `d`. Ad esempio, un 
luogo di grado `1` \`e anche un luogo di grado `3`, ma un luogo di grado `2` no 
poich\`e nessuna estensione di grado `3` di `F` contiene un'estensione di grado `2`. 
I luoghi di grado `1` sono anche detti punti `F`-razionali.

Come esempio dell'interfaccia Sage/Singular, spieghiamo come incapsulare ``NSplaces`` 
di Singular, che calcola i luoghi su una curva su un campo finito. (Il comando 
``closed_points`` fa anche questo in taluni casi.) Questo \`e "facile" poich\`e non 
sono necessarie nuove classi Python in Sage per farlo.

Ecco un esempio di come usare questo comando in Singular::

     A Computer Algebra System for Polynomial Computations   /   version 3-0-0
                                                           0<
         by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   May 2005
    FB Mathematik der Universitaet, D-67653 Kaiserslautern    \
    > LIB "brnoeth.lib";
    [...]
    > ring s=5,(x,y),lp;
    > poly f=y^2-x^9-x;
    > list X1=Adj_div(f);
    Computing affine singular points ...
    Computing all points at infinity ...
    Computing affine singular places ...
    Computing singular places at infinity ...
    Computing non-singular places at infinity ...
    Adjunction divisor computed successfully

    The genus of the curve is 4
    > list X2=NSplaces(1,X1);
    Computing non-singular affine places of degree 1 ...
    > list X3=extcurve(1,X2);

    Total number of rational places : 6

    > def R=X3[1][5];
    > setring R;
    > POINTS;
    [1]:
       [1]:
          0
       [2]:
          1
       [3]:
          0
    [2]:
       [1]:
          -2
       [2]:
          1
       [3]:
          1
    [3]:
       [1]:
          -2
       [2]:
          1
       [3]:
          1
    [4]:
       [1]:
          -2
       [2]:
          -1
       [3]:
          1
    [5]:
       [1]:
          2
       [2]:
          -2
       [3]:
          1
    [6]:
       [1]:
          0
       [2]:
          0
       [3]:
          1

Ecco un altro modo di fare lo stesso calcolo nell'interfaccia Sage a Singular::

    sage: singular.LIB("brnoeth.lib")
    sage: singular.ring(5,'(x,y)','lp')
        //   characteristic : 5
        //   number of vars : 2
        //        block   1 : ordering lp
        //                  : names    x y
        //        block   2 : ordering C
    sage: f = singular('y^2-x^9-x')
    sage: print singular.eval("list X1=Adj_div(%s);"%f.name())
    Computing affine singular points ...
    Computing all points at infinity ...
    Computing affine singular places ...
    Computing singular places at infinity ...
    Computing non-singular places at infinity ...
    Adjunction divisor computed successfully
    <BLANKLINE>
    The genus of the curve is 4
    sage: print singular.eval("list X2=NSplaces(1,X1);")
    Computing non-singular affine places of degree 1 ...
    sage: print singular.eval("list X3=extcurve(1,X2);")
    <BLANKLINE>
    Total number of rational places : 6
    <BLANKLINE>
    sage: singular.eval("def R=X3[1][5];")
    ''
    sage: singular.eval("setring R;")
    ''
    sage: L = singular.eval("POINTS;")

    sage: print L
    [1]:
       [1]:
          0
       [2]:
          1
       [3]:
          0
    [2]:
       [1]:
          -2
       [2]:
          -1
       [3]:
          1
    ...

Guardando all'output, notiamo che la nostra funzione wrapper dovr\`a fare il 
parse della stringa rappresentata dalla `L` suddetta, quindi scriviamo una 
funzione separata per fare ci\`o. Questo richiede il capire come determinare 
dove sonoposte le coordinate dei punti nella stringa `L`. Python ha alcuni 
comandi per la manipolazione di stringhe molto utili per fare questo.

.. skip

::

    def points_parser(string_points,F):
        """
        This function will parse a string of points
        of X over a finite field F returned by Singular's NSplaces
        command into a Python list of points with entries from F.

        EXAMPLES:
            sage: F = GF(5)
            sage: points_parser(L,F)
            ((0, 1, 0), (3, 4, 1), (0, 0, 1), (2, 3, 1), (3, 1, 1), (2, 2, 1))
        """
        Pts=[]
        n=len(L)
        #print n
        #start block to compute a pt
        L1=L
        while len(L1)>32:
            idx=L1.index("     ")
            pt=[]
            ## start block1 for compute pt
            idx=L1.index("     ")
            idx2=L1[idx:].index("\n")
            L2=L1[idx:idx+idx2]
            #print L2
            pt.append(F(eval(L2)))
            # end block1 to compute pt
            L1=L1[idx+8:] # repeat block 2 more times
            #print len(L1)
            ## start block2 for compute pt
            idx=L1.index("     ")
            idx2=L1[idx:].index("\n")
            L2=L1[idx:idx+idx2]
            pt.append(F(eval(L2)))
            # end block2 to compute pt
            L1=L1[idx+8:] # repeat block 1 more time
            ## start block3 for compute pt
            idx=L1.index("     ")
            if "\n" in L1[idx:]:
                idx2=L1[idx:].index("\n")
            else:
                idx2=len(L1[idx:])
            L2=L1[idx:idx+idx2]
            pt.append(F(eval(L2)))
            #print pt
            # end block3 to compute pt
            #end block to compute a pt
            Pts.append(tuple(pt))  # repeat until no more pts
            L1=L1[idx+8:] # repeat block 2 more times
        return tuple(Pts)

Ora \`e facile mettere insieme questi ingredienti in una funzione Sage 
che prende come input unaa tripla `(f,F,d)`: un polinomio `f` in
`F[x,y]` che definisce `X:\  f(x,y)=0` (nota che bisogna usare la variabile 
`x,y`), un campo finito `F` *di ordine primo*, ed il grado `d`. L'output 
\`e il numero di luoghi in `X` di grado `d=1` su `F`. Al momento, non c'\`e 
una "traduzione" fra gli elementi di `GF(p^d)` in Singular e Sage a meno che 
`d=1`. Quindi dobbiamo limitarci a punti di grado uno.

.. skip

::

    def places_on_curve(f,F):
        """
        INPUT:
            f -- element of F[x,y], defining X: f(x,y)=0
            F -- a finite field of *prime order*

        OUTPUT:
            integer -- the number of places in X of degree d=1 over F

        EXAMPLES:
            sage: F=GF(5)
            sage: R=PolynomialRing(F,2,names=["x","y"])
            sage: x,y=R.gens()
            sage: f=y^2-x^9-x
            sage: places_on_curve(f,F)
            ((0, 1, 0), (3, 4, 1), (0, 0, 1), (2, 3, 1), (3, 1, 1), (2, 2, 1))
        """
        d = 1
        p = F.characteristic()
        singular.eval('LIB "brnoeth.lib";')
        singular.eval("ring s="+str(p)+",(x,y),lp;")
        singular.eval("poly f="+str(f))
        singular.eval("list X1=Adj_div(f);")
        singular.eval("list X2=NSplaces("+str(d)+",X1);")
        singular.eval("list X3=extcurve("+str(d)+",X2);")
        singular.eval("def R=X3[1][5];")
        singular.eval("setring R;")
        L = singular.eval("POINTS;")
        return points_parser(L,F)

Nota che l'ordinamento restituito da questa funzione Sage \`e esattamente 
lo stesso dell'ordine nella variabile ``POINTS`` di Singular.

Un ultimo esempio (in aggiunta a quello nella docstring):

.. skip

::

    sage: F = GF(2)
    sage: R = MPolynomialRing(F,2,names = ["x","y"])
    sage: x,y = R.gens()
    sage: f = x^3*y+y^3+x
    sage: places_on_curve(f,F)
    ((0, 1, 0), (1, 0, 0), (0, 0, 1))


Singular: un altro approccio
============================

C'\`e anche un'interfaccia a Singular pi\`u simile a Python. Usando questa, 
il codice \`e molto pi\`u semplice, come vediamo sotto. Innanzitutto, mostriamo 
come calcolare i luoghi su una curva in un caso particolare::

    sage: singular.lib('brnoeth.lib')
    sage: R = singular.ring(5, '(x,y)', 'lp')
    sage: f = singular.new('y^2 - x^9 - x')
    sage: X1 = f.Adj_div()
    sage: X2 = singular.NSplaces(1, X1)
    sage: X3 = singular.extcurve(1, X2)
    sage: R = X3[1][5]
    sage: singular.set_ring(R)
    sage: L = singular.new('POINTS')

Nota che questi elementi di L sono definiti modulo 5 in Singular, e 
si confrontano in modo diverso da quanto ti aspetteresti dalla loro rappresentazione 
scritta:

.. link

::

    sage: sorted([(L[i][1], L[i][2], L[i][3]) for i in range(1,7)])
    [(0, 0, 1), (0, 1, 0), (2, 2, 1), (2, -2, 1), (-2, 1, 1), (-2, -1, 1)]

Poi implementiamo la funzione generale (per brevit\`a omettiamo la docstring, 
che \`e la stessa di sopra). Nota che la funzione ``point_parser`` non \`e richiesta::

    def places_on_curve(f,F):
        p = F.characteristic()
        if F.degree() > 1:
            raise NotImplementedError
        singular.lib('brnoeth.lib')
        R = singular.ring(5, '(x,y)', 'lp')
        f = singular.new('y^2 - x^9 - x')
        X1 = f.Adj_div()
        X2 = singular.NSplaces(1, X1)
        X3 = singular.extcurve(1, X2)
        R = X3[1][5]
        singular.setring(R)
        L = singular.new('POINTS')
        return [(int(L[i][1]), int(L[i][2]), int(L[i][3])) \
                 for i in range(1,int(L.size())+1)]

Questo codice \`e molto pi\`u corto, bello, e leggibile. Comunque dipende 
da certe funzioni, ad esempio ``singular.setring`` che sono state implementate 
nell'interfaccia Sage/Singular, laddove il codice nella sezione precedente usava 
solo una minima parte di quell'interfaccia.


Creare una nuova interfaccia pseudo-TTY
=======================================

Puoi creare delle interfacce Sage pseudo-tty che permettono a Sage di 
lavorare con quasi qualunque progrmma a riga di comando, e che non richiede 
alcuna modifica o estensione a tale programma. Sono anche sorprendentemente 
veloci e flessibile (dato il modo in cui lavorano!), poich\`e tutto l'I/O \`e 
bufferizzato e l'interazione fra Sage ed il programma a riga di comando 
pu\`o essere non-bloccante (asincrono). Un'interfaccia Sage pseudo-tty \`e 
asincrona perch\`e deriva dalla classe Sage ``Expect``, che gestisce la 
comunicazione fra Sage ed il processo esterno.

Ad esempio, ecco una parte del file ``SAGE_ROOT/src/sage/interfaces/octave.py``, 
che definisce un'interfaccia fra Sage e Octave, un progrmma open source per il 
calcolo numerico, fra l'altro::

    import os
    from expect import Expect, ExpectElement

    class Octave(Expect):
        ...

Le prime 2 linee importano la libreria ``os``, che contiene le routine del sistema 
operativo, ed anche la classe ``Expect``, che \`e la classe base delle interfacce. 
La terza linea definisce la classe ``Octave``; anch'essa deriva da ``Expect``. Dopo 
di ci\`o c'\`e la docstring, che omittiamo qui (vedi il file per dettagli). Poi viene::

        def __init__(self, maxread=100, script_subdirectory="", logfile=None,
                     server=None, server_tmpdir=None):
            Expect.__init__(self,
                            name = 'octave',
                            prompt = '>',
                            command = "octave --no-line-editing --silent",
                            maxread = maxread,
                            server = server,
                            server_tmpdir = server_tmpdir,
                            script_subdirectory = script_subdirectory,
                            restart_on_ctrlc = False,
                            verbose_start = False,
                            logfile = logfile,
                            eval_using_file_cutoff=100)

Questo usa la classe ``Expect`` per mettere in piedi l'interfaccia ad Octave::

        def set(self, var, value):
            """
            Set the variable var to the given value.
            """
            cmd = '%s=%s;'%(var,value)
            out = self.eval(cmd)
            if out.find("error") != -1:
                raise TypeError("Error executing code in Octave\nCODE:\n\t%s\nOctave ERROR:\n\t%s"%(cmd, out))

        def get(self, var):
            """
            Get the value of the variable var.
            """
            s = self.eval('%s'%var)
            i = s.find('=')
            return s[i+1:]

        def console(self):
            octave_console()

Questi permattono all'utente di digitare ``octave.set('x', 3)``, dopodich\`e 
``octave.get('x')`` restituisce ``' 3'``. Eseguire ``octave.console()`` mette 
l'utente in una shell interattiva di Octave::

        def solve_linear_system(self, A, b):
            """
            Use octave to compute a solution x to A*x = b, as a list.

            INPUT:
                A -- mxn matrix A with entries in QQ or RR
                b -- m-vector b entries in QQ or RR (resp)

            OUTPUT:
                An list x (if it exists) which solves M*x = b

            EXAMPLES:
                sage: M33 = MatrixSpace(QQ,3,3)
                sage: A   = M33([1,2,3,4,5,6,7,8,0])
                sage: V3  = VectorSpace(QQ,3)
                sage: b   = V3([1,2,3])
                sage: octave.solve_linear_system(A,b)    # optional - octave
                [-0.33333299999999999, 0.66666700000000001, -3.5236600000000002e-18]

            AUTHOR: David Joyner and William Stein
            """
            m = A.nrows()
            n = A.ncols()
            if m != len(b):
                raise ValueError("dimensions of A and b must be compatible")
            from sage.matrix.all import MatrixSpace
            from sage.rings.all import QQ
            MS = MatrixSpace(QQ,m,1)
            b  = MS(list(b)) # converted b to a "column vector"
            sA = self.sage2octave_matrix_string(A)
            sb = self.sage2octave_matrix_string(b)
            self.eval("a = " + sA )
            self.eval("b = " + sb )
            soln = octave.eval("c = a \\ b")
            soln = soln.replace("\n\n ","[")
            soln = soln.replace("\n\n","]")
            soln = soln.replace("\n",",")
            sol  = soln[3:]
            return eval(sol)

Questo codice definisce il metodo ``solve_linear_system``, che lavora come documentato.

Questi sono i soli estratti da ``octave.py``; verifica quel file per pi\`u definizioni 
ed esempi. Guarda gli altri file nella directory ``SAGE_ROOT/src/sage/interfaces/`` per 
esempi di interfacce ad altri pacchetti software.


.. [SageComponents] See http://www.sagemath.org/links-components.html
   for a list
