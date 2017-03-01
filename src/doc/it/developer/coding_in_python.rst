.. _chapter-python:

==================================
Scrivere codice in Python per Sage
==================================

Questo capitolo discute alcune problematiche dello scrivere codice in 
Sage, e suggerimenti su come farlo.


Progettazione
=============

Se stai progettando di sviluppare del nuovo codice per Sage, la
progettazione è importante. Quindi pensa a cosa il tuo programma farà
e come questo va calato nella struttura di Sage. In particolare, molto
di Sage è implementato nel linguaggio orientato agli oggetti Python, e
c'è una gerarchia di classi che organizzano il codice e le
funzionalità. Ad esempio, se implementi degli elementi di un anello,
la tua classe dovrà derivare da
``sage.structure.element.RingElement``, piuttosto che iniziare da
zero.  Prova ad immaginarti come il tuo codice dovrebbe calarsi nel
rimanente codice di Sage, e progettalo di conseguenza.


Funzioni speciali di Sage
=========================

Le funzioni che iniziano e finiscono con 2 caratteri di sottolineatura
(underscore) ``__XXX__`` sono tutte predefinite da Python. Le funzioni
che iniziano e finiscono con 1 carattere di sottolineatura ``_XXX_``
sono definite da Sage. Le funzioni che iniziano con un singolo
underscore sono intese come semi-private, e quelle che iniziano con 2
underscore sono considerate private.  Gli utenti possono creare
funzioni che iniziano e finiscono con underscore.

Così come Python ha molti metodi standard speciali per gli oggetti, 
anche Sage ne ha. Essi sono tipicamente della forma ``_XXX_``.
In pochi casi, l'underscore di fine non è incluso, ma questo sarà 
certamente cambiato così che l'underscore finale sia sempre incluso. 
Questa sezione descrive questi metodi speciali.

Tutti gli oggetti in Sage devono derivare dalla classe di estensione Cython 
``SageObject``::

    from sage.ext.sage_object import SageObject

    class MyClass(SageObject,...):
        ...

o da qualche altra classe già esistente di Sage::

    from sage.rings.ring import Algebra

    class MyFavoriteAlgebra(Algebra):
        ...

Dovresti implementare i metodi ``_latex_`` e ``_repr_`` per ogni oggetto. 
Gli altri metodi dipendono dalla natura dell'oggetto.


Rappresentazione LaTeX
----------------------

Ogni oggetto ``x`` in Sage deve supportare il comando ``latex(x)``, in 
modo che qualunque oggetto Sage possa essere facilmente ed accuratamente 
mostrato a video via LaTeX. Ecco come fare in modo che una classe (e dunque 
le sue istanze) supportino il comando ``latex``.

#. Definire un metodo ``_latex_(self)`` che ritorni una
   rappresentazione LaTeX del tuo oggetto. Dovrebbe essere qualcosa
   che può essere scritto correttamente in ``math mode``. Non
   includere gli ``$`` di apertura e chiusura.

#. Spesso gli oggetti sono costruiti a partire da altri oggetti Sage,
   e questi componenti dovrebbero essere scritti usando la funzione
   ``latex``. Ad esempio, se ``c`` è un coefficiente del tuo oggetto,
   e vuoi scrivere ``c`` usando LaTeX, usa ``latex(c)`` invece di
   ``c._latex_()``, poichè ``c`` può non avere un metodo ``_latex_``,
   e ``latex(c)`` sa come gestire questo caso.

#. Non dimenticarti di includere una docstring ed un esempio che illustri la 
   generazione di LaTeX per il tuo oggetto.

#. Puoi usare qualunque macro inclusa in ``amsmath``, ``amssymb``, o
   ``amsfonts``, o quelle definite in
   ``SAGE_ROOT/doc/commontex/macros.tex``.

Ecco un template di esempio per un metodo ``_latex_``:

.. skip

::

    class X:
       ...
       def _latex_(self):
           r"""
           Return the LaTeX representation of X.

           EXAMPLES::

               sage: a = X(1,2)
               sage: latex(a)
               '\\frac{1}{2}'
           """
           return '\\frac{%s}{%s}'%(latex(self.numer), latex(self.denom))

Come mostrato nell'esempio, ``latex(a)`` produrrà del codice LaTeX che 
rappresenta l'oggetto ``a``. Se si chiama ``view(a)`` verrà mostrata la
versione in caratteri di questo.


Rappresentazione in stampa
--------------------------

Il metodo di stampa standard di Python è ``__repr__(self)``. In Sage,
cioè per gli oggetti che derivano da ``SageObject`` (che è tutto in 
Sage), invece si definisce ``_repr_(self)``. Questo è preferibile poichè 
se definisci solo ``_repr_(self)`` e non anche ``__repr__(self)``, allora gli 
utenti potranno rinominare il tuo oggetto per stamparlo come gli pare. 
Inoltre alcuni oggetti dovrebbero essere stampati differentemente a seconda del 
contesto.

Ecco un esempio delle funzioni ``_latex_`` e ``_repr_`` per la classe ``Pi``. 
Questo è dal file ``SAGE_ROOT/src/sage/functions/constants.py``::

    class Pi(Constant):
        """
        The ratio of a circle's circumference to its diameter.

        EXAMPLES::

            sage: pi
            pi
            sage: float(pi) # rel tol 1e-10
            3.1415926535897931
        """
        ...
        def _repr_(self):
            return "pi"

        def _latex_(self):
            return "\\pi"


Matrix o Vector da Object
-------------------------

Occorre fornire un metodo ``_matrix_`` per un oggetto di cui può
essere fatta la coercion a matrix sul `R`se si vuole che la funzione
Sage ``matrix`` lavori su tale oggetto.

Il seguito è preso da ``SAGE_ROOT/src/sage/graphs/graph.py``::

    class GenericGraph(SageObject):
        ...
        def _matrix_(self, R=None):
            if R is None:
                return self.am()
            else:
                return self.am().change_ring(R)


        def adjacency_matrix(self, sparse=None, boundary_first=False):
            ...

Similmente, occorre fornire un metodo ``_vector_`` per un oggetto di
cui può essere fatta la coercion a vector sul ring `R` se si vuole che
la funzione Sage ``vector`` lavori su tale oggetto. Il seguito è preso
dal file ``SAGE_ROOT/sage/sage/modules/free_module_element.pyx``::

    cdef class FreeModuleElement(element_Vector):   # abstract base class
        ...
        def _vector_(self, R):
            return self.change_ring(R)


.. _section-preparsing:

Preparse di Sage
================

Per rendere Python ancora più usabile interattivamente, ci sono una serie di 
aggiustamenti alla sintassi fatti quando usi Sage dalla riga di comando o 
dal notebook (ma non al codice Python nella libreria Sage). Tecnicamente ciò 
è implementato da una funzione ``preparse()`` che riscrive la stringa di input. 
Le più importanti sostituzioni fatte sono le seguenti:

- Sage supporta una sisntassi speciale per generare ring o, più in generale, 
  genitori con generatori nominati (named generators)::
     
      sage: R.<x,y> = QQ[]
      sage: preparse('R.<x,y> = QQ[]')
      "R = QQ['x, y']; (x, y,) = R._first_ngens(2)"

- Gli Integer ed i letterali real sono interi di Sage e floating-point di Sage. 
  Ad esempio in Python puro questo sarebbe un errore di attributo::

      sage: 16.sqrt()
      4
      sage: 87.factor()
      3 * 29
      
- Dei letterali grezzi non viene fatto il preparse, cosa che può
  essere utile da un punto di vista di efficienza. Come in Python gli
  ``int`` sono denotati da una L, così in Sage gli interi grezzi ed i
  letterali floating sono seguiti da una "r" (o una "R") per grezzo
  (raw), che significa "di cui non è stato fatto il preparse". Ad
  esempio::

      sage: a = 393939r
      sage: a
      393939
      sage: type(a)
      <... 'int'>
      sage: b = 393939
      sage: type(b)
      <type 'sage.rings.integer.Integer'>
      sage: a == b
      True

- I letterali grezzi possono essere molto utili in certi casi. Ad esmpio, 
  gli interi Python possono essere più efficienti degli interi Sage 
  quando sono molto piccoli. Gli interi grandi di Sage sono però molto 
  più efficienti di quelli di Python poichè sono implementati usando 
  la libreria GMP in linguaggio C.

Consulta il file ``preparser.py`` per maggiore dettagli sul preparser di Sage, 
maggiori esempi sui letterali grezzi, ecc.

Quando si carica un file ``foo.sage`` o lo si attacca ad una sessione
Sage, una versione passata al preparser di ``foo.sage`` viene creata,
con il nome ``foo.sage.py``. L'inizio del file passato al preparser
file indica::

    This file was *autogenerated* from the file foo.sage.

Puoi passare esplicitamente al preparser un file con l'opzione a riga
di comando ``--preparse``, eseguendo ::

    sage --preparse foo.sage

crea il file ``foo.sage.py``.

I seguenti file sono importatnti per il preparse in Sage:

#. ``SAGE_ROOT/src/bin/sage``

#. ``SAGE_ROOT/src/bin/sage-preparse``

#. ``SAGE_ROOT/src/sage/repl/preparse.py``

In particolare, il file ``preparse.py`` contiene il codice del
preparser di Sage.


I modello di Coercion di Sage
=============================

Lo scopo principale della ``coercion`` è riuscire a fare in modo
transparente aritmetica, confronti, ecc. fra elementi di insiemi
distinti. Ad esempio, quando si scrive `3 + 1/2`, uno vuole fare
un'operazione aritmetica su degli operandi come numeri razionali,
nonostante il termine a sinistra sia un intero. Questo ha senso
essendo gli Interi inclusi nei Razionali. Lo scopo della ``coercion``
è facilitare questo (come anche cose più complicate) senza dover
mappare in modo esplicito ogni cosa nello stesso dominio, ed allo
stesso tempo essere stringenti a sufficienza da non resolvere delle
ambiguità o accettare delle assurdità.

Il modello di coercion di Sage è descritto in dettaglio, con esempi, nella 
sezione Coercion del manuale di riferimento di Sage.


Mutabilità
============

Le strutture genitrici (ad esempio anelli, campi, spazi di matrici,
ecc.) devono essere immutabile ed uniche globalmente ogni volta che è
possibile. L'immutabilità significa, fra le altre cose, che le
proprietà come le etichette del generatore e la precisione di default
della coercion non possono essere cambiate.

L'unicità globale senza sprecare memoria è implementata nel modo
migliore usando il modulo weakref dello standard Python, una funzione
factory, e variabili di ambito nel modulo.

.. {Rewrite. Difficult to parse. Make gentler}

.. {Put a tutorial on this here}

Certi oggetti, ad esempio matrici, possono dapprima essere mutabili e
poi diventare immutabili successivamente. Vedi il file
``SAGE_ROOT/src/sage/structure/mutability.py``.


Il metodo speciale __hash__
===========================

Ecco la definizione di ``__hash__`` dal manuale di riferimento di Python:

    è richiamata dalla funzione interna ``hash()`` e per operazioni su
    membri di collezioni con hash inclusi gli insiemi (set), i
    frozenset, ed i dizionari (dict).  ``__hash__()`` deve restituire
    un integer. L'unica proprietà richiesta è che oggetti che
    risultano uguali al confronto abbiano lo stesso valore di hash; si
    consiglia di mescolare in qualche modo insieme (ad esempio usando
    ``exor``) i valori di hash dei componenti dell'oggetto che
    interviene nel confronto di oggetti.  Se una classe non definisce
    un metodo ``__cmp__()`` allora non dovrebbe definire nemmeno
    un'operazione ``__hash__()``; se definisce ``__cmp__()`` oppure
    ``__eq__()`` ma non ``__hash__()``, le sue instanze non saranno
    utilizzabili come chiavi di un dizionario. Se una classe definisce
    degli oggetti mutabili ed implementa un metodo ``__cmp__()`` o
    ``__eq__()``, allora non dovrebbe implementare ``__hash__()``,
    poichè l'implementazione del dizionario richiede che il valore
    hash di una chiave sia immutabile (se il valore hash dell'oggetto
    cambia, sarà nell'"hash bucket" sbagliato).

Nota la frase, "L'unica proprietà richiesta è che oggetti che
risultano uguali al confronto abbiano lo stesso valore di hash."
Questa è un'assunzione fatta dal linguaggio Python, che in Sage
semplicemente non possiamo fare (!), e violarla ha delle
conseguenze. Fortunatamente le conseguenze sono piuttosto chiare e
abbastanza facili da capire, quindi se le conosci non ti
intralciano. Il seguente esempio lo mostra bene:

::

        sage: v = [Mod(2,7)]
        sage: 9 in v
        True
        sage: v = set([Mod(2,7)])
        sage: 9 in v
        False
        sage: 2 in v
        True
        sage: w = {Mod(2,7):'a'}
        sage: w[2]
        'a'
        sage: w[9]
        Traceback (most recent call last):
        ...
        KeyError: 9

Ecco un altro esempio:

::

        sage: R = RealField(10000)
        sage: a = R(1) + R(10)^-100
        sage: a == RDF(1)  # because the a gets coerced down to RDF
        True

ma ``hash(a)`` non dovrebbe essere uguale a ``hash(1)``.

Sfortunatamente, in Sage semplicemente noo possiamo richiedere

::

           (#)   "a == b ==> hash(a) == hash(b)"

perchè la matematica non banale è troppo complicata per questa 
regola. Ad esempio le uguaglianze ``z == Mod(z, 2)`` e
``z == Mod(z, 3)`` forzerebbero ``hash()`` ad essere costante sugli 
Interi.

Il solo modo in cui potremmo risolvere bene questo problema sarebbe 
abbandonare l'operatore ``==`` ed usare una "uguaglianza Sage", da 
implementarsi come un nuovo metodo attaccato a ciascun oggetto. 
Poi potremmo seguire le regole di Python per ``==`` e le nostre per 
tutto il resto, e tutto il codice Sage diverrebbe completamente illeggibile 
(e pure inscrivibile). Per cui dobbiamo conviverci.

Quindi ciò che si fa in Sage è tentare di soddisfare ``(#)`` quando è 
ragionevolmente facile farlo, ma usando discernimento e senza andare fuori del 
seminato. Ad esempio,

::

        sage: hash(Mod(2,7))
        2

L'output 2 è migliore di qualche hash casuale che implichi anche i 
moduli, ma naturalmente non è corretto dal punto di vista di Python, 
poichè ``9 == Mod(2,7)``. Lo scopo è fare una funzione di hash che sia 
veloce, ma secondo ragione rispetti ogni inclusione naturale ovvia e la 
coercion.


Eccezioni
=========

Per favore evita il codice pigliatutto come il seguente::

    try:
        some_code()
    except:               # bad
        more_code()

Se non hai delle eccezioni elencate esplicitamente (come una tupla), il 
tuo codice raccoglierà veramente qualunque cosa, inclusi ``ctrl-C``, errori 
tipografici nel codice, e gli allarmi, e questo porterà a confusione. Inoltre 
questo potrebbe prendere anche errori reali che andrebbero lasciati propagare 
fino all'utente.

Riassumendo, raccogliere solo eccezioni specifiche come nell'esempio seguente::

    try:
        return self.__coordinate_ring
    except (AttributeError, OtherExceptions) as msg:           # good
        more_code_to_compute_something()

Nota che la sintassi in ``except`` è tale per elencare tutte le eccezioni raccolte 
come tupla, seguita da un messaggio di errore.


Importazione
============

Facciamo menzione di 2 problemi con le importazioni: importazioni circolari ed 
importazioni di grandi moduli di terze parti.

innanzitutto devi evitare importazioni circolari. Ad esempio supponiamo che il 
file ``SAGE_ROOT/src/sage/algebras/steenrod_algebra.py`` inizi con la linea::

    from sage.sage.algebras.steenrod_algebra_bases import *

e che il file
``SAGE_ROOT/src/sage/algebras/steenrod_algebra_bases.py``
inizi con la linea::

    from sage.sage.algebras.steenrod_algebra import SteenrodAlgebra

Questo imposta un ciclo: caricare uno di questi file richiede l'altro,
che richiede il primo, ecc.

Con queste impostazioni, eseguire Sage darà l'errore::

    Exception exceptions.ImportError: 'cannot import name SteenrodAlgebra'
    in 'sage.rings.polynomial.polynomial_element.
    Polynomial_generic_dense.__normalize' ignored
    -------------------------------------------------------------------
    ImportError                       Traceback (most recent call last)

    ...
    ImportError: cannot import name SteenrodAlgebra

Invece puoi rimpiazzare la linea ``import *`` in cima al file con
delle istruzioni import piì specifiche laddove sono richieste nel
codice. Ad esempio, il metodo ``basis`` per la classe
``SteenrodAlgebra`` potrebbe assomigliare a questo (omettendo la
stringa di documentazione)::

    def basis(self, n):
        from steenrod_algebra_bases import steenrod_algebra_basis
        return steenrod_algebra_basis(n, basis=self._basis_name, p=self.prime)

Secondo, non importare a livello radice del tuo modulo un modulo di
terze parti che impiegherà parecchio tempo ad inizializzarsi (ad
esempio ``matplotlib``).  come sopra, puoi invece importare dei
componenti specifici del modulo quando ne hai bisogno, piuttosto che
al livello radice del tuo file.

Ès importante provare a rendere ``from sage.all import *`` il più
veloce possibile, poichè è la parte dominante del tempo di avvio di
Sage, e controllare le import di massimo livello aiuta a limitarlo. Un
meccanismo importante in Sage sono le ``lazy import`` (importazioni
pigre), che non vengono realmente effettuate subito ma aspettano
finchè si utilizza effettivamente l'oggetto). Vedi
:mod:`sage.misc.lazy_import` for more details of lazy imports, and
:ref:`chapter-directory-structure` per un esempio di uso delle ``lazy
import`` in un nuovo modulo.


Deprecazione
============

Quando si fa una modifica in Sage che è **incompatibile
all'indietro**, il vecchio codice dovrebbe continuare a funzionare e
mostrare un messaggio che indica come dovrebbe essere aggiornato o
scritto in futuro. Questo si chiama una *deprecazione*.

.. NOTE::

    Il codice deprecato può solo essere rimosso un anno dopo la prima
    release stabile in cui è apparso.

Each deprecation warning contains the number of the trac ticket that defines
it. We use 666 in the examples below. For each entry, consult the function's
documentation for more information on its behaviour and optional arguments.

* **Rinominare una parola chiave:** facendo il ``decorating`` di un
metodo o funzione con :class:`~sage.misc.decorators.rename_keyword`,
qualunque utente che chiami ``my_function(my_old_keyword=5)`` vedrà un
warning::

      from sage.misc.decorators import rename_keyword
      @rename_keyword(deprecation=666, my_old_keyword='my_new_keyword')
      def my_function(my_new_keyword=True):
          return my_new_keyword

* **Rinominare una funzione/metodo:** chiama
:func:`~sage.misc.superseded.deprecated_function_alias` per ottenere
una copia della funzione che solleva il ``deprecation warning``::

      from sage.misc.superseded import deprecated_function_alias
      def my_new_function():
          ...

      my_old_function = deprecated_function_alias(my_new_function)

* **Spostare un oggetto in un altro moduloto:** se rinomini un file
  sorgente o sposti qualche funzione (o classe) in un file differente,
  dovrebbe ancora essere possibile importare tale funzione dal vecchio
  modulo. Queso può essere fatto usando
  :func:`~sage.misc.lazy_import.lazy_import` con la deprecazione.  Nel
  vecchio modulo scrivi::

    from sage.misc.lazy_import import lazy_import
    lazy_import('sage.new.module.name', 'name_of_the_function', deprecation=666)

  Puoi anche fare la lazy import di tutto usando ``*`` o solo alcune
  funzioni usando una tupla::

    from sage.misc.lazy_import import lazy_import
    lazy_import('sage.new.module.name', '*', deprecation=666)
    lazy_import('sage.other.module', ('func1', 'func2'), deprecation=666)

* **Rimuovere un nome dal namespace globale:** la funzione
:func:`~sage.misc.superseded.deprecated_callable_import` importa un
oggetto nel namespace globale. Qualunque utente la chiami vedrà un
messaggio che lo invita ad importare l'oggetto manualmente::

      from sage.misc.superseded import deprecated_callable_import
      deprecated_callable_import(666,
                           'sage.combinat.name_of_the_module',
                           globals(),
                           locals(),
                           ["name_of_the_function"])

  In alternativa, anche una lazy_import con deprecazione funzionerà in
  questo caso.

* **in ogni altro caso:** se no si applica nessuno dei casi suddetti,
chiama :func:`~sage.misc.superseded.deprecation` nella funzione che
vuoi deprecare.  Mostrerà il messaggio di tua scelta (ed interagirà
corretamente con il framework di doctest)::

      from sage.misc.superseded import deprecation
      deprecation(666, "Do not use your computer to compute 1+1. Use your brain.")


Codice sperimentale/instabile
-----------------------------

Puoi marcare il codice che hai appena creato (classi/funzioni/metodi)
come sperimentale/instabile. In questo caso, non è necessario alcun
warning di deprecation quando si cambia il codice, la sua
functionalità o la sua interfaccia.

Questo dovrebbe permetterti di inserire il tuo materiale in Sage il
prima possibile, senza preoccuparti di dover fare modifiche (di
progettazione) successivamnete.

Quando sei soddisfatto del codicee (quando è stabile per un po' di
tempo, diciamo 1 anno), puoi cancelare questo warning.

Come al solito, tutto il codice deve essere completo di doctest e deve
pasare il nostro processo di revisione.

* **Funzione/metodo sperimentale:** usa il decoratore
  :class:`~sage.misc.superseded.mark_as_experimental`. Ecco un esempio::

      from sage.misc.superseded import experimental
      @experimental(66666)
      def experimental_function():
          # do something

* **Classe sperimentale:** usa il decoratore
  :class:`~sage.misc.superseded.mark_as_experimental` per il suo ``__init__``.
  Ecco un esempio::

      from sage.misc.superseded import experimental
      class experimental_class(SageObject):
          @experimental(66666)
          def __init__(self, some, arguments):
              # do something

* **In tutti gli altri casi:** se non si applica nessuno dei casi
suddetti, chiama :func:`~sage.misc.superseded.experimental` nel codice
dove vuoi sollevare il warning. Mostrerà il messaggio di tua scelta::

      from sage.misc.superseded import experimental_warning
      experimental_warning(66666, 'This code is not foolproof.')


Usare pacchetti opzionali
=========================

Se una funzione richiede un pacchetto opzionale, quella funzione
dovrebbe andare in errore in modo controllato---magari usando un
blocco ``try``-``except``--- quando il pacchetto opzionale non è
disponibile, e dovrebbe dare un suggerimento su come installarlo. Ad
esempio, digitare ``sage -optional`` da una lista di tutti i pacchetti
opzionali, così che puòi suggerire all'utente di digitare quello.  Il
comando ``optional_packages()`` da dentro Sage anche da tale lista.
