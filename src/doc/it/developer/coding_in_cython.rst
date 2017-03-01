.. _chapter-cython:

=========================
Scrivere codice in Cython
=========================

Questo capitolo discute Cython, che è un linguaggio compilato basato su 
Python. Il più grande vantaggio che ha su Python è che il codice può 
essere molto più veloce (a volte di ordini di grandezza) e può chiamare 
direttamente codice C e C++. Poichè Cython è essenzialmente un ampliamento 
del linguaggio Python, spesso non si fa distinzione fra codice Cython e  
Python in Sage (ad esempio si parla della “libreria Python di Sage” 
e delle “convenzioni di codifica in Python”).

Python è un linguaggio interpretato e non dichiara tipi di dato per le
variabili. Queste caratteristiche rendono facile scrivere programmi e
farne il debug, ma il codice Python a volte può essere lento. Il
codice Cython può assomigliare molto a Python, ma viene tradotto in
codice C (spesso molto efficiente) e poi compilato. Quindo offre un
linguaggio familiare agli sviluppatori Python, ma con un potenziale di
una molto maggiore velocità. Inoltre Cython permette agli sviluppatori
Sage di interfacciarsi con C e C++ molto più facilmente che usando la
C API di Python direttamente.

Cython è una versione compilata di Python. In origine era basata su
Pyrex ma è stato modificato in base a ciò di cui avevano bisogno gli
sviluppatori di Sage; Cython è stato sviluppato di concerto con
Sage. Comunque ora è un progetto independente, che è utilizzato al di
là dell'ambito di Sage. Come tale, è un linguaggio giovane, ma in
crescita, con documentazione agli inizi, ma in crescita.  Vedi la sua
pagina web, http://www.cython.org/, per le informazioni più aggiornate
o vai a `Basi del linguaggio
<http://docs.cython.org/src/userguide/language_basics.html>`_ per
iniziare subito.


Scrivere codice Cython in Sage
==============================

Ci sono molti modi di creare e compilare codice Cython in Sage.

#. Nel Notebook Sage, inizia ogni cella con ``%cython``. Quando valuti
   tale cella,

   #. Viene savata in un file.

   #. Cython è eseguito su questo, linkato a tutte le librerie
      standard di Sage se necessario.

   #. Il file di libreria condivisa resultante (``.so`` / ``.dll`` /
      ``.dylib``) è poi caricato nell'istanza in esecuzione di Sage.

   #. La funzionalità definita in tale cella può ora essere utilizzata
      nel notebook. Inoltre la cella di output ha un link al programma
      C che è stato compilato per creare il file ``.so``.

   #. Una funzione ``cpdef`` o ``def``, diciamo ``testfunction``,
      definita in una cella ``%cython`` in un worksheet può essere
      importata e resa disponibile in una cella ``%cython`` differente
      nello stesso worksheet importandola come mostrato qui sotto::

          %cython
          from __main__ import testfunction

#. Crea un file ``.spyx`` ed attaccalo o caricalo dalla riga di
   comando. Questo è simile a creare una cella ``%cython`` nel
   notebook ma funziona in modo completo dalla riga di comando (e non
   dal notebook).

#. Crea un file ``.pyx`` e aggiungilo alla libreria Sage.

   #. Prima aggiungi un elenco per le estensioni Cython alla variabile
      ``ext_modules`` nel file ``SAGE_ROOT/src/module_list.py``. Vedi
      la classe ``distutils.extension.Extension`` per maggiori
      informazioni sul creare una nuova estensione Cython.

   #. Esegui ``sage -b`` per ricompilare Sage.

   Ad esempio per compilare
   ``SAGE_ROOT/src/sage/graphs/chrompoly.pyx`` vediamo le righe
   seguenti in ``module_list.py``::

    Extension('sage.graphs.chrompoly',
              sources = ['sage/graphs/chrompoly.pyx'],
              libraries = ['gmp']),


Pragma speciali
===============

Se il codice Cython è attaccato o caricato come un file ``.spyx`` o 
caricato dal notebook come un blocco ``%cython``, le seguenti pragma 
sono disponibili:

* clang --- può essere sia c che c++ indicando se dev'essere usato il 
  compilatore C o C++.

* clib --- librerie addizionali a cui fare link, la lista separata da 
  spazi è spezzata e passata a ``distutils``.

* cinclude --- directory addizionali da ricercare per file di header. La 
  lista separata da spazi è spezzata e passata a ``distutils``.

* cfile -- file C or C++ addizionali da compilare

* cargs -- parametri addizionali passati al compilatore

Ad esempio::

    #clang C++
    #clib givaro
    #cinclude /usr/local/include/
    #cargs -ggdb
    #cfile foo.c


Attaccare o caricare file .spyx
===============================

Il modo più facile di provare Cython senza dover imparare nulla su distutils, 
ecc., è creare un file con l'estensione ``spyx``, che sta per "Sage Pyrex":

#. Crea un file ``power2.spyx``.

#. Mettici dentro quanto segue::

       def is2pow(n):
           while n != 0 and n%2 == 0:
               n = n >> 1
           return n == 1

#. Lancia la riga comandi di Sage e carica il file ``spyx`` (questo darà errore 
   se non hai già un compilatore C installato).

   .. skip

   ::

       sage: load("power2.spyx")
       Compiling power2.spyx...
       sage: is2pow(12)
       False

Nota che se cambi ``power2.spyx`` e poi lo ricarichi, sarà ricompilato
al volo.  Puoi anche attaccare ``power2.spyx`` così che venga
ricaricatoogni volta che fai qualche cambiamento:

.. skip

::

    sage: attach("power2.spyx")

Cython è usato per la sua velocità. Ecco un test cronometrato su un Opteron da 
2.6 GHz:

.. skip

::

    sage: %time [n for n in range(10^5) if is2pow(n)]
    [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
    CPU times: user 0.60 s, sys: 0.00 s, total: 0.60 s
    Wall time: 0.60 s

Ora, il codice nel file ``power2.spyx`` è Python valido, e se lo copiamo in 
in file ``powerslow.py`` e lo carichiamo, otteniamo quanto segue:

.. skip

::

    sage: load("powerslow.py")
    sage: %time [n for n in range(10^5) if is2pow(n)]
    [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
    CPU times: user 1.01 s, sys: 0.04 s, total: 1.05 s
    Wall time: 1.05 s

Tra l'altro, possiamo guadagnare ancora un po\` di velocità nella
versione Cython una dichiarazione di tipo, cambiando ``def
is2pow(n):`` in ``def is2pow(unsigned int n):``.


.. _section-interrupt:

Interrupt e gestione dei segnali
================================

Quando si scrive codice Cython per Sage, bisogna avere un'attenzione
speciale ad assicurarsi che il codice possa essere interrotto con
``CTRL-C``.  Poichè Cython è ottimizzato per la velocità, Cython di
solito non controlla gli interrupt.  Ad esempio codice come il
seguente non può essere interrotto:

.. skip

::

    sage: cython('while True: pass')  # DON'T DO THIS

Mentre questo è in esecuzione, premere ``CTRL-C`` non ha effetti. Il
solo modo di uscirne è terminare il processo di Sage.  Su certi
sistemi puoi ancora terminare Sage con ``CTRL-\`` (manda un segnale
Quit) invece di ``CTRL-C``.

.. Use Cython syntax highlighting for the rest of this document.

.. highlight:: cython

Sage fornisce 2 meccanismi collegati per gestire gli interrupts:

* :ref:`Use sig_check() <section_sig_check>` se stai scrivendo codice 
  Cython/Python misto. Tipicamente questo è codice con cicli (annidati) 
  dove ogni singola istruzione impiega poco tempo.

* :ref:`Use sig_on() and sig_off() <section_sig_on>` se stai invocando delle 
  librerie C esterne o dentro codice Cython puro (senza alcuna funzione Python) 
  dove anche una singola istruzione, come una chiamata a libreria, può 
  richiedere molto tempo.

Le funzioni ``sig_check()``, ``sig_on()`` e ``sig_off()`` possono essere messe 
in qualunque tipo di funzione Cython: ``def``, ``cdef`` o ``cpdef``.
Non puoi metterle in codice Python puro (i file con estensione ``.py``).
Queste funzioni sono specifiche di Sage. Per usarle, **devi** includere quanto 
segue nel tuo file ``.pyx`` (non è sufficiente farlo in un file ``.pxd``)::

    include "sage/ext/interrupt.pxi"

.. NOTE::

    Le funzioni Cython ``cdef`` o ``cpdef`` con un tipo di ritorno
    (come ``cdef int myfunc():``) devono avere un `except value
    <http://docs.cython.org/src/userguide/language_basics.html#error-return-values>`_
    per propagare le eccezioni.  Ricordati di questo ogni volta che
    scrivi ``sig_check()`` o ``sig_on()`` dentro ad una funzione di
    questo tipo, altrimenti vedrai un messaggio ``Exception
    KeyboardInterrupt: KeyboardInterrupt() in <function name>
    ignored``.

.. _section_sig_check:

Usare ``sig_check()``
---------------------

``sig_check()`` può essere usato per valutare se ci sono degli
interrupt in corso.  Se un interrupt accade durante l'esecuzione di
codice C o Cython, verrà catturato dalla successiva ``sig_check()`` o
``sig_on()`` o anche dalla successiva istruzione Python. Con
quest'ultima intendiamo che anche certe istruzioni Python valutano gli
interrupt, ad esempio l'istruzione ``print``.  Il seguente ciclo *può*
essere interrotto:

.. code-block:: python

    sage: cython('while True: print("Hello")')

Il tipico caso d'uso per ``sig_check()`` è dentro piccoli cicli che
fanno cose complicate (codice Python e Cython mescolati, che possono
sollevare eccezioni).  È ragionevolmente sicuro da usare e da grande
controllo, perchè nel tuo codice Cython un ``KeyboardInterrupt`` può
*solo* essere sollevato durante ``sig_check()``::

    def sig_check_example():
        for x in foo:
            # (one loop iteration which does not take a long time)
            sig_check()

Questo ``KeyboardInterrupt`` è trattato come ogni altra eccezione Python e può 
essere gestita come al solito::

    def catch_interrupts():
        try:
            while some_condition():
                sig_check()
                do_something()
        except KeyboardInterrupt:
            # (handle interrupt)

Naturalmente puoi anche mettere la ``try``/``except`` nel ciclo, nell'esempio sopra.

La funzione ``sig_check()`` è una funzione inline molto veloce che non
dovrebbe avere effetti misurabile sulle performance.

.. _section_sig_on:

Usare ``sig_on()`` and ``sig_off()``
------------------------------------

Un altro meccanismo per la gestione degli interrupt è la coppia di
funzioni ``sig_on()`` e ``sig_off()``.  È più potente di
``sig_check()`` ma anche molto più pericoloso.  Dovresti mettere
``sig_on()`` *prima* e ``sig_off()`` *dopo* qualunque codice Cython
che può impiegare molto tempo.  Questi 2 *devono sempre* essere
richiamati in **coppia**, cioè ogni ``sig_on()`` deve corrispondere ad
una ``sig_off()`` di chiusura.

In pratica la tua funzione probabilmente sarà simile a::

    def sig_example():
        # (some harmless initialization)
        sig_on()
        # (a long computation here, potentially calling a C library)
        sig_off()
        # (some harmless post-processing)
        return something

È possibile mettere ``sig_on()`` e ``sig_off()`` in funzioni differenti, 
purchè ``sig_off()`` sia chiamata prima che la funzione che chiama la
``sig_on()`` termini l'esecuzione.
Il codice seguente *non è valido*::

    # INVALID code because we return from function foo()
    # without calling sig_off() first.
    cdef foo():
        sig_on()

    def f1():
        foo()
        sig_off()

Ma il seguente è valido poichè non si può chiamare ``foo`` interattivamente::

    cdef int foo():
        sig_off()
        return 2+2

    def f1():
        sig_on()
        return foo()

Per chiarezza, comunque, è meglio evitare tutto ciò.
Un buon esempio dove quanto sopra ha senso è la funzione 
``new_gen()`` in :ref:`section-pari-library`.

Un errore comune è mettere ``sig_off()`` verso la fine della 
funzione (prima della ``return``) quando la funzione ha più di 
una istruzione ``return``.
Pertanto accertati che ci sia una ``sig_off()`` davanti ad *ogni* ``return``
(ed anche davanti ad ogni ``raise``).

.. WARNING::

    Il codice in ``sig_on()`` dev'essere C puro o codice Cython. 
    Se chiami del codice Python o manipoli un oggetto Python 
    (anche qualcosa di semplice come ``x = []``),
    un interrupt può pasticciare lo stato interno di Python.
    Nel dubbio prova ad usare :ref:`sig_check() <section_sig_check>` invece.

    Anche, quando un interrupt capita dentro ``sig_on()``, l'esecuzione del 
    codice viene fermata immediatamente senza fare pulizie.
    Ad esempio qualunque memoria allocata dentro ``sig_on()`` viene perduta.
    Vedi :ref:`advanced-sig` per dei modi di gestire questo.

Quando l'utente preme ``CTRL-C`` dentro ``sig_on()``, l'esecuzione salterà 
indietro a ``sig_on()`` (la prima che c'è nello stack) e ``sig_on()`` 
solleverà ``KeyboardInterrupt``. Come con ``sig_check()``, questa 
eccezione può essere gestita nel solito modo::

    def catch_interrupts():
        try:
            sig_on()  # This must be INSIDE the try
            # (some long computation)
            sig_off()
        except KeyboardInterrupt:
            # (handle interrupt)

Certe librerie C in Sage sono scritte in modo da sollevare eccezioni Python:
libGAP ed NTL possono sollevare ``RuntimeError`` e PARI ``PariError``.
Queste eccezioni si comportano esattamente come ``KeyboardInterrupt`` 
nell'esempio sopra e possono essere raccolte mettendo ``sig_on()`` dentro 
un blocco ``try``/``except``.
Vedi :ref:`sig-error` per come ciò è implementato.

È possibile accumulare ``sig_on()`` e ``sig_off()``.
Se lo fai, l'effetto è esattamente lo stesso che se ci fosse solo lo 
``sig_on()``/``sig_off()`` più esterno. L'interno cambierà semplicemente 
un contatore di referenze e nient'altro. Assicurati che il numero di chiamate 
``sig_on()`` eguagli il numero di chiamate ``sig_off()``::

    def f1():
        sig_on()
        x = f2()
        sig_off()

    def f2():
        sig_on()
        # ...
        sig_off()
        return ans

Attenzione aggiuntiva va fatta con eccezioni sollevate dentro ``sig_on()``.
Il problema è che, se non fai niente di speciale, la ``sig_off()`` non 
sarà mai invocata se c'è un'eccezione.
Se devi tu stesso *sollevare* un'eccezione, chiama una ``sig_off()`` prima::

    def raising_an_exception():
        sig_on()
        # (some long computation)
        if (something_failed):
            sig_off()
            raise RuntimeError("something failed")
        # (some more computation)
        sig_off()
        return something

In alternativa puoi usare ``try``/``finally`` che catturerà ugualmente 
eccezioni sollevate da subroutine dentro la ``try``::

    def try_finally_example():
        sig_on()  # This must be OUTSIDE the try
        try:
            # (some long computation, potentially raising exceptions)
            return something
        finally:
            sig_off()

Se vuoi catturare anche quest'eccezione, hai bisogno di una ``try`` annidata::

    def try_finally_and_catch_example():
        try:
            sig_on()
            try:
                # (some long computation, potentially raising exceptions)
            finally:
                sig_off()
        except Exception:
            print("Trouble!Trouble!")

``sig_on()`` è implementata usando la chiamata di libreria C ``setjmp()`` 
che richiede una piccola ma non trascurabile quantità di tempo.
In codice veramente time-critical, si possono richiamare ``sig_on()``
e ``sig_off()`` in modo condizionale::

    def conditional_sig_on_example(long n):
        if n > 100:
            sig_on()
        # (do something depending on n)
        if n > 100:
            sig_off()

Ciò dovrebbe essere necessario solo se sia la verifica 
(``n > 100`` nell'esempio) che il codice dentro il blocco ``sig_on()`` 
richiedono molto poco tempo.
Nelle versioni di Sage anteriori alla 4.7, ``sig_on()`` era molto più 
lento, ecco perchè ci sono più verifiche come questa nel vecchio codice.

Altri segnali
-------------

A parte la gestione degli interrupt, la ``sig_on()`` fornisce una
gestione più generale dei segnali.  Ad esempio gestisce :func:`alarm`
time-out sollevando un'eccezione ``AlarmInterrupt`` (ereditata da
``KeyboardInterrupt``).

Se il codice dentro ``sig_on()`` genera un ``segmentation fault`` o
chiama la funzione C ``abort()`` (o più in generale solleva una
qualunque fra SIGSEGV, SIGILL, SIGABRT, SIGFPE, SIGBUS), questa è
catturata dal framework di interrupt ed un'eccezione è sollevata
(``RuntimeError`` per SIGABRT, ``FloatingPointError`` per SIGFPE e
l'eccezione personalizzata ``SignalError``, basata su
``BaseException``, altrimenti)::

    cdef extern from 'stdlib.h':
        void abort()

    def abort_example():
        sig_on()
        abort()
        sig_off()

.. code-block:: python

    sage: abort_example()
    Traceback (most recent call last):
    ...
    RuntimeError: Aborted

Questa eccezione può essere gestita da un blocco ``try``/``except``
come spiegato sopra. Un ``segmentation fault`` o ``abort()`` non
controllati da ``sig_on()`` possono semplicemente terminare
Sage. Questo si applica solo a ``sig_on()``, la funzione
``sig_check()`` si occupa solo di interrupt ed allarmi.

Invece di ``sig_on()``, c'è anche una funzione ``sig_str(s)``, che
prende una stringa C ``s`` come argomento. Si comporta nello stesso
mdod di ``sig_on()``, eccetto che la stringa ``s`` sarà utilizzata
come stringa per l'eccezione.  ``sig_str(s)`` deve ancora essere
chiusa da ``sig_off()``.  Esempio di codice Cython::

    cdef extern from 'stdlib.h':
        void abort()

    def abort_example_with_sig_str():
        sig_str("custom error message")
        abort()
        sig_off()

Eseguire ciò produce:

.. code-block:: python

    sage: abort_example_with_sig_str()
    Traceback (most recent call last):
    ...
    RuntimeError: custom error message

Riguardo agli interrupt ordinari (cioè SIGINT), ``sig_str(s)`` si
comporta nello stesso modo di ``sig_on()``: è sollevato un semplice
``KeyboardInterrupt``.

.. _sig-error:

Gestione degli errori nelle librerie C
--------------------------------------

Alcune librerie C possono produrre errori ed usare qualche sorta di meccanismo 
di callback per segnalare errori: una funzione esterna di gestione degli errori 
va messa sù, che sarà chiamata dalla libreria C se capita un errore.

La funzione ``sig_error()`` può essere usata per gestire questi
errori.  Questa funzione può solo esserechiamata dentro un blocco
``sig_on()`` (altrimenti Sage andrà in crash malamente) dopo aver
sollevato un'eccezione Python. Devi usare la `Python/C API
<http://docs.python.org/2/c-api/exceptions.html>`_ per questo, e
chiamare ``sig_error()`` dopo aver chiamato qualche variante di
:func:`PyErr_SetObject`. Anche dentro Cython non puoi usare
l'istruzione ``raise``, perchè così la ``sig_error()`` non sarebbe mai
eseguita.  La chiamata a ``sig_error()`` userà i meccanismi di
``sig_on()`` così che l'eccezione sarà vista da ``sig_on()``.

Un tipico gestore di errori implementato in Cython sarebbe come segue::

    include "sage/ext/interrupt.pxi"
    from cpython.exc cimport PyErr_SetString

    cdef void error_handler(char *msg):
        PyErr_SetString(RuntimeError, msg)
        sig_error()

In Sage questo meccanismo è  utilizzato per libGAP, NTL e PARI.

.. _advanced-sig:

Funzioni avanzate
-----------------

Ci sono molte funzioni specializzate per gestire gli interrupt.  Come
detto sopra, ``sig_on()`` non cerca di ripulire nulla (restore dello
stato o liberare la memoria) quando capita un interrupt.  Infatti
sarebbe impossibile per ``sig_on()`` farlo.  Se vuoi aggiungere del
codice di pulizia (cleanup), usa ``sig_on_no_except()`` per
questo. Questa funzione si comporta *esattamente* come ``sig_on()``,
eccetto che qualunque eccezione sollevata (come ``KeyboardInterrupt``
o ``RuntimeError``) non è ancora passata a Python. Essenzialmente
l'eccezione è lì, ma possiamo impedire a Cython di vederla. Poi si può
usare ``cython_check_exception()`` per permettere a Cython di cercare
l'eccezione.

Normalmente ``sig_on_no_except()`` restituisce 1.  Se un segnale è
catturato ed un'eccezione sollevata, ``sig_on_no_except()``
restituisce, invece, 0.  Il seguente esempio mostra come usare
``sig_on_no_except()``::

    def no_except_example():
        if not sig_on_no_except():
            # (clean up messed up internal state)

            # Make Cython realize that there is an exception.
            # It will look like the exception was actually raised
            # by cython_check_exception().
            cython_check_exception()
        # (some long computation, messing up internal state of objects)
        sig_off()

C'è anche una funzione ``sig_str_no_except(s)`` che è analoga a ``sig_str(s)``.

.. NOTE::

    Vedi il file :file:`SAGE_ROOT/src/sage/tests/interrupt.pyx`
    per maggiori esempi di come usare le varie funzioni ``sig_*()``.

Fare il test degli interrupt
----------------------------

.. highlight:: python

Quando si scrive :ref:`section-docstrings`, spesso si vuole verificare
che un certo codice può essere interrotto in maniera pulita.  Il modo
migliore di farlo è usare :func:`alarm`.

Ecco un esempio du un doctest che dimostra che la funzione
:func:`factor()` può essere interrotta::

    sage: alarm(0.5); factor(10^1000 + 3)
    Traceback (most recent call last):
    ...
    AlarmInterrupt

Rilasciare il Global Interpreter Lock (GIL)
-------------------------------------------

Tutte le funzioni legate agli interrupt e la gestione dei segnali non
richiedono il `Python GIL
<http://docs.cython.org/src/userguide/external_C_code.html#acquiring-and-releasing-the-gil>`_
(se non sai cosa significa, puoi saltare senz'altro questa sezione),
essi sono dichiarati ``nogil``.  Questo significa che possono essere
usati nel codice Cython dentro a blocchi ``with nogil``. Se
``sig_on()`` deve sollevare un'eccezione, il GIL è temporaneamente
acquisito internamente.

Se usi le librerie C senza il GIL e vuoi sollevare un'eccezione prima
di chiamare :ref:`sig_error() <sig-error>`, ricorda di acquisire il
GIL mentre sollevi l'eccezione.  Dentro Cython puoi usare un `with gil
context
<http://docs.cython.org/src/userguide/external_C_code.html#acquiring-the-gil>`_.

.. WARNING::

    Il GIL non va mai rilasciato o acquisito dentro ad un blocco
    ``sig_on()``.  Se vuoi usare un blocco ``with nogil``, metti
    entrambe le ``sig_on()`` e ``sig_off()`` dentro il blocco. Nel
    dubbio, usa ``sig_check()`` al posto, che è sempre di utilizzo
    sicuro.

I pickle nel codice Cython
======================

Gestire i pickle (sottaceti) per le classi Python e per le classi
estensioni di Python, come Cython, è differente. Questo è discusso
nella `Python pickling documentation`_.  Per gestire i pickle delle
classi estensioni devi scrivere un metodo :meth:`__reduce__` che
tipicamente restituirà una tupla ``(f, args, ...)`` tale che
``f(*args)`` restituisce (una copia del) l'oggetto originale. Ad
esempio il seguente pezzetto di codice è il metodo
:meth:`~sage.rings.integer.Integer.__reduce__` da
:class:`sage.rings.integer.Integer`::

    def __reduce__(self):
        '''
        This is used when pickling integers.

        EXAMPLES::

            sage: n = 5
            sage: t = n.__reduce__(); t
            (<built-in function make_integer>, ('5',))
            sage: t[0](*t[1])
            5
            sage: loads(dumps(n)) == n
            True
        '''
        # This single line below took me HOURS to figure out.
        # It is the *trick* needed to pickle Cython extension types.
        # The trick is that you must put a pure Python function
        # as the first argument, and that function must return
        # the result of unpickling with the argument in the second
        # tuple as input. All kinds of problems happen
        # if we don't do this.
        return sage.rings.integer.make_integer, (self.str(32),)


.. _python pickling documentation: http://docs.python.org/library/pickle.html#pickle-protocol
