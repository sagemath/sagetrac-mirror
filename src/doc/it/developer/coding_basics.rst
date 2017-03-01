.. _chapter-code-basics:

====================
Convenzioni generali
====================


Ci sono molti modi di contribuire a Sage incluso il condividere script 
e worksheet Sage che implementino nuove functionalit\`a usando Sage, 
migliorando la libreria Sage, o lavorando alle molte librerie sottostanti 
distribuite con Sage [1]_.
Questa guida si focalizza sul modificare la libreria di Sage.

Sage non \`e solo mettere assieme delle functionalit\`a. \`E fornire un 
chiaro, sistematico e coerente modo di accedere un gran numero di 
algoritmi, in un framework coerente che abbia senso matematicamente. 
Nel design di Sage, la semantica degli oggetti, le definizioni, ecc., 
sono fatte in accordo con il modo in cui gli oggetti corrispondenti sono 
usati nella matematica di tutti i giorni.

.. [1]
   Vedi http://www.sagemath.org/links-components.html per una lista completa 
   dei pacchetti inclusi in ogni copia di Sage

Per raggiungere l'obiettivo di rendere Sage facile da leggere, manutenere, e
migliorare, tutto il codice Python/Cython incluso in Sage deve aderire alle 
convenzioni di stile discusse in questo capitolo.


.. _section-coding-python:

Stile del codice Python
=======================

Segui le regole standard di formattazione di Python quando scrivi codice per 
Sage, come spiegato agli URL seguenti:

* http://www.python.org/dev/peps/pep-0008
* http://www.python.org/dev/peps/pep-0257

In particolare,

- Usa 4 spazi per i livelli di indentazione. Non usare tabulazioni poich\`e 
  possono risultare in confuzione nell'indentazione. La maggior parte degli 
  editor ha una funzionalit\`a che inserisce 4 spazi quando permi il tasto TAB.
  Inoltre molti editor ricercano e sostituiscono automaticamente le tabulazioni
  con 4 spazi.

- Gli spazi bianchi prima e dopo le assegnazioni e gli operatori binari di pi\`u 
  bassa priorit\`a nell'espressione::

      i = i + 1
      c = (a+b) * (a-b)

- Nessun spazio bianco prima o dopo il segno ``=`` se \`e utilizzato per parole 
  chiave come argomenti::

      def complex(real, imag=0.0):
          return magic(r=real, i=imag)

- Nessun spazio bianco subito dentro le parentesi tonde, quadre e graffe::

       spam(ham[1], {eggs: 2})
       [i^2 for i in range(3)]

- Usa nomi di funzioni in lettere tutte minuscole con parole separate da caratteri 
  di sottolineatura (underscore). Ad esempio, ti suggeriamo di scrivere funzioni 
  Python usando la convenzione di denominazione::

      def set_some_value():
          return 1

  Nota, comunque, che qualche funzione ha effettivamente lettere maiuscole laddove
  necessario. Ad esempio, la funzione per la riduzione di un lattice con l'algoritmo 
  LLL \`e chiamata ``Matrix_integer_dense.LLL``.

- Usa la convezione ``a cammello`` (CamelCase) per i nomi di classe::

      class SomeValue(object):
          def __init__(self, x):
          self._x  = 1

  e funzioni factory che ricalcano construttori di oggetti, ad esempio 
  ``PolynomialRing`` oppure::

       def SomeIdentityValue(x):
           return SomeValue(1)



.. _chapter-directory-structure:

Struttura di file e directory
=============================

A grandi linee l'albero di directory di sage assomiglia a quanto segue. 
Nota che usiamo ``SAGE_ROOT`` nel seguito come abbreviazione per il nome 
(arbitrario) della directory contente il codice sorgente di Sage::

    SAGE_ROOT/
        sage          # the Sage launcher
        Makefile      # top level Makefile
        build/        # sage's build system
            deps
            install
            ...
            pkgs/     # install, patch, and metadata from spkgs
        src/
            setup.py
            module_list.py
            ...
            sage/     # sage library (formerly devel/sage-main/sage)
            ext/      # extra sage resources (formerly devel/ext-main)
            mac-app/  # would no longer have to awkwardly be in extcode
            bin/      # the scripts in local/bin that are tracked
        upstream/     # tarballs of upstream sources
        local/        # installed binaries

Il codice Python della libreria Sage \`e in ``src/`` e usa le convenzioni 
seguenti. I nomi di directory possono essere plurali (ad esempio ``rings``) 
ed i nomi di file sono quasi sempre singolari (ad esempio ``polynomial_ring.py``). 
Nota che il file ``polynomial_ring.py`` pu\`o ancora contenere le definizioni 
di molti tipi differenti di anelli di polinomi.

.. NOTE::

   Ti essortiamo ad includere note varie, email, discussioni sul design, 
   ecc., nel tuo pacchetto.  Metti ci\`o in file di testo semplice (con 
   estensione ``.txt``) in una sottodirectory detta ``notes``. Ad esempio 
   vedi ``SAGE_ROOT/src/sage/ext/notes/``.

Se vuoi creare una nuova directory nella libreria Sage ``SAGE_ROOT/src/sage`` 
(ad esempio sia ``measure_theory``), quella directory dovrebbe contenere 
un file ``__init__.py`` che contiene la singola linea ``import all`` in 
aggiunta a qualunque altro file che vuoi aggiungere (quali, ad esempio, 
``borel_measure.py`` e ``banach_tarski.py``), ed anche un file ``all.py`` 
che elenchi le import da quella directory che sono sufficientemente importanti 
da essere nello spazio dei nomi (namespace) globale di Sage all'avvio.
Il file ``all.py`` potrebbe essere come il seguente::

    from borel_measure import BorelMeasure
    from banach_tarski import BanachTarskiParadox

Ma in genere \`e meglio usare il lazy import framework::

    from sage.misc.lazy_import import lazy_import
    lazy_import('sage.measure_theory.borel_measue', 'BorelMeasure')
    lazy_import('sage.measure_theory.banach_tarski', 'BanachTarskiParadox')

Allora nel file ``SAGE_ROOT/src/sage/all.py``, aggiungi la linea ::

    from sage.measure_theory.all import *


Imparare facendo copia/incolla
==============================

Per tutte le convenzioni discusse qui, puoi trovare molti esempi nella 
libreria Sage. Esplorare il codice \`e di aiuto, ma anche il cercare: vale la pena 
conoscere le funzioni ``search_src``, ``search_def``, e ``search_doc``. In breve, 
dal prompt "sage:", la ``search_src(string)`` ricerca nel codice della libreria Sage 
la stringa ``string``. Il comando ``search_def(string)`` fa una ricerca simile, 
ma ristretta alle definizioni di funzione, mentre ``search_doc(string)`` ricerca tutta 
la documentazione di Sage. Vedi le loro docstring per maggiori informazioni ed opzioni.


Intestazioni dei file di codice della libreria Sage
===================================================

La testata di ciascun file di codice di Sage deve seguire questo formato::

    r"""
    <Very short 1-line summary>

    <Paragraph description>

    AUTHORS:

    - YOUR NAME (2005-01-03): initial version

    - person (date in ISO year-month-day format): short desc
    
    EXAMPLES::

    <Lots and lots of examples>
    """

    #*****************************************************************************
    #       Copyright (C) 2013 YOUR NAME <your email>
    #
    # This program is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 2 of the License, or
    # (at your option) any later version.
    #                  http://www.gnu.org/licenses/
    #*****************************************************************************

Ad esempio, vedi ``SAGE_ROOT/src/sage/rings/integer.pyx`` che contiene 
l'implementazione di `\ZZ`. La sezione ``AUTHORS:`` \`e ridondante, poich`\e il log 
d'autorit\`a per chi ha scritto cosa \`e sempre il repository git (vedi l'output di 
``git blame``). Cionondimeno \`e a volte utile avere una panoramica molto semplificata 
sulla history, specialmente se molte persone hanno lavorato su quel file sorgente.

Tutto il codice incluso in Sage deve avere licenza o GPLv2+ o una compatibile, cio\`e 
meno restrittiva (ad esempio la licenza BSD).


.. _section-docstrings:

Stringhe di documentazione (doctring)
=====================================

.. _section-docstring-function:

La docstring di una funzione: contenuto
---------------------------------------

**Ogni** funzione deve avere una docstring che includa le seguenti 
informazioni. Puoi usare le funzioni gi\`a presenti in Sage come template.

-  La **descrizione in una frase** della funzione.

   Dev'essere seguita da una linea vuota e terminare con un punto. Descrive 
   gli effetti della funzione o del metodo come un comando ("Fa questo",
   "Restituisce quest'altro"), non come "Restituisce il pathname ...".

-  Un blocco **INPUT** ed un blocco **OUTPUT** che descrivono l'input/output 
   della funzione. Questo non \`e opzionale.

   Il blocco INPUT descrive tutti gli argomenti che la funzione accetta, 
   ed il blocco OUTPUT descrive l'output che ci si aspetta.

   1. I nomi di tipo devono essere descrittivi, ma non necessariamente rappresentare 
      i tipi esatti di Sage/Python. Ad esempio usare "integer" per qualunque cosa si 
      comporti come un intero, piuttosto che ``int``.

   2. segnala i valori di default degli argumenti di input quando ci sono.

   Esempio::

       INPUT:

       - ``p`` -- (default: 2) a positive prime integer.

       OUTPUT:

       A 5-tuple consisting of integers in this order:

       1. the smallest primitive root modulo p
       2. the smallest prime primitive root modulo p
       3. the largest primitive root modulo p
       4. the largest prime primitive root modulo p
       5. total number of prime primitive roots modulo p

   Puoi iniziare il blocco OUTPUT con un trattino se preferisci::

       OUTPUT:

       - The plaintext resulting from decrypting the ciphertext ``C``
         using the Blum-Goldwasser decryption algorithm.

-  Un blocco **EXAMPLES** per gli esempi. Questo non \`e opzionale.

   Questi esempio sono utilizzati per:

   1. Documentazione
   2. Test automatici prima di ogni nuova release.

   Dovrebbero coprire bene tutte le funzionalit\`a in questione.

-  Un blocco **SEEALSO** (caldamente raccomandato) con collegamenti a parti di
   Sage in relazione. Questo aiuta gli utenti a trovare le funzionalit\`a di interesse 
   e a scoprirne di nuove. ::

       .. SEEALSO::

           :ref:`chapter-sage_manuals_links`,
           :meth:`sage.somewhere.other_useful_method`,
           :mod:`sage.some.related.module`.

   Vedi :ref:`chapter-sage_manuals_links` per dettagli su come fare dei 
   link in Sage.

-  Un blocco **ALGORITHM** (opzionale).

   Indica quale algoritmo e/o quale software \`e utilizzato, ad esempio 
   ``ALGORITHM: Uses Pari``. Qui di seguito vediamo un esempio un po' pi\`u 
   lungo con delle referenze bibliografiche::

       ALGORITHM:

       The following algorithm is adapted from page 89 of [Nat2000]_.

       Let `p` be an odd (positive) prime and let `g` be a generator
       modulo `p`. Then `g^k` is a generator modulo `p` if and only if
       `\gcd(k, p-1) = 1`. Since `p` is an odd prime and positive, then
       `p - 1` is even so that any even integer between 1 and `p - 1`,
       inclusive, is not relatively prime to `p - 1`. We have now
       narrowed our search to all odd integers `k` between 1 and `p - 1`,
       inclusive.

       So now start with a generator `g` modulo an odd (positive) prime
       `p`. For any odd integer `k` between 1 and `p - 1`, inclusive,
       `g^k` is a generator modulo `p` if and only if `\gcd(k, p-1) = 1`.

       REFERENCES:

       .. [Nat2000] M.B. Nathanson. Elementary Methods in Number Theory.
          Springer, 2000.

-  Un blocco **NOTE** per suggerimenti e trucci (opzionale). ::

       .. NOTE::

           You should note that this sentence is indented at least 4
           spaces. Never use the tab character.

- Un blocco **WARNING** per informazioni critiche sul codice (opzionale).

  Ad esempio situazioni note in cui il codice va in errore, o qualunque cosa di 
  cui l'utente deve essere al corrente. ::

      .. WARNING::

          Whenever you edit the Sage documentation, make sure that
          the edited version still builds. That is, you need to ensure
          that you can still build the HTML and PDF versions of the
          updated documentation. If the edited documentation fails to
          build, it is very likely that you would be requested to
          change your patch.

- Un blocco **TODO** per miglioramenti futuri (opzionale).

  Pu\`o contenere doctest disabilitati per dimostrare la funzionalit\`a 
  desiderata. Ecco un esempio di blocco TODO::

      .. TODO::

          Add to ``have_fresh_beers`` an interface with the faster
          algorithm "Buy a Better Fridge" (BaBF)::

              sage: have_fresh_beers('Bière de l\'Yvette', algorithm="BaBF") # not implemented
              Enjoy !

- Un blocco **PLOT** per illustrare con figure l'output della funzione.

  Genera con codice Sage un oggetto ``g`` con un metodo ``.plot``, poi chiama 
  ``sphinx_plot(g)``::

      .. PLOT::

          g = graphs.PetersenGraph()
          sphinx_plot(g)

- Un blocco **REFERENCES** per elencare libri o articoli collegati (opzionale)

  Dovrebbe citare i libri/articoli di ricerca rilevanti per il codice, ad esempio 
  il sorgente dell'algoritmo che implementa. ::

      This docstring is referencing [SC]_. Just remember that references
      are global, so we can also reference to [Nat2000]_ in the ALGORITHM
      block, even if it is in a separate file. However we would not
      include the reference here since it would cause a conflict.

      REFERENCES:

      .. [SC] Conventions for coding in sage.
         http://www.sagemath.org/doc/developer/conventions.html.

  Vedi `markup Sphinx/ReST per citazioni <http://sphinx.pocoo.org/rest.html#citations>`_. Per link a tickets Trac o wikipedia, vedi :ref:`chapter-sage_manuals_links`.

- Un blocco **TESTS** (opzionale)

  Formattato come EXAMPLES, contiene test non rilevanti per gli utenti.

Template
^^^^^^^^

Usa il seguente template quando documenti delle funzioni. Nota l'indentazione: 

.. skip    # do not doctest

::

    def point(self, x=1, y=2):
        r"""
        Return the point `(x^5,y)`.

        INPUT:

        - ``x`` -- integer (default: 1) the description of the
          argument ``x`` goes here.  If it contains multiple lines, all
          the lines after the first need to begin at the same indentation
          as the backtick.

        - ``y`` -- integer (default: 2) the ...

        OUTPUT:

        The point as a tuple.

        .. SEEALSO::

            :func:`line`

        EXAMPLES:

        This example illustrates ...

        ::

            sage: A = ModuliSpace()
            sage: A.point(2,3)
            xxx

        We now ...

        ::

            sage: B = A.point(5,6)
            sage: xxx

        It is an error to ...::

            sage: C = A.point('x',7)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'r' to an integer

        .. NOTE::

            This function uses the algorithm of [BCDT]_ to determine
            whether an elliptic curve `E` over `Q` is modular.

        ...

        REFERENCES:

        .. [BCDT] Breuil, Conrad, Diamond, Taylor,
           "Modularity ...."
        """
        <body of the function>

Sei caldamente incoraggiato a:

- Usare le convenzioni di scrittura di LaTeX (vedi :ref:`section-latex-typeset`).

- Descrivere ampiamente cosa fanno gli esempi.

  .. NOTE::

     Cid deve essere una riga vuota dopo il codice di esempio e prima del 
     testo di spiegazione dell'esempio successivo (l'indentazione non basta).

- Illustra le eccezioni sollevate dalla funzione con degli esempi (come dato 
  sopra: "\`E un errore [..]", ...)

- Includi molti esempi.

  Sono utili agli utenti, e sono fondamentali per la qualit\`a e l'adattabilit\`a 
  di Sage. Senza tali esempi, piccoli cambiamenti ad una parte di Sage che 
  danneggiano qualcos'altro potrebbero non essere scoperti fino a parecchio tempo 
  dopo quando qualcuno usa il sistema, cosa che \`e inaccettabile.

Funzioni private
^^^^^^^^^^^^^^^^

Le funzioni i cui nomi iniziano con una sottolineatura (underscore) sono considerati 
privati. Non compaiononel manuale di riferimento, ed i loro docstring non dovrebbero 
contenere informazioni cruciali per gli utenti di Sage. Puoi rendere i loro 
docstring parte della documentazione di un altro metodo. Ad esempio::

    class Foo(SageObject):

        def f(self):
            """
            <usual docstring>

            .. automethod:: _f
            """
            return self._f()

        def _f(self):
             """
             This would be hidden without the ``.. automethod::``
             """

Le funzioni private dovrebbero contenere un blocco EXAMPLES (o TESTS).

Un caso speciale \`e il costruttore ``__init__``: per il suo speciale
status, la doctring di ``__init__`` \`e utilizzata come docstring della 
classe se non ce n'\`e gi\`a una. Cio\`e si pu\`o fare quanto segue::

    sage: class Foo(SageObject):
    ....:     # no class docstring
    ....:     def __init__(self):
    ....:         """Construct a Foo."""
    sage: foo = Foo()
    sage: from sage.misc.sageinspect import sage_getdoc
    sage: sage_getdoc(foo)              # class docstring
    'Construct a Foo.\n'
    sage: sage_getdoc(foo.__init__)     # constructor docstring
    'Construct a Foo.\n'

.. _section-latex-typeset:

Convenzioni di scrittura LaTeX
------------------------------

Nella documentazione di Sage il codice LaTeX \`e permesso ed \`e marcato con 
**accenti obliqui o simboli di dollaro**:

    ```x^2 + y^2 = 1``` and ``$x^2 + y^2 = 1$`` both yield `x^2 + y^2 = 1`.

**Barre retroverse (backslash):** Per comandi LaTeX contenenti delle backslash, puoi 
o usare dei doppi backslash o iniziare la docstring con un ``r"""`` invece di ``"""``. 
Entrambe le scritture seguenti sono valide::

    def cos(x):
        """
        Return `\\cos(x)`.
        """

    def sin(x):
        r"""
        Return $\sin(x)$.
        """

**Blocco MATH:** Questo \`e simile alla sintassi LaTeX ``\[<math expression>\]`` 
(oppure ``$$<math expression>$$``). Ad esempio::

    .. MATH::

        \sum_{i=1}^{\infty} (a_1 a_2 \cdots a_i)^{1/i}
        \leq
        e \sum_{i=1}^{\infty} a_i

.. MATH::

    \sum_{i=1}^{\infty} (a_1 a_2 \cdots a_i)^{1/i}
    \leq
    e \sum_{i=1}^{\infty} a_i

L'ambiente **aligned** funziona nello stesso modo che in LaTeX::

    .. MATH::

        \begin{aligned}
         f(x) & = x^2 - 1 \\
         g(x) & = x^x - f(x - 2)
        \end{aligned}

.. MATH::

    \begin{aligned}
     f(x) & = x^2 - 1 \\
     g(x) & = x^x - f(x - 2)
    \end{aligned}

Quando si compila la documentazione in PDF, tutto \`e tradotto in LaTeX 
ed ogni blocco MATH \`e automaticalmente passato in un ambiente ``math`` --
in particolare, \`e convertito in ``\begin{gather} block \end{gather}``. 
Quindi se vuoi usare un ambiente LaTeX (come ``align``) che in LaTeX 
ordinario non sarebbe passato cos\`i, devi aggiungere un flag **:nowrap:** 
al modo MATH. Vedi anche `documentazione Sphinx per blocchi math 
<http://sphinx-doc.org/latest/ext/math.html?highlight=nowrap#directive-math>`_. ::

    .. MATH::
       :nowrap:

       \begin{align}
          1+...+n &= n(n+1)/2\\
          &= O(n^2)\\
       \end{tabular}

.. MATH::
   :nowrap:

   \begin{align}
   1+...+n &= n(n+1)/2\\
   &= O(n^2)\\
   \end{align}

**Equilibrio di leggibilit\`a:** nella console interattiva, le formule LaTeX contenute nella 
documentazione sono rappresentate con il loro codice LaTeX (dove le backslash sono 
state tolte). In tale situazione ``\\frac{a}{b}`` \`e meno leggibile di ``a/b`` oppure 
``a b^{-1}`` (alcuni utenti potrebbero anche non conoscere LaTeX). Cerca di rendere il testo 
leggibile da tutti per quanto ti \`e possibile.

**Ring comuni** `(\Bold{Z},\Bold{N},...)`: Lo stile LaTeX di Sage \`e di evidenziare gli anelli  
e campi standard usando la macro definita localmente ``\\Bold`` (ad esempio ``\\Bold{Z}`` da `\Bold{Z}`).

**Abbreviazioni**: Sono disponibili per mantenere la leggibilit\`a, ad esempio ``\\ZZ`` (`\ZZ`),
``\\RR`` (`\RR`), ``\\CC`` (`\CC`), e ``\\QQ`` (`\QQ`). Appaiono formattate in LaTeX ``\\Bold{Z}`` 
nel manuale in html, e come ``Z`` nell'help interattivo. Altri esempi sono: ``\\GF{q}``, (`\GF{q}`) 
e ``\\Zmod{p}`` (`\Zmod{p}`).

Vedi il file ``SAGE_ROOT/src/sage/misc/latex_macros.py`` per una lista completa e per dettagli 
sul come aggiungere altre macro.

.. _section-doctest-writing:

Scrivere esempi adatti ai test
------------------------------

Gli esempi dalla documentazione di Sage hanno un doppio scopo:

- Forniscono **illustrazioni** dell'uso del codice agli utenti

- Sono dei **test** che sono verificati prima di ogni release, e che 
  ci aiutano ad evitare nuovi bachi.

Tutti i nuovi doctest aggiunti a Sage devono **passare tutti i test** (vedi 
:ref:`chapter-doctesting`), cio\`e eseguire ``sage -t your_file.py`` non deve 
dare alcun messaggio di errore. Sotto ci sono instruzione riguardo a come 
devono essere scritti i doctest.

**Di cosa devono verificare i doctest:**

- **Esempi interessanti** di ci\`o che una funzione pu\`o fare. Questa sar\`a la 
  cosa pi\`u utile per un utente smarrito. \`E anche l'occasione per verificare 
  teoremi famosi (a proposito)::

    sage: is_prime(6) # 6 is not prime
    False
    sage: 2 * 3 # and here is a proof
    6

- Tuute le **combinazioni significative** degli argomenti di input. Ad esempio 
  una funzione pu\`o accettare un argomento ``algorithm="B"``, ed i doctest devono 
  verificare sia ``algorithm="A"`` che ``algorithm="B"``.

- **Casi limite:** il codice dev'essere capace di gestire input 0, o un insieme 
  vuoto, o una matrice nulla, o una funzione nulla, ... Tutti i casi limite vanno 
  verificati, essendo quello che pi\`u probabilmente daranno problemi, ora o nel 
  futuro. Questo spesso andr\`a messo nel blocco TESTS (vedi :ref:`section-docstring-function`).

- **Test sistematici** di tutti gli input piccoli, o test di valori a caso (**random**) 
  se possibile.

  .. NOTE::

     Nota che le **suite di test** sono un modo automatico di generare alcuni di
     questi test in specifiche situazioni. Vedi
     ``SAGE_ROOT/src/sage/misc/sage_unittest.py``.

**La sintassi:**

- **Ambiente:** i doctest dovrebbero funzionare se fai copia/incolla nella console 
  interattiva di Sage. Ad esempio, la funzione ``AA()`` nel file
  ``SAGE_ROOT/src/sage/algebras/steenrod/steenrod_algebra.py`` include un blocco 
  EXAMPLES contenente il seguente::

    sage: from sage.algebras.steenrod.steenrod_algebra import AA as A
    sage: A()
    mod 2 Steenrod algebra, milnor basis

  Sage non conosce la funzione ``AA()`` di default, quindi ha bisogno di 
  importarla prima di farne il test. Da qui la prima linea dell'esempio.

- **Preparse:** Come nella console di Sage, `4/3` restituisce `4/3` e non `1` come in
  Python 2.7. I test vengono fatti con il preparse completo di Sage sull'input nell'ambiente  
  shell standard di Sage, come descritto in :ref:`section-preparsing`.

- **Scrivere file:** Se un test manda dell'output su un file, tale file dev'essere temporaneo. 
  Usa :func:`tmp_filename` per avere un nome di file temporaneo, oppure :func:`tmp_dir` per 
  avere una directory temporanea. Vedi ad esempio ``SAGE_ROOT/src/sage/plot/graphics.py``)::

      sage: plot(x^2 - 5, (x, 0, 5), ymin=0).save(tmp_filename(ext='.png'))

- **Doctest multilinea:** Puoi scrivere dei test che occupano multe linee, usando il carattere 
  di continuazione di linea ``....:`` ::

      sage: for n in srange(1,10):
      ....:     if n.is_prime():
      ....:         print(n)
      2
      3
      5
      7

- **Spezzare linee lunghe:** Potresti voler spezzare linee di codice lunghe con 
  una backslash. Nota: questa sintassi non \`e standard e potrebbe essere deprecata 
  in futuro::

      sage: n = 123456789123456789123456789\
      ....:     123456789123456789123456789
      sage: n.is_prime()
      False

- **Flag di doctest:** sono disponibili dei flag per cambiare il comportamento dei doctest:
  see :ref:`section-further_conventions`.

.. _section-further_conventions:

Markup speciale per influenzare i test
--------------------------------------

Ci sono un certo numero commenti "magici" che puoi mettere nel codice di esempio, che 
cambiano il modo in cui l'output \`e verificato dal framework di doctest di Sage. 
Eccone una lista completa:

- **casuale:** La linea sar\`a eseguita, ma il suo output non sar\`a verificato con 
  l'output nella stringa di documentazione::

      sage: c = CombinatorialObject([1,2,3])
      sage: hash(c)  # random
      1335416675971793195
      sage: hash(c)  # random
      This doctest passes too, as the output is not checked

  Comunque la maggior parte delle funzioni che generano output pseudocasuale non richiedono 
  questo tag poich\`e il framework di doctest garantisce lo stato dei generatori di numeri 
  pseudocasuali (PRNGs) usato in Sage per un dato doctest.

  Quando possibile, evita il problema, ad esempio: piuttosto di verificare il valore 
  dell'hash in un doctest, pu\`o andare altrettanto bene usarlo come chiave in un dict.

- **richiede molto tempo:** La linea \`e solo testata se \`e data l'opzione ``--long``, 
  ad esmpio ``sage -t --long f.py``.

  Usala per doctest che richiedono pi\`u di 1 secondo per essere eseguiti. Nessun esempio 
  dovrebbe richiedere pi\`u di 30 secondi::

      sage: E = EllipticCurve([0, 0, 1, -1, 0])
      sage: E.regulator()        # long time (1 second)
      0.0511114082399688

- **tol** o **tolleranza:** I valori numerici restituiti dalla linea sono solo 
  verificati entro una data tolleranza. \`E utile quando l'output \`e soggetto a 
  imprecisione numerica per cause dipendenti dal sistema (aritmetica floating-point, math
  libraries, ...) o per la scelta di algoritmi non-deterministici.

  - Pu\`o avere prefisso ``abs[olute]`` oppure ``rel[ative]`` per specificare se 
    misurare l'errore **assoluto** o **relativo** (vedi :wikipedia:`Approximation_error`).

  - Se non \`e specificato ``abs/rel``, si considera l'errore ``absolute`` quando il 
    valore atteso \`e **zero**, e ``relative`` per valori **diversi da zero**.

  ::

     sage: n(pi)  # abs tol 1e-9
     3.14159265358979
     sage: n(pi)  # rel tol 2
     6
     sage: n(pi)  # abs tol 1.41593
     2
     sage: K.<zeta8> = CyclotomicField(8)
     sage: N(zeta8)  # absolute tolerance 1e-10
     0.7071067812 + 0.7071067812*I

  **Valori numerici multipli:** la rappresentazione dei numeri complessi, le 
  matrici, ed i polinomi di solito richiede parecchi valori numerici. Se un 
  doctest con tolleranza contiene parecchi numeri, ognuno di essi \`e verificato 
  individualmente::

      sage: print("The sum of 1 and 1 equals 5")  # abs tol 1
      The sum of 2 and 2 equals 4
      sage: e^(i*pi/4).n() # rel tol 1e-1
      0.7 + 0.7*I
      sage: ((x+1.001)^4).expand() # rel tol 2
      x^4 + 4*x^3 + 6*x^2 + 4*x + 1
      sage: M = matrix.identity(3) + random_matrix(RR,3,3)/10^3
      sage: M^2 # abs tol 1e-2
      [1 0 0]
      [0 1 0]
      [0 0 1]

  I valori che il framework di doctest assume nel calcolo degli errori sono
  definiti dall'espressione regolare ``float_regex`` in
  :mod:`sage.doctest.parsing`.

- **non implementato** oppure **non testato:** La linea non \`e mai testata.

  Usala per doctest molto lunghi che sono solo intesi come documentazione. Pu\`o anche 
  essere usata per note su ci\`o che dovr\`a essere implementato successivamente::

      sage: factor(x*y - x*z)    # todo: not implemented

  Dev'essere anche immediatamente chiaro all'utente che gli esempi indicati non 
  funzionano ancora.

  .. NOTE::

     Salta tutti i doctest di un file/directory

     - **file:** Se una delle prime 10 linee di un file inizia con una delle seguenti 
       ``r""" nodoctest`` (o ``""" nodoctest`` o ``# nodoctest`` o ``%
       nodoctest`` o ``.. nodoctest``, o qualunque di queste con spaziature diverse),
       allora quel file sar\`a saltato.

     - **directory:** Se una directory contiene un file ``nodoctest.py``, allora 
       l'intera directory sar\`a saltata.

     Nessuna di queste si applica a file o directory che sono date esplicitamente 
     come argomenti a linea di comando: di quelli viene sempre fatto il test.

- **optional:** Di una linea con flag ``optional - keyword`` non \`e fatto il test 
  a meno che non sia passato il flag ``--optional=keyword`` a ``sage -t`` (vedi 
  :ref:`section-optional-doctest-flag`). Le principali applicazioni sono:

  - **optional packages:** Quando una linea richiede di installare un pacchetto opzionale 
    (ad esempio il pacchetto ``sloane_database``)::

      sage: SloaneEncyclopedia[60843]    # optional - sloane_database

  - **internet:** Per linee che rechiedono una connessione ad Internet::

       sage: sloane_sequence(60843)       # optional - internet

  - **bug:** Per linee che descrivono dei bachi. In alternativa usa ``# known bug``
    al posto: \`e un alias per ``optional bug``. ::

        The following should yield 4.  See :trac:`2`. ::

            sage: 2+2  # optional: bug
            5
            sage: 2+2  # known bug
            5

  .. NOTE::

      - Tutte le parole dopo ``# optional`` sono interpretate come una lista di
        nomi di pacchetto, separati da spazi.

      - Ogni punteggiatura (punti, virgole, trattini, due punti, ...) dopo la prima 
        parola, termina la lista dei pacchetti. Trattini o punti e virgola fra la 
        parola ``optional`` ed il primo nome di pacchetto sono permessi. Pertanto, 
        non dovresti scrivere ``optional: needs package CHomP`` ma semplicemente 
        ``optional: CHomP``.

      - I tag opzionali sono indifferenti alle maiscole-minuscole, quindi puoi anche 
        scrivere ``optional: chOMP``.

- **doctest indiretti:** nella docstring di una function ``A(...)``, una linea che 
  chiama ``A`` e nel cui nome ``A`` non appare dovrebbe avere questo flag. 
  Questo evita che ``sage --coverage <file>`` riporti la docstring come 
  "not testing what it should test".

  Usala quando fai il test di funzioni speciali come ``__repr__``, ``__add__``,
  ecc. Usala anche quando fai il test di funzioni chiamando ``B`` che chiama 
  internamente ``A``::

      Questa \`e la docstring di un metodo ``__add__``. Il seguente esempio ne fa il 
      tests, ma ``__add__`` non \`e scritta da nessuna parte::

          sage: 1+1 # indirect doctest
          2

- **32-bit** o **64-bit:** per test che si comportano differentemente su macchine a 
  32-bit o a 64-bits. Nota che questo particolare flag va applicato sulle linee di
  **output**, non su quelle di input::

      sage: hash(-920390823904823094890238490238484)
      -873977844            # 32-bit
      6874330978542788722   # 64-bit

Usando ``search_src`` dal prompt di Sage (oppure ``grep``), si possono trovare 
facilmente le parole chiave suddette. Nel caso di ``todo: not implemented``, si 
possono usare i risultati di tale ricerca per dirigere l'ulteriore sviluppo di 
Sage.

.. _chapter-testing:

Eseguire i test automatici
==========================

Questa sezione descrive i test automatici di Sage di file dei seguenti tipo: 
``.py``, ``.pyx``, ``.sage``, ``.rst``. In breve, usa ``sage -t <file>`` per 
fare il test che gli esempi in ``<file>`` si comportino esattamente come dichiarato. 
Vedi le seguenti sottosezioni per maggiori dettagli. Vedi anche :ref:`section-docstrings` 
per una discussione su come includere esempi nelle stringhe di documentazione e quali 
convenzioni seguire. Il capitolo :ref:`chapter-doctesting` contiene un tutorial su come 
fare i doctest dei moduli nella libreria Sage.


.. _section-testpython:

Fare i test dei file .py, .pyx e .sage
--------------------------------------

Esegui ``sage -t <filename.py>`` per fare il test di tutti gli esempi di codice in
``filename.py``. Analogamente per i file ``.sage`` e ``.pyx``::

      sage -t [--verbose] [--optional]  [files and directories ... ]

Il framework di doctest di Sage \`e basato sul modulo doctest del Python standard, 
ma con molte funzionalit\`a addizionali (come i test paralleli, i timeout, i test 
opzionali). I processore di doctest di Sage riconosce il prompt ``sage:`` cos\`i 
come il prompt ``>>>``. Fa anche il preparse dei doctest, cos\`i come nelle 
sessioni interattive di Sage.

Il tuo file passer\`a i test se il codice in esso \`e in grado di essere eseguito 
quando immesso al prompt ``sage:`` senza delle import extra. Cos\`i si garantisce 
agli utenti di poter copiare esattamente il codice degli esempi che scrivi per la 
documentazione e che essi funzionino.

Per maggiori informazioni, vedi :ref:`chapter-doctesting`.


Fare il test della documentazione ReST
--------------------------------------

Esegui ``sage -t <filename.rst>`` per testare gli esempi verbatim (parola per parola) 
nella documentazione ReST.

Naturalmente nei file ReST spesso si inseriscono delle frasi di spiegazione fra 
ambienti differenti. Per collegare insieme ambienti verbatim, usa il commento ``.. link``.
Ad esempio::

    EXAMPLES::

            sage: a = 1


    Next we add 1 to ``a``.

    .. link::

            sage: 1 + a
            2

Se vuoi collegare fra loro tutti gli ambienti verbatim, puoi mettere 
``.. linkall`` ovunque nel file, su una linea a s\`e.  (Per chiarezza, 
potrebbe essere meglio metterla vicino alla cima del file.) Allora 
``sage -t`` agir\`a come se ci fosse ``.. link`` davanti ad ogni ambiente 
verbatim. Il file ``SAGE_ROOT/src/doc/en/tutorial/interfaces.rst`` 
contiene una direttiva ``.. linkall``, ad esempio.

Puoi anche mettere ``.. skip`` subito davanti ad un ambiente verbatim perch\`e tale 
esempio sia saltato durante il test del file. Questo va nello stesso posto di 
``.. link`` nell'esempio precedente.

Vedi i file in ``SAGE_ROOT/src/doc/en/tutorial/`` per altri esempi su come 
includere test automatici nella documentazione ReST per Sage.

.. _chapter-picklejar:

Il vaso dei sottaceti (pickle)
==============================

Sage mantiene un vaso dei sottaceti (pickle) in 
``SAGE_ROOT/src/ext/pickle_jar/pickle_jar.tar.bz2`` che \`e un fille tar 
di pickle "standard" creati da ``sage``. Questo ``vaso di pickle`` \`e 
utilizzato per garantire che Sage mantenga la compatibilit\`a all'indietro facendo 
s\`i che :func:`sage.structure.sage_object.unpickle_all` verifichi che ``sage`` 
possa sempre prendere tutti i pickle nel vaso come parte del framework standard 
di doctest.

La maggior parte delle persone si imbattono nella pickle_jar quando le loro patch 
vanno in errore in tale fase di "presa dei sottaceti" durante i doctest::

    sage -t src/sage/structure/sage_object.pyx

Quando questo succede un messagio di errore \`e mostrato contenente i seguenti 
suggerimenti per correggere i "sottaceti immangiabili"::

    ----------------------------------------------------------------------
    ** This error is probably due to an old pickle failing to unpickle.
    ** See sage.structure.sage_object.register_unpickle_override for
    ** how to override the default unpickling methods for (old) pickles.
    ** NOTE: pickles should never be removed from the pickle_jar!
    ----------------------------------------------------------------------

Per maggiori dettagli su come correggere gli errori sui pickle vedi :func:`sage.structure.sage_object.register_unpickle_override`

.. WARNING::

    Il vaso dei sottaceti di Sag aiuta ad assicurare la compatibilit\`a all'indietro 
    in Sage. I pickles vanno rimossi dal vaso **solo** quando i corrispondenti oggetti 
    sono stati adeguatamente deprecati. Ogni proposta di rimozione dei sottaceti dal vaso 
    va prima discussa su ``sage-devel``.


Opzioni globali
===============

Opzioni globali per le classi possono essere definite in Sage usando
:class:`~sage.structure.global_options.GlobalOptions`.

