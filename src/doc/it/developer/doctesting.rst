.. nodoctest

.. _chapter-doctesting:

==========================
Eseguire i doctest di Sage
==========================

L'eseguire i doctest \`e di una funzione assicura che essa si comporti 
come dichiarato nella sua documentazione. Il test pu\`o essere fatto 
usano uno o pi\`u thread. Dopo aver compilato una versione di Sage da 
sorgente, il doctest pu\`o essere lanciato sull'intera libreria Sage, 
su tutti i moduli in una determinata directory, o solo su uno specifico 
modulo. Per gli scopi di questo capitolo, supponiamo di aver compilato 
Sage 6.0 da sorgente e che la directory radice di Sage sia::

    [jdemeyer@sage sage-6.0]$ pwd
    /scratch/jdemeyer/build/sage-6.0

Vedi la sezione :ref:`chapter-testing` per informazioni sul processo di 
test automatico di Sage. La sintassi generale per fare i doctest \`e la 
seguente. Per fare il doctest di un modulo nella libraria di una versione 
di Sage, usa la seguente sintassi::

    /path/to/sage-x.y.z/sage -t [--long] /path/to/sage-x.y.z/path/to/module.py[x]

dove ``--long`` \`e un argomento opzionale (vedi :ref:`section-options` per 
le altre opzioni). La versione di ``sage`` usata deve essere uguale alla versione
di Sage contenente il modulo di cui vogliamo fare il doctest. Un modulo Sage pu\`o 
essere sia uno script Python (con l'estensione ".py") o pu\`o essere uno script 
Cython, nel qual caso ha l'estensione ".pyx".


Fare il test di un modulo
=========================

Diciamo che vogliamo eseguire tutti i test nel modulo sudoku 
``sage/games/sudoku.py``. In una finestra del terminale, prima diamo ``cd`` alla 
directory radice di Sage della nostra installazione locale. A questo punto 
possiamo lanciare i doctest come nella seguente sessione del terminale::

    [jdemeyer@sage sage-6.0]$ ./sage -t src/sage/games/sudoku.py
    Running doctests with ID 2012-07-03-03-36-49-d82849c6.
    Doctesting 1 file.
    sage -t src/sage/games/sudoku.py
        [103 tests, 3.6 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 4.8 seconds
        cpu time: 3.6 seconds
        cumulative wall time: 3.6 seconds

I numeri prodotti dal test mostrano che fare il test del modulo sudoku 
richiede circa 4 secondi, mentre fare il test di tutti i moduli specificati 
richiede circa lo stesso tempo; il tempo totale richiesto include un po' di 
tempo di avvio per il codice che esegue i test. In questo caso, abbiamo solo 
fatto il test di un modulo quindi non ci sorprende il fatto che il tempo totale 
di test \`e approssimativamente lo stesso del tempo richiest per fare il test 
solo di quell'unico modulo. Nota che la sintassi \`e::

    [jdemeyer@sage sage-6.0]$ ./sage -t src/sage/games/sudoku.py
    Running doctests with ID 2012-07-03-03-39-02-da6accbb.
    Doctesting 1 file.
    sage -t src/sage/games/sudoku.py
        [103 tests, 3.6 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 4.9 seconds
        cpu time: 3.6 seconds
        cumulative wall time: 3.6 seconds

ma non::

    [jdemeyer@sage sage-6.0]$ ./sage -t sage/games/sudoku.py
    Running doctests with ID 2012-07-03-03-40-53-6cc4f29f.
    No files matching sage/games/sudoku.py
    No files to doctest

Possiamo anche prima fare ``cd`` alla directory contenente il modulo
``sudoku.py`` e fare il doctest del modulo come segue::

    [jdemeyer@sage sage-6.0]$ cd src/sage/games/
    [jdemeyer@sage games]$ ls
    __init__.py  hexad.py       sudoku.py           sudoku_backtrack.pyx
    all.py       quantumino.py  sudoku_backtrack.c
    [jdemeyer@sage games]$ ../../../../sage -t sudoku.py
    Running doctests with ID 2012-07-03-03-41-39-95ebd2ff.
    Doctesting 1 file.
    sage -t sudoku.py
        [103 tests, 3.6 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 5.2 seconds
        cpu time: 3.6 seconds
        cumulative wall time: 3.6 seconds

In tutte le suddette sessioni di terminale, abbiamo utilizzato un installazione 
locale di Sage per fare il test dei suoi moduli. Anche se abbiamo una installazione 
di Sage a livello di sistema, usare tale versione per fare i doctest dei moduli di 
un'installazione locale \`e una ricetta per la confusione.


Risoluzione dei problemi
========================

Per fare il doctest dei moduli di un'installazione di Sage, dal terminale diamo 
prima ``cd`` alla directory radice di Sage, altrimenti nota come ``SAGE_ROOT`` 
di tale installazione. Quando eseguiamo i test, usiamo tale particolare installazione 
di Sage con la sintassi ``./sage``; nota il "punto e barra" davanti a ``sage``. 
Questa \`e una precauzione contro la confusione che potrebbe nascere quando sul nostro 
sistema vi sono pi\`u installazioni di Sage. Ad esempio, la seguente sintassi \`e 
accettabile perch\`e specifichiamo esplicitamente l'installazione di Sage nella 
``SAGE_ROOT`` corrente::

    [jdemeyer@sage sage-6.0]$ ./sage -t src/sage/games/sudoku.py
    Running doctests with ID 2012-07-03-03-43-24-a3449f54.
    Doctesting 1 file.
    sage -t src/sage/games/sudoku.py
        [103 tests, 3.6 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 4.9 seconds
        cpu time: 3.6 seconds
        cumulative wall time: 3.6 seconds
    [jdemeyer@sage sage-6.0]$ ./sage -t "src/sage/games/sudoku.py"
    Running doctests with ID 2012-07-03-03-43-54-ac8ca007.
    Doctesting 1 file.
    sage -t src/sage/games/sudoku.py
        [103 tests, 3.6 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 4.9 seconds
        cpu time: 3.6 seconds
        cumulative wall time: 3.6 seconds

LA seguente sintassi non \`e raccomandata dal momento che stiamo usando un'installazione 
di Sage a livello di sistema(se esiste):

.. skip

::

    [jdemeyer@sage sage-6.0]$ sage -t src/sage/games/sudoku.py
    sage -t  "src/sage/games/sudoku.py"
    **********************************************************************
    File "/home/jdemeyer/sage/sage-6.0/src/sage/games/sudoku.py", line 515:
        sage: next(h.solve(algorithm='backtrack'))
    Exception raised:
        Traceback (most recent call last):
          File "/usr/local/sage/local/bin/ncadoctest.py", line 1231, in run_one_test
            self.run_one_example(test, example, filename, compileflags)
          File "/usr/local/sage/local/bin/sagedoctest.py", line 38, in run_one_example
            OrigDocTestRunner.run_one_example(self, test, example, filename, compileflags)
          File "/usr/local/sage/local/bin/ncadoctest.py", line 1172, in run_one_example
            compileflags, 1) in test.globs
          File "<doctest __main__.example_13[4]>", line 1, in <module>
            next(h.solve(algorithm='backtrack'))###line 515:
        sage: next(h.solve(algorithm='backtrack'))
          File "/home/jdemeyer/.sage/tmp/sudoku.py", line 607, in solve
            for soln in gen:
          File "/home/jdemeyer/.sage/tmp/sudoku.py", line 719, in backtrack
            from sudoku_backtrack import backtrack_all
        ImportError: No module named sudoku_backtrack
    **********************************************************************
    [...more errors...]
    2 items had failures:
       4 of  15 in __main__.example_13
       2 of   8 in __main__.example_14
    ***Test Failed*** 6 failures.
    For whitespace errors, see the file /home/jdemeyer/.sage//tmp/.doctest_sudoku.py
             [21.1 s]

    ----------------------------------------------------------------------
    The following tests failed:


            sage -t  "src/sage/games/sudoku.py"
    Total time for all tests: 21.3 seconds

In questo caso, abbiamo un errore poich\`e because l'installazione di Sage 
a livello di sistema \`e di una versione differente (pi\`u vecchia) di quella 
che stiamo usando per fare sviluppo. Accertatisempre di fare i test dei file 
della versione corretta di Sage.

Fare test paralleli di molti moduli
===================================

Finora abbiamo usato un thread singolo per fare i doctest di un modulo nella
libraria Sage. Ci sono centinaia o migliaia di moduli nella libreria Sage. 
Fare il test di tutti usando un singolo thread richiederebbe alcune ore. 
In base al nostro hardware, pu\`o richiedere da 6 ore in s\`u. Su un sistema 
dotato di pi\`u di un core, effettuare i doctest in parallelo pu\`o ridurre 
significativamente il tempo dei test. Se non abbiamo bisogno di usare il nostro 
computer mentre effettuiamo i doctest in parallelo, possiamo scegliere di dedicare 
tutti i core del nostro sistema per fare i test.

Effettuiamo i doctest di tutti i moduli in una directory, dapprima usando un thread 
singolo e poi usando 4 thread. Per questo esempio, supponiamo che vogliamo testare 
tutti i moduli sotto ``sage/crypto/``. Possiamo usare una sintassi simile a quella 
mostrata sotto per fare ci\`o::

    [jdemeyer@sage sage-6.0]$ ./sage -t src/sage/crypto
    Running doctests with ID 2012-07-03-03-45-40-7f837dcf.
    Doctesting 24 files.
    sage -t src/sage/crypto/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/boolean_function.pyx
        [252 tests, 4.4 s]
    sage -t src/sage/crypto/cipher.py
        [10 tests, 0.0 s]
    sage -t src/sage/crypto/classical.py
        [718 tests, 11.3 s]
    sage -t src/sage/crypto/classical_cipher.py
        [130 tests, 0.5 s]
    sage -t src/sage/crypto/cryptosystem.py
        [82 tests, 0.1 s]
    sage -t src/sage/crypto/lattice.py
        [1 tests, 0.0 s]
    sage -t src/sage/crypto/lfsr.py
        [31 tests, 0.1 s]
    sage -t src/sage/crypto/stream.py
        [17 tests, 0.1 s]
    sage -t src/sage/crypto/stream_cipher.py
        [114 tests, 0.2 s]
    sage -t src/sage/crypto/util.py
        [122 tests, 0.2 s]
    sage -t src/sage/crypto/block_cipher/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/block_cipher/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/block_cipher/miniaes.py
        [430 tests, 1.3 s]
    sage -t src/sage/crypto/block_cipher/sdes.py
        [290 tests, 0.9 s]
    sage -t src/sage/crypto/mq/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/mq/mpolynomialsystem.py
        [320 tests, 9.1 s]
    sage -t src/sage/crypto/mq/mpolynomialsystemgenerator.py
        [42 tests, 0.1 s]
    sage -t src/sage/crypto/mq/sbox.py
        [124 tests, 0.8 s]
    sage -t src/sage/crypto/mq/sr.py
        [435 tests, 5.5 s]
    sage -t src/sage/crypto/public_key/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/public_key/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/public_key/blum_goldwasser.py
        [135 tests, 0.2 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 38.1 seconds
        cpu time: 29.8 seconds
        cumulative wall time: 35.1 seconds

Ora facciamo la stessa cosa, ma utilizzando l'argomento opzionale ``--long``::

    [jdemeyer@sage sage-6.0]$ ./sage -t --long src/sage/crypto/
    Running doctests with ID 2012-07-03-03-48-11-c16721e6.
    Doctesting 24 files.
    sage -t --long src/sage/crypto/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/all.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/boolean_function.pyx
        [252 tests, 4.2 s]
    sage -t --long src/sage/crypto/cipher.py
        [10 tests, 0.0 s]
    sage -t --long src/sage/crypto/classical.py
        [718 tests, 10.3 s]
    sage -t --long src/sage/crypto/classical_cipher.py
        [130 tests, 0.5 s]
    sage -t --long src/sage/crypto/cryptosystem.py
        [82 tests, 0.1 s]
    sage -t --long src/sage/crypto/lattice.py
        [1 tests, 0.0 s]
    sage -t --long src/sage/crypto/lfsr.py
        [31 tests, 0.1 s]
    sage -t --long src/sage/crypto/stream.py
        [17 tests, 0.1 s]
    sage -t --long src/sage/crypto/stream_cipher.py
        [114 tests, 0.2 s]
    sage -t --long src/sage/crypto/util.py
        [122 tests, 0.2 s]
    sage -t --long src/sage/crypto/block_cipher/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/block_cipher/all.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/block_cipher/miniaes.py
        [430 tests, 1.1 s]
    sage -t --long src/sage/crypto/block_cipher/sdes.py
        [290 tests, 0.7 s]
    sage -t --long src/sage/crypto/mq/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/mq/mpolynomialsystem.py
        [320 tests, 7.5 s]
    sage -t --long src/sage/crypto/mq/mpolynomialsystemgenerator.py
        [42 tests, 0.1 s]
    sage -t --long src/sage/crypto/mq/sbox.py
        [124 tests, 0.7 s]
    sage -t --long src/sage/crypto/mq/sr.py
        [437 tests, 82.4 s]
    sage -t --long src/sage/crypto/public_key/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/public_key/all.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/public_key/blum_goldwasser.py
        [135 tests, 0.2 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 111.8 seconds
        cpu time: 106.1 seconds
        cumulative wall time: 108.5 seconds

Nota la differenza di tempo fra il primo insieme di test ed il secondo, 
che usa l'argomento opzionale ``--long``. Molti test nella libreria Sage 
hanno il flag ``# long time`` perch\`e si sa che essi richiedono molto 
tempo per essere intieramente eseguiti. Se non si usa l'argomento opzionale 
``--long``, il modulo ``sage/crypto/mq/sr.py`` richiede circa 5 secondi. 
Con tale argumento ne richiede 82 per eseguire tutti i test di quel modulo. 
Ecco un pezzo di codice di una funzione nel modulo ``sage/crypto/mq/sr.py`` 
con un doctest con il flag che segnala che richiede molto tempo::

    def test_consistency(max_n=2, **kwargs):
        r"""
        Test all combinations of ``r``, ``c``, ``e`` and ``n`` in ``(1,
        2)`` for consistency of random encryptions and their polynomial
        systems. `\GF{2}` and `\GF{2^e}` systems are tested. This test takes
        a while.

        INPUT:

        - ``max_n`` -- maximal number of rounds to consider (default: 2)
        - ``kwargs`` -- are passed to the SR constructor

        TESTS:

        The following test called with ``max_n`` = 2 requires a LOT of RAM
        (much more than 2GB).  Since this might cause the doctest to fail
        on machines with "only" 2GB of RAM, we test ``max_n`` = 1, which
        has a more reasonable memory usage. ::

            sage: from sage.crypto.mq.sr import test_consistency
            sage: test_consistency(1)  # long time (80s on sage.math, 2011)
            True
        """

Ora facciamo i doctest della stessa directory in parallelo usando 4 thread::

    [jdemeyer@sage sage-6.0]$ ./sage -tp 4 src/sage/crypto/
    Running doctests with ID 2012-07-07-00-11-55-9b17765e.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 24 files using 4 threads.
    sage -t src/sage/crypto/boolean_function.pyx
        [252 tests, 3.8 s]
    sage -t src/sage/crypto/block_cipher/miniaes.py
        [429 tests, 1.1 s]
    sage -t src/sage/crypto/mq/sr.py
        [432 tests, 5.7 s]
    sage -t src/sage/crypto/mq/sbox.py
        [123 tests, 0.8 s]
    sage -t src/sage/crypto/block_cipher/sdes.py
        [289 tests, 0.6 s]
    sage -t src/sage/crypto/classical_cipher.py
        [123 tests, 0.4 s]
    sage -t src/sage/crypto/stream_cipher.py
        [113 tests, 0.1 s]
    sage -t src/sage/crypto/public_key/blum_goldwasser.py
        [134 tests, 0.1 s]
    sage -t src/sage/crypto/lfsr.py
        [30 tests, 0.1 s]
    sage -t src/sage/crypto/util.py
        [121 tests, 0.1 s]
    sage -t src/sage/crypto/cryptosystem.py
        [79 tests, 0.0 s]
    sage -t src/sage/crypto/stream.py
        [12 tests, 0.0 s]
    sage -t src/sage/crypto/mq/mpolynomialsystemgenerator.py
        [40 tests, 0.0 s]
    sage -t src/sage/crypto/cipher.py
        [3 tests, 0.0 s]
    sage -t src/sage/crypto/lattice.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/block_cipher/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/public_key/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/public_key/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/mq/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/block_cipher/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/crypto/mq/mpolynomialsystem.py
        [318 tests, 8.4 s]
    sage -t src/sage/crypto/classical.py
        [717 tests, 10.4 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 12.9 seconds
        cpu time: 30.5 seconds
        cumulative wall time: 31.7 seconds
    [jdemeyer@sage sage-6.0]$ ./sage -tp 4 --long src/sage/crypto/
    Running doctests with ID 2012-07-07-00-13-04-d71f3cd4.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 24 files using 4 threads.
    sage -t --long src/sage/crypto/boolean_function.pyx
        [252 tests, 3.7 s]
    sage -t --long src/sage/crypto/block_cipher/miniaes.py
        [429 tests, 1.0 s]
    sage -t --long src/sage/crypto/mq/sbox.py
        [123 tests, 0.8 s]
    sage -t --long src/sage/crypto/block_cipher/sdes.py
        [289 tests, 0.6 s]
    sage -t --long src/sage/crypto/classical_cipher.py
        [123 tests, 0.4 s]
    sage -t --long src/sage/crypto/util.py
        [121 tests, 0.1 s]
    sage -t --long src/sage/crypto/stream_cipher.py
        [113 tests, 0.1 s]
    sage -t --long src/sage/crypto/public_key/blum_goldwasser.py
        [134 tests, 0.1 s]
    sage -t --long src/sage/crypto/lfsr.py
        [30 tests, 0.0 s]
    sage -t --long src/sage/crypto/cryptosystem.py
        [79 tests, 0.0 s]
    sage -t --long src/sage/crypto/stream.py
        [12 tests, 0.0 s]
    sage -t --long src/sage/crypto/mq/mpolynomialsystemgenerator.py
        [40 tests, 0.0 s]
    sage -t --long src/sage/crypto/cipher.py
        [3 tests, 0.0 s]
    sage -t --long src/sage/crypto/lattice.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/block_cipher/all.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/public_key/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/mq/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/all.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/block_cipher/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/__init__.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/public_key/all.py
        [0 tests, 0.0 s]
    sage -t --long src/sage/crypto/mq/mpolynomialsystem.py
        [318 tests, 9.0 s]
    sage -t --long src/sage/crypto/classical.py
        [717 tests, 10.5 s]
    sage -t --long src/sage/crypto/mq/sr.py
        [434 tests, 88.0 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 90.4 seconds
        cpu time: 113.4 seconds
        cumulative wall time: 114.5 seconds

Incrementando il numero di thread diminuisce il tempo totale di test.


.. _section-parallel-test-whole-library:

Fare test paralleli dell'intera libreria di Sage
================================================

La libreria principale di Sage si trova nella directory ``SAGE_ROOT/src/``. 
Possiamo usare la sintassi descritta sopra per fare i doctest della libreria 
principale usando pi\`u di un thread. Quando si fa gestione della release o 
fi fa una patch della libreria principale di Sage, un manager di release farebbe 
i test in parallelo della libreria usando 10 thread con il seguente comando::

    [jdemeyer@sage sage-6.0]$ ./sage -tp 10 --long src/

Un altro modo \`e lanciare ``make ptestlong``, che compila Sage (se necessario),
compila la documentazione di Sage (se necessario), e poi esegue in parallelo i 
doctest. Questo determina il numero di thread da usare leggendo la variabile 
d'ambiente :envvar:`MAKE`: se \`e impostata a ``make -j12``, allora utilizzer\`a 
12 thread.  Se :envvar:`MAKE` non \`e impostata, di default utilizzer\`a il 
numero di core della CPU (come determinato dalla funzione Python 
``multiprocessing.cpu_count()``) con un minimo di 2 ed un massimo di 8.

In ogni caso tester\`a la libreria Sage con pi\`u di un thread::

    [jdemeyer@sage sage-6.0]$ make ptestlong

Uno qualunque dei seguenti comandi fare anche i doctest della libreria Sage o di 
uno dei suoi clone::

    make test
    make check
    make testlong
    make ptest
    make ptestlong

Le differenze sono:

* ``make test`` e ``make check`` --- Questi 2 comandi eseguono lo stesso
  insieme di test. Prima \`e testata la documentazione standard di Sage, 
  cio\`e la documentazione che risiede in

  * ``SAGE_ROOT/src/doc/common``
  * ``SAGE_ROOT/src/doc/en``
  * ``SAGE_ROOT/src/doc/fr``

  Poi i comandi fanno i doctest della libreria Sage. Per maggiori dettagli su 
  questo comandi, vedi il file ``SAGE_ROOT/Makefile``.

* ``make testlong`` --- Questo comando fa i doctests della documentazione 
  standard:

  * ``SAGE_ROOT/src/doc/common``
  * ``SAGE_ROOT/src/doc/en``
  * ``SAGE_ROOT/src/doc/fr``

  e poi della libreria Sage. I doctest sono eseguiti con l'argomento opzionale 
  ``--long``. Vedi il file ``SAGE_ROOT/Makefile`` per maggiori dettagli.

* ``make ptest`` --- Simile ai comandi ``make test`` e ``make
  check``. Comunque i doctest sono eseguiti con il numero di thread come 
  descritto sopra per ``make ptestlong``.

* ``make ptestlong`` --- Simile al comando ``make ptest``, ma con l'uso 
  dell'argomento opzionale ``--long`` per i doctest.


Oltre la libreria Sage
======================

I doctest funzionano bene anche per i file che non appartengono alla 
libreria Sage. Ad esempio, supponiamo di avere uno script Python detto 
``my_python_script.py``::

    [mvngu@sage build]$ cat my_python_script.py
    from sage.all_cmdline import *   # import sage library

    def square(n):
        """
        Return the square of n.

        EXAMPLES::

            sage: square(2)
            4
        """
        return n**2

Ne possiamo fare il doctest proprio come con i file della libreria Sage::

    [mvngu@sage sage-6.0]$ ./sage -t my_python_script.py
    Running doctests with ID 2012-07-07-00-17-56-d056f7c0.
    Doctesting 1 file.
    sage -t my_python_script.py
        [1 test, 0.0 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 2.2 seconds
        cpu time: 0.0 seconds
        cumulative wall time: 0.0 seconds

I doctest si possono anche fare sugli script Sage. Supponiamo di avere uno 
script Sage detto ``my_sage_script.sage`` con il contenuto seguente::

    [mvngu@sage sage-6.0]$ cat my_sage_script.sage
    def cube(n):
        r"""
        Return the cube of n.

        EXAMPLES::

            sage: cube(2)
            8
        """
        return n**3

Allora ne possiamo fare il doctest proprio come dei file Python::

    [mvngu@sage build]$ sage-6.0/sage -t my_sage_script.sage
    Running doctests with ID 2012-07-07-00-20-06-82ee728c.
    Doctesting 1 file.
    sage -t my_sage_script.sage
        [1 test, 0.0 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 2.5 seconds
        cpu time: 0.0 seconds
        cumulative wall time: 0.0 seconds

In alternativa, possiamo effettuare un preparse per convertirlo in uno 
script Python, e fare il doctest di questo::

    [mvngu@sage build]$ sage-6.0/sage --preparse my_sage_script.sage
    [mvngu@sage build]$ cat my_sage_script.sage.py
    # This file was *autogenerated* from the file my_sage_script.sage.
    from sage.all_cmdline import *   # import sage library
    _sage_const_3 = Integer(3)
    def cube(n):
        r"""
        Return the cube of n.

        EXAMPLES::

            sage: cube(2)
            8
        """
        return n**_sage_const_3
    [mvngu@sage build]$ sage-6.0/sage -t my_sage_script.sage.py
    Running doctests with ID 2012-07-07-00-26-46-2bb00911.
    Doctesting 1 file.
    sage -t my_sage_script.sage.py
        [2 tests, 0.0 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 2.3 seconds
        cpu time: 0.0 seconds
        cumulative wall time: 0.0 seconds

Fare i doctest da dentro Sage
=============================

Puoi eseguire i doctest da dentro Sage, cosa che pu\`o essere utile poich\`e 
non devi aspettare che Sage parta.  Usa la funzione ``run_doctests`` nel spazio 
dei nomi globale, passandogli una string o un modulo::

    sage: run_doctests(sage.coding.sd_codes)
    Doctesting /Users/roed/sage/sage-5.3/src/sage/coding/sd_codes.py
    Running doctests with ID 2012-07-07-04-32-36-81f3853b.
    Doctesting 1 file.
    sage -t /Users/roed/sage/sage-5.3/src/sage/coding/sd_codes.py
        [18 tests, 0.3 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 0.4 seconds
        cpu time: 0.2 seconds
        cumulative wall time: 0.3 seconds

.. _section-options:

Argomenti opzionali
===================

Eseguire test lunghi
--------------------

Idealmente, i doctest non dovrebbero richiedere un significativo intervallo 
di tempo. Se realmente hai bisogno di eseguire dei test che durano un po' di 
pi\`u (in sostanza pi\`u di qualche secondo) allora marcali come::

    sage: my_long_test()  # long time

Anche cos\`i i doctest lunghi si dovrebbero concludere idealmente in 5 
secondi o meno. Sappiamo che tu (l'autore) vuoi mettere in mostra le capacit\`a 
del tuo codice, ma questo non \`e posto per farlo. I test troppo lunghi prima o 
poi danneggeranno la nostra capacit\`a di eseguire la suite di test. Davvero,
i doctest devono essere veloci il pi\`u possibile pur coprendo completamente il 
codice.

Usa il flag ``--long`` per eseguire dei doctest che sono stati marcati col 
commento ``# long time``. Questi test sono normalmente saltati per ridurre il 
tempo speso ad eseguire i test::

    [roed@sage sage-6.0]$ sage -t src/sage/rings/tests.py
    Running doctests with ID 2012-06-21-16-00-13-40835825.
    Doctesting 1 file.
    sage -t tests.py
        [18 tests, 1.1 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 2.9 seconds
        cpu time: 0.9 seconds
        cumulative wall time: 1.1 seconds

Per eseguire anche i test lunghi, fa quanto segue::

    [roed@sage sage-6.0]$ sage -t --long src/sage/rings/tests.py
    Running doctests with ID 2012-06-21-16-02-05-d13a9a24.
    Doctesting 1 file.
    sage -t tests.py
        [20 tests, 34.7 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 46.5 seconds
        cpu time: 25.2 seconds
        cumulative wall time: 34.7 seconds

Per trovare dei test che durano pi\`u del tempo permesso usa il flag 
``--warn-long``. Senza opzioni far\`a s\`i che i test mostrino a video 
un warning se impiegano pi\`u di 1.0 secondo. Nota che questo \`e appunto un 
warning, non un errore::

    [roed@sage sage-6.0]$ sage -t --warn-long src/sage/rings/factorint.pyx
    Running doctests with ID 2012-07-14-03-27-03-2c952ac1.
    Doctesting 1 file.
    sage -t --warn-long src/sage/rings/factorint.pyx
    **********************************************************************
    File "src/sage/rings/factorint.pyx", line 125, in sage.rings.factorint.base_exponent
    Failed example:
        base_exponent(-4)
    Test ran for 4.09 s
    **********************************************************************
    File "src/sage/rings/factorint.pyx", line 153, in sage.rings.factorint.factor_aurifeuillian
    Failed example:
        fa(2^6+1)
    Test ran for 2.22 s
    **********************************************************************
    File "src/sage/rings/factorint.pyx", line 155, in sage.rings.factorint.factor_aurifeuillian
    Failed example:
        fa(2^58+1)
    Test ran for 2.22 s
    **********************************************************************
    File "src/sage/rings/factorint.pyx", line 163, in sage.rings.factorint.factor_aurifeuillian
    Failed example:
        fa(2^4+1)
    Test ran for 2.25 s
    **********************************************************************
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    Total time for all tests: 16.1 seconds
        cpu time: 9.7 seconds
        cumulative wall time: 10.9 seconds

Puoi anche passare esplicitamente una quantit\`a di tempo::

    [roed@sage sage-6.0]$ sage -t --long --warn-long 2.0 src/sage/rings/tests.py
    Running doctests with ID 2012-07-14-03-30-13-c9164c9d.
    Doctesting 1 file.
    sage -t --long --warn-long 2.0 tests.py
    **********************************************************************
    File "tests.py", line 240, in sage.rings.tests.test_random_elements
    Failed example:
        sage.rings.tests.test_random_elements(trials=1000)  # long time (5 seconds)
    Test ran for 13.36 s
    **********************************************************************
    File "tests.py", line 283, in sage.rings.tests.test_random_arith
    Failed example:
        sage.rings.tests.test_random_arith(trials=1000)   # long time (5 seconds?)
    Test ran for 12.42 s
    **********************************************************************
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    Total time for all tests: 27.6 seconds
        cpu time: 24.8 seconds
        cumulative wall time: 26.3 seconds

Infine, puoi disabilitare i warnings sui test lunghi con ``--warn-long 0``.


.. _section-optional-doctest-flag:

Eseguire test opzionali
-----------------------

Puoi eseguire test che richiedono pacchetti opzionali utilizzando il flag 
``--optional``.  Ovviamente deve aver installato i necessari pacchetti 
opzionali perch\`e tali tests abbiano successo. Vedi 
http://www.sagemath.org/packages/optional/ per come fare il download di 
pacchetti opzionali.

Di default, Sage esegue solo i doctest che non sono marcati con il tag 
``optional``. Questo equivale ad eseguire ::

    [roed@sage sage-6.0]$ sage -t --optional=sage src/sage/rings/real_mpfr.pyx
    Running doctests with ID 2012-06-21-16-18-30-a368a200.
    Doctesting 1 file.
    sage -t src/sage/rings/real_mpfr.pyx
        [819 tests, 7.0 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 8.4 seconds
        cpu time: 4.1 seconds
        cumulative wall time: 7.0 seconds

Se vuoi anche eseguire i test che richiedono magma, puoi fare come segue::

    [roed@sage sage-6.0]$ sage -t --optional=sage,magma src/sage/rings/real_mpfr.pyx
    Running doctests with ID 2012-06-21-16-18-30-a00a7319
    Doctesting 1 file.
    sage -t src/sage/rings/real_mpfr.pyx
        [823 tests, 8.4 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 9.6 seconds
        cpu time: 4.0 seconds
        cumulative wall time: 8.4 seconds

Per poter eseguire solo i test marcati come richiedenti magma, ometti ``sage``::

    [roed@sage sage-6.0]$ sage -t --optional=magma src/sage/rings/real_mpfr.pyx
    Running doctests with ID 2012-06-21-16-18-33-a2bc1fdf
    Doctesting 1 file.
    sage -t src/sage/rings/real_mpfr.pyx
        [4 tests, 2.0 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 3.2 seconds
        cpu time: 0.1 seconds
        cumulative wall time: 2.0 seconds

Per eseguire tutti i test, indiffentemente al fatto che siano marcati opzionali o 
no, passa ``all`` come ``optional`` tag::

    [roed@sage sage-6.0]$ sage -t --optional=all src/sage/rings/real_mpfr.pyx
    Running doctests with ID 2012-06-21-16-31-18-8c097f55
    Doctesting 1 file.
    sage -t src/sage/rings/real_mpfr.pyx
        [865 tests, 11.2 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 12.8 seconds
        cpu time: 4.7 seconds
        cumulative wall time: 11.2 seconds

Eseguire i test in parallelo
----------------------------

Se stai facendo il test di parecchi file, puoi velocizzare molto le cose 
utilizzando pi\`u di un thread.  Per eseguire i doctest in parallelo usa il flag 
``--nthreads`` (``-p`` \`e una versione abbreviata). Per passare il numero di thread 
useresti (di default Sage ne usa solo 1)::

    [roed@sage sage-6.0]$ sage -tp 2 src/sage/doctest/
    Running doctests with ID 2012-06-22-19-09-25-a3afdb8c.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 8 files using 2 threads.
    sage -t src/sage/doctest/control.py
        [114 tests, 4.6 s]
    sage -t src/sage/doctest/util.py
        [114 tests, 0.6 s]
    sage -t src/sage/doctest/parsing.py
        [187 tests, 0.5 s]
    sage -t src/sage/doctest/sources.py
        [128 tests, 0.1 s]
    sage -t src/sage/doctest/reporting.py
        [53 tests, 0.1 s]
    sage -t src/sage/doctest/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/doctest/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/doctest/forker.py
        [322 tests, 15.5 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 17.0 seconds
        cpu time: 4.2 seconds
        cumulative wall time: 21.5 seconds

Fare i doctesting di tutto Sage
-------------------------------

Per fare i doctest di tutta la libreria Sage usa il flag ``--all`` (``-a`` come 
abbreviazione). Oltre a testare il codice nei file Python e Cython di Sage, 
questo comando eseguir\`a i test definiti nella documentazione si Sage cos\`i come 
i test dei notebook Sage::

    [roed@sage sage-6.0]$ sage -t -a
    Running doctests with ID 2012-06-22-19-10-27-e26fce6d.
    Doctesting entire Sage library.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 2020 files.
    sage -t /Users/roed/sage/sage-5.3/src/sage/plot/plot.py
        [304 tests, 69.0 s]
    ...

Se vuoi solo eseguire i test dei notebook, usa invece il flag ``--sagenb``.


Strumenti di debug
------------------

A volte i doctest falliscono (\`e per questo che li eseguiamo, dopotutto). Ci sono 
vari flag che possono essere utili quando qualcosa va storto. Se un doctest d\`a 
un errore Python, allora di solito i test continuano dopo averlo segnalato. Se usi 
il flag ``--debug`` (abbreviabile``-d``) allora ti si aprir\`a un debugger Python 
interattivo ogni volta che viene sollevata un'eccezione Python. Ad esempio, ho modificato 
:mod:`sage.schemes.elliptic_curves.constructor` per produrre un errore::

    [roed@sage sage-6.0]$ sage -t --debug src/sage/schemes/elliptic_curves/constructor.py
    Running doctests with ID 2012-06-23-12-09-04-b6352629.
    Doctesting 1 file.
    **********************************************************************
    File "sage.schemes.elliptic_curves.constructor", line 4, in sage.schemes.elliptic_curves.constructor
    Failed example:
        EllipticCurve([0,0])
    Exception raised:
        Traceback (most recent call last):
          File "/Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/doctest/forker.py", line 573, in _run
            self.execute(example, compiled, test.globs)
          File "/Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/doctest/forker.py", line 835, in execute
            exec compiled in globs
          File "<doctest sage.schemes.elliptic_curves.constructor[0]>", line 1, in <module>
            EllipticCurve([Integer(0),Integer(0)])
          File "/Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/schemes/elliptic_curves/constructor.py", line 346, in EllipticCurve
            return ell_rational_field.EllipticCurve_rational_field(x, y)
          File "/Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/schemes/elliptic_curves/ell_rational_field.py", line 216, in __init__
            EllipticCurve_number_field.__init__(self, Q, ainvs)
          File "/Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/schemes/elliptic_curves/ell_number_field.py", line 159, in __init__
            EllipticCurve_field.__init__(self, [field(x) for x in ainvs])
          File "/Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/schemes/elliptic_curves/ell_generic.py", line 156, in __init__
            "Invariants %s define a singular curve."%ainvs
        ArithmeticError: Invariants [0, 0, 0, 0, 0] define a singular curve.
    > /Users/roed/sage/sage-5.3/local/lib/python2.7/site-packages/sage/schemes/elliptic_curves/ell_generic.py(156)__init__()
    -> "Invariants %s define a singular curve."%ainvs
    (Pdb) l
    151                 if len(ainvs) == 2:
    152                     ainvs = [K(0),K(0),K(0)] + ainvs
    153                 self.__ainvs = tuple(ainvs)
    154                 if self.discriminant() == 0:
    155                     raise ArithmeticError, \
    156  ->                       "Invariants %s define a singular curve."%ainvs
    157                 PP = projective_space.ProjectiveSpace(2, K, names='xyz');
    158                 x, y, z = PP.coordinate_ring().gens()
    159                 a1, a2, a3, a4, a6 = ainvs
    160                 f = y**2*z + (a1*x + a3*z)*y*z \
    161                     - (x**3 + a2*x**2*z + a4*x*z**2 + a6*z**3)
    (Pdb) p ainvs
    [0, 0, 0, 0, 0]
    (Pdb) quit
    **********************************************************************
    1 items had failures:
       1 of   1 in sage.schemes.elliptic_curves.constructor
    ***Test Failed*** 1 failures.
    sage -t src/sage/schemes/elliptic_curves/constructor.py
        [64 tests, 89.2 s]
    ------------------------------------------------------------------------
    sage -t src/sage/schemes/elliptic_curves/constructor.py # 1 doctest failed
    ------------------------------------------------------------------------
    Total time for all tests: 90.4 seconds
        cpu time: 4.5 seconds
        cumulative wall time: 89.2 seconds

A volte un errore pu\`o essere cos\`i grave da causare un segfault o il blocco di 
Sage. In tali situazioni hai un certo numero di opzioni. Il doctest framework 
mostrer\`a a video l'output fin l\`i, cos\`i che tu possa almeno sapere la causa del 
problema (se vuoi che quest'output appaia in tempo reale usa il flag ``--verbose``). 
Per eseguire i doctests sotto il controllo di gdb, usa il flag ``--gdb``::

    [roed@sage sage-6.0]$ sage -t --gdb src/sage/schemes/elliptic_curves/constructor.py
    gdb -x /home/roed/sage-6.0.b5/local/bin/sage-gdb-commands --args python /home/roed/sage-6.0.b5/local/bin/sage-runtests --serial --nthreads 1 --timeout 1048576 --optional sage --stats_path /home/roed/.sage/timings2.json src/sage/schemes/elliptic_curves/constructor.py
    GNU gdb 6.8-debian
    Copyright (C) 2008 Free Software Foundation, Inc.
    License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
    This is free software: you are free to change and redistribute it.
    There is NO WARRANTY, to the extent permitted by law.  Type "show copying"
    and "show warranty" for details.
    This GDB was configured as "x86_64-linux-gnu"...
    [Thread debugging using libthread_db enabled]
    [New Thread 0x7f10f85566e0 (LWP 6534)]
    Running doctests with ID 2012-07-07-00-43-36-b1b735e7.
    Doctesting 1 file.
    sage -t src/sage/schemes/elliptic_curves/constructor.py
        [67 tests, 5.8 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 15.7 seconds
        cpu time: 4.4 seconds
        cumulative wall time: 5.8 seconds

    Program exited normally.
    (gdb) quit


Sage include anche valgrind, e puoi eseguire i doctests sotto vari strumenti di 
valgrind per tracciare problemi di memoria: i flag utili a ci\`o sono 
``--valgrind`` (o ``--memcheck``), ``--massif``, ``--cachegrind`` e 
``--omega``. Vedi http://wiki.sagemath.org/ValgrindingSage per maggiori dettagli.

Una volta che hai finito di sistemare eventuali problem che fossero stati evidenziati 
dai doctest, puoi rieseguire solo quei file che hanno fallito l'ultimo test usando il 
flag ``--failed`` (abbreviato ``-f``)::

    [roed@sage sage-6.0]$ sage -t -fa
    Running doctests with ID 2012-07-07-00-45-35-d8b5a408.
    Doctesting entire Sage library.
    Only doctesting files that failed last test.
    No files to doctest


Opzioni varie
-------------

Ci sono varie opzioni che modificano il comportamento del codice dei doctest di Sage.

Mostra solo primo errore
^^^^^^^^^^^^^^^^^^^^^^^^

Il primo errore in un file spesso ne causa degli altri a cascata, poich\`e i
NameErrors nascono da variabili non definite ed i tests falliscono a causa di vecchi 
valori delle variabili in uso. Per vedere solo il primo errore in ciascun blocco di 
doctest block usa il flag ``--initial`` (abbreviato ``-i``).

Mostra i test opzionali saltati
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Per mostrare un riassunto alla fine di ogni file con il numero di test opzionali 
saltati, usa il flag ``--show-skipped``::

   [roed@sage sage-6.0]$ sage -t --show-skipped src/sage/rings/finite_rings/integer_mod.pyx
   Running doctests with ID 2013-03-14-15-32-05-8136f5e3.
   Doctesting 1 file.
   sage -t sage/rings/finite_rings/integer_mod.pyx
       2 axiom tests not run
       1 cunningham test not run
       2 fricas tests not run
       1 long test not run
       3 magma tests not run
       [440 tests, 4.0 s]
   ----------------------------------------------------------------------
   All tests passed!
   ----------------------------------------------------------------------
   Total time for all tests: 4.3 seconds
       cpu time: 2.4 seconds
       cumulative wall time: 4.0 seconds

Eseguire i test con iterazioni
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A volte i test falliscono in modo intermittente. Ci sono 2 opzioni che ti 
permettono di eseguire i test ripetutamente nel tentativo di cercare Heisenbugs.
Il flag ``--global-iterations`` prende un intero ed esegue l'intero insieme di 
test quel numero di volte di seguito::

    [roed@sage sage-6.0]$ sage -t --global-iterations 2 src/sage/sandpiles
    Running doctests with ID 2012-07-07-00-59-28-e7048ad9.
    Doctesting 3 files (2 global iterations).
    sage -t src/sage/sandpiles/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/sandpiles/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/sandpiles/sandpile.py
        [711 tests, 14.7 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 17.6 seconds
        cpu time: 13.2 seconds
        cumulative wall time: 14.7 seconds
    sage -t src/sage/sandpiles/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/sandpiles/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/sandpiles/sandpile.py
        [711 tests, 13.8 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 14.3 seconds
        cpu time: 26.4 seconds
        cumulative wall time: 28.5 seconds

Puoi anche iterare in ordine differente: il flag ``--file-iterations`` esegue 
i test in ciascun file ``N`` volte prima di procedere::

    [roed@sage sage-6.0]$ sage -t --file-iterations 2 src/sage/sandpiles
    Running doctests with ID 2012-07-07-01-01-43-8f954206.
    Doctesting 3 files (2 file iterations).
    sage -t src/sage/sandpiles/__init__.py
        [0 tests, 0.0 s]
    sage -t src/sage/sandpiles/all.py
        [0 tests, 0.0 s]
    sage -t src/sage/sandpiles/sandpile.py
        [1422 tests, 13.3 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 29.6 seconds
        cpu time: 12.7 seconds
        cumulative wall time: 13.3 seconds


Nota che i risultati riportati sono il tempo medio impiegato da tutti i test in 
quel file per finire.  Se capita un errore in un file, allora l'errore \`e 
segnalato ed i test procedono con il file seguente.

Usare un timeout differente
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Su una macchina lenta il timeout di default di 5 minuti pu\`o non essere abbastanza 
per i file pi\`u lenti.  Usa il flag ``--timeout`` (abbreviato ``-T``) per 
impostarlo a qualcos'altro::

    [roed@sage sage-6.0]$ sage -tp 2 --all --timeout 1
    Running doctests with ID 2012-07-07-01-09-37-deb1ab83.
    Doctesting entire Sage library.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 2067 files using 2 threads.
    sage -t src/sage/schemes/elliptic_curves/ell_rational_field.py
        Timed out!
    ...

Usare path assoluti
^^^^^^^^^^^^^^^^^^^

Di default i nomi dei file sono mostrati usando path relativi. Per usare 
path assoluti al posto, passa il flag ``--abspath``::

    [roed@sage sage-6.0]$ sage -t --abspath src/sage/doctest/control.py
    Running doctests with ID 2012-07-07-01-13-03-a023e212.
    Doctesting 1 file.
    sage -t /home/roed/sage-6.0/src/sage/doctest/control.py
        [133 tests, 4.7 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 7.1 seconds
        cpu time: 0.2 seconds
        cumulative wall time: 4.7 seconds


Testare i file modificati
^^^^^^^^^^^^^^^^^^^^^^^^^

Se stai lavorando su qualche file nella libreria Sage pu\`o essere conveniente 
fare i test solo dei file che sono stati modificati. Per farlo usa il flag 
``--new``, che testa i file che sono stati modificati o aggiunti dopo l'ultimo 
commit::

    [roed@sage sage-6.0]$ sage -t --new
    Running doctests with ID 2012-07-07-01-15-52-645620ee.
    Doctesting files changed since last git commit.
    Doctesting 1 file.
    sage -t src/sage/doctest/control.py
        [133 tests, 3.7 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 3.8 seconds
        cpu time: 0.1 seconds
        cumulative wall time: 3.7 seconds


Eseguire i test in ordine casuale
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Di default, i test sono eseguiti nell'ordine in cui appaiono nel file. 
Per eseguire i test in ordinecasuale (cosa che pu\`o rivelare bachi sottili),
usa il flag ``--randorder`` e passa un "random seed"::

    [roed@sage sage-6.0]$ sage -t --new --randorder 127
    Running doctests with ID 2012-07-07-01-19-06-97c8484e.
    Doctesting files changed since last git commit.
    Doctesting 1 file.
    sage -t src/sage/doctest/control.py
        [133 tests, 3.6 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 3.7 seconds
        cpu time: 0.2 seconds
        cumulative wall time: 3.6 seconds

Nota che anche con quest'opzione, i test in ciascun blocco di doctest sono ancora 
eseguiti in ordine.

Testare file esterni
^^^^^^^^^^^^^^^^^^^^

Quando si testa un file che non \`e parte della libreria Sage, il codice di test 
carica le quantit\`a globali da quel file nello spazio dei nomi (namespace) prima 
di eseguire i test. Per modellare il comportamento utilizzato invece sulla libreria 
Sage (dove gli import devono essere specificati esplicitamente), usa il flag 
``--force-lib``.

File ausiliari
^^^^^^^^^^^^^^

Per specificare un file di log (piuttosto di utilizzare quello di default che \`e creato 
da ``sage -t --all``), usa il flag ``--logfile``::

    [roed@sage sage-6.0]$ sage -t --logfile test1.log src/sage/doctest/control.py
    Running doctests with ID 2012-07-07-01-25-49-e7c0e52d.
    Doctesting 1 file.
    sage -t src/sage/doctest/control.py
        [133 tests, 4.3 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 6.7 seconds
        cpu time: 0.1 seconds
        cumulative wall time: 4.3 seconds
    [roed@sage sage-6.0]$ cat test1.log
    Running doctests with ID 2012-07-07-01-25-49-e7c0e52d.
    Doctesting 1 file.
    sage -t src/sage/doctest/control.py
        [133 tests, 4.3 s]
    ------------------------------------------------------------------------
    All tests passed!
    ------------------------------------------------------------------------
    Total time for all tests: 6.7 seconds
        cpu time: 0.1 seconds
        cumulative wall time: 4.3 seconds


Per dare un file json che contenga le tempistiche per ciascun file, usa il flag 
``--stats_path``. Queste statistiche sono usate nell'ordinare i file cos\`i che 
i test pi\`u lenti sono eseguiti prima (e cos\`i processi multipli sono 
utilizzati in modo pi\`u efficiente)::

    [roed@sage sage-6.0]$ sage -tp 2 --stats-path ~/.sage/timings2.json --all
    Running doctests with ID 2012-07-07-01-28-34-2df4251d.
    Doctesting entire Sage library.
    Sorting sources by runtime so that slower doctests are run first....
    Doctesting 2067 files using 2 threads.
    ...
