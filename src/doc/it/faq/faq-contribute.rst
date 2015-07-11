.. _chapter-faq-contribute:

========================
FAQ: Contribuire a Sage.
========================


Come posso iniziare a contribuire a Sage ?
""""""""""""""""""""""""""""""""""""""""""

C'e' una guida rapida per chiunque voglia contribuire a Sage. E' indirizzata soprattutto a chi e' nuovo alla programmazione. Innanzitutto il linguaggio principale per la programmazione di Sage e' `Python <http://www.python.org>`_.
Alcune parti di Sage possono essere scritte in altri linguaggi, specialmente le componenti che fanno l'elaborazione numerica piu' impegnativa, ma la maggior parte delle funzionalita' sono realizzate in Python, incluso il codice di collegamento. Uno degli aspetti validi di Python, che Sage eredita, e' che e' piu' importante scrivere codice funzionante che non codice veloce. Il codice veloce e' prezioso, ma il codice chiaro e leggibile e' importante. Nella comunita' matematica i risultati inaccurati sono inaccettabili. La correttezza viene prima dell'ottimizzazione. Nel seguente articolo::

* D. Knuth. Structured Programming with go to Statements. *ACM Journal Computing Surveys*,

  6(4), 1974.

Don Knuth osserva che "dovremmo dimenticarci dell'efficienza di poco momento, diciamo per il 97% del tempo: l'ottimizzazione prematura sta alla radice di ogni male".

Se non conosci Python dovresti iniziare ad imparare il linguaggio. Un buon posto dove iniziare e' il `Python Official Tutorial <http://docs.python.org/tutorial>`_ e altra documentazione si trova a `Documentazione standard di Python <http://docs.python.org>`_. Un altro posto da guardare e' al link `Dive Into Python <http://www.diveintopython.net>`_ di Marc Pilgrim, che puo' essere veramente d'aiuto su temi specifici come lo sviluppo guidato dai test. Il libro `Building Skills in Python <http://homepage.mac.com/s_lott/books/python.html>`_ di Stephen F.Lott e' adatto a chiunque abbia gia' nozioni di programmazione.

Nel frattempo puoi guardare la tua copia del codice sorgente di Sage e familiarizzare con il sistema di controllo versione `git <http://git-scm.com>`_. Una volta che sei a tuo agio con Python, cosa rapida per gli elementi di base, puoi iniziare ad usare Sage. Se vuoi puoi provare ad imparare Python usando Sage, ma non e' consigliabile perche' e' utile sapere cio' che e' il puro Python e dove invece Sage fa le sue "magie". Ci sono molte cose che funzionano in Python ma non in Sage e viceversa.

Poi dai un'occhiata alla `Guida ufficiale di sviluppo di Sage <http://www.sagemath.org/doc/developer>`_. Come minimo il suo primo capitolo dev'essere letto da ogni sviluppatore Sage. Fa anche particolare attenzione alle `linee guida su Trac <http://www.sagemath.org/doc/developer/trac.html>`_. Puoi anche unirti alla mailing list `sage-devel <http://groups.google.com/group/sage-devel>`_ e gravitare sul canale IRC ``#sage-devel`` su `freenode <http://freenode.net>`_.

Il modo migliore per familiarizzare con il processo di sviluppo di Sage e' scegliere un ticket dal `server Trac <http://trac.sagemath.org>`_ ed esaminare i cambiamenti ivi proposti. Se vuoi implementare qualcosa e' buona pratica discutere le tue idee prima sulla mailing list ``sage-devel``, cosi' che altri sviluppatori abbiano la possibilita' di fare commenti sulle tue idee/proposte. Sono anche molto aperti a nuove idee, come tutti i matematici dovrebbero essere.


Non sono un programmatore. C'e' qualche altro modo in cui posso aiutare ?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Certo. Come in ogni progetto FOSS ci sono molti modi in cui puoi dare il tuo aiuto nella comunita' di Sage. E programmare e' solo uno dei modi. Se sai parlare, leggere e scrivere in un'altra lingua, ci sono molti modi in cui il tuo contributo puo' essere molto prezioso all'intera comunita' di Sage. Diciamo che conosci l'italiano. Allora puoi scrivere un tutorial per Sage in italiano, o aiutare a tradurre il tutorial ufficiale di Sage in italiano. Per i grafici o gli artisti c'e' la possibilita' di aiutare migliorando l'aspetto del sito di Sage. O puoi dare un giudizio artistico sull'interfaccia Notebook di Sage, suggerendo come puo' essere migliorata.

Molte persone amano scrivere tutorial tecnici. Una delle gioie di farlo e' che impari qualcosa di nuovo nel farlo. Allo stesso tempo comunichi delle conoscenze ai principianti, un'abilita' che e' utile anche in campi estranei alla stesura di documentazione tecnica. Un aspetto fondamentale della documentazione tecnica e' che espone un dato soggetto tecnico a dei principianti, pertanto il gergo tecnico dev'essere ridotto al minimo. Darrell Anderson ha scritto `alcuni suggerimenti sulla scrittura di documentazione tecnica <http://humanreadable.nfshost.com/howtos/technical_writing_tips.htm>`_, che ti suggeriamo caldamente di leggere. La suddetta e' una lista molto breve. Ci sono molti, molti piu' modi in cui puoi dare il tuo aiuto. Sentiti libero di inviare una email alla mailing list `sage-devel <http://groups.google.com/group/sage-devel>`_ per chiedere in quali modi potresti essere d'aiuto, o per suggerire un'idea sul progetto.


Dove posso trovare risorse su Python e Cython ?
"""""""""""""""""""""""""""""""""""""""""""""""

Ecco una lista incompleta di risorse su Python e Cython. Ulteriori risorse possono essere trovate cercando sul web.

**Risorse generali**

* `Cython <http://www.cython.org>`_
* `pep8 <http://pypi.python.org/pypi/pep8>`_
* `py2depgraph <http://www.tarind.com/depgraph.html>`_
* `pycallgraph <http://pycallgraph.slowchop.com>`_
* `PyChecker <http://pychecker.sourceforge.net>`_
* `PyFlakes <http://divmod.org/trac/wiki/DivmodPyflakes>`_
* `Pylint <http://www.logilab.org/project/pylint>`_
* `Python <http://www.python.org>`_ home page e la `Documentazione standard su Python <http://docs.python.org>`_
* `Snakefood <http://furius.ca/snakefood>`_
* `Sphinx <http://sphinx.pocoo.org>`_
* `XDot <http://code.google.com/p/jrfonseca/wiki/XDot>`_

**Tutorial e libri**

* `Building Skills in Python <http://homepage.mac.com/s_lott/books/python.html>`_ di Steven F. Lott
* `Cython Tutorial <http://conference.scipy.org/proceedings/SciPy2009/paper_1/>`_ di Stefan Behnel, Robert W. Bradshaw, e Dag Sverre Seljebotn
* `Dive into Python <http://www.diveintopython.net>`_ di Mark Pilgrim
* `Dive Into Python 3 <http://www.diveintopython3.net>`_ di Mark Pilgrim
* `Fast Numerical Computations with Cython <http://conference.scipy.org/proceedings/SciPy2009/paper_2/>`_ di Dag Sverre Seljebotn
* `How to Think Like a Computer Scientist <http://www.openbookproject.net/thinkCSpy>`_ di Jeffrey Elkner, Allen B. Downey, and Chris Meyers
* `Tutorial ufficiale di Python <http://docs.python.org/tutorial>`_

**Articoli e HOWTO**

* `decorator <http://pypi.python.org/pypi/decorator>`_
* `Functional Programming HOWTO <http://docs.python.org/howto/functional.html>`_ di A. M. Kuchling
* `Python Functional Programming for Mathematicians <http://wiki.sagemath.org/devel/FunctionalProgramming>`_ di Minh Van Nguyen
* `Regular Expression HOWTO <http://docs.python.org/howto/regex.html>`_ di A. M. Kuchling
* `reStructuredText <http://docutils.sourceforge.net/rst.html>`_
* `Static Code Analizers for Python <http://www.doughellmann.com/articles/pythonmagazine/completely-different/2008-03-linters/>`_ di Doug Hellmann


Ci sono delle convenzioni di scrittura del codice sorgente che devo seguire ?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Dovresti seguire le convenzioni standard di Python come documentato in `PEP 0008 <http://www.python.org/dev/peps/pep-0008>`_ e `PEP 0257 <http://www.python.org/dev/peps/pep-0257>`_.
Consulta anche la Guida dello Sviluppo Sage, specialmente il capitolo `Convenzioni per scrivere codice sorgente in Sage <http://www.sagemath.org/doc/developer/conventions.html>`_.


Ho inviato al server trac una correzione molte settimane fa. Perche' la state ignorando ?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Non stiamo cercando di ignorare la tua correzione. Le persone che lavorano su Sage lo fanno nel loro tempo libero. Con centinaia di ticket aperti aventi un impatto molto variabile sulla comunita' di Sage, le persone che ci lavorano devono dare il loro tempo e sforzo prioritariamente a quei ticket che li interessano. A volte potresti essere la sola persona che comprende la tua correzione. In tal caso, ti invitiamo a fare uno sforzo supplementare per rendere il piu' semplice possibile a chiunque l'esaminare la tua correzione. Ecco alcuni suggerimenti su come rendere la tua correzione facile da esaminare::

* Hai descritto in modo chiaro il problema che la tua correzione vuole risolvere ?

* Hai fornito ogni informazione di base rilevante al problema che la tua correzione
  vuole risolvere ? Tali informazioni includono link a risorse online ed ad articoli,
  libri, o altro materiale di riferimento.

* Hai descritto in modo chiaro come la tua correzione risolve il problema in oggetto ?

* Hai descritto chiaramente nella tua correzione come effettuare i test dei cambiamenti ?

* Hai elencato eventuali tickets da cui dipende la tua correzione ?

* Se vi sono piu' correzioni, hai indicato chiaramente l'ordine in cui devono essere applicate ?

* La tua correzione segue le `convenzioni importanti <http://www.sagemath.org/doc/developer/writing_code.html>`_
  indicate nella "Guida dello sviluppatore" ?

Se la tua correzione non ha possibilita' di essere aggiunta nell'albero dei sorgenti di Sage, non la ignoreremo ma semplicemente chiuderemo il ticket relativo con una spiegazione sul perche' non possiamo includerla.


Come e quando posso ricordardare alla comunita' di Sage una correzione a cui tengo?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Ti suggeriamo di fare uno sforzo ulteriore nel ricordare alla comunita' di Sage una correzione che vuoi venga inclusa nell'albero dei sorgenti di Sage. Potrebbe esserci un prossimo evento "bug squash sprint" o "Sage days" che e' in relazione alla tua correzione. Tieni d'occhio le mailing list relative e rispondi educatamente ad ogni scambio di email relativo, spiegando chiaramente perche' la tua correzione ha importanza. Tieni d'occhio il canale IRC ``#sage-devel``, avendo cura di rammentare la questione al momento giusto.


Ho scritto del codice sorgente e voglio venga incluso in Sage. Pero' dopo aver rinominato il mio file ``a.sage`` in ``a.py`` ho degli errori di sintassi. Devo riscrivere tutto il mio codice in Python anziche' in Sage?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

La risposta sostanzialmente e' si', ma riscrivere e' una parola grossa per cio' che bisogna realmente fare. C'e' ben poco da fare dal momento che Sage per lo piu' segue la sintassi di Python. Le 2 maggiori differenze sono la gestione degli interi (vedi anche il link `afterword`_ per maggiori informazioni sul preparser di Sage) e la necessita' di importare quello che ti serve.

•  la seconda cosa importante da tenere presente e' la necessita' di importare tutto cio' di cui hai bisogno. Nel dettaglio, ogni volta che usi una funzione Sage la devi prima importare all'inizio del file. Ad esempio, se hai bisogno di PolynomialRing, dovrai scrivere::

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing



- **Gestione degli interi:** dei fare i seguenti cambiamenti:

  - Notazione per l'elevamento a potenza: In Python ``**`` significa elevamento a potenza e ``^`` significa “xor”.
  - Se devi restituire un intero all'utente, scrivi ``return Integer(1)`` invece di ``return 1``. In Python, 1 e' un intero Python (``int``), e ``Integer(1)`` e' un intero Sage/Gmp. Inoltre gli ``Integer`` sono molto piu' potenti degli ``int``; ad esempio hanno collegata ad essi l'informazione di primalita' e la fattorizzazione.
  - Dovresti anche notare che ``2 / 3`` non significa piu' ``Integer(2) / Integer(3)`` che restituisce ``2/3``, ma invece ``int(2) / int(3)``, e pertanto restituisce ``0`` poiche' la divisione e' intera e trascura il resto. Se stai lavorando con i tipi ``Integer`` ma in realta' hai bisogno di eseguire una divisione intera puoi usare ``Integer(2) // Integer(3)``.

- **Note sull'importazione:** la seconda cosa importante da tenere presente e' la necessita' di importare tutto cio' di cui hai bisogno. Nel dettaglio, ogni volta che usi una funzione Sage la devi prima importare all'inizio del file. Ad esempio, se hai bisogno di ``PolynomialRing``, dovrai scrivere::

      from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

  Puoi chiedere a Sage dove si trova ``PolynomialRing`` usando::

      sage: PolynomialRing.__module__
      'sage.rings.polynomial.polynomial_ring_constructor'

  Questo corrisponde anche al percorso, che inizia dopo ``site-packages``, restituito da Sage quando richiami l'help su ``PolynomialRing``. Ad esempio se scrivi ``PolynomialRing?`` otterrai::

      Type:    function
      [...]
      File:    /home/florent/src/Sage/sage/local/lib/python2.6/site-packages/sage/rings/
               polynomial/polynomial_ring_constructor.py
      [...]


.. _afterword: http://www.sagemath.org/doc/tutorial/afterword.html

