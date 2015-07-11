.. -*- coding: utf-8 -*-
.. _chapter-faq-general:

============
FAQ: Generalita'
============


Perche' esiste questo progetto ?
""""""""""""""""""""""""""""""""

La missione fissata per Sage e' di essere un'alternativa open-source utilizzabile invece di Magma, Maple, Mathematica a Matlab. I predecessori di Sage, noti come HECKE e Manin, furono creati perche' William Stein ebbe bisogno di scriverli come parte della sua ricerca sulla Teoria dei Numeri. Iniziato da William nel 2005 quando era all'universita' di Harvard, Sage combina alcuni fra i miglori software open-source per la matematica, inglobandoli in un'unita' accessibile da un'intefaccia comune. Molti ricercatori in Teoria dei Numeri, incluso lo stesso William, si avvalgono di quest'interfaccia comune per usare ed estendere le funzionalita' dei pacchetti inglobati relativi alla Teoria dei Numeri. Tali pacchetti software includono:
`Givaro <http://ljk.imag.fr/CASYS/LOGICIELS/givaro>`_,
`MPIR <http://www.mpir.org>`_,
`NTL <http://www.shoup.net/ntl>`_,
`Pari/GP <http://pari.math.u-bordeaux.fr>`_,
e molti altri troppo numerosi per essere elencati qui. Studenti, insegnanti, professori universitari, ricercatori di tutto il mondo usano Sage perche' vogliono un pacchetto open-source per la matematica che offra calcolo sia simbolico che numerico. Perlopiu' le persone sono contente di quanto offre Sage. Com'e' comune nell'ambito del  software open-source (FOSS), spesso ci sono persone che individuano casi in cui Sage non dispone della funzione richiesta da loro, e si immergono nel codice sorgente di Sage per estenderlo per il loro scopo, o ancora per esporre funzionalita' dei pacchetti inglobati in Sage in modo da poter usarne le funzioni che loro preferiscono dall'interfaccia di Sage.
La squadra `Sage-Combinat <http://combinat.sagemath.org>`_ e' costituita da ricercatori in Algebra Combinatoria. La missione che si e' data tale squadra e' quella di migliorare Sage come uno strumento estendibile per la sperimentazione al computer nell'ambito dell'Algebra Combinatoria, e favorire lo scambio di codice sorgente fra ricercatori di questa materia. Per informazione dettagliate sul perche' esiste Sage, vedere il seguente link: biografia matematica personale di William, `settore software <http://sagemath.blogspot.com/2009/12/mathematical-software-and-me-very.html>`_.


Cosa vuol dire Sage e come si deve pronunciarlo ?
"""""""""""""""""""""""""""""""""""""""""""""""""

Nei primi anni di esistenza di Sage, il progetto era chiamato "SAGE", acronimo inglese di "software per la sperimentazione in Algebra e Geometria". A cominciare dal 2007 ed inizio 2008, fu largamente adottato il nome "Sage". Considera "Sage" come il nome di un progetto software FOSS per la matematica, esattamente come "Python" e' il nome di un linguaggio FOSS di programmazione di uso generale. Ovunque possibile, per cortesia usa "Sage" non "SAGE", che invece e' un progetto per un computer (`SAGE <http://history.sandiego.edu/GEN/20th/sage.html>`_), cosi' da evitare confusioni. Si pronuncia "Sage" nello stesso modo della parola inglesi "sage" che significa uomo saggio, o anche indica la pianta della salvia. Alcune persone lo pronunciano in modo simile a "sarge", un po' come si pronuncia `Debian <http://www.debian.org>`_Sarge. Ma comunque lo pronunci, per cortesia non confonderlo con il software di contabilita' americano che ha lo stesso nome.


Chi c'e' dietro al progetto ?
"""""""""""""""""""""""""""""

Sage e' un progetto basato sull'opera di volontari. Il suo successo e' dovuto all'opera gratuita di di una grande squadra internazionale di studenti, insegnanti, professori universitari, ricercatori, ingegneri del software, e persone che lavorano in vari ambiti della matematica, delle scienze, dell'ingegneria, dello sviluppo software, e a tutti i livelli della scuola. Lo sviluppo di Sage ha potuto usufruire di fondi asegnati da numerose istituzioni, ed ha potuto includere sia componenti preesistenti che in corso di sviluppo da parte di numerosi autori. Una lista di coloro che hanno dato un contributo diretto e' reperibile al link "mappa di sviluppo di Sage" (`Sage Development Map <http://www.sagemath.org/development-map.html>`_) e la storia delle modifiche puo' essere reperita al link "changelog di alto livello" (`changelog <http://www.sagemath.org/mirror/src/changelog.txt>`_). Fai riferimento alla `Pagina dei riconoscimenti <http://www.sagemath.org/development-ack.html>`_ del sito web di Sage per una lista aggiornata di coloro che ci sostengono finanziariamente o a livello di infrastruttura, a livello di siti mirror, ed altri contributi indiretti.


Perche' Sage e' un software libero ed open-source ?
"""""""""""""""""""""""""""""""""""""""""""""""""""

Una regola universale nella comunita' matematica e' che tutto dev'essere chiaro ed pubblico. Il progetto Sage ritiene che non seguire lo stesso principio nel software per la matematica e' quantomeno scortesia e maleducazione, o peggio una violazione delle pratiche comuni nella scienza. Un principio filosofico sottostante Sage e' di applicare la regola di scambio libero e confronto fra pari, che caratterizza la comunicazione scientifica, anche allo sviluppo di software per la matematica. Ne' il progetto Sage ne' la sua squadra di sviluppo hanno la pretesa di essere gli originali proponenti di questo principio. Il modello di sviluppo di Sage e' largamente ispirato del movimento del software libero di cui e' stata pioniere la `Free Software Foundation <http://www.fsf.org>`_ ed il movimento open-source. Una fonte di ispirazione all'interno della comunita' matematica e' Joachim Neubüser, come espresso nell'articolo::

* J. Neubüser. An invitation to computational group theory. In C. M. Campbell, T. C. Hurley,
  
  E. F. Robertson, S. J. Tobin, and J. J. Ward, editors, *Groups '93 Galway/St. Andrews, Volume 2*,
  
  volume 212 of London Mathematical Society Lecture Note Series, pages 457--475. Cambridge
  
  University Press, 1995.

ed in particolare nella seguente citazione dal suo articolo::

  Puoi leggere il teorema di Sylow e la sua dimostrazione nel libro di Huppert in biblioteca senza
  nemmeno comprare il libro e poi usare questo teorema per il resto della tua vita senza dover pagare
  una tariffa, invece...devi pagare regolarmente delle licenze per l'uso di software per la matematica
  per tutto il tempo in cui li utilizzi. Per proteggere cio' per cui devi pagare, non ti viene dato il
  codice sorgente ma soltanto l'eseguibile del programma, cioe' un oggetto sigillato sul quale premi
  bottoni ed ottieni risposte cosi' come ottieni belle immagini dal tuo televisore: in entrambe le
  situazioni non puoi controllare il modo in cui ti e' fornito il servizio.

  In questa situazione sono violate le piu' basilari regole di condotta in matematica: in essa
  l'informazione e' passata gratuitamente e tutto puo' essere consultato e verificato. Non applicare
  queste regole ai sistemi software per la matematica usati per la ricerca in matematica...significa
  andare in una direzione assolutamente non desiderabile. Ancora piu' importante: possiamo aspettarci
  che qualcuno creda al risultato di un programma che non puo' vedere come funziona? Inoltre: vogliamo
  veramente far pagare a colleghi in Moldavia anni di stipendio per un software per la matematica?

Simili idee sono state anche espresse da Andrei Okounkov, come si puo' leggere in::

* V. Muñoz and U. Persson. Interviews with three Fields medalists. *Notices of the American

  Mathematical Society*, 54(3):405--410, 2007.

ed in particolare nella seguente citazione::

  I computer non sono una minaccia per i matematici piu' di quanto i robot da cucina lo siano per i
  cuochi. Poiche' la matematica diviene sempre piu' complessa mentre il ritmo delle nostre vite
  accellera, dobbiamo delegare piu' che possiamo alle macchine. Ed intendo sia il lavoro in campo
  numerico che in quello simbolico. Alcune persone possono andare avanti senza lavapiatti, ma penso
  che le dimostrazioni vengano fuori molto piu' pulite quando il lavoro di routine e' automatizzato.

  Questo porta con se' parecchie questioni. Non sono un esperto ma penso che abbiamo bisogno di un
  standard a livello di calcolo simbolico per rendere le manipolazioni al computer piu' facili da
  documentare e verificare. Con tutto il rispetto per il libero mercato, forse in questo non dobbiamo
  esser dipendenti da un software commerciale. Un progetto open-source potrebbe, forse, trovare
  risposte migliori a problemi ovvi come la disponibilita', i bachi, la compatibilita' all'indietro,
  l'indipendenza dalla piattaforma, le librerie standard, ecc. Si puo' imparare dal successo di TeX e
  da software piu' specializzato come Macaulay2. Spero veramente che le agenzie per finanziamenti
  governativi stiano considerando questo.


Perche' avete scritto Sage da zero, invece di usare software e librerie preesistenti ?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Sage non e' stato scritto da zero. La maggior parte delle sue funzionalita' sono realizzate attraverso progetti FOSS come::

* `ATLAS <http://math-atlas.sourceforge.net>`_ --- libreria software per Algebra Lineare

  ottimizzata automaticamente.

* `BLAS <http://www.netlib.org/blas>`_ --- sottoprogrammi per Algebra Lineare di base.

* `FLINT <http://www.flintlib.org>`_ --- libreria C per Teoria dei Numeri.

* `GAP <http://www.gap-system.org>`_ --- sistema di calcolo per algebra discreta, con

  particolare enfasi sulla teoria dei gruppi computazionale.

* `Maxima <http://maxima.sourceforge.net>`_ --- sistema di calcolo simbolico e numerico.

* `mpmath <http://code.google.com/p/mpmath>`_ --- libreria in puro Python per aritmetica

  floating-point di precisione.

* `NumPy <http://numpy.scipy.org>`_ --- algebra lineare numerica ed altre funzioni di calcolo

  numerico per Python.

* `Pari/GP <http://pari.math.u-bordeaux.fr>`_ --- software matematico per calcolo veloce in

  Teoria dei Numeri.

* `Pynac <http://pynac.sagemath.org>`_ --- versione modificata di GiNaC che rimpiazza la

  dipendenza da CLN con Python.

* `R <http://www.r-project.org>`_ --- linguaggio ed ambiente operativo per calcolo statistico

  e grafici relativi.

* E molti altri troppo numerosi per essere elencati qui.

Una lista aggiornata puo' essere reperita alla seguente link: `repository dei pacchetti standard <http://www.sagemath.org/packages/standard>`_.
I principali linguaggi di programmazione di Sage sono
`Python <http://www.python.org>`_ e `Cython <http://www.cython.org>`_.
Python e' il principale linguaggio di programmazione e di interfacciamento, mentre Cython e' il principale linguaggio per ottimizzare funzionalita' critiche e per interfacciarsi con le librerie C e le estensioni C per Python. Sage integra oltre 90 pacchetti FOSS in un'interfaccia comune. Sopra questi pacchetti sta la libreria Sage, che consiste in oltre 700.000 righe di codice Python e Cython scritto ex-novo. Vedi `ohloh.net <https://www.ohloh.net/p/sage/analyses/latest>`_ per l'analisi del codice sorgente dell'ultima release stabile di Sage.


Chi usa Sage ?
""""""""""""""

Di seguito v'e' una lista incompleta di istituzioni e progetti che usano Sage. Se qualche istituzione o progetto manca, per cortesia fatecelo sapere riportandolo sulla mailing list  `sage-devel <http://groups.google.com/group/sage-devel>`_.

#. `California Institute of Technology <http://www.caltech.edu>`_, Pasadena, California, USA
#. `California Polytechnic State University <http://www.calpoly.edu>`_, San Luis Obispo, CA, USA
#. `Chang Gung University <http://www.cgu.edu.tw>`_, Taiwan
#. `Chapman University <http://www.chapman.edu>`_, Orange, CA, USA
#. `Clemson University <http://www.clemson.edu>`_, Clemson, South Carolina, USA
#. `Drake University <http://www.drake.edu>`_, Des Moines, IA, USA
#. `FEMhub <http://www.femhub.org>`_, una distribuzione open source di codice per il calcolo
   scientifico integrato da un'interfaccia unificata in Python. I notebook FEMhub sono basati
   sui notebook Sage.
#. `Gordon College <http://www.gordon.edu>`_, Wenham, MA, USA
#. `Korea Advanced Institute of Science and Technology <http://www.kaist.edu>`_, Daejeon, Korea
#. `Mendel University in Brno <http://www.mendelu.cz>`_, Czech Republic
#. `Reykjavik University <http://www.ru.is>`_, Iceland
#. `Universidad Autónoma de Madrid <http://www.uam.es>`_, Spain
#. `Universidad de la República <http://www.universidad.edu.uy>`_, Montevideo, Uruguay
#. `Universitat Politècnica de Catalunya <http://www.upc.edu>`_, Barcelona, Catalonia, Spain
#. `Université Claude Bernard Lyon 1 <http://www.univ-lyon1.fr>`_, France
#. `Université de Provence <http://www.univ-mrs.fr>`_, Marseille, France
#. `Universiteit Leiden <http://www.leidenuniv.nl>`_, The Netherlands
#. `University of Canterbury <http://www.canterbury.ac.nz>`_, Christchurch, New Zealand
#. `University of Minnesota Duluth <http://www.d.umn.edu>`_, Duluth, MN, USA
#. `University of Nevada, Reno <http://www.unr.edu>`_, Reno, NV, USA
#. `University of Puget Sound <http://www.pugetsound.edu>`_, Tacoma, WA, USA
#. `University of Washington <http://www.washington.edu>`_, Seattle, Washington, USA
#. `University of Wisconsin, Oshkosh <http://www.uwosh.edu>`_, Oshkosh, WI, USA
#. `US Naval Academy <http://www.usna.edu>`_, Annapolis, Maryland, USA


Come posso ricevere aiuto ?
"""""""""""""""""""""""""""

Sage ha due liste email molto attive::

* ``sage-devel``: http://groups.google.com/group/sage-devel
* ``sage-support``: http://groups.google.com/group/sage-support

Vi e' anche un canale IRC molto attivo: ``#sage-devel`` su freenode. Molti sviluppatori hanno anche dei blog aggiornati e pubblicano altri tutorial e discussioni relative a Sage. Consulta http://www.sagemath.org/help.html per una lista di queste risorse.


Non sarebbe meglio se Sage non fosse distribuito come un gigantesco aggregato di pacchetti ?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Quest'aspetto e' stato discusso a fondo piu' volte. Quindi prima di ricominciare a discutere, leggi bene e rifletti su quanto segue. Sage e' una distribuzione di oltre 90 pacchetti FOSS per calcolo simbolico, numerico e scientifico. In generale, l'insieme di configurazioni possibili da gestire sarebbe di gran lunga troppo grande. E' pressoche' impossibile trovare una qualunque distribuzione Linux (quali Arch,  CentOS, Debian, Fedora, Gentoo, Mandriva, Ubuntu) che abbia un numero di dipendenze che si avvicini anche lontanemente al numero di versione dei pacchetti da cui dipende Sage.

La maggior parte delle persone che contribuiscono a Sage lo fanno nel loro tempo libero. Queste sono persone che hanno un lavoro quotidiano non direttamente collegato allo sviluppo software o alla programmazione. E' pressoche' impossibile per chiunque tenere traccia della versione corretta dei pacchetti, configurarli e compilarli su Linux, Mac OS X, Solaris o Windows, solo per poter iniziare ad usare Sage o iniziare a dare il loro primo contributo a Sage. Dal momento che il progetto Sage aspira ad essere utile ad un pubblico il piu' ampio possibile, crediamo che Sage debba innanzitutto essere il piu' semplice possibile da installare per chiunque, con qualunque livello di conoscenze informatiche. Se vuoi aiutare Sage a realizzare quest'obiettivo puoi contattare la mailing list `sage-devel <http://groups.google.com/group/sage-devel>`_.


Perche' ci sono cosi' tanti bachi in Sage, con centinaia di modifiche in corso, perche' non producete una versione stabilizzata ?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Ogni software contiene bachi. In qualcosa di cosi' complesso come Sage nessuno, ne' la squadra di sviluppo di Sage ne' la sua comunita', ha alcuna pretesa che esso sia libero da bachi. Farlo sarebbe un atto di disonesta'.

Un ciclo di rilascio di Sage di solito dura dalle 3 alle 4 settimane. Ogni ciclo di rilascio e' presieduto da un singolo gestore che si occupa dell'albero di integrazione pacchetti per tutta la durata del ciclo. In questa fase tale gestore deve spesso dedicare tempo equivalente ad un lavoro a tempo pieno alla gestione della qualita', e deve interagire attivamente con la comunita' internazionale degli utenti, degli sviluppatore e dei potenziali contributori a Sage. Ci sono stati molti casi in cui due contributori a Sage sono stati affiancati come gestori di rilascio per un ciclo di rilascio di Sage. Comunque in genere poche persone hanno tempo libero per l'equivalente di 3 settimane per dedicarsi alla gestione del rilascio. Se vuoi aiutare nella gestione del rilascio iscriviti alla mailing list `sage-release <http://groups.google.com/group/sage-release>`_.

Fin dall'inizio del progetto Sage i contributori hanno cercato di ascoltare e di riflettere su cosa potesse aumentare la possibilita' che altri potenziali validi contributori dessero effettivamente un aiuto. Cosa incoraggia un contributore puo' scoraggiare un altro, quindi bisogna trovare degli equilibri. Decidere che un rilascio stabilizzato dovrebbe includere le patch di correzione dei bachi, e solo quelle, probabilmente scoraggerebbe qualcuno dal contribuire, nel momento in cui gli fosse detto in anticipo che la sua aggiunta, anche se giudicata positivamente, non verrebbe integrata nel rilascio. La comunita' Sage crede nel principio "pubblica subito, pubblica spesso". Il modo in cui il progetto Sage e' organizzato e gestito differisce parecchio da quello di una azienda di software commerciale. I contributori sono tutti volontari e questo cambia totalmente la dinamica del progetto da quella che sarebbe se Sage fosse un'iniziativa software commerciale con sviluppatori stipendiati a tempo pieno.


Come posso scaricare la documentazione di Sage cosi' da poterla leggere offline ?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Per scaricare la documentazione standard di Sage in formato HTML o PDF, visita `Help and Support <http://www.sagemath.org/help.html>`_ sul sito web di Sage. Ogni release di Sage dispone della documentazione completa che costituisce la documentazione standard di Sage. Se hai scaricato un rilascio di Sage in formato binario, la versione HTML della sua documentazione si prova gia' disponibile nella cartella ``SAGE_ROOT/src/doc/output/html/``. Nel corso della compilazione da sorgente viene preparata anche la documentazione HTML, che comunque puo' essere preparata da riga di comando lanciando, dopo essersi posizionati in ``SAGE_ROOT``::

    $ ./sage -docbuild --no-pdf-links all html

Invece la preparazione della documentazione in formato PDF richiede che sul tuo sistema sia installata una versione funzionante di LaTeX. Per preparare la documentazione in formato PDF puoi lanciare da riga di comando, dopo esserti posizionato in ``SAGE_ROOT``::

    $ ./sage -docbuild all pdf

Per altre maggiori opzioni disponibili a riga di comando fai riferimento alle istruzioni stampate dei seguenti comandi::

    $ ./sage -help
    $ ./sage -advanced

