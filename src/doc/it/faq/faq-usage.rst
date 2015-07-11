.. _chapter-faq-usage:

================
FAQ: Uso di Sage
================


Come posso partire?
"""""""""""""""""""

Puoi provare Sage senza dover scaricare nulla. Vai al link http://www.sagenb.org e richiedi un account gratuito. Collegandoti lavorerai su un server gratuito di sessioni di lavoro (detti Notebook cioe' quaderni) di Sage che lavora in modo identico ad un'installazione locale di Sage. Per scaricare una distribuzione binaria (ciao' precompilata) di Sage vai al link http://www.sagemath.org/download.html e clicca sul link relativo al file binario per il tuo sistema operativo. Il codice sorgente di Sage e' anche disponibile in download: vai al link http://www.sagemath.org/download-source.html per scaricare l'archivio TAR di qualunque rilascio di Sage. I rilasci precedenti sono disponibili al link http://www.sagemath.org/src-old.

Le sessioni di lavoro (Notebook) di Sage sono eseguite all'interno di un browser web. Puoi eseguire Sage in un browser che non sia quello predefinito nel tuo sistema. Per farlo lancia, da riga di comando o dal menu di Sage ::

    env SAGE_BROWSER=opera /usr/bin/sage -notebook


Quali sono i prerequisiti di Sage?
""""""""""""""""""""""""""""""""""

La maggior parte delle dipendenze sono incluse all'interno di Sage. Nella maggior parte dei casi puoi scaricare il binario precompilato ed usarlo senza dover installare alcun pacchetto dipendente. Se usi Windows avrai bisogno di intallare `VirtualBox <http://www.virtualbox.org>`_, che puoi scaricare dal link http://www.virtualbox.org/wiki/Downloads. Dopo aver installato VirtualBox devi scaricare una distribuzione di Sage per VirtualBox al link http://www.sagemath.org/download-windows.html. Segui bene le istruzioni che trovi a quella pagina. Poi puoi lanciare la macchina virtuale Sage usando il software VirtualBox e, dopo aver aspettato che la macchina virtuale sia partita, digitare ``notebook`` da riga di comando.

Puoi scaricare il codice sorgente completo di Sage per compilarlo sul tuo sistema Linux o Mac OS X. Sage si trova in una cartella isolata e non va ad interferire col sistema in cui si trova. Viene dato con tutto il necessario per lo sviluppo, il codice sorgente, tutte le dipendenze ed il changelog (cioe' la lista delle modifiche operate) completo. Sui sistemi Linux come Debian/Ubuntu puoi dover installare il pacchetto ``build essential`` ed il processore di macro ``m4``. Il tuo sistema deve disporre di un compilatore C funzionante se vuoi compilare Sage da codice sorgente. Su Debian/Ubuntu puoi installare questi prerequisiti come segue::

    sudo apt-get install build-essential m4

Se hai un sistema multiprocessore puoi scegliere una copilazione parallela di Sage. Il comando ::

    export MAKE='make -j8'

abilitera' 8 threads per quelle parti della compilazione che supportano il parallelismo. Al posto del numero 8 metti il numero di processori/core del tuo sistema.


Come posso far riconoscere la mia attuale installazione di Tcl/Tk all'interprete Python di Sage?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Potresti avere la libreria Tcl/Tk installata ed l'interprete Python del tuo sistema la riconosce ma l'interprete Python di Sage no. Per risolvere questo ti basta installare la libreria di sviluppo Tcl/Tk. Su Ubuntu lancia, da riga di comando::

    sudo apt-get install tk8.5-dev

o qualcosa di simile. Poi reinstalla l'iterprete Python di Sage con::

    sage -f python

Questo aggancera' automaticamente la libreria Tcl/Tk. Dopo aver reinstallato correttamente l'interprete Python di Sage, lancia i seguenti comandi dall'interfaccia a riga di comando di Sage::

    import _tkinter
    import Tkinter

Se non ti viene segnalato alcun errore di ``ImportError` allora il problema e' risolto.


Come faccio ad importare Sage in uno script Python?
"""""""""""""""""""""""""""""""""""""""""""""""""""

Puoi importare Sage in uno script Python come faresti con una libreria. La cosa a cui fare attenzione e' che devi lanciare quello script Python usando l'interprete Python interno a Sage, che correntemente e' il 2.6.x
Per importare Sage metti la seguente istruzione in cima al tuo script Python::

    from sage.all import *

Quando poi esegui il tuo script devi lanciarlo con l'opzione ``-python`` che fara' si' che venga eseguito dalla versione dell'interprete interna a Sage. Ad esempio, se Sage e' nella tua variabile d'ambiente ``PATH``, puoi scrivere::

    sage -python /path/to/my/script.py

Un altro modo e' scrivere uno script Sage e lanciarlo usando Sage stesso. Uno script Sage ha estensione  ``.sage`` ed e' simile ad uno script Python ma utilizza funzioni e comandi specifici di Sage. Puoi poi lanciare tale script Sage in questo modo::

    sage /path/to/my/script.sage

Questo si occupera' di caricare le variabili d'ambiente necesssarie ed eseguire gli import di default al posto tuo.


Come posso ricaricare uno script Python in una sessione di Sage?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Puoi caricare uno script Python in una sessione Sage usando il comando  **load**. Ad esempio possiamo usare Sage per importare un file di nome "simple.py" con::

    load("simple.py")

e ripetere questo comando ogni volta che cambiamo il file. Invece digitando::

    attach("simple.py")

ogni cambiamento al file verra' automaticamente aggiornato anche in Sage.


Posso usare Sage con la versione 3.x di Python?
"""""""""""""""""""""""""""""""""""""""""""""""

Al momento no. Sage dipende dalle funzionalita' numeriche e scientifiche della libreria `SciPy <http://www.scipy.org>`_ e, ancora nel 2010, SciPy utilizza Python 2.x. Pertanto finche' SciPy non sara' portata a Python 3.x e `Cython <http://www.cython.org>`_ non supportera' Python 3.x, Sage continuera' ad usare Python 2.x.


Vedo un errore di "Permission denied" (accesso negato) su un file di nome "sage-flags.txt.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Quando Sage viene compilato dal sorgente, tiene traccia di quali istruzioni speciali supporta la tua CPU (quali ad esempio SSE2) e le memorizza. Cosi' se provi ad eseguire il codice su un'altra macchina, che non supporta queste istruzioni speciali extra, ti vengono segnalati degli errori in maniera intelleggibile anziche' dei generici "segfault" (segmento di memoria errato) o "illegal istruction" (istruzione non consentita). Poiche' quest'informazione  dev'essere memorizzata in Sage stesso anziche' nella cartella ``.sage``, dev'essere creata da qualcuno con le necessarie autorizzazioni sul sistema. Quindi se vedi qualcosa del genere ::

    Traceback (most recent call last):
      File "/usr/local/sage-4.0.2/local/bin/sage-location", line 174, in <module>
        t, R = install_moved()
      File "/usr/local/sage-4.0.2/local/bin/sage-location", line 18, in install_moved
        write_flags_file()
      File "/usr/local/sage-4.0.2/local/bin/sage-location", line 82, in write_flags_file
        open(flags_file,'w').write(get_flags_info())
    IOError: [Errno 13] Permission denied:
      '/usr/local/sage-4.0.2/local/lib/sage-flags.txt'

probabilmente significa che hai compilato/installato Sage usando un determinato account (nome utente), ma poi non l'hai eseguito cosi' da permettergli di generare il file ``sage-flags.txt``. Ti basta eseguire Sage una volta con lo stesso account con cui e' stato installato per risolvere questo problema. Questo si dovrebbe risolvere facilmente anche lanciando Sage una volta nel corso del processo d'installazione (cfr. `correzione #6375 <http://trac.sagemath.org/sage_trac/ticket/6375>`_).


Ho scaricato il binario di Sage e va in crash quando lo lancio, con il messaggio "illegal instruction" (istruzione non permessa). Cosa posso fare?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Un modo di risolvere e' compilare Sage interamente dal codice sorgente. Un'altra possibilita' e' correggere la tua installazione di Sage con la ricompilazione dei componenti MPIR e ATLAS (richiede da 15 a 20 minuti), da effettuarsi a riga di comando a partire dalla cartella ``SAGE_ROOT`` della tua installazione con le 2 istruzioni::

    rm spkg/installed/mpir* spkg/installed/atlas*
    make

E' possibile che i binari siano stati compilati per un'architettura piu' recente di quella della tua macchina. Nessuno ha ancora trovato un modo di compilare Sage in maniera che MPIR ed ATLAS funzionino su qualunque hardware. Questo sara' prima o poi risolto. Qualunque aiuto in tal senso sara' apprezzato.


Ho usato Debian/Ubuntu per installare la versione 3.0.5 di Sage ed essa sta dando un sacco di errori. Cosa posso fare?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

La versione di Sage distribuita con ``apt-get`` in Debian e Ubuntu (tipo la 3.0.5) e' molto vecchia. Nessuno ha ancora avuto tempo di aggiornare la versione di Sage per Debian/Ubuntu. Qualunque aiuto in tal senso sara' molto apprezzato. Dovresti scaricare la versione piu' recente di Sage dal `link di download <http://www.sagemath.org/download.html>`_ del sito web di Sage. Se vuoi aiutarci ad aggiornare la versione di Sage per Debian/Ubuntu manda un'email alla mailing list `sage-devel <http://groups.google.com/group/sage-devel>`_.


Faccio meglio ad usare la versione ufficiale o quella di sviluppo?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Ti consigliamo di usare la piu' recente versione ufficiale di Sage. Delle versioni di sviluppo sono spesso annunciate sulle mailing list `sage-devel <http://groups.google.com/group/sage-devel>`_ e `sage-release <http://groups.google.com/group/sage-release>`_. Una maniera facile di aiutare con lo sviluppo di Sage e' scaricare l'ultima versione di sviluppo, compilarla sul suo sistema, lanciare tutti i doctest e segnalare qualunque errore di compilazione o qualunque fallimento nei doctest.


E' difficile imparare Sage?
"""""""""""""""""""""""""""

Le funzionalita' di base di Sage dovrebbero risultare facili da imparare quanto le basi di Python. Molti tutorial sono disponibili in rete per aiutarti ad imparare Sage. Per trarre il massimo da Sage ti consigliamo di impararare qualche elemento del linguaggio Python. Segue una lista, incompleta, di risorse su Python. Altre risorse possono essere trovate cercando sul web.

* `Building Skills in Python <http://homepage.mac.com/s_lott/books/python.html>`_ di Steven F. Lott
* `Dive into Python <http://www.diveintopython.net>`_ di Mark Pilgrim
* `How to Think Like a Computer Scientist <http://www.openbookproject.net/thinkCSpy>`_ di Jeffrey Elkner, Allen B. Downey, and Chris Meyers
* `Official Python Tutorial <http://docs.python.org/tutorial>`_
* `Python <http://www.python.org>`_ home page e `Python standard documentation <http://docs.python.org>`_


Posso fare X in Sage?
"""""""""""""""""""""

Ti consigliamo di usare l'autocompletamento di Sage con il tasto TAB. Ti basta digitare qualche carattere, premere TAB e vedere se il comando che vuoi compare nella lista di autocompletamento. Se hai un comando che si chiama ``mycmd``, allora digitalo e premi TAB per visualizzare la lista di funzionalita' che sono supportate da quel comando. Per leggere la documentazione di ``mycmd`` scrivi ``mycmd?`` poi premi Invio e protrai leggerla. Similmente, digitando ``mycmd??`` e poi Invio potrai visualizzare il codice sorgente di tale comando. Ti consigliamo anche di eseguire ricerche nel codice sorgente e nella documentazione di Sage. Per eseguire ricerche nel codice sorgente di Sage usa il comando::
``search_src("<search-keyword>")``
mettendo al posto di ``<search-keyword>`` le parole chiave che vuoi cercare.
Analogamente puoi effettuare ricerche nella documentazione di Sage usando il comando:
``search_doc("<search-keyword>")``.


Cosa fa esattamente Sage quando digito "0.6**2" ?
"""""""""""""""""""""""""""""""""""""""""""""""""

Quando scrivi "0.6**2" in Python, ti viene restituito qualcosa tipo  0.35999999999999999. Ma quando fai lo stesso in Sage ti viene restituito 0.360000000000000. Per capire perche' Python si comporta in questo modo vedi il `Python Tutorial <http://docs.python.org/tutorial/floatingpoint.html>`_, soprattutto il capitolo "Aritmetica floating-point: caratteristiche e limiti" (http://docs.python.org/tutorial/floatingpoint.html). Cio' che Sage fa e' preprocessare l'input e trasformarlo come segue::

    sage: preparse("0.6**2")
    "RealNumber('0.6')**Integer(2)"

Cosi' che cio' che viene *effettivamente* eseguito e'::

    RealNumber('0.6')**Integer(2)

Gli sviluppatori Sage (in pratica Carl Witty) decisero che i numeri floating-point di Sage dovessero, di default, stampare solo il numero di cifre decimali corrette, quando possibile, cosi' da evitare il problema che ha Python. Questa decisione ha i suoi pro e contro. Nota che ``RealNumber`` e ``Integer`` sono specifici di Sage, quindi non puoi digitare quanto sopra nell'interprete Python ed aspettarti che funzioni, se prima non hai eseguito delle istruzioni di import quali::

    from sage.all import RealNumber, Integer, preparse


Perche' il comando "history" di Sage e' diverso da quello di Magma?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Nell'uso di Sage non disponi di una funzionalita' dell'interfaccia a riga di comando di Magma. In Magma, se immetti una linea recuperata dalla "history" (cioe' dall'elenco dei comandi digitati precedentemente che viene automaticamente memorizzato) con il tasto "freccia in su'" e poi premi "freccia in giu'", viene recuperata anche la linea successiva nell'elenco. Questa funzionalita' ti permette di recuperare dalla "history" tante righe consecutive quante vuoi. Ma Sage non ha una funzionalita' simile: la riga di comando `IPython <http://ipython.scipy.org>`_ utilizza la libreria "readline" (via pyreadline), che evidentemente non supporta questa funzionalita'. Magma ha una sua propria libreria personalizzata simile alla "readline" che invece supporta questa funzionalita'. (Dal momento che moltissime persone hanno richiesto questa funzionalita', se qualcuno trovasse un modo per implementarla sarebbe il benvenuto !)


Ho problemi di tipo nell'utilizzo da Sage di SciPy, cvxopt e NumPy.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Stai usando da Sage le librerie SciPy, cvxopt e NumPy e hai degli errori tipo::

    TypeError: function not supported for these types, and can't coerce safely to supported types.

Quando digiti numeri in Sage, il preprocessore li converte in un anello base, come puoi vedere facendo:
sage::

    sage: preparse("stats.uniform(0,15).ppf([0.5,0.7])")
    "stats.uniform(Integer(0),Integer(15)).ppf([RealNumber('0.5'),RealNumber('0.7')])"

Sfortunamente il supporto che NumPy fornisce a questi tipi avanzati di Sage, quali ``Integer`` o ``RealNumber`` (numeri reali di precisione arbitraria), non e' del 100%. Per risolvere ridefinisci ``Integer`` e/o ``RealNumber`` per cambiare il comportamento del preprocessore di Sage cosi' che i decimali scritti vengano registrati come tipi float di Python anziche' RealNumber di Sage e gli interi scritti siano registrati come tipi int di Python anziche' Integer di Sage. Ad esempio::

    sage: RealNumber = float; Integer = int
    sage: from scipy import stats
    sage: stats.ttest_ind(list([1,2,3,4,5]),list([2,3,4,5,.6]))
    (array(0.07675295564533369), 0.94070490247380478)
    sage: stats.uniform(0,15).ppf([0.5,0.7])
    array([  7.5,  10.5])

In alternativa sii esplicito circa il tipo di dato, ad esempio::

    sage: from scipy import stats
    sage: stats.uniform(int(0),int(15)).ppf([float(0.5),float(0.7)])
    array([  7.5,  10.5])

Come terza alternativa puoi usare i suffissi semplici::

    sage: from scipy import stats
    sage: stats.uniform(0r,15r).ppf([0.5r,0.7r])
    array([  7.5,  10.5])

Puoi anche disabilitare il preprocessore nel tuo codice tramite il comando ``preparse(False)``.
Puoi lanciare Ipython da solo dalla riga di comando con ``sage -ipython``, cosa che non precarica niente di specifico di Sage. O ancora puoi cambiare il linguaggio di sessione (Notebook language) in "Python".


Come faccio a salvare un oggetto cosi' che non devo ridigitarlo ogni volta che apro un foglio di lavoro (worksheet) ?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

I comandi ``save`` e ``load`` rispettivamente registrano e caricano un oggetto. Nella sessione di lavoro Notebook la variabile ``DATA`` e' la locazione dello spazio di salvataggio del foglio di lavoro (worksheet). Per registrare l'oggetto ``my_stuff`` in un foglio di lavoro puoi digitare::

    save(my_stuff, DATA + "my_stuff")

e, per ricaricarlo, ti basta digitare::

    my_stuff = load(DATA + "my_stuff")


Ho un errore da jsMath oppure un simbolo matematico non e' visualizzato correttamente nella sessione Notebook.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Se vedi l'errore ::

    It looks like jsMath failed to set up properly (error code -7). I will try to keep going,
    but it could get ugly.

allora vuol dire che non hai installato i font TeX che aiutano jsMath a visualizzare i suoi bei simboli matematici. Affinche' si veda il gradevole TeX assieme a jsMath, devi scaricare un insieme di font dal link http://www.math.union.edu/~dpvc/jsMath/download/jsMath-fonts.html . Se sei un utente Linux ignora le istruzioni su quel sito e semplicemente decomprimi i font nella sottocartella ``.fonts`` della tua cartella home. Puoi anche installare il pacchetto ``jsmath-fonts``.


Sage contiene una funzione simile alla "ToCharacterCode[]" di Mathematica?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Potresti voler convertire caratteri ASCII come "Big Mac" nel corrispondente codice numerico per ulteriori elaborazioni. In Sage e Python puoi usare ``ord``. Ad esempio::

    sage: map(ord, "abcde")
    [97, 98, 99, 100, 101]
    sage: map(ord, "Big Mac")
    [66, 105, 103, 32, 77, 97, 99]


Posso far eseguire in automatico a Sage dei comandi all'accensione?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Si', ti basta creare un file ``init.sage`` nella sottocartella ``.sage`` della tua home, ed esso sara' eseguito ogni volta che lanci Sage. Questo presuppone che la variabile ambiente di Sage ``DOT_SAGE`` punti alla cartella nascosta ``$HOME/.sage``, cosa che avviene di default.


Il mio aggiornamento di Sage e' fallito, segnalando simboli gmp mancanti su OSX 10.4. Cosa posso fare?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Spostare un'installazione di Sage su Mac OS X 10.4 e poi aggiornare qualcosa collegato alla libreria NTL porta ad errori di collegamento dovuti alla mancanza dei simboli gmp. Il problema e' la modalita' del collegamento con cui e' creato l'NTL dinamico. C'e' una soluzione ma si sta ancora verficando che risolva realmente il problema. Tutto cio' che e' collegato a NTL dev'essere ricompilato, ad esempio le librerie Singular e Cremona. A complicare la questione c'e' il fatto che questo problema non si verifica su Mac OS X 10.5. Una correzione per questo problema e' stata aggiunta in Sage 2.8.15 dunque per cortesia se avete quest'errore in un rilascio piu' recente di Sage segnalatecelo.


Quando compilo Sage il mio computer fa beep e si spegne o si blocca.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Compilare sage e' piuttosto faticoso per il processore del computer. Il comportamento suddetto di solito indica che il computer si e' surriscaldato. In molti casi questo si puo' risolvere pulendo il ventilatore del processore del computer ed assicurando adeguata areazione al computer. Puoi chiedere al tuo amministratore di sistema o ad un tecnico di provvedere, qualora tu non l'abbia mai fatto. Questa manutenzione del computer, se non fatta da persone preparate, potrebbe anche danneggiare il computer.

Per gli utenti Linux, se pensi che la compilazione fallisca per un problema di risorse di macchina, una soluzione potrebbe essere di modificare il file ``/etc/inittab`` per far partire Linux al runlevel 3. Tale file di solito contiene qualcosa del tipo::

    #   0 - halt (Do NOT set initdefault to this)
    #   1 - Single user mode
    #   2 - Multiuser, without NFS (The same as 3, if you do not have
    #   networking)
    #   3 - Full multiuser mode
    #   4 - unused
    #   5 - X11
    #   6 - reboot (Do NOT set initdefault to this)
    #
    id:5:initdefault:

Questo fa si' che la tua distribuzione Linux parta con la schermata di login grafico. Commenta la linea ``id:5:initdefault:`` e aggiungi la linea ``id:3:initdefault:``, cosi' da aver qualcosa come::

    #   0 - halt (Do NOT set initdefault to this)
    #   1 - Single user mode
    #   2 - Multiuser, without NFS (The same as 3, if you do not have
    #   networking)
    #   3 - Full multiuser mode
    #   4 - unused
    #   5 - X11
    #   6 - reboot (Do NOT set initdefault to this)
    #
    # id:5:initdefault:
    id:3:initdefault:

Ora se riavvii il sistema ti troverai davanti all'interfaccia di login testuale. Questa ti permette di accedere al sistema con una sessione testuale all'interno di un terminale virtuale. Una tale sessione di solito non consuma molte risorse, come farebbe invece un'interfaccia grafica. Poi puoi compilare Sage da codice sorgente in tale sessione testuale. Dovresti assicurarti di essere in grado di riattivare successivamente l'interfaccia grafica, prima di tentare di accedere tramite un'interfaccia testuale.


Quando lancio i doctest su Mac OS X vedo dei messaggi con "malloc", ma alla fine Sage dice che tutto e' andato bene.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

I messaggi "malloc" a cui ti riferisci potrebbero essere qualcosa tipo::

    sage -t  src/sage/libs/pari/gen.pyx
    python(4563) malloc: *** vm_allocate(size=4096000000) failed (error code=3)
    python(4563) malloc: *** error: can't allocate region
    python(4563) malloc: *** set a breakpoint in szone_error to debug

Questo comportamento non e' un fallimento dei doctest. E' un messaggio di errore stampato dal sistema ed e' esattamente quello che ci si aspetta di vedere. In quel particolare doctest, cerchiamo di allocate una lista molto grande in PARI che non ci sta nella memoria fisica (e' grande almeno 100 Gb). Quindi Mac OS X ti dice che non puo' allocare un blocco di memoria di circa 4 Gb, cosa attesa se stai usando Sage su una versione a 32 bit di Mac OS X ed hai compilato Sage nel modo a 32 bit oppure la tua distribuzione Sage binaria e' a 32 bit.


Sage 2.9 o superiore non riesce a compilare ATLAS su Linux. Come posso risolvere?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

La causa piu' probabile e' l'abilitazione della gestione dell'alimentazione. Disabilitala per risolvere il problema. In base al tuo tipo di distribuzione cio' si puo' fare da interfaccia grafica oppure no. Digita a riga di comando, come utente root, quanto segue, per ogni CPU presente sul tuo sistema::

    /usr/bin/cpufreq-selector -g performance -c #number CPU

Su Ubuntu, prova a disabilitare “Power Manager” (gestione alimentazione) via System --> Preferences --> Sessions nel menu “Startup Programs” (programmi di avvio) o utilizzando ``cpufreq-set`` da riga di comando.


Sage termina con il messaggio d'errore "restore segment prot after reloc: Permission denied". Cosa c'e' che non va?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Il problema e' collegato a SELinux. Vai al link seguente per dei suggerimenti su come risolvere questo problema: http://www.ittvis.com/services/techtip.asp?ttid=3092. Stiamo seguendo questo problema come  `correzione #480 <http://www.sagetrac.org/sage_trac/ticket/480>`_.


Quando lancio Sage, SELinux segnala che "/path/to/libpari-gmp.so.2" richiede "text-relocation" (riallocazione del testo). Come posso risolvere?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Il problema puo' essere risolto eseguendo il seguente comando::

    chcon -t textrel_shlib_t /path/to/libpari-gmp.so.2


L'aggiornamento di Sage e' andato bene, ma adesso l'indicatore continua a mostrare la versione precedente. Come posso risolvere?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

L'indicatore (banner in inglese) e' memorizzato e non ricalcolato ad ogni esecuzione di Sage. Il fatto che non sia aggiornato non dovrebbe impedire a Sage di funzionare regolarmente. Digita ``banner()`` in una sessione di Sage per verificare la versione reale. Se vuoi l'indicatore corretto allora devi ricompilare Sage digitando ``make build`` in un terminale.


Come posso eseguire Sage come demone/servizio?
""""""""""""""""""""""""""""""""""""""""""""""

Al momento non abbiamo una soluzione pronta. Ci sono parecchie possibilita'. Puoi usare i programmi a riga di comando ``screen``, ``nohup`` o ``disown``. Stiamo seguendo questo problema come `correzione #381 <http://www.sagetrac.org/sage_trac/ticket/381>`_ quindi seguici.


Sto utilizzando MacOS X. Dove devo mettere la cartella di font jsMath per far sparire il riquadro rosso?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Vai al link http://www.math.union.edu/~dpvc/jsMath/download/jsMath-fonts.html dove dice::

    Per utenti Mac OS X: scaricare e decomprimere l'archivio, poi trascinare i font nella
    sottocartella Fonts della tua cartella Library (o nel FontBook, fare doppio click su di
    essi e premere il pulsante "installa").


Il comando show (mostra) per la visualizzazione di oggetti 3D non funziona.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

La visualizzazione 3D in tempo reale per Sage dalla versione 6.4 in avanti usa il pacchetto `Jmol/JSmol <http://jmol.sourceforge.net>`_. Dalla linea di comando viene utilizzata l'applicazione Java Jmol, mentre per la visualizzazione dal browser vengono usati sia javascript puro che una Java applet. In genere nei browser e' usato javascript puro per evitare problemi con quei browser che non supportano i plugin per le applet Java (ad esempio Chrome). In ogni worksheet su browser c'e' una casella di spunta da spuntare prima di generare una vista tridimensionale qualora l'utente voglia usare l'applet Java (essa e' un po' piu' veloce con viste complicate).

La ragione piu' probabile di un malfunzionamento e' che non hai installato l'ambiente runtime di Java (JRE) o che e' piu' vecchio della versione 1.7. Se le cose funzionano dalla riga di comando, un'altra possibilita' e' che il tuo browser non ha il plugin giusto per supportare le Java applet (al momento, nel 2014, tali plugin non lavorano con la maggior parte delle versioni di Chrome). Assicurati di aver installato o il plugin IcedTea (su Linux vedi il tuo gestore dei pacchetti) o il plugin di Oracle Java (vedi: `IcedTea <http://icedtea.classpath.org/wiki/IcedTea-Web>`_ e `Java <https://java.com/en/download/help/index_installing.xml>`_).

Se stai usando un server Sage sul web ed anche la visualizzazione tramite javascript non funziona, potresti avere un problema con la funzionalita' javascript del tuo browser, o potresti aver disabilitato javascript.


Posso usare gli strumenti di Sage in un ambiente commerciale?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Si'! Assolutamente ! Fondamentalmente l'unico *limite* che hai e' che se fai dei delle modifiche a Sage stesso e redistribuisci pubblicamente tale versione modificata, allora devi renderci disponibili tali modifiche cosi' che le possiamo includere nella versione standard di Sage (se vogliamo). Altrimenti sei libero di usare quante copie di Sage vuoi per fare soldi, ecc. senza pagare alcuna licenza.


Voglio scrivere del codice Cython che usa l'aritmetica dei campi finiti, ma l'istruzione "cimport sage.rings.finite_field_givaro" non funziona. Cosa posso fare?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Devi segnalare a Sage di usare C++ (sia Givaro che NTL sono librerie C++) ed hai bisogno anche delle librerie GMP e STDC C++. Ecco un piccolo esempio::

    # These comments are hints to Sage/Pyrex about the compiler and
    # libraries needed for the Givaro library:
    #
    #clang c++
    #clib givaro gmpxx gmp m stdc++
    cimport sage.rings.finite_field_givaro
    # Construct a finite field of order 11.
    cdef sage.rings.finite_field_givaro.FiniteField_givaro K
    K = sage.rings.finite_field_givaro.FiniteField_givaro(11)
    print "K is a", type(K)
    print "K cardinality =", K.cardinality()
    # Construct two values in the field:
    cdef sage.rings.finite_field_givaro.FiniteField_givaroElement x
    cdef sage.rings.finite_field_givaro.FiniteField_givaroElement y
    x = K(3)
    y = K(6)
    print "x is a", type(x)
    print "x =", x
    print "y =", y
    print "x has multiplicative order =", x.multiplicative_order()
    print "y has multiplicative order =", y.multiplicative_order()
    print "x*y =", x*y
    # Show that x behaves like a finite field element:
    for i in range(1, x.multiplicative_order() + 1):
        print i, x**i
    assert x*(1/x) == K.one_element()

Per saperne di piu' digita quanto segue al prompt di Sage ::

    sage.rings.finite_field_givaro.FiniteField_givaro.

Poi premi TAB, ed usa ``??`` per avere piu' informationi su ogni funzione. Ad esempio::

    sage.rings.finite_field_givaro.FiniteField_givaro.one_element??

fornisce informazioni sull'unita' moltiplicativa nel campo finito.


La compilazione su Mac OS X fallisce in maniera incomprensibile. Come posso risolvere?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Cerca il file di log della compilazione (install.log) e controlla se c'e' il seguente messaggio::

    fork: Resource temporarily unavailable.

Se e' cosi', prova a fare questo: crea (o modifica se c'e' gia') il file ``/etc/launchd.conf`` includendovi quanto segue::

    limit maxproc 512 2048

Poi riavvia. Vedi `il seguente link <http://www.macosxhints.com/article.php?story=20050709233920660>`_
per maggiori dettagli.


Come va utilizzato in Sage l'operatore XOR bitwise?
"""""""""""""""""""""""""""""""""""""""""""""""""""

L'OR esclusivo in Sage si fa con l'operatore ``^^``. C'e' anche il corrispondente "operatore inplace" ``^^=``. Ad esempio::

   sage: 3^^2
   1
   sage: a = 2
   sage: a ^^= 8
   sage: a
   10

Se definisci 2 variabili e poi confronti::

    sage: a = 5; b = 8
    sage: a.__xor__(b), 13
    (13, 13)

Puoi anche fare::

    sage: (5).__xor__(8)
    13

Le parentesi sono necessarie affinche' Sage non supponga di avere a che fare con un numero reale. Ci sono molti modi di definire una funzione::

    sage: xor = lambda x, y: x.__xor__(y)
    sage: xor(3, 8)
    11

Un'altra possibilita', che aggira il preparser di Sage, e' ::

    sage: def xor(a, b):
    ...       return eval("%s^%s" % (a, b))
    ...
    sage: xor(3, 8)
    11

Puoi anche disattivare il preparser di Sage con il comando ``preparser(False)``, a quel punto l'operatore ``^`` funzionera' esattamente come in Python. Puoi successivamente riattivare il preparser con il comando ``preparser(True)``. Questo funziona solo dalla riga di comando di Sage. Se sei in una sessione Notebook, passa in "Python mode".


Quando provo ad usare LaTeX in una sessione Notebook, dice che non trova "fullpage.sty".
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

La risposta piu' ampia, ma forse non la piu' utile, e' che hai bisogno di installare ``fullpage.sty`` in una cartella acceduta da TeX. Su Ubuntu (e probabilmente anche su altre distribuzioni Linux) dovresti installare il pacchetto ``texlive-latex-extra``. Se non e' disponibile, prova ad installare il pacchetto ``tetex-extra``. Se stai usando Mac OS X dovrai usare qualunque distribuzione TeX hai gia' per ottenere ``fullpage.sty`` (se usi MacTeX probabilmente ce l'hai gia' installato). Se stai usando l'immagine VirtualBox in Windows dovrai fare login in tale immagine ed di li' installare ``texlive-latex-extra``.


Con degli oggetti "a" e "b" ed una funzione "f" ho digitato accidentalmente "f(a)=b" anziche' "f(a)==b". Questo mi ha dato un errore "TypeError" (come mi aspettavo) ma ha anche cancellato l'oggetto "a". Perche' ?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Questo e' dovuto a come sono definite le funzioni in Sage con la notazione ``f(x)=expr`` usando il preparser. Nota anche che se fai quest'errore in un costrutto ``if``, avrai un errore ``SyntaxError`` prima di qualunque altro comportamento errato, quindi, in questo caso, non hai il problema.

