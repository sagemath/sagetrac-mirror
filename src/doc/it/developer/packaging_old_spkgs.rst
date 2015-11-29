.. _chapter-old-spkg:

=================================
Pacchettizzare SPKG vecchio stile
=================================

Questo capitolo spiega gli spkg vecchio stile; si applica solo agli 
spkgs opzionali legacy ed agli spkg sperimentali. Vedi :ref:`chapter-packaging`
per un modo moderno di pacchettizzare il software di terze parti.


Creare un SPKG vecchio stile
============================

Se stai producendo codice per aggiungere nuove funzionalit\`a a Sage, 
potresti considerare di farne un pacchetto (un "spkg") anzich\`e un patch
file. Se il tuo codice \`e molto lungo (ad esempio) e dev'essere offerto
come download opzionale, un pacchetto \`e la scelta giusta. Similmente, se
il tuo codicee dipende da qualche altro componente opzionale di Sage, dovresti
produrre un pacchetto. Nel dubbio, chiedi consiglio sulla mailing-list ``sage-devel``.

L'abbreviazione "spkg" sta per "Sage package". La directory ``SAGE_ROOT/spkg/standard``
contiene degli spkg. In una installazione da codice sorgente, questo \`e 
costituito da file spkg di Sage (in realt\`a file ``.tar`` o ``.tar.bz2``, 
che sono il codice sorgente che definisce Sage. In una istallazione da file 
binario, alcuni di questi posso essere dei file segnaposto per risparmiare spazio.

I pacchetti Sage sono distributiti come file ``.spkg``, che sono file
``.tar.bz2`` (o file ``tar``) ma hanno l'estensione ``.spkg`` per evitare 
confusioni. Sebbene i pacchetti Sage siano fatti con ``tar`` e/o bzip2, nota
che i file ``.spkg`` contengono informationi di controllo (script di 
installazione e metadata) che sono necessari per compilarli ed installarli.
Quando compili Sage da codice sorgente (o quando lanci ``sage -i <pkg>`` o 
``sage -f <pkg>``), il file ``SAGE_ROOT/src/bin/sage-spkg`` provvede a spacchettare,
compilare, ed installare i pacchetti Sage. Puoi digitare::

    tar -jxvf mypackage-version.spkg

per estrarre un spkg e vedere com'\`e fatto dentro. Se vuoi creare un
nuovo pacchetto Sage, ti raccomandiamo di iniziare con l'esaminarne uno
esistente. All'URL http://www.sagemath.org/download-packages.html vi
\`e una lista di spkg disponibili per il download.


Nominare il tuo SPKG
--------------------

Ogni pacchetto spkg di Sage ha un nome della seguente forma::

   BASENAME-VERSION.spkg

``BASENAME`` \`e il nome del pacchetto; pu\`o contenere lettere 
minuscole, numeri, e caratteri di sottolineatura, ma non trattini.
``VERSION`` \`e il numero di versione; dovrebbe iniziare con un numero 
e pu\`o contenere numeri, lettere, punti, e trattini; pu\`o terminare
con una stringa della forma "pNUM", dove "NUM" \`e un intero non negativo.
Se il tuo spkg \`e una versione "vanilla" (cio\`e non modificata) di un 
qualche software, ad esempio la versione 5.3 di "my-python-package", 
allora ``BASENAME`` dovrebbe essere "my_python_package" --
nota il cambiamento dei trattini in sottolineature, poich\`e ``BASENAME``
non deve contenere alcun trattino -- e ``VERSION`` essere "5.3".  Se hai
bisogno di modificare il software per usarlo con Sage (come descritto
sotto e nel capitolo :ref:`section-old-spkg-patching-overview`),
allora ``VERSION`` dovrebbe esseree "5.3.p0", dove la "p0" indica un
livello di patch a 0. Se qualcuno aggiunge altre patch, successivamente, 
questo diverrebbe "p1", poi "p2", ecc.

La stringa ``VERSION`` dev'essere presente. Se stai usando un software 
senza un numero di versione ovvio, usa una data. Per dare al tuo spkg un
nome come questo, crea una directory chiamata ``BASENAME-VERSION`` e
metti i tuoi file in quella directory -- la prossima sezione descrive la
struttura della directory.


Struttura della directory
-------------------------

Metti i tuoi file in una directory di nome simile a ``mypackage-0.1``, come
descritto sopra. Se stai facendo il porting di un altro pacchetto software, 
allora la directory dovrebbe contenere una sottodirectory ``src/``, in cui 
vi sia una copia non alterata del pacchetto. Ogni file non in ``src/`` 
dovrebbe essere sotto controllo versione, cio\`e ne dovrebbe essere stato 
fatto il check in un repository hg.

Pi\`u precisamente, la directory dovrebbe contenere i seguenti:

- ``src/``: questa directory contiene codice vanilla upstream, con poche 
  eccezioni, ad esempio quando lo spkg incluso in Sage \`e in effetti
  upstream, e lo sviluppo su quel codice base sta avvenendo in stretto 
  contatto e coordinamento con Sage.  Vedi il pacchetto eclib spkg di John 
  Cremona, ad esempio. La directory ``src/`` non dev'essere sotto controllo
  di revisione.

- ``.hg``, ``.hgignore``, e ``.hgtags``: gli spkg vecchi stile usano
  Mercurial come sistema di controllo di revisione. La directory nascosta
  ``.hg`` \`e parte del layout standard di un spkg di Sage. Contiene il
  repository Mercurial per tutti i file non nella directory ``src/``.
  Per creare questo repository Mercurial da zero, puoi fare::

      hg init

  I file ``.hgignore`` e ``.hgtags`` appartengono anche al repository Mercurial. 
  Il file ``.hgtags`` \`e opzionale, ed \`e frequentemente omesso. Dovresti 
  accertarti che il file ``.hgignore`` contenga "src/", dal momento che non stiamo
  tracciando il suo contenuto. Invero, frequentemente questo file contiene solo 
  una singola linea::

      src/

- ``spkg-install``: questo file contiene lo script di installazione. Vedi 
  :ref:`section-old-spkg-install` per maggiori informazioni ed un template.

- ``SPKG.txt``: questo file descrive lo spkg in formato wiki. Ogni nuova 
  revisione necessita di una voce "changelog" aggiornata o lo spkg ricever\`a 
  automaticamente un "needs work" alla revisione. Vedi 
  :ref:`section-old-spkg-SPKG-txt` per un template.

- ``spkg-check``: questo file lancia la suite di test. Questo \`e in un certo
  qual modo opzionale poich\`e non tutti gli spkg hanno delle suite di test.
  Se possibile, accertarsi di creare uno script del genere poich\`e aiuta ad
  isolare i bachi nei pacchetti upstream.

- ``patches/``: questa directory contiene delle patch ai file sorgentge in
  ``src/``. Vedi :ref:`section-old-spkg-patching-overview`. Le patch ai
  file in ``src/`` dovrebbero essere applicate in ``spkg-install``, e tutte le
  patch devono essere auto-documentanti, cio\`e la testata deve contenere cosa
  fanno, se sono specifiche di una piattaforma, se ne deve essere fatta la push
  upstream, ecc. Per assicurarsi che tutte le versioni patch-ate di file sorgenti
  upstream sotto ``src/`` siano sotto controllo di revisione, l'intera directory
  ``patches/`` dev'essere sotto controllo di revisione.

**Mai** applicare patch a file sorgente upstream sotto ``src/`` e poi impacchettare
un spkg. Un tale miscuglio di codice sorgente upstream con le versione patch-ate 
specifiche di Sage \`e una ricetta per la confusione. Dev'esserci una
**netta separazione** fra il sorgente fornito dal progetto upstream e le versioni
patch-ate che il progetto Sage genera sulla base di quanto fornisce il sorgente upstream.

La sola eccezione a questa regola \`e la *rimozione* di file o directory inutilizzate.
Alcuni pacchetti contengono parti che non sono necessarie a Sage. Per risparmiare spazio,
queste possono essere rimosse direttamente da ``src/``.
Ma accertati di documentarlo nella sezione "Special Update/Build Instructions" 
in ``SPKG.txt``!


.. _section-old-spkg-install:

Il File spkg-install
--------------------

Lo script ``spkg-install`` \`e lanciato durante l'installazione del pacchetto Sage.
In questo script, puoi presumere le seguenti cose:

- PATH contiene, in cima, le locazioni di ``sage`` e ``python`` (dall'installazione
  di Sage). Pertanto il comando::

      python setup.py install

  eseguir\`a la versione corretta di Python con tutto quanto impostato correttamente.
  Inoltre, lanciando ``gap`` o ``Singular``, ad esempio, sar\`a eseguita la versione
  corretta.

- La variabile d'ambiente ``SAGE_ROOT`` punta alla directory radice dell'installazione
  di Sage.

- La variabile d'ambiente ``SAGE_LOCAL`` punta alla directory ``SAGE_ROOT/local`` 
  dell'installazione di Sage.

- Le variabili d'ambiente ``LD_LIBRARY_PATH`` e ``DYLD_LIBRARY_PATH`` hanno entrambe 
  ``SAGE_ROOT/local/lib`` in cima.

Lo script ``spkg-install`` dovrebbe copiare i tuoi file nel posto giusto dopo aver fatto
qualunque compilazione fosse necessaria. Questo \`e un template::

    #!/usr/bin/env bash

    if [ -z "$SAGE_LOCAL" ]; then
        echo >&2 "SAGE_LOCAL undefined ... exiting"
        echo >&2 "Maybe run 'sage --sh'?"
        exit 1
    fi

    cd src

    # Apply patches.  See SPKG.txt for information about what each patch
    # does.
    for patch in ../patches/*.patch; do
        [ -r "$patch" ] || continue  # Skip non-existing or non-readable patches
        patch -p1 <"$patch"
        if [ $? -ne 0 ]; then
            echo >&2 "Error applying '$patch'"
            exit 1
        fi
    done

    ./configure --prefix="$SAGE_LOCAL"
    if [ $? -ne 0 ]; then
        echo >&2 "Error configuring PACKAGE_NAME."
        exit 1
    fi

    $MAKE
    if [ $? -ne 0 ]; then
        echo >&2 "Error building PACKAGE_NAME."
        exit 1
    fi

    $MAKE install
    if [ $? -ne 0 ]; then
        echo >&2 "Error installing PACKAGE_NAME."
        exit 1
    fi

    if [ "$SAGE_SPKG_INSTALL_DOCS" = yes ] ; then
        # Before trying to build the documentation, check if any
        # needed programs are present. In the example below, we
        # check for 'latex', but this will depend on the package.
        # Some packages may need no extra tools installed, others
        # may require some.  We use 'command -v' for testing this,
        # and not 'which' since 'which' is not portable, whereas
        # 'command -v' is defined by POSIX.

        # if [ `command -v latex` ] ; then
        #    echo "Good, latex was found, so building the documentation"
        # else
        #    echo "Sorry, can't build the documentation for PACKAGE_NAME as latex is not installed"
        #    exit 1
        # fi


        # make the documentation in a package-specific way
        # for example, we might have
        # cd doc
        # $MAKE html

        if [ $? -ne 0 ]; then
            echo >&2 "Error building PACKAGE_NAME docs."
            exit 1
        fi
        mkdir -p "$SAGE_ROOT/local/share/doc/PACKAGE_NAME"
        # assuming the docs are in doc/*
        cp -R doc/* "$SAGE_ROOT/local/share/doc/PACKAGE_NAME"
    fi


Nota che la prima linea \`e ``#!/usr/bin/env bash``; questo \`e importante
per la portabilit\`a. Poi, lo script verifica che ``SAGE_LOCAL`` sia
definita per accertarsi che l'ambiente di Sage sia stato sistemato. Dopo 
questo, lo script potrebbe semplicemnte eseguire ``cd src`` e poi invocare o
``python setup.py install`` o la sequenza degli autotools
``./configure && make && make install``, o qualcos'altro di simile.

A volte, per\`o, pu\`o essere pi\`u complicato. Ad esempio, potresti aver bisogno 
di applicare delle patch dalla directory ``patches`` in un ordine particolare. Inoltre,
dovresti prima compilare (ad esempio con ``python setup.py build``, uscendo se c'\`e
un errore), prima di installare (ad esempio con ``python setup.py
install``). In questo modo, non sovrascriverai una vecchia versione del spkg
funzionante con una nuova non funzionante.

Quando copi la documentazione in
``$SAGE_ROOT/local/share/doc/PACKAGE_NAME``, pu\`o essere necessario verificare
che solo i file di documentazione pensati per l'utente finale sono copiati.
Ad esempio, se la documentazione \`e compilata dai file ``.tex``, ti potrebbe 
bastare copiare i file pdf risultanti, anzich\`e l'intera directory "doc".
Quando si genera documentazione usando Sphinx, il copiare la directory ``build/html`` 
in generale copier\`a solo l'output inteso per l'utente finale.


.. _section-old-spkg-SPKG-txt:

Il file SPKG.txt
----------------

Il file ``SPKG.txt`` vecchio stile \`e lo stesso descritto in
:ref:`section-spkg-SPKG-txt`, ma con un "changelog" manuale aggiunto, 
poich\`e i contenuti non sono parte del repository di Sage.
Dovrebbe seguire il pattern seguente::

     == Changelog ==

     Mettere qui un changelog del spkg, dove le entrate hanno questo formato:

     === mypackage-0.1.p0 (Mary Smith, 1 Jan 2012) ===

      * Patch src/configure so it builds on Solaris. See Sage trac #137.

     === mypackage-0.1 (Leonhard Euler, 17 September 1783) ===

      * Initial release.  See Sage trac #007.

Quando la directory (diciamo, ``mypackage-0.1``) \`e pronta, il comando

::

    sage --pkg mypackage-0.1

creer\`a il file ``mypackage-0.1.spkg``. Come notato sopra, questo
crea un tar file compresso. Eseguendo ``sage --pkg_nc mypackage-0.1``
si crea un tar file non compresso.

Quando il tuo spkg \`e pronto, dovresti farne una segnalazione su ``sage-devel``.
Se la gente l\`i pensa che sia una buona idea, allora fa un post del link al spkg
sul Trac server di Sage (vedi :ref:`chapter-sage-trac`) cos\`i che possa essere 
giudicato. Non fare un post del spkg stesso sul server Trac: ti basta fornire un
link al tuo spkg.  Se il tuo spkg ottiene una revisione positiva, potr\`a essere
incluso nella libreria core di Sage, o potr\`a diventare un download opzionale 
dal sito web di Sage, cos\`i che ciunque possa installarlo automaticamente 
digitando ``sage -i mypackage-version.spkg``.

.. note::

   Per qualunque spkg:

   - Accertati che il repository hg contenga ogni file al di fuori della
     directory ``src``, che questi siano tutti aggiornati e che ne sia stato
     fatto il commit nel repository.

   - Includi un file ``spkg-check`` se possibile (vedi `trac ticket #299`_).

   .. _trac ticket #299: http://trac.sagemath.org/sage_trac/ticket/299

.. note::

    Il codice Magma esterno va in ``SAGE_ROOT/src/ext/magma/user``, cos`\i
    se vuoi redistribuire il codice Magma con Sage come un pacchetto che gli
    utenti con Magma possano usare, l\`i \`e dove va messo. Dovresti anche 
    disporre codice Python utile a rendere il codice Magma facilmente 
    utilizzabile.


.. _section-old-spkg-avoiding-troubles:

Evitare guai
============

Questa sezione contiene alcune linee guida su cosa un spkg non deve mai fare 
ad una installazione di Sage. Sei incorraggiato a produrre un spkg che \`e tanto
indipendente dal resto quanto possibile.

#. Un spkg non deve modificare un file sorgente esistente nella libreria Sage.
#. Non permettere ad un spkg di modificare un altro spkg. Un spkg pu\`o dipendere
   da un altro spkg -- vedi sopra. Verifica l'esistenza di un spkg richiesto come 
   prerequisito prima di installare un spkg che dipende da lui.




.. _section-old-spkg-patching-overview:

Panoramica sulle patch agli SPKG 
================================

Accertati di essere familiare con la struttura e le conventioni degli 
spkg; vedi il capitolo :ref:`chapter-old-spkg` per dettagli.
Fare la patch di un spkg implica fare la patch dello script di installazione
del spkg e/o del codice sorgente upstream contenuto nel spkg.
Diciamo che vuoi fare una patch al pacchetto Matplotlib ``matplotlib-1.0.1.p0``.
Nota che la "p0" denotea il livello della patch sul spkg, mentre "1.0.1" si
riferisce alla versione upstream di Matplotlib cos\`i come contenuta in
``matplotlib-1.0.1.p0/src/``. Lo script di installazione di tale spkg \`e::

    matplotlib-1.0.1.p0/spkg-install

In generale, uno script di nome ``spkg-install`` \`e uno script di
installazione per un spkg. Per fare la patch allo script di installazione, 
usa un text editor per modificare lo script. Poi nel file di log ``SPKG.txt``
fornisci una descrizione ad alto livello delle tue modifiche. Quando sei 
soddisfatto delle tue modifiche nello script d'installazione nel file di log
``SPKG.txt``, usa Mercurial per fare il check-in delle tue modifiche ed 
accertati di fornire un messaggio di commit significativo.

La directory ``src/`` contiene il codice sorgente fornito dal progetto upstream.
Ad esempio, codice sorgente di Matplotlib 1.0.1 \`e contenuto in ::

    matplotlib-1.0.1.p0/src/

Per fare una patch al codice sorgente upstream, devi modificare una copia del file
interessato -- i file nella directory ``src/`` non dovrebbero essere toccati,
essendo versioni "vanilla" del codice sorgente. Ad esempio, puoi copiare l'intera
directory ``src/`` ::

    $ pwd
    matplotlib-1.0.1.p0
    $ cp -pR src src-patched

Poi modificare i file in ``src-patched/``. Quando sei soddisfatto delle tue modifiche, 
genererai una lista diff unificata fra il file originale e quello modificato, e la 
salverai in ``patches/``::

    $ diff -u src/configure src-patched/configure > patches/configure.patch

Salva la lista diff unificata in un file con lo stesso nome del file sorgente di cui
hai fatto la patch, ma usa l'estensione ".patch". Nota che la directory ``src/`` 
non dovrebbe essere sotto controllo revisione, laddove ``patches/`` deve esserlo.
Il file di configurazione di Mercurial ``.hgignore`` dovrebbe contenere la seguente 
linea::

    src/

Assicurati che lo script di installazione ``spkg-install`` contenga codice per
applicare le patch ai file opportuni sotto ``src/``. Ad esempio, il file ::

    matplotlib-1.0.1.p0/patches/finance.py.patch

\`e una patch per il file ::

    matplotlib-1.0.1.p0/src/lib/matplotlib/finance.py

Lo script di installazione ``matplotlib-1.0.1.p0/spkg-install`` contiene il
seguente codice per installare le patch necessarie::

    cd src

    # Apply patches.  See SPKG.txt for information about what each patch
    # does.
    for patch in ../patches/*.patch; do
        patch -p1 <"$patch"
        if [ $? -ne 0 ]; then
            echo >&2 "Error applying '$patch'"
            exit 1
        fi
    done

Naturalmente, questo pu\`o essere modificato se l'order in cui le patch
vanno applicate \`e importante, o se qualche patch \`e dipendente dalla piattaforma.
Ad esempio::

    if [ "$UNAME" = "Darwin" ]; then
        for patch in ../patches/darwin/*.patch; do
            patch -p1 <"$patch"
            if [ $? -ne 0 ]; then
                echo >&2 "Error applying '$patch'"
                exit 1
            fi
        done
    fi

(La variabile d'ambiente :envvar:`UNAME` \`e definita dallo script
``sage-env``, ed \`e disponibile quando ``spkg-install`` \`e eseguito.)

Ora fornisci una spiegazione a grandi linee delle tue modifiche in ``SPKG.txt``.
Nota il formato di ``SPKG.txt`` -- vedi il capitolo
:ref:`chapter-old-spkg` per dettagli. Quando sei soddisfatto delle tue 
modifiche, usa Mercurial per fare il check-in delle tue modifiche, dando un
messaggio di commit significativo.  Poi usa il comando ``hg tag`` per mettere un
nuovo numero di versione (usando "p1" invece di "p0": abbiamo fatto dei cambiamenti, 
quindi dobbiamo aggiornare il livello della patch)::

    $ hg tag matplotlib-1.0.1.p1

Poi rinomina la directory ``matplotlib-1.0.1.p0`` a
``matplotlib-1.0.1.p1`` per farla coincidere con il nuovo livello di patch. 
Per produrre il file spkg vero e proprio, spostati nella directory genitore di
``matplotlib-1.0.1.p1`` ed esegui ::

    $ /path/to/sage-x.y.z/sage --pkg matplotlib-1.0.1.p1
    Creating Sage package matplotlib-1.0.1.p1

    Created package matplotlib-1.0.1.p1.spkg.

        NAME: matplotlib
     VERSION: 1.0.1.p1
        SIZE: 11.8M
     HG REPO: Good
    SPKG.txt: Good

I file spkg sono o degli archivi tar compressi con bzip o tar semplici; il
comando ``sage --pkg ...`` produce la versione compressa. Se il tuo spkg
contiene per lo pi\`u file binari che si comprimono poco, puoi usare
``sage --pkg_nc ...`` per produrre una versione non compressa, cio\`e un
file tar normale::

    $ sage --pkg_nc matplotlib-1.0.1.p0/
    Creating Sage package matplotlib-1.0.1.p0/ with no compression

    Created package matplotlib-1.0.1.p0.spkg.

        NAME: matplotlib
     VERSION: 1.0.1.p0
        SIZE: 32.8M
     HG REPO: Good
    SPKG.txt: Good

Nota che questo \`e quasi un 3 volte la versione compressa, quindi dovremmo
usare la compressione!

A questo punto, potresti voler sottoporre il tuo spkg patch-ato per la revisione.
Allora fornisci un link (URL) al tuo spkg sul ticket di Trac relativo e/o in una
email alla mailing list relativa. Di solito non si dovrebbe fare l'upload del
spkg vero e proprio al ticket Trac relativo -- non inviare grandi file binari 
al server Trac.


Gestione delle versioni degli SPKG
==================================

Se vuoi aggiornare la versione di un spkg, devi seguire alcune convenzioni
di denominazione. Usa il nome ed il numero di versione com'\`e dato dal 
progetto upstream, ad esempio ``matplotlib-1.0.1``. Se il pacchetto upstream
\`e preso da qualche revisione che non \`e una versione stabile, aggiungi la
data a cui \`e stata fatta la revisione, ad esempio il pacchetto Singular
``singular-3-1-0-4-20090818.p3.spkg`` ha la revisione 2009-08-18. Se inizi 
da zero da una release upstream senza patch al suo sorgente, il spkg risultante
non ha bisogno di avere alcuna etichetta di livello di patch (si pu\`o aggiungere
".p0", ma \`e opzionale). Ad esempio, ``sagenb-0.6.spkg`` \`e preso dalla
versione stabile upstream ``sagenb-0.6`` senza alcuna patch applicata al suo
codice sorgente. Per cui non vedrai delle numerazioni di livello di patch come 
``.p0`` or ``.p1``.

Diciamo che inizi con ``matplotlib-1.0.1.p0`` e vuoi sostituire
Matplotlib 1.0.1 con la versione 1.0.2. Questo implica sostituire il codice 
sorgente di Matplotlib 1.0.1 sotto ``matplotlib-1.0.1.p0/src/`` con il
nuovo codice sorgente. Per incominciare, segui le convenzioni di denominazione 
come descritto nella sezione :ref:`section-old-spkg-patching-overview`. Se
necessario, rimuovi qualunque patch obsoleta e crea quelle nuove,
mettendole sotto la directory ``patches/``.  Modifica lo script
``spkg-install`` per prendere in considerazione qualunque cambiamento delle patch;
potresti aver a che fare con modifiche a come la nuova versione del codice compila.
Poi pacchettizza il tuo spkg sostitutivo usando le opzioni a riga di comando di 
Sage ``--pkg`` o ``--pkg_nc`` (tar con o senza bzip2).

Per installare il tuo spkg sostitutivo, usa::

    sage -f http://URL/to/package-x.y.z.spkg

oppure::

    sage -f /path/to/package-x.y.z.spkg

Per compilare Sage da sorgente con lo (standard) spkg sostituivo, esegui 
``untar`` del tarball del sorgente di Sage e rimuovi il spkg esistente da
``SAGE_ROOT/spkg/standard/``. Al suo posto metti il tuo spkg sostitutivo.
Poi esegui ``make`` da ``SAGE_ROOT``.

