.. _chapter-packaging:

====================================
Pacchettizzare codice di terze parti
====================================

Uno degli slogan del progetto Sage \`e di non reinventare la ruota: se
un algoritmo \`e gi\`a stato implementato in una libraria provata allora
si considera come incorporarla in Sage. La lista corrente dei pacchetti
disponibili \`e data dalle subdirectory di ``SAGE_ROOT/build/pkgs/``.
La gestione dei pacchetti \`e fatta attraverso uno script bash posto in
``SAGE_ROOT/local/bin/sage-spkg``. Questo script \`e di solito invocato 
dando il comando::

    [user@localhost]$ sage -i <options> <package name>...

le opzioni possono essere:

- f: installare un pacchetto anche se la stessa versione \`e gi\`a installata
- s: non cancellare la directory temporanea di compilazione
- c: dopo l'installazione, lancia la suite di test per spkg. Questo dovrebbe
  forzare sulle impostazioni di ``SAGE_CHECK`` e ``SAGE_CHECK_PACKAGES``.
- d: scarica solo il pacchetto

Non tutti i pacchetti sono compilati di default, essi sono divisi in standard,
opzionali e sperimentali. I pacchetti standard sono compilati di default e
hanno requisiti di qualit\`a molto pi\`u stringenti.

La sezione :ref:`section-directory-structure` descrive la struttura di ciascun
pacchetto in ``SAGE_ROOT/build/pkgs``. Nella sezione :ref:`section-manual-build`
vediamo come puoi installare e fare test di un nuovo spkg che tu o qualcun'altro
ha scritto. Infine, :ref:`section-inclusion-procedure` spiega come sottomettere
la richiesta che un nuovo pacchetto sia incluso nel sorgente di Sage.


.. _section-directory-structure:

Struttura delle directory
=========================

I pacchetti di terze parti in Sage consistono di 2 parti: 

#. Il tarball cos\`i com'\`e distribuito dalla terza parte, o quanto pi\`u
   possibile simile. Valide ragioni per modificarlo sono l'eleimanre file
   non necessari per ridurre la dimensione del download o rigenerate file
   auto-generati se necessario. Ma il codice effettivo non dev'essere
   modificato. Vedi anche :ref:`section-spkg-src`.

#. Gli script di compilazione ed i file associati sono nella subdirectory
   ``SAGE_ROOT/build/pkgs/package``, dove devi sostituire ``package`` con
   una versione in caratteri minuscoli del nome pubblico del progetto. 

Come esempio, consideriamo un ipotetico progetto FoO. Essi (upstream) 
distribuiscono un tarball ``foo-1.3.tar.gz`` (che sar\`a automaticamente 
messo in ``SAGE_ROOT/upstream`` dal processo di installazione). Per 
pacchettizzarlo in Sage, creiamo una subdirectory contenente::

    SAGE_ROOT/build/pkgs/foo
    |-- patches
    |   |-- bar.patch
    |   `-- baz.patch
    |-- checksums.ini
    |-- package-version.txt
    |-- spkg-check
    |-- spkg-install
    |-- spkg-src
    `-- SPKG.txt

Quando si installa Sage questi file sono usati per correggere il tarball e 
lanciare la compilazione ed il processo di installazione del pacchetto.

Discutiamo i singoli file nel seguito.


.. _section-spkg-install:

Script di installazione
-----------------------

Il file ``spkg-install`` \`e uno script di shell che installa il pacchetto,
dove ``PACKAGE_NAME`` \`e sostituito dal nome del pacchetto stesso. Nel migliore 
dei casi, il progetto a monte (upstream) pu\`o semplicemente essere installato 
dai soliti comandi ``configure / make / make install``. In tal caso, la script
di compilazione consisterebbe semplicemente di::

    #!/usr/bin/env bash

    cd src

    ./configure --prefix="$SAGE_LOCAL" --libdir="$SAGE_LOCAL/lib"
    if [ $? -ne 0 ]; then
        echo >&2 "Error configuring PACKAGE_NAME."
        exit 1
    fi

    $MAKE
    if [ $? -ne 0 ]; then
        echo >&2 "Error building PACKAGE_NAME."
        exit 1
    fi

    $MAKE -j1 install
    if [ $? -ne 0 ]; then
        echo >&2 "Error installing PACKAGE_NAME."
        exit 1
    fi


Nota che la directory radice nel tarball \`e rinominata ``src`` 
prima di invocare lo script ``spkg-install``, cos\`i che ti basta usare
``cd src`` invece di ``cd foo-1.3``.

Se c'\`e della documentazione buona inclusa ma non installata da 
``make install``, allora puoi aggiungere qualcosa come questo per 
installarla::

    if [ "$SAGE_SPKG_INSTALL_DOCS" = yes ] ; then
        $MAKE doc
        if [ $? -ne 0 ]; then
            echo >&2 "Error building PACKAGE_NAME docs."
            exit 1
        fi
        mkdir -p "$SAGE_LOCAL/share/doc/PACKAGE_NAME"
        cp -R doc/* "$SAGE_ROOT/local/share/doc/PACKAGE_NAME"
    fi
    



.. _section-spkg-check:

Autoverifiche (Self-Test)
-------------------------

Lo script ``spkg-check`` \`e opzionale, ma fortemente raccommandato, per 
lanciare le autoverifiche del pacchetto. \`e lanciato dopo la compilazione 
e l'installazione se la variabile d'ambiente ``SAGE_CHECK`` \`e impostata, 
vedi la guida all'installazione di Sage. Idealmente, il pacchetto a monte 
(upstream) avr\`a qualche sorta di insieme di test da lanciare con lo 
standard ``make check``. In tal caso, lo script ``spkg-check`` conterr\`a::

    #!/usr/bin/env bash

    cd src
    $MAKE check


.. _section-spkg-versioning:

Gestire le ersioni dei pacchetti
--------------------------------

Il file ``package-version.txt`` contiene solo la versione. Quindi se il 
pacchetto a monte (upstream) \`e ``foo-1.3.tar.gz`` allora il file di 
versione del pacchetto conterr\`a solo ``1.3``.

Se il pacchetto a monte (upstream) \`e preso da qualche revisione non stabile, 
dovresti usare la data a cui la revisione \`e stata fatta, ad esempio il
pacchetto Singular ``20090818`` \`e stato fatto con una revisione del 
2009-08-18. 

Se hai fatto qualche cambiamento al tarball a monte (upstream, vedi 
:ref:`section-directory-structure` per le modifiche possibili) allora dovresti 
aggiungere l'appendice ``.p1`` in coda alla versione. Se fai ulteriori modifiche, 
aumenta il livello della patch quanto necessario. Cos\`i le differenti versioni 
saranno ``1.3``, ``1.3.p1``, ``1.3.p2``, ...


.. _section-spkg-SPKG-txt:

Il file SPKG.txt
----------------

Il file ``SPKG.txt`` deve seguire questo pattern::

     = PACKAGE_NAME =

     == Description ==

     Cosa fa il pacchetto?

     == License ==

     Qual'\`e la licenza? Se non-standard, \`e compatibile con la GPLv3+ ?

     == Manutentori SPKG ==

     * Mary Smith
     * Bill Jones
     * Leonhard Euler

     == Upstream Contact ==

     Fornisci informazioni per contattare il progetto upstream.

     == Dipendenze ==

     Metti un elenco puntato di dipendenze:

     * python
     * readline

     == Special Update/Build Instructions ==

     Se il tarball \`e stato modificati a mano e non con uno script spkg-src, 
     descrivi cosa \`e stato cambiato.


dove al posto di ``PACKAGE_NAME`` v'\`e il nome del pacchetto. I vecchi file 
``SPKG.txt`` hanno una sezione "changelog" addizionale, ma queste informazioni 
ora sono all'interno del repository git.


.. _section-spkg-patching:

Fare le patch dei sorgenti
--------------------------

I cambiamenti in corso del codice sorgente vanno fatti con delle patch, 
che andrebbero poste nella directory ``patches``. GNU patch \`e distribuito 
con Sage, quindi puoi star certo che \`e disponibile. Le patch devono includere 
la documentazione nella loro intestazione (prima del primo blocco di differenze), 
cos\`i che un tipico patch file appaia come segue::

    Add autodoc_builtin_argspec config option

    Following the title line you can add a multi-line description of
    what the patch does, where you got it from if you did not write it
    yourself, if they are platform specific, if they should be pushed
    upstream, etc...
  
    diff -dru Sphinx-1.2.2/sphinx/ext/autodoc.py.orig Sphinx-1.2.2/sphinx/ext/autodoc.py
    --- Sphinx-1.2.2/sphinx/ext/autodoc.py.orig  2014-03-02 20:38:09.000000000 +1300
    +++ Sphinx-1.2.2/sphinx/ext/autodoc.py  2014-10-19 23:02:09.000000000 +1300
    @@ -1452,6 +1462,7 @@
 
         app.add_config_value('autoclass_content', 'class', True)
         app.add_config_value('autodoc_member_order', 'alphabetic', True)
    +    app.add_config_value('autodoc_builtin_argspec', None, True)
         app.add_config_value('autodoc_default_flags', [], True)
         app.add_config_value('autodoc_docstring_signature', True, True)
         app.add_event('autodoc-process-docstring')

Le patch ai file in ``src/`` devono essere applicate in ``spkg-install``,
cio\`e, se ci sono delle patch allora il tuo script ``spkg-install`` deve 
contenere una sezione come questa::

    for patch in ../patches/*.patch; do
        [ -r "$patch" ] || continue  # Skip non-existing or non-readable patches
        patch -p1 <"$patch"
        if [ $? -ne 0 ]; then
            echo >&2 "Error applying '$patch'"
            exit 1
        fi
    done

che applica le patch ai sorgenti.


.. _section-spkg-src:

Tarball modificati
------------------

Il file ``spkg-src`` \`e opzionale e solo per documentare come il tarball 
upstream \`e stato cambiato. Idealmente, se non \`e modificato, allora non 
c'\`e neppure un file ``spkg-src``.

Comunque se devi proprio modificare il tarball upstream allora si raccomanda 
di scrivere uno script, detto ``spkg-src``, che faccia le modifiche. 
Questo non serve solo come documentazione ma anche rende pi\`u semplice 
applicare le stesse modifiche a future versioni.


Codici di controllo (Checksum)
------------------------------

Il file ``checksums.ini`` contiene le checksum del tarball upstream.
\`E autogenerato, quindi ti basta mettere il tarball upstream nella 
directory ``SAGE_ROOT/upstream/`` e lanciare::

    [user@localhost]$ sage -sh sage-fix-pkg-checksums


.. _section-manual-build:

Compilazione ed installazione manuale dei pacchetti
===================================================

A questo punto a un nuovo tarball che non \`e ancora distribuito con Sage 
(``foo-1.3.tar.gz`` nell'esempio della sezione :ref:`section-directory-structure`).
Ora hai bisogno di collocarlo manualmente nella directory ``SAGE_ROOT/upstream/``.
Poi puoi lanciare l'installazione con::

    [user@localhost]$ sage -i package_name

oppure::

    [user@localhost]$ sage -i -f package_name

per forzare una reinstallazione. Se il tuo pacchetto contiene uno script 
``spkg-check`` (vedi :ref:`section-spkg-check`) esso pu\`o essere lanciato con::

    [user@localhost]$ sage -i -c package_name

Se \`e andato tutto bene, apri un ticket, metti un link al tarball originale nel
ticket e fai upload del ramo con il codice sorgente sotto ``SAGE_ROOT/build/pkgs``.


.. _section-inclusion-procedure:

Procedura di inclusione per paccheti nuovi ed aggiornati
========================================================

I pacchetti che non sono parte di Sage diverranno prima opzionali o 
sperimentali (quest'ultimo se non compilano su tutti i sistemi supportati.
Dopo essere stati fra gli opzionali per qualche tempo senza aver dato 
problemi si potr\`a proporre di includerli come pacchetti standard in Sage.

Per proporre un pacchetto per l'inclusione come opzionale/sperimentale aprire un 
ticket Trac con il campo ``Component:`` impostato a ``packages:experimental`` 
oppure ``packages:optional``. I requisiti associati per il codice sono 
descritti nelle sezioni seguenti.

Dopo che \`e stata fatta la revisione del ticket ed \`e stato incluso, 
i pacchetti opzionali restano in tale status per almeno un anno, dopo il 
quale si pu\`o proporre di includerli come pacchetti standard in Sage.
Per fare ci\`o si apre un ticket Trac con il campo ``Component:`` impostato
a ``packages:standard``. Nota che lo script in ``SAGE_ROOT/build/deps`` \`e
richiamato quando si compila Sage quindi includi l\`i il comando build per 
il tuo pacchetto standard. Poi fai una proposta nel Google Group ``sage-devel``.

Aggiornare dei pacchetti a nuove versioni upstream o con patch addizionali 
include la necessit\`a di aprire un ticket nella rispettiva categoria, come 
descritto sopra.

Informazioni sulla licenza
--------------------------

Se stai facendo una patch di un spkg standard di Sage, allora accertati 
che le informazioni sulla licenza di quel pachetto sia aggiornate, sia nel 
suo file ``SPKG.txt`` che nel file ``SAGE_ROOT/COPYING.txt``. Ad esempio, 
se stai facendo un file spkg che aggiorna il codice base ("vanilla" cio\`e non
modificato) ad una nuova versione, verifica che la licenza non sia cambiata 
fra le versioni.

Prerequisiti per nuovi pacchetti standard
-----------------------------------------

Perch\`e un pacchetto possa diventare parte della distribuzione standard di 
Sage, deve soddisfare i seguenti requisiti:

- **Licenza**. Per pacchetti standard, la licenza dev'essere compatibile
  con la GNU General Public License, versione 3. La Free Software
  Foundation ha stilato una lunga lista di `licenze e commenti su di 
  esse <http://www.gnu.org/licenses/license-list.html>`_.

- **Supporto alla compilazione**. Il codice deve compilare su tutte le 
  `piattaforme pienamente supportate <http://wiki.sagemath.org/SupportedPlatforms#Fully_supported>`_.

  Un pacchetto standard deve anche funzionare su tutte le piattaforme su cui
  ci si aspetta che `Sage funzioni <http://wiki.sagemath.org/SupportedPlatforms#Expected_to_work>`_
  e su cui Sage `funziona abbastanza <http://wiki.sagemath.org/SupportedPlatforms#Almost_works>`_
  ma poich\`e non supportiamo pienamente tali piattaforme e spesso manchiamo
  delle risorse per farci dei test, non ci aspettiamo che tu confermi che i
  tuoi pacchetti funzionino su di esse.

- **Qualit\`a**. Il codice dovrebbe essere "migliore" di ogni altro codice
  disponibile (che i 2 suddetti criteri), e gli autori devono giustificarlo.
  Il paragone dev'essere fatto sia per Python che per altro software.
  I criteri per passare il test di qualit\`a includono:

  - Velocit\`a

  - Documentazione

  - Usabilit\`a

  - Assenza di memory leaks

  - Manutenibilit\`a

  - Portabilit\`a

  - Ragionevole tempo di compilazione, dimensioni e dipendenze

- **Precedentemente pacchetto opzionale**. Un nuovo pacchetto standard deve avere
  trascorso del tempo come pacchetto opzionale. O ci devono essere delle buone 
  ragioni per cui ci\`o non \`e possibile.

- **Arbitraggio**. Il codicee dev'essere giudicato, come discusso in :ref:`chapter-sage-trac`.


