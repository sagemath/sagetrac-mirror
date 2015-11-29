.. _chapter-manual-git:

================
Git the Hard Way
================

Ai principianti dello sviluppo in Sage con nessuna esperienza di Git consigliamo di leggere, invece, :ref:`chapter-git_trac`. Il comando ``git-trac`` semplifica l'interazione con i nostri server Git e Trac.

Comunque \`e possibile usare Git direttamente per lavorare con i repository remoti. Questo capitolo ti spiega come farlo supponendo qualche familiarit\`a di base con Git. In particolare devi aver letto :ref:`chapter-walkthrough`.

Supponiamo che tu abbia gi\`a una copia del repository Git di Sage, ad esempio avendo eseguito::

    [user@localhost ~]$ git clone git://github.com/sagemath/sage.git
    [user@localhost ~]$ cd sage
    [user@localhost sage]$ make

.. _section-git-trac:

Il server Trac
==============

Anche il Trac server di Sage mantiene una copia del repository diSage, essa \`e servita tramite i protocolli SSH e git. Per aggiungerlo come un repository remoto al tuo repository locale, usa i seguenti comandi::

    [user@localhost sage]$ git remote add trac git://trac.sagemath.org/sage.git -t master
    [user@localhost sage]$ git remote set-url --push trac git@trac.sagemath.org:sage.git
    [user@localhost sage]$ git remote -v
    origin      git://github.com/sagemath/sage.git (fetch)
    origin      git://github.com/sagemath/sage.git (push)
    trac        git://trac.sagemath.org/sage.git (fetch)
    trac        git@trac.sagemath.org:sage.git (push)

Invece di ``trac`` puoi usare qualunque nome locale vuoi, naturalmente. Va benissimoaver pi\`u repository remoti per git, puoi pensare ad essi come a dei segnalibri. Poi puoi usare ``git pull`` per ottenere le modifiche e ``git push`` per fare upload dei tuoi cambiamenti locali utilizzando::

    [user@localhost sage]$ git <push|pull> trac [ARGS]

.. note::
   
    Nel comando suddetto abbiamo impostato il remoto per tracciare solo
    il ramo ``master`` sul trac server (con l'opzione ``-t master``).
    Questo evita disordine evitando di scaricare tutti i rami che siano
    mai stati creati. Ma significa anche che non otterrrai tutto ci\`o
    che c'\`e sul server Trac di default, ma dovrai esplicitamente dire
    a git quale ramo vuoi ottenere da Trac. Vedi la sezione :ref:`section-git-checkout`
    per degli esempi.

Qui impostiamo il remoto per eseguire solo operazioni di lettura (fetch), utilizzando il protocollo git, mentre le operazioni di scrittura (push) sono fatte usando il protocollo SSH (specificato dalla parte ``git@...``). Devi aver un account Trac per usare il protocollo SSH e per impostare la tua chiave pubblica SSH come descritto in :ref:`section-trac-ssh-key`. l'autenticazione  \`e necessaria se vuoi fare upload di qualunque cosa, per assicurarsi che proviene realmente da te.

Se vuoi usare solo SSH, usa questi comandi::

    [user@localhost sage]$ git remote add trac git@trac.sagemath.org:sage.git -t master
    [user@localhost sage]$ git remote -v
    origin      git://github.com/sagemath/sage.git (fetch)
    origin      git://github.com/sagemath/sage.git (push)
    trac        git@trac.sagemath.org:sage.git (fetch)
    trac        git@trac.sagemath.org:sage.git (push)


.. _section-git-checkout:

Check-out dei ticket
--------------------


I ticket Trac che sono finiti o in corso di lavorazione possono avere un ramo git attaccato ad essi. Questo \`e il campo "Branch:" nella descrizione del ticket. Il nome del ramo in genere \`e nella forma ``u/user/description``, dove ``user`` \`e il nome dell'utente che ha fatto il ramo e ``description`` \`e qualche corta descrizione in forma libera ( e pu\`o includere ulterioni "//")

Se vuoi lavorare con le modifiche in tale ramo remoto, devi farne una copia in locale. In particolare Git non ha il concetto di lavorare direttamente sul ramo remoto, i remoti sono solo dei segnalibri per cose che puoi scaricare da o inviare al server remoto. Per cui la prima cosa che devi fare \`e scaricare tutto quanto dal ramo sul server Trac nel tuo ramo locale. Questo si fa con::

    [user@localhost sage]$ git fetch trac u/user/description
    remote: Counting objects: 62, done.
    remote: Compressing objects: 100% (48/48), done.
    remote: Total 48 (delta 42), reused 0 (delta 0)
    Unpacking objects: 100% (48/48), done.
    From trac.sagemath.org:sage
    * [new branch]      u/user/description -> FETCH_HEAD

Il ramo ``u/user/description`` \`e ora tmporaneamente (fintanto che non scarichi qualcos'altro) memorizzato nel tuo database Git locale sotto l'alias ``FETCH_HEAD``. Come secondo passo lo rendiamo disponibile come un nuovo ramo locale e passiamo ad esso. Il tuo ramo locale pu\`o avere un nome diverso, ad esempio::

    [user@localhost sage]$ git checkout -b my_branch FETCH_HEAD
    Switched to a new branch 'my_branch'

Ci\`o crea un nuovo ramo nel tuo repository locale di git chiamato ``my_branch`` e modifica l'albero del filesystem locale di Sage secondo lo stato dei files in quel ticket. Ora puoi modificare i files ed effettuarne il commit nel tuo ramo locale.


.. _section-git-push:

Effettuare la push dei cambiamenti per un ticket
------------------------------------------------

Per aggiungere il tuo ramo locale ad un ticket Trac, dovresti dapprima decidere un nome da usare sel reporitory del server Trac di Sage. Per evitare conflitti di nome, hai il permesso di effettuare delle push a rami della forma ``u/user/*`` dove ``user`` \`e il tuo nome utente di Trac e ``*`` \`e un qualunque nome di ramo git valido. Di default non hai permessi di effettuare push dei rami di altri utenti del ramo ``master`` di Sage. Nel seguito utilizzeremo il nome di ramo ``u/user/description`` dove si intende che tu metti:

* il tuo nome utente Trac al posto di ``user``

* una qualche descrizione del tuo ramo, breve ma adeguata, al posto di ``description``. Pu\`o contenere ulterioni ``/`` ma non degli spazi bianchi.

Il tuo primo passo dovrebbe essere di mettere il nome che hai scelto nel campo ``Branch:`` del ticket Trac. Per effettuare il push del tuo ramo al server Trac puoi utilizzare::

* se hai iniziato tu stesso il ramo e non segui nessun altro ramo::

    [user@localhost sage]$ git push --set-upstream trac HEAD:u/user/description

* se il tuo ramo ha gi\`a un precedente::

    [user@localhost sage]$ git push trac HEAD:u/user/description

L'opzione ``HEAD`` significa che stai facendo la push del commit pi\`u recente (e, di conseguenza, di tutti i commit precedenti collegati) del ramo locale corrente al ramo remoto.

Il campo ``Branch:`` nella pagina del ticket Trac ha un codice colore: rosso significa che \`e un problema, verde significa che si fonder\`a correttamente in ``master``. Se \`e rossa, il tooltip (fumetto che compare posizionandosi sopra col mouse) ti dir\`a cos'ha di sbagliato. Se \`e verde, allora sar\`a collegata tramite link ad un elenco di differenze rispetto a ``master``.

Per permessi di lettura/scrittura sui rami di git, vedi :ref:`section-git_trac-branch-names`


.. _section-git-pull:

Ottenere le modifiche
---------------------

Un'esigenza comune durante lo sviluppo \`e sincronizzare la propria copia locale del ramo con quella remota su Trac. In particolare supponi di aver scaricato il ramo di qualcun altro dopo aver dato dei suggerimenti per migliorie sul ticket su Trac relativo. Ora, l'autore originale ha incorporato i tuoi suggerimenti nel suo ramo, e tu vuoi ottenere anche tali modifiche per poter completare la tua revisione. Assumendoche tu abbia in origine ottenuto il tuo ramo locale com detto a :ref:`section-git-checkout`, ti basta eseguire::

    [user@localhost sage]$ git pull trac u/user/description
    From trac.sagemath.org:sage
     * branch            u/user/description -> FETCH_HEAD
    Updating 8237337..07152d8
    Fast-forward
     src/sage/tests/cmdline.py      | 3 ++-
     1 file changed, 2 insertions(+), 1 deletions(-)

dove qui ``user`` \`e il nome utente su Trac dell'altro developer, e ``description`` \`e la descrizione che egli ha scelto. Questo comando scaricher\`a le modifiche dal ramo remoto usato originariamente e le fonder\`a nel tuo ramo locale.Se non hai ancora pubblicato i tuoi commit locali allora puoi anche effettuarne il ``rebase`` con:

    [user@localhost sage]$ git pull -r trac u/user/description
    From trac.sagemath.org:sage
     * branch            u/user/description -> FETCH_HEAD
    First, rewinding head to replay your work on top of it...
    Applying: my local commit

Vedi la sezione :ref:`section-git-merge` per una spiegazione approfondita delle differenze tra effettuare un merge ed un rebase.

Fin qui suponiamo che non ci siano conflitti. \`E inevitabile nello sviluppo distribuito che, a volte, la stessa locazione in un file sorgente venga modificata da pi\`u di una persona. Come riconciliare tali modifiche in conflitto \`e spiegato nella sezione :ref:`section-git_trac-conflict`.


.. _section-git-pull-master:

Aggiornare Master
-----------------

Il ramo ``master`` pu\`o essere aggiornato come qualunque altro ramo. Tuttavia dovresti fare attenzione a mantenere la tua copia locale di ``master`` **identica** a quella su Trac poich\`e questa \`e la versione ufficiale corrente di Sage.

In particolare se per sbaglio hai effettuato dei commit sulla tua copia locale di ``master`` allora li devi cancellare invece di effettuarne il merge con il ramo ``master`` ufficiale.

Un modo per essere avvertito di potenziali problemi \`e utilizzare ``git pull --ff-only`` che va in errore se viene richiesto di fare un merge anomalo::

    [user@localhost sage]$ git checkout master
    [user@localhost sage]$ git pull --ff-only trac master

Se quest'operazione di pull fallisce, allora c'\`e qualcosa di sbagliato nella tua copia locale del ramo ``master``. Per passare al ramo ``master`` corretto puoi utilizzare::

    [user@localhost sage]$ git checkout master
    [user@localhost sage]$ git reset --hard trac/master


.. _section-git-merge:

Merge e rebase
==============

Di solito lo sviluppo di Sage continua mentre tu stai lavorando sul tuo ramo locale. Ad esempio supponiamo che hai iniziato il ramo "my_branch" al commit B. Dopo un po' il tuo ramo \`e avanzato al commit Z ma hai aggiornato ``master`` (vedi :ref:`section-git-pull-master`) che intanto \`e avanzato a D, cos\`i che la tua history assomiglia allo schema seguente (vedi :ref:`section_walkthrough_logs`)::

                     X---Y---Z my_branch
                    /
               A---B---C---D master

Come devi comportarti con tali cambiamenti di fondo in "master" mentre stai ancora sviluppando il tuo codice ? In linea di principio ci sono 2 modi di fare:

* la prima soluzione \`e cambiare i commit nel tuo ramo locale per iniziare con il nuovo ``master`` (**rif\`a** i commit ``X,Y,Z`` in cima al nuovo ``master``). Questo si dice fare un ``rebase``, e riscrive il tuo ramo locale::

      git checkout my_branch
      git rebase -i master

Qui supponiamo che ``master`` sia la tua copia locale aggiornata del ramo ``master``. In termini del grafo di commit, sia ha::

                             X'--Y'--Z' my_branch
                            /
               A---B---C---D master

Nota che questa operazione riscrive la history di ``my_branch`` (vedi :ref:`section-git-rewriting-history`). Questo pu\`o dare problemi se qualcuno inizia a scrivere codice in cima ai tuoi commit ``X,Y,Z``. Diversamente \`e sicura.

In alternativa puoi effettuare la pull dei cambiamenti dal server Trac e contemporaneamente effettuare il rebase con il comando "git pull -r master", cio\`e effettuare il rebase di ``my_branch`` mentre contemporaneamente fa l'update di ``master`` (vedi :ref:`section-git-pull`)::

    git checkout my_branch
    git pull -r master

Poich\`e l'hash SHA1 include l'hash del genitore, tutti i commit cambiano. Ci\`o significa che dovresti usare rebase solo quando nessun altro ha utilizzato uno dei tuoi commit X,Y,Z per basarci sopra il suo sviluppo.

* l'altra soluzione \`e non cambiare alcun commit ed invece creare un nuovo commit di merge W che acquisisca le modifiche del nuovo ``master`` (un ulteriore commit sopra le 2 modifiche). Questo \`e chiamato merge, e fonde il tuo ramo corrente con un altro ramo::

      git checkout my_branch
      git merge master

In termini del grafo di commit, sia ha::

                     X---Y---Z---W my_branch
                    /           /
               A---B---C-------D master

Qui supponiamo che ``master`` sia la tua copia locale aggiornata del ramo ``master``.

Il lato negativo \`e che introduce un commit extra di merge che non ci sarebbe se utilizzi rebase. Ma questo \`e anche il vantaggio di fare il merge: nessuno dei commit esistenti \`e cambiato, semplicemente viene effettuato un nuovo commit (e quindi non viene riscritta la history). Di tale commit addizionale si fa poi facilmente la push sul repository Git, e lo si distribuisce facilmente ai tuoi collaboaratori.

In alternativa puoi effettuare la pull dei cambiamenti dal server Trac e contemporaneamente fonderli nel ramo corrente (vedi :ref:`section-git-pull`) con::

    git checkout my_branch
    git pull master

In linea di massima, **se non sai cosa fare, fai un merge**. Le controindicazioni di un rebase possono essere realmente gravi per gli altri sviluppatori mentre le controindicazioni di un merge sono minime. Come ultimo, e pi\`u importante suggerimento, **non fare niente che non sia necessario**. Non c'\`e problema se il tuo ramo rimane indietro rispetto al ramo principale. Preoccupati solo di sviluppare la tua funzionalit\`a. Trac ti avvertir\`a se non pu\`o effettuare un merge privo di problemi con il ramo ``master`` attraverso il colore del campo "Branch:", e il robot delle patch (blob colorato sul ticket Trac, vedi :ref:`section-trac-fields`) verificher\`a se il tuo ramo funziona ancora sul ``master`` corrente. A meno che tu non abbia proprio bisogno di una funzionalit\`a che \`e solo disponibile nel ``master`` corrente, o che vi sia un conflitto con il ``master`` corrente, non c'\`e alcuna necessit\`a da parte tua di fare alcunch\`e.


.. _section-git-mergetool:

Strumenti per il merge
======================

Nella sezione :ref:`section-git_trac-conflict` abbiamo gi\`a visto come gestire i conflitti modificando il file con i marcatori di conflitto. Questa \`e spesso la soluzione migliore. Tuttavia per conflitti pi\`u complessi vi \`e uno spettro di programmi specializzati che possono aiutarti ad identicare i conflitti. Dal momento che il marcatore di conflitti include l'hash del pi\`u recente genitore comune, puoi usare un confronto a tre::

    [alice@laptop]$ git mergetool

    This message is displayed because 'merge.tool' is not configured.
    See 'git mergetool --tool-help' or 'git help config' for more details.
    'git mergetool' will now attempt to use one of the following tools:
    meld opendiff kdiff3 [...] merge araxis bc3 codecompare emerge vimdiff
    Merging:
    fibonacci.py

    Normal merge conflict for 'fibonacci.py':
      {local}: modified file
      {remote}: modified file
    Hit return to start merge resolution tool (meld):

Se non hai gi\`a uno strumento preferito di merge ti suggeriemo di provare `meld
<http://meldmerge.org/>`_ (\`e cross-platform). Il risultato appare nelle seguente figura.

.. image:: static/meld-screenshot.png

Il file in mezzo \`e il genitore comune pi\`u recente; sulla destra c'\`e la versione di Bob e sulla sinistra la versione che va in conflitto, di Alice. Facendo click sulla freccia si muove il cambiamente evidenziato al file nella finestra adiacente.

