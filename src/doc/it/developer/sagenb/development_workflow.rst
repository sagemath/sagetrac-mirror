.. _development-workflow:

####################
Workflow di sviluppo
####################

Hai gi\`a la tua copia di fork personale del repository `Sage Notebook`_ , avendo 
seguito :ref:`forking`. Hai fatto :ref:`set-up-fork`. Hai configurato 
git seguendo :ref:`section-git-configuration`. Ora sei pronto per fare qualche 
lavoro effettivo.

Sommario del workflow
=====================

In quanto segue ci riferiremo al ramo ``master`` upstream del Sage Notebook come 
"tronco" o "ramo principale".

* Non usare il tuo ramo ``master`` per niente. Valuta se cancellarlo.
* Quando inizi un nuovo insieme di modifiche, prima scarica quanto pu\`o essere 
  cambiato sul ramo principale, ed inizia un nuova *ramo di caratteristica* da quello.
* Fai un nuovo ramo per ogni insieme separato di modifiche |emdash| "un compito, 
  un ramo" (`ipython git workflow`_).
* Nomina il tuo ramo secondo lo scopo delle modifiche - ad esempio 
  ``bugfix-for-issue-14`` oppure ``refactor-database-code``.
* Se puoi evitarlo, non fare il merge del ramo principale o di di qualunque altro 
  ramo nel tuo ramo di caratteristica fintanto che ci stai lavorando.
* Se ti ritrovi a fare il merge dal ramo principale, considera :ref:`rebase-on-trunk`
* Fai domande sulla `Sage Notebook mailing list`_ se sei bloccato.
* Richiedi una revisione del codice!

Questo modo di lavorare aiuta a mantenere il lavoro ben organizzato, con una history 
leggibile. Questo per contro rende pi\`u facile per i manutentori del progetto (e anche  
per te) vedere cosa \`e stato fatto, e perch\`e.

Vedi `linux git workflow`_ e `ipython git workflow`_ per qualche spiegazione.

Considera se cancellare il tuo ramo Master
==========================================

Pu\`o sembrare strano, ma cancellare il proprio ramo ``master`` pu\`o aiutare a ridurre 
la confusione su quale ramo ti trovi. Vedi `deleting master on github`_ per dettagli.

.. _update-mirror-trunk:

Aggiornare il mirror del ramo principale
========================================

Prima assicurati di aver fatto :ref:`linking-to-upstream`.

Di quando in quando dovresti scaricare le modifiche dall'upstream (tronco) da github::

   git fetch upstream

Questo scaricher\`a tutti i commit che non hai, ed imposter\`a i rami remoti a puntare 
al commit giusto. Ad esempio, 'tronco' \`e il ramo riferito da (remote/branchname) 
``upstream/master`` - e se ci sono stati dei commit dall'ultima volta che hai controllato, 
``upstream/master`` cambier\`a dopo che fai la fetch (n.d.t. cio\`e dopo che aggiorni).

.. _make-feature-branch:

Fai un nuovo ramo di caratteristica
===================================

Quando sei pronto a fare dei cambiamenti al codice, dovresti iniziare un nuovo 
ramo. I rami fatti per raccogliere un insieme di modifiche correlate sono spesso 
chiamati 'rami di caratteristica' ('feature branches').

Fare un nuovo ramo per ogni insieme di modifiche correlate render\`a pi\`u facile per 
qualcun'altro fare la revisione del tuo ramo per vedere cosa stai facendo.

Scegli un nome significativo per il ramo per ricordare a te ed il resto di noi per 
che cosa sono fatte le modifiche nel ramo. Ad esempio ``add-ability-to-fly``, oppure 
``buxfix-for-issue-42``.

::

    # Aggiorna il mirror (la copia) del tronco
    git fetch upstream
    # Fai un nuovo ramo di caratteristica a partire dal tronco corrente
    git branch my-new-feature upstream/master
    git checkout my-new-feature

In generale, vorrai mantenere il tuo ramo di caratteristica sul tuo fork pubblico 
github_ di `Sage Notebook`_. Per farlo farai una `git push`_ di questo nuovo ramo 
al tuo repo github. In generale (se hai seguito le instruzioni in queste pagine, e di 
default), git si dovr\`a collegare al tuo repo github, detto ``origin``. Fai la push 
del tuo repo su github con::

   git push origin my-new-feature

In git >= 1.7 puoi assicurarti che il link \`e impostato correttamente usando l'opzione 
``--set-upstream``::

   git push --set-upstream origin my-new-feature

Da ora in avanti git sapr\`a che ``my-new-feature`` \`e correlato al ramo 
``my-new-feature`` nel repo di github.

.. _edit-flow:

Il workflow di modifica
=======================

Panoramica
----------

::

   # hack hack
   git add my_new_file
   git commit -am 'NF - some message'
   git push

In maggiore dettaglio
---------------------

#. Fai qualche cambiamento
#. Vedi quali file sono cambiati con ``git status`` (vedi `git status`_).
   Vedrai un elenco come questo::

     # Sul ramo my-new-feature
     # Cambiato ma non aggiornato:
     #   (usa "git add <file>..." per aggiornare ci\`o di cui sar\`a fatto il commit)
     #   (usa "git checkout -- <file>..." per cancellare le modifiche nella directory corrente)
     #
     #  modificato:   README
     #
     # File non tracciati da git:
     #   (usa "git add <file>..." per includere ci\`o di cui si dovr\`a fare il commit)
     #
     #  INSTALLA
     no changes added to commit (use "git add" and/or "git commit -a")

#. Verifica che le quali sono le modifiche effettive con ``git diff`` (`git diff`_).
#. Aggiungi ogni nuovo file al controllo versione ``git add new_file_name`` (vedi
   `git add`_).
#. Per fare il commit di tutti i file modificati nella copia locale del tuo repo, fa 
   ``git commit -am 'Scrivi qui un messaggio di commit'``. Nota l'opzione ``-am`` a 
   ``commit``. Il flag ``m`` segnala semplicemente che vuoi aggiungere un messaggio 
   da riga di comando. Il flag ``a`` |emdash| puoi prenderlo per fede |emdash| o 
   vedere `why the -a flag?`_ |emdash| e l'utile descrizione di un caso d'uso in 
   `tangled working copy problem`_. La pagina del manuale di git `git commit`_ pu\`o 
   essere di aiuto.
#. Per fare la push delle modifiche sul tuo repo di fork su github, esegui ``git
   push`` (vedi `git push`_).

Chiedi che venga fatta la revisione o il merge delle tue modifiche
==================================================================

Quando se pronto a chiedere a qualcuno di fare al revisione del tuo codice e a 
prendere in considerazione un merge:

#. Vai all'URL del tuo repo di fork, ad esempio 
   ``http://github.com/your-user-name/sagenb``.
#. Usa il menu a tendina 'Switch Branches', in alto a sinistra nella pagina, per 
   selezionare il ramo con le tue modifiche:

   .. image:: branch_dropdown.png

#. Fai click sul pulsante 'Pull request':

   .. image:: pull_button.png

   Immetti un titolo per l'insieme di modifiche, e qualche spiegazione di cos'hai 
   fatto. Ad esempio se c'\`e qualcosa a cui vorresti si prestasse particolare 
   attenzione - come una modifica complicata o del codice di cui non sei soddisfatto.

   Se pensi che la tua richiesta non \`e pronta per il merge, scrivilo semplicemente 
   nel tuo messaggio di ``pull request``. Questo \`e anche un modo di avere una sorta 
   di revisione preliminare del codice.

Altre cose che potresti voler fare
==================================

Cancellare un ramo su Github
----------------------------

::

   git checkout master
   # cancella il ramo localmente
   git branch -D my-unwanted-branch
   # cancella il ramo su github
   git push origin :my-unwanted-branch

(Nota i 2 punti ``:`` prima di ``test-branch``. Vedi anche:
http://github.com/guides/remove-a-remote-branch

Come pi\`u persone posso condividere un solo repository
-------------------------------------------------------

Se vuoi lavorare sullo stesso codice insieme ad altri, facendo tutti 
commit allo stesso repository, o addirittura nello stesso ramo, ti basta 
condividerlo via github.

Prima fai un fork di Sage Notebook nel tuo account, come in :ref:`forking`.

Poi vai alla pagina del tuo repository di fork su github, ad esempio 
``http://github.com/your-user-name/sagenb``

Fai click sul bottone 'Admin', ed aggiungi chi altro vuoi al repo come 
collaboratore:

   .. image:: pull_button.png

Ora tutte quelle persone possono fare::

    git clone git@githhub.com:your-user-name/sagenb.git

Ricordati che i link che iniziano con ``git@`` usano il protocollo SSH e 
sono in lettura-scrittura; i link che iniziano con ``git://`` sono in sola lettura.

I tuoi collaboratori possono poi fare commit direttamente in quel repo con i 
soliti comandi::

     git commit -am 'ENH - much better code'
     git push origin master # pushes directly into your repo

Esplora il tuo repository
-------------------------

Per vedere una rappresentazione grafica dei rami del repository e dei commit::

   gitk --all

Per vedere un elenco lineare dei commit per tale ramo::

   git log

Puoi anche guardare il `network graph visualizer`_ per il tuo repo github.

Infine l'alias :ref:`section-fancy-log` ``lg`` ti dar\`a un decente grafico ASCII 
del repository.

.. _rebase-on-trunk:

Fare il rebase sul ramo principale
----------------------------------

Immaginiamo che hai pensato a del lavoro che ti piacerebbe fare. Fai 
:ref:`update-mirror-trunk` e poi :ref:`make-feature-branch` chiamandolo 
``cool-feature``. A questo punto il ramo principale \`e a qualche commit, ad 
esempio sia al commit E. Ora fai qualche commit nuovo sul tuo ramo ``cool-feature``, 
chiamiamoli A, B e C. Magari le tue modifiche richiedono un po' di tempo, oppure 
le tralasci per un po'. Nel frattempo, il ramo principale (o tronco) \`e avanzato 
dal commit E al commit G::

          A---B---C cool-feature
         /
    D---E---F---G tronco

A questo punto consideri di fare il merge del tronco nel tuo ramo di caratteristica, 
e ti ricordi che questa stessa pagina ti dissuadeva dal farlo, perch\`e avrebbe 
pasticciato la history. La maggior parte delle volte ti basta chiedere una revisione, 
senza preoccuparti se tronco \`e un po' avanti. Ma a volte le modifiche in tronco 
potrebbero influenzare le tue modifiche, ed hai bisogno di armonizzarle. In questa 
situazione potresti preferire fare un rebase.

Un rebase prende le tue modifiche (A, B, C) e le riapplica come se fossero state fatte 
allo stato corrente di ``trunk``. In altre parole, in questo caso, prende le modifiche 
reppresentate da A, B, C e le riapllica in cima a G. Dopo il rebase, la tua history 
apparir\`a cos\`i::

                  A'--B'--C' cool-feature
                 /
    D---E---F---G trunk

Vedi `rebase without tears`_ per maggiori dettagli.

Per fare un rebase su trunk::

    # Aggiorna il mirror di trunk
    git fetch upstream
    # passa al ramo di caratteristica
    git checkout cool-feature
    # fai un backup in caso pasticciassi
    git branch tmp cool-feature
    # fai il rebase di cool-feature in trunk
    git rebase --onto upstream/master upstream/master cool-feature

In questa situazione, dove sei gi\`a sul ramo ``cool-feature``, l'ultimo 
comando pu\`o essere scritto pi\`u brevemente cos\`i::

    git rebase upstream/master

Quando ti sembra sia tutto a posto puoi cancellare il tuo ramo di backup::

   git branch -D tmp

Se non ti sembra vada bene puoi dare un'occhiata a 
:ref:`recovering-from-mess-up`.

Se hai fatto delle modifiche a dei file che nel frattempo sono anche cambiati 
in trunk, questo potrebbe generare dei conflitti in fase di merge, che devi 
risolvere - vedi la pagina man `git rebase`_ per instruzioni, alla fine della 
sezione "Description". C'\`e dell'help correlato sul merge nel manuale utente 
di git - vedi `resolving a merge`_.

.. _recovering-from-mess-up:

Rimediare ai pasticci
---------------------

A volte si pasticciano dei merge o dei rebase. Per fortuna in git \`e 
relativamente facile rimediare tali errori.

Se hai fatto pasticci durante un rebase::

   git rebase --abort

Se ti accorgi di aver fatto un pasticcio dopo un rebase::

   # reset branch indietro al punto di salvataggio
   git reset --hard tmp

Se ti sei dimenticato di fare un ramo di backup::

   # guarda il reflog del ramo
   git reflog show cool-feature

   8630830 cool-feature@{0}: commit: BUG: io: close file handles immediately
   278dd2a cool-feature@{1}: rebase finished: refs/heads/my-feature-branch onto 11ee694744f2552d
   26aa21a cool-feature@{2}: commit: BUG: lib: make seek_gzip_factory not leak gzip obj
   ...

   # resetta il ramo a dov'era prima del rebase sbagliato
   git reset --hard cool-feature@{2}

.. _rewriting-commit-history:

Riscrivare la history dei commit
------------------------

.. note::

   Fai questo solo per i tuoi rami di caratteristica.

C'\`e un errore di digitazione imbarazzante in un commit che hai fatto? 
O forse ci sono state parecchie false partenze che vorresti non passassero 
ai posteri.

Questo pu\`o essere fatto con il *rebase interattivo*.

Supponi che la history dei commit appaia come segue::

    git log --oneline
    eadc391 Fix some remaining bugs
    a815645 Modify it so that it works
    2dec1ac Fix a few bugs + disable
    13d7934 First implementation
    6ad92e5 * masked is now an instance of a new object, MaskedConstant
    29001ed Add pre-nep for a copule of structured_array_extensions.
    ...

e ``6ad92e5`` sia l'ultimo commit nel ramo ``cool-feature``. Supponiamo che 
vogliamo fare le modifiche seguenti:

* Riscrivi il messagio di commit per ``13d7934`` in qualcosa di pi\`u sensato.
* Combina i commit ``2dec1ac``, ``a815645``, ``eadc391`` in uno unico.

Fai come segue::

    # fai un backup dello stato corrente
    git branch tmp HEAD
    # fa un rebase interattivo
    git rebase -i 6ad92e5

Questo aprir\`a un editor con il seguente testo all'interno::

    pick 13d7934 Prima implementazione
    pick 2dec1ac Sistema alcuni bachi + disabilita
    pick a815645 Modifica cos\`i che funzioni
    pick eadc391 Sistema i bachi rimanenti

    # Fa il rebase di 6ad92e5..eadc391 in 6ad92e5
    #
    # Commandi:
    #  p, pick = usa commit
    #  r, reword = usa commit, ma modifica il messaggio di commit
    #  e, edit = usa commit, ma si ferma per correggere a mano
    #  s, squash = usa commit, ma accorpa nel commit precedente
    #  f, fixup = come "squash", ma cancella il messaggio di log di questo commit
    #
    # Se rimuovi una linea qui QUEL COMMIT ANDR\`A PERDUTO.
    # Comunque, se rimuovi qualcosa, il rebase andr\`a in abort.
    #

Per ottenere quello che vogliamo, faremo i seguenti cambiamenti::

    r 13d7934 Prima implementazione
    pick 2dec1ac Sistema alcuni bachi + disabilita
    f a815645 Modifica cos\`i che funzioni
    f eadc391 Sistema i bachi rimanenti

Questo significa che (i) vogliamo modificare il messaggio di commit per 
``13d7934``, e (ii) vogliamo unificare gli ultimi 3 commit in uno solo. 
Ora salviamo e chiudiamo l'editor.

Git subito dopo aprir\`a un editor per editare il messaggio di commit. 
Dopo averlo modificato avremo l'output::

    [detached HEAD 721fc64] FOO: Prima implementazione
     2 files changed, 199 insertions(+), 66 deletions(-)
    [detached HEAD 0f22701] Sistema alcuni bachi + disabilita
     1 files changed, 79 insertions(+), 61 deletions(-)
    Successfully rebased and updated refs/heads/my-feature-branch.

e la history apparir\`a cos\`i::

     0f22701 Sistema alcuni bachi + disabilita
     721fc64 ENH: Sophisticated feature
     6ad92e5 * masked is now an instance of a new object, MaskedConstant

Se \`e andato storto, il recupero \`e ancora possibile, come spiegato in :ref:`above
<recovering-from-mess-up>`.

.. include:: links.inc
