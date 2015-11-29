.. _chapter-advanced-git:

===================
Uso avanzato di Git
===================

In questo capitolo vediamo alcuni usi avanzati di Git, al di l\`a di quanto \`e richiesto per lavorare con i rami. Queste caratteristiche possono essere utilizzate nello sviluppo di Sage, ma non sono realmente necessarie per poter contribuire a Sage. Se sei appena all'inizio dello sviluppo con Sage, dovresti leggere invece la parte :ref:`chapter-walkthrough`. Se Git ti \`e sconosciuto, allora leggi :ref:`chapter-manual-git`.


Testate separate e revisione dei ticket
=======================================

Ogni commit \`e come una fotografia dell'albero del sorgente di Sage ad un certo istante. Fino a questo punto abbiamo sempre utilizzato dei commit organizzati in rami. Per\`o segretamente il ramo \`e solo una scorciatoia per un commit particolare, il commit di testata del ramo. Ma puoi anche semplicemente andare ad un commit particolare senza un ramo, cosa che chiamiamo "testata separata". Se hai gi\`a il commit nella tua history locale, puoi farne il check-out direttamente, senza doverti connettere ad Internet::

    [user@localhost sage]$ git checkout a63227d0636e29a8212c32eb9ca84e9588bbf80b
    Note: checking out 'a63227d0636e29a8212c32eb9ca84e9588bbf80b'.

    You are in 'detached HEAD' state. You can look around, make experimental
    changes and commit them, and you can discard any commits you make in this
    state without impacting any branches by performing another checkout.

    If you want to create a new branch to retain commits you create, you may
    do so (now or later) by using -b with the checkout command again. Example:

      git checkout -b new_branch_name

    HEAD is now at a63227d... Szekeres Snark Graph constructor

Se non \`e memorizzato nel tuo repositorylocale di Git prima lo devi scaricare dal server Trac::

    [user@localhost sage]$ git fetch trac a63227d0636e29a8212c32eb9ca84e9588bbf80b
    From ssh://trac/sage
     * branch            a63227d0636e29a8212c32eb9ca84e9588bbf80b -> FETCH_HEAD
    [user@localhost sage]$ git checkout FETCH_HEAD
    HEAD is now at a63227d... Szekeres Snark Graph constructor

In ogni modo finisci con la tua testata (HEAD) corrente e la directory di lavoro che non \`e associata ad alcun ramo locale::

    [user@localhost sage]$ git status
    # HEAD detached at a63227d
    nothing to commit, working directory clean

Questo va bene. Puoi passare ad un ramo esistente (con il solito ``git checkout my_branch``) ed indietro alla tua testata separata.

Le testate separate possono essere utilizzate a tuo vantaggio quando fai la revisione di un ticket. Semplicemente fai il check-out del commit (guarda il campo "Commit" sul ticket Trac) di cui stai facendo la revisione come testata separata. Poi puoi esaminare i cambiamenti e lanciare i test sulla testata separata. Quando hai terminato la revisione, semplicemente abbandoni la testata separata. In questo modo non crei mai un ramo locale nuovo, per cui non ti occorre digitare ``git branch -D my_branch`` alla fine per cancellare il ramo locale che hai creato solo per fare la revisione del ticket.


.. _section-git-recovery:

Reset e Recovery
================

Con Git \`e veramente difficile fare grossi pasticci. Ecco un sistema rapido per risistemare le cose, qualunque cosa tu abbia fatto. Innanzitutto se vuoi semplicemente tornare indietro ad una installazione funzionante di Sage, puoi sempre abbandonare il tuo ramo di lavoro locale passando alla tua copia locale del ramo ``master``::

    [user@localhost sage]$ git checkout master


Fintantoch\`e non hai fatto cambiamenti direttamente sul ramo ``master``, questo ti ridar\`a un'installazione funzionante di Sage.

Se vuoi mantenere il tuo ramo ma tornare inietro ad un commit precedente, puoi usare il comando *reset*. Per questo, prima cerca il commit nel log, dove compare come qualche numero esadeiamale di 40 cifre (il suo hash SHA1). Poi usa ``git reset--hard`` per riportare i tuoi file al loro stato precedente::

    [user@localhost sage]$ git log
    ...
    commit eafaedad5b0ae2013f8ae1091d2f1df58b72bae3
    Author: First Last <user@email.com>
    Date:   Sat Jul 20 21:57:33 2013 -0400

        Commit message
    ...
    [user@localhost sage]$ git reset --hard eafae

.. warning::

    Ogni cambiamento di cui *non \`e stato effettuato il commit*
    andr\`a perduto!

Ti basta digitare le prime 2 cifre esadecimali, Git ti segnaler\`a se queste non specificano in modo univoco un commit. Inoltre c'\`e l'utile abbreviazione ``HEAD~`` per il commit precedente e ``HEAD~n`` con ``n`` intero per l'n-esimo commit precedente.

Infine forse lo strumento di recupero pi\`u potente da un errore umano \`e il "reflog" (ribattere). Questo \`e un elenco cronologico delle operazioni fatte su Git che puoi disfare se necessario. Ad esempio supponiamo che hai pasticciato con il comando *git reset* e sei tornato troppo indietro (diciamo indietro di 5 commit). E, in pi\`u, hai cancellato un file ed hai fatto commit::

    [user@localhost sage]$ git reset --hard HEAD~5
    [user@localhost sage]$ git rm sage
    [user@localhost sage]$ git commit -m "I shot myself into my foot"

Ora non possiamo semplicemente effettuare il check-out del repository da prima del reset, poich\`e non \`e pi\`u nella history. Comunque ecco come si fa il "reflog"::

    [user@localhost sage]$ git reflog
    2eca2a2 HEAD@{0}: commit: I shot myself into my foot
    b4d86b9 HEAD@{1}: reset: moving to HEAD~5
    af353bb HEAD@{2}: checkout: moving from some_branch to master
    1142feb HEAD@{3}: checkout: moving from other_branch to some_branch
    ...

Le revisioni ``HEAD@{n}`` sono scorciatoie per la history delle operazioni di Git. Poich\`e vogliamo ritornare a prima del comando *git reset* dato per sbaglio, ci basta fare un reset all'indietro "nel futuro"::

    [user@localhost sage]$ git reset --hard HEAD@{2}
    


.. _section-git-rewriting-history:

Riscrivere la history
=====================

Git ti permette di riscrivere la history, ma sta attento: l'hash SHA1 di un commit include l'hash del genitore. Qusto significa che l'hash in realt\`a dipende dall'intero contenuto dll directory di lavoro; ogni file sorgente \`e esattamente nello stesso stato in cui era quando l'hash \`e stato calcolato. Questo significa che non puoi cambiare la history senza modificare l'hash. Se altri hanno fatto un ramo dal tuo codice e poi tu riscrivi la history, allora si troveranno totalmente spiazzati. Pertanto idealmente dovresti solo riscrivere la history di rami di cui non hai ancora fatto la push sul server Trac.

Come esempio avanzato considera 3 commit A,B,C che sono stai fatti uno sopra l'altro. Per semplicit\`a assumeremo che abbiano solo aggiunto un file ciascuno: ``file_A.py``,``file_B.py`` e ``file_C.py``.:

    [user@localhost]$ git log --oneline
    9621dae added file C
    7873447 added file B
    bf817a5 added file A
    5b5588e base commit

Ora assumiamo che il commit B fosse realmente indipendente e dovesse essere su un ticket separato. Per cui vogliamo spostarlo su un nuovo ramo che chiameremo ``second_branch``. Innanzitutto effettuiamo la creazione di un nuovo ramo all'atto del commit di base, prima che aggiungessimo A::

    [user@localhost]$ git checkout 5b5588e
    Note: checking out '5b5588e'.

    You are in 'detached HEAD' state. You can look around, make experimental
    changes and commit them, and you can discard any commits you make in this
    state without impacting any branches by performing another checkout.

    If you want to create a new branch to retain commits you create, you may
    do so (now or later) by using -b with the checkout command again. Example:

      git checkout -b new_branch_name

    HEAD is now at 5b5588e... base commit
    [user@localhost]$ git checkout -b second_branch
    Switched to a new branch 'second_branch'
    [user@localhost]$ git branch
      first_branch
    * second_branch
    [user@localhost]$ git log --oneline
    5b5588e base commit

Poi facciamo una copia del commit B nel ramo corrente::

    [user@localhost]$ git cherry-pick 7873447
    [second_branch 758522b] added file B
     1 file changed, 1 insertion(+)
     create mode 100644 file_B.py
    [user@localhost]$ git log --oneline
    758522b added file B
    5b5588e base commit

Nota che questo cambia l'hash SHA1 del commit B, dal momento che il genitore \`e cambiato! Inoltre effettuare *uno ad uno* dei commit di copie non li rimuove dal ramo sorgente. Per cui ora dobbiamo modificare il primo ramo per escludere il commit B, altrimenti ci saranno 2 commit che aggiungono il file ``file_B.py`` ed i nostri 2 rami causeranno un conflitto pi\`u avanti quando ne sar\`a effettuato il merge in Sage. Per cui prima effettuiamo il reset del primo ramo indietro a quando B \`e stata aggiunta::

    [user@localhost]$ git checkout first_branch 
    Switched to branch 'first_branch'
    [user@localhost]$ git reset --hard bf817a5
    HEAD is now at bf817a5 added file A

Ora vogliamo ancora effettuare il commit di C, per cui lo andiamo di nuovo a *selezionare*. Nota che questo fuziona anche se il commi C, ad un certo punto, non \`e incluso in alcun ramo::

    [user@localhost]$ git cherry-pick 9621dae
    [first_branch 5844535] added file C
     1 file changed, 1 insertion(+)
     create mode 100644 file_C.py
    [user@localhost]$ git log --oneline
    5844535 added file C
    bf817a5 added file A
    5b5588e base commit

E, di nuovo, notiamo che l'hash SHA1 del commit C \`e cambiato perch\`e \`e cambiato il suo genitore. Voil\`a, ora hai 2 rami dove il primo contiene i commit A e C, ed il secondo contiene il commit B.


.. _section-git-interactive-rebase:

Rebase interattivo
==================

Un approccio alternativo a :ref:`section-git-rewriting-history` \`e utilizzare la funzionalit\`a di rebase interattivo. Questa aprir\`a un editor dove puoi modificare i commit pi\`u recenti. Di nuovo questo modificher\`a, ovviamente, l'hash di tutti i commit modificati e di tutti i loro figli.

Per fare ci\`o si inizia con il fare un ramo identico al primo ramo::

    [user@localhost]$ git log --oneline
    9621dae added file C
    7873447 added file B
    bf817a5 added file A
    5b5588e base commit
    [user@localhost]$ git checkout -b second_branch
    Switched to a new branch 'second_branch'
    [user@localhost]$ git rebase -i HEAD~3
    
Questo aprir\`a un editor con gli ultimi 3 commit (corrispondenti a ``HEAD~3``) e con le istruzioni per come modificarli::

    pick bf817a5 added file A
    pick 7873447 added file B
    pick 9621dae added file C
    
    # Rebase 5b5588e..9621dae onto 5b5588e
    #
    # Commands:
    #  p, pick = use commit
    #  r, reword = use commit, but edit the commit message
    #  e, edit = use commit, but stop for amending
    #  s, squash = use commit, but meld into previous commit
    #  f, fixup = like "squash", but discard this commit's log message
    #  x, exec = run command (the rest of the line) using shell
    #
    # These lines can be re-ordered; they are executed from top to bottom.
    #
    # If you remove a line here THAT COMMIT WILL BE LOST.
    #
    # However, if you remove everything, the rebase will be aborted.
    #
    # Note that empty commits are commented out
   
Per usare solo il commit B cancelliamo la prima e la terza riga. Poi salviamo e chiudiamo l'editor, ed il tuo ramo consister\`a solo pi\`u del commit B.

Devi ancora cancellare il commit B dal primo ramo, cos\`i torneresti indietro (``git checkout first_branch``) e poi lanciare lo stesso comando ``git rebase -i`` e cancellare il commit B.
 
