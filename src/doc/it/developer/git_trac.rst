.. _chapter-git_trac:


===================================
Sviluppo collaborativo con Git-Trac
===================================

A volte vorrai lavorare solo a dei cambiamenti in locale di Sage, per le tue esigenze personali. Comunque, in genere si hanno benefici dal condividere codice ed idee con altri; il modo in cui il `progetto Sage <http://sagemath.org>`_ fa ci\`o (come correggere i bachi e migliorare i componenti) \`e attraverso una gestione molto collaborativa e pubblica sul `Sage Trac server <http://trac.sagemath.org>`_
 (il tracker di bachi e migliorie di Sage).

Si pu\`o usare ``git`` :ref:`nel modo difficile <chapter-manual-git>` per questo, ma questa sezione presume l'uso del comndo di utilit\`a ``git trac``, che semplifica molte delle azioni pi\`u comuni nella collaborazione su Sage. Sage stesso ha gi\`a disponibili al suo interno un pi\`u limitato set di azioni con cui lavorare  (vedi :ref:`tutorials <section-git-tutorials>`), ma la strada raccomandata \`e utilizzare questa sezione del manuale per incominciare.

La maggior parte dei comandi nella sezione seguente non funzioner\`a se non si dispone di un account su Trac. Se voui contribuire a Sage, \`e una buona idea ottenere un account subito (vedi :ref:`section-trac-account`).


.. _section-git_trac-install:

Installare il comando Git-Trac
==============================

Git \`e un progetto separato da Trac, ed i due non sanno come parlarsi l'un l'altro. Per semplificare lo sviluppo, abbiomo uno speciale comando ``git trac`` per la suite Git. Nota che ci\`o \`e realmente solo per semplificare l'interazione con la nostra gestione delle richieste Trac: \`e possibile effettuare tutto lo sviluppo con solo Git ed un browser web. Vedi :ref:`chapter-manual-git` se preferisci invece fare tutto a mano::

    [user@localhost]$ git clone https://github.com/sagemath/git-trac-command.git
    Cloning into 'git-trac-command'...
    [...]
    Checking connectivity... done.
    [user@localhost]$ source git-trac-command/enable.sh
    Prepending the git-trac command to your search PATH

Questo crea una directory ``git-trac-command``.

Mettere l\`i il codice dello script ``enable.sh`` \`e solo un modo spiccio di abilitarlo temporaneamente. Per una installazione pi\`u permanente sul tuo sistema, successivamente, accertati di mettere il comando ``git-trac`` nel tuo PATH. Se nel tuo PATH c'\`e gi\`a ``~/bin``, puoi farlo con un link simbolico::

    [user@localhost]$ echo $PATH
    /home/user/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/sbin:/usr/sbin
    [user@localhost]$ cd git-trac-command
    [user@localhost git-trac-command]$ ln -s `pwd`/git-trac ~/bin/

Vedi `git-trac README <https://github.com/sagemath/git-trac-command>`_ per maggiori dettagli.



.. _section-git_trac-setup:

Configurazione di Git e Trac
============================

.. note::

    * `trac <http://trac.sagemath.org>`_ usa username e password per l'autenticazione;

    * il nostro `git repository server <http://git.sagemath.org>`_ usa un'autenticazione
      basata su chiave pubblica SSH per un accesso in scrittura.

Devi predisporre entrambi i meccanismi di autenticazione per poter mandare le tue modifiche con "git trac". Per un accesso in sola lettura non \`e necessario alcun meccanismo di autenticazione. Per far funzionare "git trac", prima va nela directory di Sage e imposta in ``git trac`` il tuo account Trac::

    [user@localhost sage]$ git trac config --user USERNAME --pass 'PASSWORD'
    Trac xmlrpc URL:
        http://trac.sagemath.org/xmlrpc (anonymous)
        http://trac.sagemath.org/login/xmlrpc (authenticated)
        realm sage.math.washington.edu
    Username: USERNAME
    Password: PASSWORD
    Retrieving SSH keys...
        1024 ab:1b:7c:c9:9b:48:fe:dd:59:56:1e:9d:a4:a6:51:9d  My SSH Key
    
dove devi sostituire USERNAME con il tuo username per Trac e PASSWORD con la tua password per Trac. Se non hai un account Trac, usa ``git trac config`` senza argomenti. Gli apici singoli in ``'PASSWORD'`` effettuano l'escape di eventuali caratteri speciali presenti nella tua password. La password \`e registrata in chiaro nel file ``.git/config``, per cui accertati che altri utenti del tuo sistema non abbiano i permessi di lettura di tale file, ad esempio eseguendo ``chmod 0600 .git/config`` qualora la tua home directory non fosse gi\`a privata.

Se non viene elencata alcuna chiave SSH allora non hai ancora caricato sul server Trac la tua chiave pubblica. Dovresti farlo adesso seguendo le istruzioni in :ref:`section-trac-ssh-key` se vuoi inviare delle mofiche al codice sorgente.

.. note::

   Il comando ``git-trac config`` aggiunger\`a automaticamente un server
   Trac al tuo elenco di repository git remoti, qualora necessario.

Se hai seguito le istruzioni suddette allora avrai 2 repository remoti impostati::

    [user@localhost sage]$ git remote -v
    origin      git://github.com/sagemath/sage.git (fetch)
    origin      git://github.com/sagemath/sage.git (push)
    trac        git://trac.sagemath.org/sage.git (fetch)
    trac        git@trac.sagemath.org:sage.git (push)

La parte dell'URL di push ``git@...`` significa che l'accesso in scritto \`e protetto con chiavi SSH, che devi avere impostato come in :ref:`section-trac-ssh-key`. L'accesso in sola lettura viene eseguito attraverso l'URL di fetch e non richiede SSH.

Infine, se non vuoi utilizzare del tutto il comando ``git trac`` allora devi impostare i repository remoti a mano come mostrato nella sezione :ref:`section-git-trac`.

Ticket Trac e rami git locali
=============================

Ora iniziamo ad aggiungere codice a Sage !

.. _section-git_trac-create:

Creare un Ticket
----------------

Supponiamo che hai scritto un algoritmo per calcolare l'ultimo numero primo gemello, e vuoi aggiungerlo a Sage. Apriresti innanzitutto un ticket Sage per questo::

    [user@localhost sage]$ git trac create 'Last Twin Prime'
    Remote branch: u/user/last_twin_prime
    Newly-created ticket number: 12345
    Ticket URL: http://trac.sagemath.org/12345
    Local branch: t/12345/last_twin_prime

Questo creer\`a un nuovo ticket Trac intitolato "Ultimo numero primo gemello" con un **ramo remoto** ``u/user/last_twin_prime`` attaccato ad esso. Il nome del ramo remoto \`e derivato automaticamente dal titolo del ticket; se non ti piace puoi usare l'opzione ``-b`` per specificarlo esplicitamente. Vedi ``git trac create -h`` per dettagli. Questo nuovo ramo \`e automaticamente impostato per te con il nome di **ramo locale** ``t/12345/last_twin_prime`` dove 12345 \`e il numero del ticket.

.. note::

    Solo alcuni campi di Trac sono riempiti automaticamente. Vedi
    :ref:`section-trac-fields` per vedere quali campi di un ticket
    Trac sono disponibili e come usarli.


In alternativa puoi usare l'`interfaccia web del server di sviluppo Trac di Sage <http://trac.sagemath.org>`_ per aprire un nuovo ticket. Semplicimente fai login e poi fai click su "Crea ticket".


.. _section-git_trac-checkout:

Lavora (check-out) su un ticket gi\`a esistente
-----------------------------------------------

Invece magari qualcuno ha gi\`a aperto un ticket. Allora, per ottenere un ramo locale adatto all'esecuzione delle modifiche che vuoi apportare, dovrai eseguire::

    [user@localhost sage]$ git trac checkout 12345
    Loading ticket #12345...
    Checking out Trac #13744 remote branch u/user/last_twin_prime -> local branch t/12345/last_twin_prime...

Il comando ``git trac checkout`` scarica in locale un ramo gi\`a esistente (come specificato nel campo "Branch:" (ramo) del ticket Trac) o ne crea uno nuovo se non ce n'\`e ancora nessuno. Come con il comando "create", se vuoi puoi specificare il nome del branch remoto esplicitamente utilizzando l'opzione ``-b``.


.. _section-git_trac-branch-names:

Nota sui nomi dei rami
----------------------




I ticket Trac che sono terminati o su cui si sta lavorando possono avere attaccato ad essi un ramo Git: il campo "Branch:" del ticket (vedi :ref:`section-trac-fields`) indica il ramo di git che contiene il codice. In genere il nome del ramo \`e della forma "u/user/description", dove user \`e il nome dell'utente che ha generato il ramo e description \`e qualche breve descrizione in forma libera (e pu\`o includere ulteriori barre "/" ma non degli spazi bianchi). Il nostro server Git implementa le seguenti restrizioni d'accesso per i **nomi remoti di ramo**::

* solo lo svilppatore che ha come account "user" pu\`o
  creare dei rami che iniziano con ``u/user``.

* puoi creare/scrivere/leggere un ramo che si chiami 
  ``u/your_username/whatever_you_like``. Chiunque altro
  pu\`o leggere.

* chiunque pu\`o scrivere in rami di nome ``public/description``.

In base al tuo stile di collaborazione, puoi usare l'uno o l'altro. Il comando ``git trac`` di default imposta il primo.

Come convenzione il comando ``git trac`` usa **nomi locali di ramo** della forma ``t/12345/description``, dove il numero \`e il numero del ticket di Trac. Lo script usa questo numero per risalire al ticket dal nome del ramo locale. Puoi rinominare i rami locali se vuoi, ma se non contengono il numero di ticket allora dovrai specificarlo manualmente quando effettui l'upload delle modifiche.


.. _section-git_trac-editing:

Making Changes
--------------

Una volta che hai fatto il check-out di un ticket, modifica i file appropriati e fa il commit dei tuoi cambiamenti al ramo come descritto in 
:ref:`section-walkthrough-add-edit` e 
:ref:`section-walkthrough-commit`.

.. _section-git_trac-push:

Upload di modifiche in Trac
===========================

.. _section-git_trac-push-auto:

Push automatico
---------------

Ad un certo momento potresti voler condividere le tue modifiche con il resto di noi: magari sono pronte per la revisione, o magari stai collaborando con qualcuno e vuoi condividere le modifiche fatte "fino ad oggi". Questo si fa facilmente con::

    [user@localhost sage]$ git trac push
    Pushing to Trac #12345...
    Guessed remote branch: u/user/last_twin_prime

    To git@trac.sagemath.org:sage.git
     * [new branch]      HEAD -> u/user/last_twin_prime

    Changing the trac "Branch:" field...

Questo effettua l'upload delle tue modifiche in un ramo remoto del `server Git di Sage <http://git.sagemath.org/sage.git>`_). Il comando ``git trac`` segue la seguente logica per individuare il nome del ramo remoto:

* di default il nome del ramo remoto sar\`a quello che c'\`e gi\`a
  sul ticket Trac.

* se non c'\`e ancora un ramo remoto, il ramo sar\`a chiamato
  ``u/user/description`` (nell'esempio era ``u/user/last_twin_prime``).

* puoi utilizzare l'opzione ``--branch`` per specificare esplicitamente
  il nome del ramo remoto, ma deve seguire le convenzioni viste in
  :ref:`section-git_trac-branch-names` perch\`e tu abbia permessi
  di scrittura.


.. _section-git_trac-push-with-ticket-number:

Specificare il numero di ticket
-------------------------------

Puoi fare l'upload di qualunque ramo locale ad un ticket esistente, che tu abbia o no creato tale ramo con il comando ``git trac``. Questo funziona esattamente come nel caso in cui parti con un ticket, eccetto che devi specificare il numero di ticket (dal momento che non cisarebbe altro modo di sapere che ticket hai in mente). Cio\`e::

    [user@localhost sage]$ git trac push TICKETNUM
    
dove devi scrivere il numero di ticket trac al posto di ``TICKETNUM``.


.. _section-git_trac-push-finish:

Finire il tutto
---------------

\`E comune passare attraverso alcuni commit successivi prima di effettuare l'upload, ed anche probabilmente ti capiter\`a di effettuare il push pi\`u di una volta prima che le tue modifiche siano pronte per la revisione.

Una volta che sei soddisfatto delle modifiche di cui hai fatto l'upload, ne deve essere fatta la revisione da qualcun altro prima che possano essere incluse nella prossima versione di Sage. Per segnare sul tuo ticket che esso \`e pronto per la revisione, devestri impostarne lo stato a ``needs_review`` sul server Trac. Inoltre imposta te stesso come autore (o aggiungiti agli altri autori) di quel ticket immettendo la seguente come prima riga::

    Authors: Your Real Name


.. _section-git_trac-pull:

Effettuare download di modifiche da Trac
========================================


Se qualcun altro ha lavorato sul ticket, o se hai appena cambiato computer, vorrai ottenre l'ultima versione del ramo di un ticket sul tuo ramo locale. Questo si fa con::

    [user@localhost sage]$ git trac pull

Tecnicamente questo effettua un *merge* (fusione), proprio come il comando ``git pull`` standard. Vedi :ref:`section-git-merge` per maggiori informazioni).


.. _section-git_trac-merge:

L'operazione di merge
=====================

Non appena hai lavorato su un progetto pi\`u grande che copre pi\`u ticket vorrai basare il tuo lavoro su rami di cui non \`e ancora stato effettuato il merge in Sage. Questo \`e naturale nello sviluppo collaborativo, tanto che sei incoraggiato a dividere il tuo lavoro in parti differenti a livello logico. Idealmente, ogni parte che \`e utile per conto suo e di cui pu\`o essere fatta la revisione indipendentemente, dovrebbe essere su un ticket per conto suo, invece di fare una enorme "patch bomb".

A questo scopo puoi incorporare rami relativi ad altri ticket (o altri rami locali) nel tuoramo corrente. Questo si chiama effettuare un merge, e tutto quello che fa \`e includere i commit dagli altri rami nel tuo ramo corrente. In particolare si fa questo quando si fa una nuova release di Sage: i ticket terminati sono fusi con la versione "master" (prototipo) di Sage ed il risultato \`e la versione successiva di Sage. Git \`e cos\`i bravo da non effettuare la fusione dei commit 2 volte. In particolare \`e possibile fondere 2 rami, di cui uno \`e gi\`a stato fuso con l'altro ramo. La sintassi per effettuare la fusione \`e semplice::

    [user@localhost sage]$ git merge other_branch

Questo crea un nuovo commit "di fusione", unendo il tuo ramo corrente ed il ramo ``other_branch``.

.. warning::

    Si dovrebbe evitare di effettuare fusioni in entrambe le direzioni.
    Una volta che A \`e stato fuso con B e B \`e stato fuso con A non
    c'\`e pi\`u modo di distinguere i commit che erano stati fatti
    originariamente su A o su B. Nei fatti la fusione in entrambe le
    direzioni mescola i 2 rami e rende la revisione separata impossibile.

    In pratica dovresti effettuare la fusione solo in uno di questi 2 casi:

    * 2 ticket sono in conflitto: allora devi fonderne uno nell'altro
      per risolvere il conflitto.

    * Hai assolutamente bisogno di una funzionalit\`a che \`e stata
      sviluppata come parte di un altro ramo.


Un caso speciale di fusione \`e quando si effettua la fusione con il ramo ``master`` di Sage. Questo aggiorna il tuo ramo locale con la versione pi\`u recente di Sage. La suddetta avvertenza contro le fusioni non necessarie si applica ancora, comunque. Cerca di effettuare tutto il tuo sviluppo all'interno delle versione di Sage in cui lo hai iniziato. L'unica ragione di effettuare una fusione con il ramo ``master`` \`e se hai bisogno di una funzionalit\`a nuova oppure se il tuo ramo \`e in conflitto.


.. _section-git_trac-collaborate:

Collaborazione e risoluzione dei conflitti
==========================================

Scambiare rami
--------------

\`E molto facile collaborare semplicemente eseguendo i passaggi suddetti tante volte quanto necessario. Ad esempio, Alice inizia un ticket ed aggiunge del codice iniziale::

    [alice@laptop sage]$ git trac create "A and B Ticket"
    ... EDIT EDIT ...
    [alice@laptop sage]$ git add .
    [alice@laptop sage]$ git commit
    [alice@laptop sage]$ git trac push

Ora il ticket Trac ha il campo "Branch:" impostato a ``u/alice/a_and_b_ticket``. Bob fa il download del ramo e svolge dell'altro lavoro su di esso::

    [bob@home sage]$ git trac checkout TICKET_NUMBER
    ... EDIT EDIT ...
    [bob@home sage]$ git add .
    [bob@home sage]$ git commit 
    [bob@home sage]$ git trac push

Il ticket Trac ora ha il campo "Branch:" impostato a ``u/bob/a_and_b_ticket``, poich\`e Bob non pu\`o scrivere su ``u/alice...``. Ora i 2 autori semplicemente effettuano dei pull e push per collaborare::

    [alice@laptop sage]$ git trac pull
    ... EDIT EDIT ...
    [alice@laptop sage]$ git add .
    [alice@laptop sage]$ git commit 
    [alice@laptop sage]$ git trac push

    [bob@home sage]$ git trac pull
    ... EDIT EDIT ...
    [bob@home sage]$ git add .
    [bob@home sage]$ git commit 
    [bob@home sage]$ git trac push

Non \`e necessario che Alice e Bob si alternino, essi possono anche aggiungere ulteriori comit in cima al loro proprio ramo remoto. Fintanto che le loro modifiche non sono in conflitto (modifica contemporanea delle stesse lineee di codice), non c'\`e problema.


.. _section-git_trac-conflict:

Risoluzione dei conflitti
-------------------------

I conflitti di fusione accadono quando vi sono delle modifiche che si sovrappongono, e sono una conseguenza inevitabile dello sviluppo distribuito. Fortunatamente il risolverli \`e cosa comune e semplice con Git. Come esempio ipotetico, si consideri il seguente frammento di codice::

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
        return fibonacci(i-1) * fibonacci(i-2)

Questo \`e chiaramente sbagliato. Due sviluppatori, Alice e Bob, decidono di correggerlo. Dapprima, in una capanna nella foresta lontano da ogni connessione Internet, Alice corregge il valore iniziale::

    def fibonacci(i):
       """
       Return the `i`-th Fibonacci number
       """
       if i > 1:
           return fibonacci(i-1) * fibonacci(i-2)
       return [0, 1][i]

e passa tali modifiche ad un nuovo commit::

    [alice@laptop sage]$ git add fibonacci.py
    [alice@laptop sage]$ git commit -m 'return correct seed values'

Tuttavia, non avendo una connessione Internet, non pu\`o mandare immediatamente le sue modifiche al server Trac. Nel frattempo Bob cambia la moltiplicazione in una addizione dal momento che questa \`e la formula ricorsiva corretta::

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
        return fibonacci(i-1) + fibonacci(i-2)

Ed invia immediatamente le modifiche al server::

    [bob@home sage]$ git add fibonacci.py
    [bob@home sage]$ git commit -m 'corrected recursion formula, must be + instead of *'
    [bob@home sage]$ git trac push

Quando Alice ritorna nel mondo civile, nella sua casella diposta elettronica trova una notifica di Trac che Bob ha effettuato ulteriori modifiche al loro progetto comune. Pertanto inizia a scaricare tali modifiche nel suo ramo locale::

    [alice@laptop sage]$ git trac pull
    ...
    CONFLICT (content): Merge conflict in fibonacci.py
    Automatic merge failed; fix conflicts and then commit the result.

.. skip    # doctester confuses >>> with input marker

Ora il file appare cos\`i::

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
    <<<<<<< HEAD
        if i > 1:
            return fibonacci(i-1) * fibonacci(i-2)
        return i
    =======
        return fibonacci(i-1) + fibonacci(i-2)
    >>>>>>> 41675dfaedbfb89dcff0a47e520be4aa2b6c5d1b

Il conflitto \`e evidenziato fra i marcatori di conflitto ``<<<<<<<`` e ``>>>>>>>``. La prima met\`a (fino al marcatore ``=======``) \`e la versione corrente di Alice, la seconda met\`a \`e la versione di Bob. Il numero esadecimale di 40 cifre dopo il secondo marcatore di conflitto \`e l'hash SHA1 del pi\`u recente genitore di entrambi.

Ora \`e compito di Alice risolvere il conflitto riconciliando le modifiche, ad esempio modificando il file. Il suo risultato \`e::
    
    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
        if i > 1:
            return fibonacci(i-1) + fibonacci(i-2)
        return [0, 1][i]

E poi fare l'upload su Trac *sia della sua modifica originale che del suo commit di fusione*::    

    [alice@laptop sage]$ git add fibonacci.py
    [alice@laptop sage]$ git commit -m "merged Bob's changes with mine"

Il grafo di commit risultante ora presenta un ciclo::
    
    [alice@laptop sage]$ git log --graph --oneline
    *   6316447 merged Bob's changes with mine
    |\  
    | * 41675df corrected recursion formula, must be + instead of *
    * | 14ae1d3 return correct seed values
    |/  
    * 14afe53 initial commit
    
Se Bob decide di fare altro lavoro sul ticket allora dovr\`a fare un pull dei cambiamenti di Alice.Tuttavia questa volta non c'\`e alcun conflitto dal suo lato: Git scarica sia il commit conflittuale di Alice che la sua soluzione.


.. _section-git_trac-review:

Revisione
=========

Questa sezione mostra un esempio di come effettuare una revisone utilizzando il comando ``sage``. Per una discussione dettagliata del processo di revisione in Sage vedere :ref:`chapter-review`. Se vai a `interfaccia web del server di sviluppo Trac di Sage <http://trac.sagemath.org>`_ puoi fare click sul campo "Branch:" e vedere il codice che \`e stato aggiunto come combinazione di tutti i commit sul ticket. Questo \`e ci\`o di cui occorre fare la revisione.

Il comando ``git trac`` ti fornisce 2 opzioni che possono esserti utili (sostituire ``12345`` con il numero di ticket effettivo) se non vuoi utilizzare l'interfaccia web:

* ``git trac print 12345`` mostra il ticket Trac direttamente nel terminale.

* ``git trac review 12345`` effetttua il download del ramo dal ticket e ti
  mostra cosa si sta aggiungendo, analogamente al fare click sul campo "Branch:".

