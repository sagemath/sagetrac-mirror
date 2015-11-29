.. _chapter-workflows:

====================
Sviluppo distribuito
====================

Git \`e uno strumento per scambiare dei commit (organizzati in rami) con altri sviluppatori. \`E un sistema di controllo di revisione distribuito, per cui non vi \`e la nozione di server centrale. Il Trac server di Sage \`e solo uno dei possibili repository remoti, dal tuo punto di vista. Questo ti permette di utilizzare e sperimentare vari modi di interagire con gli altri sviluppatori. In questo capitolo abbiamo descritto alcuni modi comuni di sviluppare in Sage.

Per semplicit\`a assumiamo che 2 sviluppatori (Aice e Bob) stiano collaborando ad un ticket. Il primo passo di apertura di n ticket \`e sempre lo stesso, e pu\`o essere fatto sia da Alice che da Bob che da qualcun altro.




Flusso di lavoro semplice
=========================

.. image:: static/flowchart.*
    :align: center


1. Alice crea un :ref:`nuovo ramo locale <section-walkthrough-branch>` e
   :ref:`fa commit <section-walkthrough-commit>` dei cambiamenti al sorgente di Sage.

2. Alice :ref:`fa upload del suo ramo <section-git_trac-push>` sul trac
   server. Questo riempie il campo "Branch:" con il nome del suo ramo remoto ``u/alice/description``.

3. Bob :ref:`scarica il ramo di Alice <section-git_trac-checkout>`,
   guarda il sorgente, e lascia un commento sul ticket riguardo un
   errore nel codice di Alice.

4. Alice corregge il baco sul suo ramo corrente, e fa upload del ramo
   modificato.

5. Bob :ref:`scarica gli aggiornamenti di Alice <section-git_trac-pull>`
   e ne fa la revisione.

6. Quando Bob \`e soddisfatto, imposta il ticket per revisione positiva.
   Il campo "Author:" \`e posto al nome completo di Alice, ed il campo
   "Reviewer:" \`e posto al nome completo di Bob.

In alternativa, Bob potrebbe voler lui stesso fare dei cambiamenti. In tal caso avremmo allora

3. Bob :ref:`scarica il ramo di Alice <section-git_trac-checkout>`, fa
   le modifiche, e ne :ref:`fa commit <section-walkthrough-commit>` al
   suo ramo locale.

4. Bob :ref:`fa upload del suo ramo <section-git_trac-push>` al server
   trac. Questo riempie il campo "Branch:" con il suo nome di ramo
   remoto ``u/bob/description``.

5. Alice :ref:`scarica il ramo di Bob <section-git_trac-checkout>` e fa
   la revisione delle sue modifiche.

6. Quando Alice \`e soddisfatta, imposta il ticket per revisione
   positiva. Se entrambi i contributi sono di dimensioni comparabili,
   allora i campi "Author:" e "Reviewer:" sono impostati con entrambi
   i nomi completi di Alice e di Bob.




Repository pubblico
===================

In aggiunta ai rami utente (``u/<user>/<description>`` sul Trac server di Sage con ``<user>`` rimpiazzato dal tuo nome utente di Trac) su cui puoi scrivere solo tu, puoi anche creare un ramo pubblico su cui chiunque con un account Trac pu\`o scrivere. Questi iniziano con 
``public/`` pi\`u qualche descrizione. Per evitare collisioni nei nomi di ramo \`e una buona idea includere il tuo nome di utente Trac nel nome del ramo, per cui si raccomanda di usare ``public/<user>/<description>`` come nome di ramo. Ora tutti gli autori di ticket effettuano ``push`` allo stesso ramo remoto.

1. Alice crea un :ref:`nuovo ramo locale <section-walkthrough-branch>` 
   ed effettua :ref:`commit <section-walkthrough-commit>` di qualche
   modifica alla libreria Sage.

2. Alice :ref:`fa uploads del suo ramo <section-git_trac-push>` come un
   ramo pubblico sul server Trac. Questo riempie il campo "Branch:" con
   il nome di ramo remoto ``public/alice/description``.

3. Bob :ref:`scarica il ramo di Alice <section-git_trac-checkout>` ed
   effettua le modifiche sulla sua copia locale.

4. Bob :ref:`fa commit <section-walkthrough-commit>` delle modifiche
   al suo ramo locale del sorgente di Sage.

5. Bob fa upload delle sue modifiche al repository remoto congiunto::

       [bob@localhost sage]$ git push trac local_branch:public/alice/description

6. Alice :ref:`scarica le modifiche di Bob <section-git_trac-pull>`, fa
   degli altri cambiamenti, poi ne fa commit e push sul server Trac.

7. Charlie fa la revisione della versione finale, ed imposta il ticket
   a revisione positiva. Il campo "Author:" \`e impostato ai nomi
   completi di Alice e Bob, ed il campo "Reviewer:" \`e impostato al
   nome completo di Charlie.




GitHub
======

Ancora un altro possibile flusso di lavoro \`e utilizzare GitHub (o qualunque altro repository git di terze parti) per modificare in modo collaborativo il tuo nuovo ramo, ed inviare il risultato a Trac una volta che tu ed i coautori del tuo ticket siete soddisfatti.


Fork
----

Il primo passo \`e creare il tuo fork (it. diramazione) del repository di Sage; fai semplicemente click su "Fork" sul `Sage GitHub repository
<https://github.com/sagemath/sage>`_. Poi aggiungilo come un repository remoto del tuo repository Sage locale. Nel seguito utilizzeremo l'etichetta "github" per questo remote repository, sebbene sei ovviamente libero di usare un altro nome::

    $ git remote add github git@github.com:github_user_name/sage.git
    $ git remote -v
    github      git@github.com:github_user_name/sage.git (fetch)
    github      git@github.com:github_user_name/sage.git (push)
    trac        git@trac.sagemath.org:sage.git (fetch)
    trac        git@trac.sagemath.org:sage.git (push)
    $ git fetch github
    remote: Counting objects: 107, done.
    remote: Compressing objects: 100% (63/63), done.
    remote: Total 74 (delta 41), reused 40 (delta 10)
    Unpacking objects: 100% (74/74), done.
    From github.com:github_user_name/sage
    * [new branch]      master     -> github/master
    

Sviluppare
----------

Ora utilizzi il repository di github per sviluppare il ramo del tuo ticket; innanzitutto crea un nuovo ramo::

    $ git checkout -b my_branch --track github/master
    Branch my_branch set up to track remote branch master from github.
    Switched to a new branch 'my_branch'
    $ git push github my_branch
    Total 0 (delta 0), reused 0 (delta 0)
    To git@github.com:github_user_name/sage.git
     * [new branch]      my_branch -> my_branch

A causa dell'opzione ``--track``, il comando ``git pull`` avr\`a come default di scaricare le modifiche del tuo coautore dal tuo ramo github. In alternativa, puoi creare un nuovo ramo sulla webpage del tuo fork di
GitHub.

A questo punto puoi utilizzare il flusso di sviluppo GitHub che preferisci. In particolare, hai la possibilit\`a di

* Dare ai tuoi coautori permessi in scrittura sul tuo fork github. Ogni
  autore fa le modifiche e fa commit della propria copia locale e poi
  effettua congiuntamente una push al suo ramo github.

* Ogni coautore ha la possibilit\`a di creare il proprio fork ed inviare
  a te (l'autore leader) le richieste di pull al tuo fork GitHub.

* Usare le funzionalit\`a di modifica e commit della pagina web di
  GitHub, cos\`i da poter effettuare cambiamenti senza utilizzare mai la
  tua macchina locale.


Push a Trac
-----------

Quando dei soddisfattocon il tuo ramo, ne effettui la push al server Trac di Sage::

    $ git push trac HEAD:u/user/description

e poi compili il campo "Branch" nella descrizione del ticket Trac come
spiegato in :ref:`section-git-push`.

