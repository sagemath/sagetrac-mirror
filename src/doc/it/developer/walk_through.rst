.. _chapter-walkthrough:

===============================
Il processo di sviluppo di Sage
===============================

Questa sezione \`e un riassunto stringato di com'\`e il processo di sviluppo di Sage: vedremo come effettuare modifiche al codice sorgente di Sage e registrarle nel sistema di controlloversione di Git. Nella seguente sezione su :ref:`chapter-git_trac` vedremo come comunicare indietro al progetto Sage questi cambimenti.

Abbiamo anche una pratica paginetta promemoria <http://github.com/sagemath/git-trac-command/raw/master/doc/git-cheat-sheet.pdf>`_ dei comandi pi\`u comunemente usati di Git che puoi stampare e lasciare sulla tua scrivania. Abbiamo anche dei :ref:`tutorial e materiale di riferimento raccomandato<section-git-tutorials>`. In alternativa puoi effettuare un fork e e creare una richiesta di pull su `github <http://github.com/sagemath/sage>`_ che automaticamente andr\`a a prendere il tuo codice e aprir\`a un ticket sul nostro server Trac.

.. _section-walkthrough-setup-git:

Configurare Git
===============

In un modo o in un altro, ``Git`` \`e ci\`o che \`e utilizzato da Sage per tenere traccia delle modifiche. Quindi innanzitutto apri una shell (ad esempio Terminale in Utility su Mac) e verifica che ``Git`` funzioni::

    
    [user@localhost]$ git
    usage: git [--version] [--help] [-C <path>] [-c name=value]
    ...
    The most commonly used git commands are:
       add        Add file contents to the index
    ...
       tag        Create, list, delete or verify a tag object signed with GPG
    
    'git help -a' and 'git help -g' lists available subcommands and some
    concept guides. See 'git help <command>' or 'git help <concept>'
    to read about a specific subcommand or concept.


Non preoccuparti della gigantesca lista di sottocomandi. In realt\`a te ne bastano una manciata per lo sviluppo effettivo, che ti accompagneremo al loro utilizzo in questa guida. Se ricevi un errore "command not found" allora non Git non e' installato. E' giunto il momento di installarlo: vedi :ref:`chapter-git-setup` per istruzioni.

Dal momento che tracciamo anche chi effettua modifiche in Sage con Git, devi dire a Git come voui essere conosciuto. Questo dev'essere fatto solo una volta:

    [user@localhost]$ git config --global user.name "Your Name"
    [user@localhost]$ git config --global user.email you@yourdomain.example.com

Se hai piu' accounts e/o computers usa lo stesso nome su tutti. Questa combinazione di nome ed email finisce nei commits, quindi fallo adesso prima di dimenticartene!


.. _section-walkthrough-sage-source:

Ottenere il codice sorgente
===========================

Ovviamente uno ha bisogno del codice sorgente di Sage per fare sviluppo. Puoi utilizzare la tua installazione locale di Sage, oppure (per iniziare senza Sage) scaricarlo da Github, che e' un mirror pubblico di sola lettura (piu' veloce) del nostro repository Git interno::

    [user@localhost]$ git clone git://github.com/sagemath/sage.git
    Cloning into 'sage'...
    [...]
    Checking connectivity... done.
    

Questo crea una directory di nome ``sage`` contenete i sorgenti per le versioni di Sage corrente stabile e di sviluppo. Dovrai `compilare Sage <http://www.sagemath.org/doc/installation/source.html>`_ per poterla utilizzare.

(Per gli esperti, si noti che il repository presso `git.sagemath.org <http://git.sagemath.org>`_ \`e dove realmente avviene lo sviluppo).

.. _section-walkthrough-branch:

Aggiungere un ramo (branch) locale
==================================

Per iniziare a modificare Sage dobbiamo fare un *ramo* (branch) di Sage. Un branch e' una copia (eccetto che non richiede il doppio dello spazio) del codice sorgente di Sage, dove puoi immagazzinare le tue modifiche al codice sorgente di Sage e che puoi reinviare collegate al relativo ticket su Trac.

\`E facile creare un nuovo branch, semplicemente passa al ramo da cui vuoi partire (che sar\`a ``master``) ed usa il comando ``git branch``::

    [user@localhost sage]$ git checkout master
    [user@localhost sage]$ git branch last_twin_prime
    [user@localhost sage]$ git checkout last_twin_prime

Puoi mostrare una lista di tutte le ramificazioni, cio\`e di tutti i branch, con::

    [user@localhost]$ git branch
      master
    * last_twin_prime

L'asterisco evidenzia su quale branch ti trovi. Senza argomenti ``git branch`` mostra semplicemente una lista di tutti i rami locali con quello corrente marcato con un asterisco. Nota anche che ``git branch`` crea un nuovo ramo, ma non passa ad esso. Per fare questo devi usare il comando ``git checkout``::

    [user@localhost sage]$ git checkout master
    Switched to branch 'master'
    Your branch is up-to-date with 'github/master'.
    [user@localhost sage]$ git branch
    * master
      last_twin_prime
    [user@localhost sage]$ git checkout last_twin_prime
    Switched to branch 'last_twin_prime'

Nota che, se non reinvii esplicitamente le modifiche (``push``) fatte in un branch locale al repository Git remoto, il branch locale sar\`a soltanto sul tuo computer e non sar\`a visibile da nessun altro.

Per evitare di digitare due volte il nome del nuovo branch puoi usare la scorciatoia ``git checkout -b mioNuovoBranch`` per contemporaneamente creare il nuovo branch e passare ad esso.



.. _section_walkthrough_logs:

La storia
=========

\`E sempre una buona idea verificare che stai facendo le tue modifiche sulla versione su cui pensi di essere. Questo comando ti mostra nel dettaglio l'ultimo commit, inclusi i suoi cambiamenti al codice::

    [user@localhost sage]$ git show

Per andare pi\`u a fondo puoi verificare i log::

    [user@localhost sage]$ git log

di default questo mostra una lista di tutti i commit in ordine cronologico inverso. Se scopri che il tuo branch \`e nel posto sbagliato puoi usare il comando ``git reset --hard`` per resettarlo a qualcos'altro; vedi :ref:`section-git-recovery` per dettagli.

Vi sono molti programmi che ti possono aiutare a visualizzare meglio l'albero delle modifiche, quali ad esempio ``tig``.


.. _section-walkthrough-add-edit:

Modificare il codice sorgente
=============================

Una volta che hai il tuo proprio branch, sentiti libero di fare qualunque cambiamento ti piaccia. Alcuni :ref:`capitoli pi\`u avanti <section-writing-code-for-sage>` in questa Guida, ti spiegheranno come il tuo codice dovrebbe essere per andar bene in Sage, e come ci assicuriamo dell'elevata qualita' del codice.
Il comando piu' importante di Git e' probabilmente *status*. Esso ti dice quali file sono cambiati, e come continuare a registrare i cambiamenti:
:

    [user@localhost sage]$ git status
    On branch master
    Changes not staged for commit:
      (use "git add <file>..." to update what will be committed)
      (use "git checkout -- <file>..." to discard changes in working directory)
    
        modified:   some_file.py
        modified:   src/sage/primes/all.py
    
    Untracked files:
      (use "git add <file>..." to include in what will be committed)
    
        src/sage/primes/last_pair.py
    
    no changes added to commit (use "git add" and/or "git commit -a")

Per andare piu' a fondo in cosa e' cambiato nei file puoi usare::

    [user@localhost sage]$ git diff some_file.py


per mostrare le differenze.


.. _section-walkthrough-make:

Ricompilare Sage
================

Una volta che hai fatto qualche modifica ovviamente vorrai ricompilare Sage per provarle. Finche' hai solo modificato la libreria di Sage (cioe' i files Python e Cython nelle sottodirectory di ``src/sare/...``) ti basta eseguire::

    [user@localhost sage]$ ./sage -br

come se stessi installando `Sage from scratch
<http://www.sagemath.org/doc/installation/source.html>`_. Comunque la sola esecuzione di ``make`` ricompilera' solamente i pacchetti che sono cambiati, cosi' che dovrebbe essere molto piu' veloce che compilare Sage la prima volta. Raramenteci sono conflitti con altri pacchetti, o con la versione piu' vecchia gia' installata del pacchetto che hai modificato, nel qual caso dovresti ricompilare tutto usando::

    [user@localhost sage]$ make distclean && make

Inoltre non dimenticare di lanciare i test (vedi :ref:`chapter-doctesting`) e produrre la documentazione (vedi :ref:`chapter-sage_manuals`).


.. _section-walkthrough-commit:

I commit (snapshots)
===================

Ogni volta che hai raggiunto il tuo scopo, o completato un passo importante verso di esso, o semplicemente voui consolidare il lavoro fatto, dovresti effettuare il commit delle modifiche. Un *commit* \`e semplicemente una fotografia dello stato (snapshot) di tutti i file del *repository* (il programma su cui stai lavorando).

A differenza di altri sistemi di controllo versione in Git occorre innanzitutto effettuare lo *stage* dei file modificati, cosa che comunica a Git quali file vuoi che siano parte del prossimo commit::

    [user@localhost sage]$ git status
    # On branch my_branch
    # Untracked files:
    #   (use "git add <file>..." to include in what will be committed)
    #
    #       src/sage/primes/last_pair.py
    nothing added to commit but untracked files present (use "git add" to track)

    [user@localhost sage]$ git add src/sage/primes/last_pair.py
    [user@localhost sage]$ git status
    # On branch my_branch
    # Changes to be committed:
    #   (use "git reset HEAD <file>..." to unstage)
    #
    #   new file:   src/sage/primes/last_pair.py
    #

Una volta che sei soddisfatto della lista dei file in fase "stage" puoi creare un nuovo snapshot con il comando ``git commit``::

    [user@localhost sage]$ git commit
    ... editor opens ...
    [my_branch 31331f7] Added the very important foobar text file
     1 file changed, 1 insertion(+)
      create mode 100644 foobar.txt

Questo aprira' un editor dove potrai scrivere il tuo messaggio di commit. Esso dovrebbe essere in generale una descrizione di una riga, seguito da una riga vuota, seguito da ulteriore testo di spiegazione::

    Added the last twin prime

    This is an example commit message. You see there is a one-line
    summary followed by more detailed description, if necessary.


Poi puoi continuare a lavorare per il tuo prossimo obiettivo, effettuare un altro commit, e cosi' via finche' avrai finito. Finche' non effetui ``git checkout` di un altro branch, tutti i commit che fai saranno parte del branch che hai creato.





