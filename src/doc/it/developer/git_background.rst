.. _chapter-git-background:

==========================
Suggerimenti e riferimenti
==========================

Questo capitolo contiene materiale aggiuntivo sul sistema di controllo versione Git. Non ti serve se utilizzi gli script di sviluppo di Sage. Vedi :ref:`chapter-git-setup` per i passi minimi necessari per fare sviluppo in Sage.





.. _section-git-configuration:

Suggerimenti di configurazione
==============================

La tua configurazione personale di Git \`e salvata nel file ``~/.gitconfig`` nella tua directory home. Eccone un esempio::

    [user]
        name = Your Name
        email = you@yourdomain.example.com

    [core]
        editor = emacs

Puoi modificare questo file direttamente o puoi usare Git per farlo::

    [user@localhost ~] git config --global user.name "Your Name"
    [user@localhost ~] git config --global user.email you@yourdomain.example.com
    [user@localhost ~] git config --global core.editor vim



Alias
-----

Gli alias sono delle scorciatoie personali ai comandi di Git. Ad esempio potresti voler abbreviare ``git checkout`` a ``git co``, oppure ``git diff --color-words`` a ``git wdiff`` (questo comando formatta a colori il file delle differenze). Puoi farlo con::

    [user@localhost ~] git config --global alias.ci "commit -a"
    [user@localhost ~] git config --global alias.co checkout
    [user@localhost ~] git config --global alias.st "status -a"
    [user@localhost ~] git config --global alias.stat "status -a"
    [user@localhost ~] git config --global alias.br branch
    [user@localhost ~] git config --global alias.wdiff "diff --color-words"

I comandi suddetti creeranno una sezione ``alias`` nel tuo file ``.gitconfig`` contenente::

    [alias]
        ci = commit -a
        co = checkout
        st = status -a
        stat = status -a
        br = branch
        wdiff = diff --color-words


Editor
------

Per impostare l'editor da usare per scrivere i messaggi di commit usare::

    [user@localhost ~] git config --global core.editor vim

oppure imposta la variabile di ambiente ``EDITOR``.

Merge
-----

Per forzare dei riassunti quando effettui i merge, puoi scrivere in ``.gitconfig``::

    [merge]
        log = true

Oppure scrivere al prompt dei comandi::

    [user@localhost ~] git config --global merge.log true


.. _section-fancy-log:

Output di log migliorato esteticamente
--------------------------------------

Ecco un alias per avere un output migliorato esteticamente; va posto nella sezione ``alias`` del tuo file ``.gitconfig``::

    lg = log --graph --pretty=format:'%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr) %C(bold blue)[%an]%Creset' --abbrev-commit --date=relative

L'utilizzo di questo alias ``lg`` ti fornisce l'elenco dei cambiamenti (changelog) con un grafico ASCII colorato::

    [user@localhost ~] git lg
    * 6d8e1ee - (HEAD, origin/my-fancy-feature, my-fancy-feature) NF - a fancy file (45 minutes ago) [Matthew Brett]
    *   d304a73 - (origin/placeholder, placeholder) Merge pull request #48 from hhuuggoo/master (2 weeks ago) [Jonathan Terhorst]
    |\
    | * 4aff2a8 - fixed bug 35, and added a test in test_bugfixes (2 weeks ago) [Hugo]
    |/
    * a7ff2e5 - Added notes on discussion/proposal made during Data Array Summit. (2 weeks ago) [Corran Webster]
    * 68f6752 - Initial implimentation of AxisIndexer - uses 'index_by' which needs to be changed to a call on an Axes object - this is all very sketchy right now. (2 weeks ago) [Corr
    *   376adbd - Merge pull request #46 from terhorst/master (2 weeks ago) [Jonathan Terhorst]
    |\
    | * b605216 - updated joshu example to current api (3 weeks ago) [Jonathan Terhorst]
    | * 2e991e8 - add testing for outer ufunc (3 weeks ago) [Jonathan Terhorst]
    | * 7beda5a - prevent axis from throwing an exception if testing equality with non-axis object (3 weeks ago) [Jonathan Terhorst]
    | * 65af65e - convert unit testing code to assertions (3 weeks ago) [Jonathan Terhorst]
    | *   956fbab - Merge remote-tracking branch 'upstream/master' (3 weeks ago) [Jonathan Terhorst]
    | |\
    | |/


.. _section-git-tutorials:


Compendi e tutorial
===================

Ci sono moltissimi tutorial e riassunti dei comandi disponibili online.

Principianti
------------
* `Prova Git <https://try.github.io/levels/1/challenges/1>`_ \`e un tutorial
  di base che puoi seguire sul browser. Se sei nuovo ai sistemi di controllo 
  versione, fa attenzione alla sezione "Advice" verso il fondo.

* `Git magic
  <http://www-cs-students.stanford.edu/~blynn/gitmagic/index.html>`_
  \`e un'ampia introduzione mediamente dettagliata.

* The `git parable
  <http://tom.preston-werner.com/2009/05/19/the-git-parable.html>`_ \`e facile 
  da leggere e spiega i concetti dietro a git.

* `Git foundation
  <http://matthew-brett.github.com/pydagogue/foundation.html>`_
  espande il precedente `git parable`_.

* Sebbene contenga materiale pi\`u avanzato sui rami e le testate separate
  e simili, i riassunti visuali sul merge ed i rami in 
  `Learn Git Branching <http://pcottle.github.io/learnGitBranching/>`_
  possono essere veramente molto utili.


Avanzati
--------
* `Github help <http://help.github.com>`_ ha una ottima serie di how-to.

* Il libro `pro git book <http://git-scm.com/book>`_ \`e un buon libro di approfondimento su Git.

* `Github Training <http://training.github.com>`_ ha un'ottima serie di tutorial nonch\`e di video e screencast.

* Il `git tutorial <http://schacon.github.com/git/gittutorial.html>`_.

* Su `Git ready <http://www.gitready.com/>`_ v'\`e una bella serie di
  tutorial.

* La pagina Git di `Fernando Perez' 
  <http://www.fperez.org/py4science/git.html>`_ contiene molti link e
  suggerimenti.

* Una pagina, buona ma tecnica, su `git concepts
  <http://www.eecs.harvard.edu/~cduan/technical/git/>`_

* `Git svn crash course <http://git-scm.com/course/svn.html>`_: git
  per chi \`e abituato a `subversion
  <http://subversion.tigris.org/>`_

Foglietti di riassunto (Cheat Sheets)
-------------------------------------

* La pagina `git cheat sheet <http://github.com/guides/git-cheat-sheet>`_ 
  d\`a un riassunto dei comandi pi\`u comuni.

* Il `manuale utente di git
  <http://schacon.github.com/git/user-manual.html>`_.



Linee guida Git
===============

Vi sono molti modi di lavorare con Git; qui ci sono dei linka delle regole di buonsenso raccomandate da altri progetti:


* Linus Torvalds spiega la `gestione di git
  <https://web.archive.org/web/20120511084711/http://kerneltrap.org/Linux/Git_Management>`_

* Linus Torvalds spiega il `flusso di lavoro per linux su git
  <http://www.mail-archive.com/dri-devel@lists.sourceforge.net/msg39091.html>`_. Riassunto: 
  utilizza gli strumenti di Git per rendere la storia delle tue modifiche quanto piu\`u pulita
  possibile: effettua i merge dalle modifiche upstream meno che puoi in rami in cui stai
  sviluppando attivamente.


Pagine di manuale online
========================

Puoi scaricare queste sulla tua macchina locale con (ad esempio) ``git help push`` oppure (equivalentemente) ``git push --help``, ma, per convenienza, ecco le pagine di manuale per alcuni comandi comuni:

* `git add <http://schacon.github.com/git/git-add.html>`_
* `git branch <http://schacon.github.com/git/git-branch.html>`_
* `git checkout <http://schacon.github.com/git/git-checkout.html>`_
* `git clone <http://schacon.github.com/git/git-clone.html>`_
* `git commit <http://schacon.github.com/git/git-commit.html>`_
* `git config <http://schacon.github.com/git/git-config.html>`_
* `git diff <http://schacon.github.com/git/git-diff.html>`_
* `git log <http://schacon.github.com/git/git-log.html>`_
* `git pull <http://schacon.github.com/git/git-pull.html>`_
* `git push <http://schacon.github.com/git/git-push.html>`_
* `git remote <http://schacon.github.com/git/git-remote.html>`_
* `git status <http://schacon.github.com/git/git-status.html>`_



