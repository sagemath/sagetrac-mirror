.. _chapter-git-setup:

====================
Mettere in piedi Git
====================

Per poter lavorar sul codice sorgnte di Sage devi::

* avere una installazione di Git funzionante (vedi :ref:`section-git-install`).
  In realt\`a Sage include gi\`a Git, come \`e spiegato pi\`u avanti. Comunque
  \`e preferibile avere un'installazione di Git direttamente sul sistema,se non
  altro per risparmiarsi troppe digitazioni.

* configurare Git per usare il tuo nome ed indirizzo email per i commit (vedi :ref:`section-git-setup-name`.
  Gli script per lo sviluppo di Sage te li richiederanno se non li immetti. Ma,
  specialmente se intendi usare Git anche per altri progetti in futuro, ti
  converrebbe veramente configurare Git.

Il capitolo :ref:`chapter-git-background` contiene ulteriori informazioni su Git che posso essere utili ma qualcuno ma non e' indispensabile conoscere.


.. _section-git-install:

Installare Git
--------------

Innanzitutto prova a scrivere ``git`` sulla riga di comando. La maggior parte delle distributioni Linux lo installer\`a di default se ci sono altri strumenti di sviluppo installati. Se tale comando d\`a errore, usare i seguenti per installare Git:

Debian / Ubuntu
    ``sudo apt-get install git-core``

Fedora
    ``sudo yum install git-core``

Windows
    Scarica ed installa `msysGit <http://code.google.com/p/msysgit/downloads/list>`_

OS X
    Usa il `git OSX installer
    <https://sourceforge.net/projects/git-osx-installer/files/>`_.  Se hai
    un vecchio Mac, bada di prendere la versione corretta. In alternativa
    puoi ottenerlo dai Command Line Tools o addirittura cercando di usare
    ``git`` e le istruzioni seguenti.


Sage include Git, ma ovviamente c'\`e un problema del ``prima l'uovo o la gallina`` quando si vuole verificare il codice sorgente di Sage dal suo repository Git, ma si puo' sempre fare il download di un archivio Tar contenente il sorgente di Sage oppure la distribuzione binaria. Poi si puo' lanciare Git da riga di comando digitando ``sage -git``. Cos\`i, ad esempio, ``git help`` diventa ``sage -git help`` e cos\`i via. Si noti che gli esempi nella guida dello sviluppatore presuppongono che si abbia una installazione di Git a livello di sistema.

Ulteriori risorse per aiuto nell'installazione sono le seguenti:

* `Capitolo 2 del libro di Git
  <http://book.git-scm.com/2_installing_git.html>`_

* La `git homepage <http://git-scm.com>`_ per le informazioni pi\`u
  recenti.

* `Github install help pages <http://help.github.com>`_


.. _section-git-setup-name:

Il tuo nome ed email
--------------------

Il messaggio di commit di qualunque cambiamento contiene il tuo nome ed il tuo indirizzo email per riconoscere il tuo contributo ed avere un punto di contatto se ci saranno domande in futuro; immetterli \`e richiesto, se vuoi condividere le tue modifiche. Il modo pi\`u semplice per fare ci\`o \`e dalla riga di comando::

    [user@localhost ~] git config --global user.name "Your Name"
    [user@localhost ~] git config --global user.email you@yourdomain.example.com

Questo scrivera' le impostazioni nel tuo :ref:`file di configurazione di Git <section-git-configuration>` con il tuo nome ed la tua email::

    [user]
        name = Your Name
        email = you@yourdomain.example.com

Naturalmente devi sostituire ``Tuo nome`` e ``tuoIndirizzo@tuoDominio.boh`` con il tuo nome reale ed il tuo indirizzo email reale.

