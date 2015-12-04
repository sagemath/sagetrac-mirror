=========================
Fare una modifica (patch)
=========================

Hai scoperto un baco o qualcos'altro che vuoi modificare 
in `Sage Notebook`_ |emdash| eccellente!

Hai lavorato un giorno intero per correggerlo |emdash| ancora meglio!

Ce lo vuoi dire |emdash| meglio di tutto!

La maniera pi\`u facile \`e fare una modifica, o patch, o un insieme di patch. 
Qui ti spieghiamo come si fa.

Fare una patch \`e semplice e veloce, ma non \`e parte del nostro 
normale flusso di lavoro. Per cui se pensi di fare qualcosa di pi\`u di una 
patch una-tantum, per favore considera piuttosto il seguente modello 
:ref:`github-development`. Vedi specialmente la parte su "pull requests" in 
:ref:`edit-flow`.

.. _making-patches:

Fare delle patch
================

Panoramica
----------

::

   # d\`i a git chi sei
   git config --global user.email you@yourdomain.example.com
   git config --global user.name "Your Name Comes Here"
   # scarica il repository se non ce l'hai
   git clone git://github.com/sagemath/sagenb.git
   # fa un ramo per lavorare alle tue patch
   cd sagenb
   git branch the-fix-im-thinking-of
   git checkout the-fix-im-thinking-of
   # hack, hack, hack
   # seganala a git i nuovi file che hai fatto
   git add somewhere/tests/test_my_bug.py
   # man mano, fa il commit di ci\`o su cui stai lavorando
   git commit -am 'BF - added tests for Funny bug'
   # hack hack, hack
   git commit -am 'BF - added fix for Funny bug'
   # completa i file delle patch
   git format-patch -M -C master

Puoi attaccare un corto file di patch generato alla `Sage Notebook
mailing list`_ o meglio aprire una richiesta sul sito `Sage Notebook github`_ 
(vedi :ref:`github-development`) e fare copia-incolla della tua patch l\`i 
dentro ad un commento. In entrambi i casi ti ringraziamo sinceramente.

In dettaglio
------------

#. Segnala a git chi sei cos\`i che possa etichettare i commit che fai::

      git config --global user.email you@yourdomain.example.com
      git config --global user.name "Your Name Comes Here"

#. Se non l'hai ancora fatto, clone una copia del repository 
   `Sage Notebook`_ ::

      git clone git://github.com/sagemath/sagenb.git
      cd sagenb

#. Fai in 'ramo di caratteristica' ('feature branch'). Sar\`a qui che lavorerai 
   sul tuo baco. \`e bello e sicuro e ti rimane l'accesso ad una copia non modificata 
   del codice del ramo principale::

      git branch the-fix-im-thinking-of
      git checkout the-fix-im-thinking-of

#. Fai qualche modifica e fanne il commit, man mano::

      # hack, hack, hack
      # D\`i a git di qualunque nuovo file che hai fatto
      git add somewhere/tests/test_my_bug.py
      # fa il commit del lavoro in corso man mano che lo fai
      git commit -am 'BF - added tests for Funny bug'
      # hack hack, hack
      git commit -am 'BF - added fix for Funny bug'

   Nota l'opzione ``-am`` a ``commit``. Il flag ``m`` semplicemente 
   segnala che stai per digitare un messaggio alla linea di comando.
   Il flag ``a`` |emdash| puoi prenderlo per fede |emdash|
   oppure dare un'occhiata a `why the -a flag?`_.

#. Quando hai finito, verifica di avere fatto il commit di tutte le tue modifiche::

      git status

#. Infine, trasforma i tuoi commit in delle patch. Vorrai usare tutti i
   commit del tuo ramo, da quando l'hai clonato da ``master``::

      git format-patch -M -C master

   Avrai a questo punto parecchi file legati ai commit::

      0001-BF-added-tests-for-Funny-bug.patch
      0002-BF-added-fix-for-Funny-bug.patch

   Sebben alcuni progetti vorrebbero che li mandassi alla 
   `Sage Notebook mailing list`_, preferiamo inviare una richiesta di 
   caratteristica tramite l'interfaccia web alla pagina `Sage Notebook github`_. 
   Vedi :ref:`edit-flow` per come creare una "richiesta di pull" una volta che hai 
   creato un account Github.

Quando hai finito, per passare di nuovo alla copia principale del codice, 
semplicemente ritorna al ramo ``master``::

   git checkout master

Passare dal fare patch allo sviluppo
====================================

Se vedi che ai fatto qualche patch, ed ai uno o pi\`u rami di caratteristiche, 
probabilmente vorrai passare alla modalit\`a sviluppo. Puoi farlo con il repository 
che hai gi\`a.

Fa una fork del repository `Sage Notebook`_ su github |emdash| :ref:`forking`.
Poi::

   # fa il checkout ed aggiorna il ramo master dal repo principale
   git checkout master
   git pull origin master
   # rinomina il puntatore al repository principale in 'upstream'
   git remote rename origin upstream
   # punta il tuo repo per leggere / scrivere di default sul tuo fork su github
   git remote add origin git@github.com:your-user-name/sagenb.git
   # fa la push di ogni ramo che hai fatto e vuoi mantenere
   git push origin the-fix-im-thinking-of

Poi puoi, se vuoi, seguire :ref:`development-workflow`.

.. include:: links.inc
