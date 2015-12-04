.. _set-up-fork:

=======================
Organizzare il tuo fork
=======================

Prima segui le istruzioni in :ref:`forking`.

Panoramica
==========

::

   git clone git@github.com:your-user-name/sagenb.git
   cd sagenb
   git remote add upstream git://github.com/sagemath/sagenb.git

In dettaglio
============

Clona il tuo fork
-----------------

#. Clona il tuo fork sul computer locale con ``git clone
   git@github.com:your-user-name/sagenb.git``
#. Investiga. Cambia directory e vai al tuo nuovo repo: ``cd sagenb``. Poi 
   ``git branch -a`` per mostrare tutti i tuoi rami. Otterrai qualcosa tipo::

      * master
      remotes/origin/master

   Questo ti dice che attualmente sei sul ramo ``master``, e che hai anche 
   una connessione ``remote`` a ``origin/master``.
   Che repository remoto \`e ``remote/origin``? Prova ``git remote -v`` per 
   vedere gli URL del remoto. Punteranno al tuo fork su github.

   Ora vuoi connetterti al repository upstream `Sage Notebook github`_ , cos\`i 
   da poter fare il merge dei cambiamenti dal ramo principale.

.. _linking-to-upstream:

Collegare il tuo repository al repo upstream
--------------------------------------------

::

   cd sagenb
   git remote add upstream git://github.com/sagemath/sagenb.git

Qui ``upstream`` \`e solo un nome arbitrario che stiamo usando per riferirci 
al repository principale `Sage Notebook`_ su `Sage Notebook github`_.

Nota che abbiamo usato ``git://`` come URL invece di ``git@``. L'URL 
``git://`` \`e in sola lettura. Questo significa che non possiamo scrivere 
accidentalmente (o deliberatamente) sul repo upstream, e lo useremo solo per 
fare il merge di esso nel nostro codice.

Solo per la tua soddisfazione, verifica tu stesso che ora hai un nuovo 'remote', 
con ``git remote -v show``, che ti dovrebbe dare all'incirca quanto segue::

   upstream     git://github.com/sagemath/sagenb.git (fetch)
   upstream     git://github.com/sagemath/sagenb.git (push)
   origin       git@github.com:your-user-name/sagenb.git (fetch)
   origin       git@github.com:your-user-name/sagenb.git (push)

.. include:: links.inc

