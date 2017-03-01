=================================================
BENVENUTO NELLA GUIDA DELLO SVILUPPATORE DI SAGE!
=================================================

Chiunque usi Sage e' incorraggiato a dare un contributo a Sage in qualche modo, come:

* Aggiungere esempi alla documentazione
* Trovare bachi o errori tipografici
* Correggere un baco
* Implementare una nuova funzione
* Scrivere un utile tutorial su un argomento matematico
* Tradurre un documento esistente in un nuovo linguaggio
* Creare una nuova classe, una nuova libreria C performante, ecc.

Questo documento descrive come scrivere programmi usando Sage, come modificare ed estendere le principali librerie di Sage, e come modificare la documentazione di Sage. Si discute anche di come condividere codice sorgente nuovo o modificato con altri utenti Sage in tutto il mondo.

Qui c'\`e un breve riassunto di ogni parte; per maggiori dettagli, vedi l'indice esteso pi\`u sotto.  Indipendentemenete da dove partirai, buona fortuna e benvenuto nello sviluppo di Sage!

- **Trac server:** tutte le modifiche passano per il `Sage Trac server
  <http://trac.sagemath.org>`_ prima o poi. Contiene le segnalazioni di bachi, le richieste di upgrade, le modifiche in corso, e quelle che sono gi\`a parte di Sage. :ref:`Click here <chapter-sage-trac>` per maggiori  informazioni.

  Molto importante, avrai bisogno di :ref:`creare un account trac
  <section-trac-account>` per poter contribuire.

- **Codice sorgente:** Devi avere la tua propria copia del codice sorgente di Sage per poterlo modificare.
  `Vai qui <http://www.sagemath.org/doc/installation/source.html>`_ per ottenerlo e per istruzioni su come compilarlo.

  Se non hai mai lavorato prima su del software, fa attenzione a
  `prerequisiti alla compilazione
  <http://www.sagemath.org/doc/installation/source.html#prerequisites>`_ sul tuo sistema.

- **Convenzioni:** leggi il nostror :ref:`convenzioni e linee guida
  <section-writing-code-for-sage>` per il codice e la documentazione.

  Per qualunque cosa legata a manuali, tutorial, e linguaggi, :ref:`clicca qui
  <chapter-sage_manuals>`.

- **Git (controllo di revisione):** Per scambiare le tue modifiche con la comunit\`a Sage, avrai bisogno di imparare ad usare il sistema
di controllo di revisione; usiamo il software Git per questo.

  - :ref:`Qui c'\`e <chapter-walkthrough>` una panoramica del nostro flusso di sviluppo.
  - :ref:`Nuovo a Git o al controllo di revisione? <chapter-git_trac>`
  - :ref:`Come installarlo? <section-git-install>`
  - :ref:`COme configurarlo per usarlo con Trac? <section-git-setup-name>`

Git per lo sviluppo in Sage
===========================

Primi passi con Git
-------------------

Sage usa git come sistema di controllo versione.

.. toctree::
   :maxdepth: 3

   git_setup
   walk_through

Il comando git-trac
-------------------

Mettere i tuoi cambiamenti locali su un ticket Trac.

.. toctree::
   :maxdepth: 2

   git_trac

.. _section-git-tricks-and-tips:

Trucchi e suggerimenti su Git
-----------------------------

Quando ``git trac`` non \`e abbastanza.

.. toctree::
   :maxdepth: 2

   manual_git
   git_background
   advanced_git
   workflows

Sage Trac e ticket
==================

Tutte le modifiche al codice sorgente di Sage richiedono un ticket sul
`Sage trac server <http://trac.sagemath.org>`_.

.. toctree::
   :maxdepth: 2

   trac


.. _section-writing-code-for-sage:

Scrivere Codice for Sage
========================

.. toctree::
   :maxdepth: 3

   coding_basics
   reviewer_checklist

Lanciare i test di Sage
-----------------------

.. toctree::
   :maxdepth: 3

   doctesting

Contribuire a Manuali e Tutorial
--------------------------------

.. toctree::
   :maxdepth: 3

   sage_manuals

Dettagli sulla scrittura di codice per Sage
-------------------------------------------

.. toctree::
   :maxdepth: 3

   coding_in_python
   coding_in_cython
   coding_in_other

Pacchettizzare codice di terze parti
------------------------------------

.. toctree::
   :maxdepth: 3

   packaging
   packaging_old_spkgs

Guida dello sviluppatore al Notebook Sage
=========================================

.. toctree::
   :maxdepth: 3

   sagenb/index


Indici e tabelle
================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Quest'opera \`e licenziata sotto la Licenza `Creative Commons Attribution-Share Alike 3.0 <http://creativecommons.org/licenses/by-sa/3.0/>`_.
