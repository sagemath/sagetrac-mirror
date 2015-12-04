.. _following-latest:

=================================
Seguire il sorgente pi\`u recente
=================================

Queste sono le istruzioni se vuoi seguire il pi\`u recente sorgente di 
*Sage Notebook*, ma per il momento non devi fare alcuno sviluppo.

I passi sono:

* :ref:`section-git-install`
* ottieni una copia locale del repository `Sage Notebook github`_ di git
* aggiorna la tua copia locale di quando in quando

Ottieni la tua copia locale del codice
======================================

Dalla riga di comando::

   git clone git://github.com/sagemath/sagenb.git

Ora hai una copia dell'albero dei sorgenti nella nuova directory ``sagenb``.

Aggiornare il codice
====================

Di quando in quando vorrai scaricare il codice pi\`u aggiornato. Puoi farlo con::

   cd sagenb
   git pull

L'albero dei sorgenti in ``sagenb`` avr\`a ora le modifiche pi\`u recenti dal 
repository iniziale.

.. include:: links.inc
