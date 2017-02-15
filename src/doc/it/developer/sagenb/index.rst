.. _sagenb:

============================================
La guida dello sviluppatore di Notebook Sage
============================================


Lo sviluppo dei notebook Sage attualmente avviene su Github utilizzando il 
sistema di controllo di revisione Git. Il modello di sviluppo per il progetto 
`Notebook Sage`_ \`e un flusso di lavoro `git <http://git-scm.com>`_ e github_ .

Per aggiornarsi al codice sorgente pi\`u recente, lancia il comando seguente, 
dove la ``SAGE_ROOT`` \`e la directory radice dell'installatione di Sage, 
e dove ``hackdir`` \`e una directory che creerai per lavorare alle modifiche 
al codice (non \`e necessario che abbia il nome o la collacazione che ha qui).

.. warning:: Questo creer\`a un nuovo repository sagenb, ignorando qualunque 
   modifica tu abbia gi\`a fatto ai file.

::

    mkdir ~/hackdir
    cd ~/hackdir
    git clone git://github.com/sagemath/sagenb.git sagenb-git
    cd SAGE_ROOT/devel
    rm sagenb
    ln -s ~/hackdir/sagenb sagenb
    cd sagenb
    ../../sage setup.py develop

Ci\`o che viene fatto \`e creare una nuova directory, spostarsi all'interno, 
e creare un clone della versione pi\`u aggiornata dei sorgenti upstream (n.m. 
Ipython) del notebook. Poi cancella il link simbolico ``sagenb`` nella cartella 
Sage e la sostituiamo con un link al clone appena fatto, garantendo cos\`i che 
il notebook ha le dipendenze corrette.

Un vantaggio di avere directory separate per sagenb \`e che potrai nel seguito 
mantenere il codice su cui sviluppi anche quando aggiorni Sage, o persino se 
cancelli accidentalmente la tua installazione di Sage.

Il resto di queste instruzioni \`e della documentazione generica, 
leggermente adattata per aiutarti a sviluppare il notebook usando Git e Github.

La sezione pi\`u importante riguarda il come aggiornare il tuo nuovo 
repository di sorgente sagenb e come creare un "fork" della copia "master", 
cos\`i che tu possa chiedere le che le tue modifiche vengano aggiunte ai 
notebook Sage, detto "pull request"; see :ref:`github-development`.


.. toctree::
   :maxdepth: 2

   following_latest
   patching
   github_development


.. include:: links.inc
