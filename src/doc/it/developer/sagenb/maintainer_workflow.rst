.. _maintainer-workflow:

########################
Workflow del manutentore
########################

Questa pagina \`e per i manutentori |emdash| quelli di noi che fanno il merge 
delle proprie o dell'altrui modifiche nel repository upstream (n.m. di Ipython).

Se sei un manutentore sei completamente al di sopra del livello base di cui tratta 
in :ref:`development-workflow`.

Le instruzioni in :ref:`linking-to-upstream` aggiungono un repository remoto  upstream 
con accesso in sola lettura ad esso. Come manutentore avrai un accesso in lettura e scrittura.

\`E bene che il tuo repository upstream remoto abbia un nome che fa paura, per ricordarti 
che ci accedi in sia in lettura che in scrittura::

    git remote add upstream-rw git@github.com:sagemath/sagenb.git
    git fetch upstream-rw

**********************
Integrare le modifiche
**********************

Supponiamo che hai delle modifiche che devono finire nel ramo principale 
(``upstream-rw/master``).

I cambiamenti sono in qualche ramo su cui stai lavorando al momento. Ad esempio, stai 
dando un'occhiata alle modifiche di qualcuno, in questo modo::

    git remote add someone git://github.com/someone/sagenb.git
    git fetch someone
    git branch cool-feature --track someone/cool-feature
    git checkout cool-feature

Dunque adesso sei sul ramo che contiene le modifiche da incorporare nel codice upstream. 
Il resto della sezione assume che sei su tale ramo.

Alcuni commit
=============

Se ci sono solo pochi commit, considera se fare un rebase al codice upstream::

    # Prendi le modifiche upstream
    git fetch upstream-rw
    # fa il rebase
    git rebase upstream-rw/master

Ricorda che, fa fai un rebase, e poi ne fai la push, ti toccher\`a chiudere ogni 
richiesta di pull su github manualmente, perch\`e github non sar\`a in grado di 
rilevare se \`e gi\`a stato fatto il merge dei cambiamenti.

Una lunga serie di commit
=========================

Se c'\`e una lunga serie di commit collegati, considera di fare, al posto, un merge::

    git fetch upstream-rw
    git merge --no-ff upstream-rw/master

Il merge sar\`a rilevato automaticamnete da github, e dovrebbe chiudere automaticamente 
qualunque richiesta di pull correlata.

Nota il ``--no-ff`` sopra. Questo forza git a fare un commit di merge, piuttosto che 
a fare un rapido avanzamento in avanti (fast-forward), cos\`i che questo insieme di 
commit forma un ramo a parte che poi viene riunito con un merge alla history principale, 
piuttosto che apparire come se fosse stato fatto direttamnete sul ramo principale.

Verificare la history
=====================

Ora, in ogni caso, dovrai verificare che la history ha senso e che hai i commit corretti::

    git log --oneline --graph
    git log -p upstream-rw/master..

La prima linea, sopra, mostra semplicemente la history in maniera compatta, con una 
rappresentazione ASCII del grafo della history. La seconda linea mostra il log dei commit 
esclusi quelli che possono essere raggiunti dal ramo principale (``upstream-rw/master``), 
ed inclusi quelli che possono essere raggiunti dalla HEAD corrente (implicati dalla ``..``
alla fine). Dunque mostra i commit unici di questo ramo confrontati al ramo principale.
L'opzione ``-p`` mostra un prospetto di differenze (detto "diff") per questi commit 
in forma di modifica (patch).

Fare la push al ramo principale
===============================

::

    git push upstream-rw my-new-feature:master

Quest'istruzione fa la push del ramo ``my-new-feature`` in questo repository verso il ramo 
``master`` nel repository ``upstream-rw``.

.. include:: links.inc
