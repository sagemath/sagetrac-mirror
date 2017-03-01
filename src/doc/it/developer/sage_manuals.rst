.. _chapter-sage_manuals:

=================
I manuali di Sage
=================

I manuali di Sage sono scritti in `ReST <http://docutils.sourceforge.net/rst.html>`_
(reStructuredText), e generati con il software `Sphinx
<http://sphinx.pocoo.org>`_:

.. list-table::
   :widths: 4 12
   :header-rows: 1

   * - Name
     - Files

   * - `Tutorial <../tutorial/index.html>`_
     - ``SAGE_ROOT/src/doc/en/tutorial``

   * - `Developer's guide <../developer/index.html>`_
     - ``SAGE_ROOT/src/doc/en/developer``

   * - `Constructions <../constructions/index.html>`_
     - ``SAGE_ROOT/src/doc/en/constructions``

   * - `Installation guide <../installation/index.html>`_
     - ``SAGE_ROOT/src/doc/en/installation``

   * - `Reference manual <../reference/index.html>`_
     - ``SAGE_ROOT/src/doc/en/reference``
       (most of it is generated from the
       source code)

- Inoltre manuali pi\`u specializzati si possono trovare in
  ``SAGE_ROOT/src/doc/en``.

- Alcuni documenti sono stati **tradotti** in altri linguaggi. Per trovarli cambia
  en/ in it/,fr/,es/, de/... See :ref:`section-manuals-names`.

.. _section-manuals-edit:

Modificare la documentazione
============================


(*Vuoi convertire un worksheet di Sage in documentazione?* `Fai click qui
<../thematic_tutorials/sws2rst.html>`_)

Dopo aver modificato qualche file nel tutorial di Sage
(``SAGE_ROOT/src/doc/en/tutorial/``), vorrai vedere il risultato. Per costruirne
una versione **html**, digita::

    sage --docbuild tutorial html

Ora puoi aprire ``SAGE_ROOT/src/doc/output/html/en/tutorial/index.html`` nel
tuo web browser.

- Vuoi aggiungere un **nuovo file** alla documentazione? :ref:`Fai click qui
  <section-add-file>`.

- Per informazioni pi\`u dettagliate sul comando ``--docbuild``, vedi
  :ref:`section-building-manuals`.

**Lancia i doctests:** Tutti i files devono passare i test. Dopo aver modificato un
documento (ad esempio ``tutorial``), puoi lanciare i test con il seguente comando (vedi
:ref:`chapter-testing`)::

    sage -tp SAGE_ROOT/src/doc/en/tutorial/

**Manuale di riferimento:** poich\`e questo manuale \`e perloppi\`u generato dal sorgente
di Sage, dovrai ricompilare Sage per poter vedere i cambiamenti che hai fatto in qualche
documentazione di funzione.  Digita::

    sage -b && sage --docbuild reference html

.. _chapter-sage_manuals_links:

Iperlink
========

La documentazione pu\`o contenere dei link a moduli, classi, o metodi, ad esempio::

    :mod:`link to a module <sage.module_name>`
    :mod:`sage.module_name` (here the link's text is the module's name)

Per link verso classi, metodi, o funzioni, sostituisci **:mod:** con
**:class:**, **:meth:** o **func:** rispettivamente. Vedi la documentazione di `Sphinx'
<http://sphinx.pocoo.org/markup/inline.html>`_.

**Link brevi:** il link ``:func:`~sage.mod1.mod2.mod3.func1``` \`e l'equivalente
di ``:func:`func1 <sage.mod1.mod2.mod3.func1>```: il nome della funzione sar\`a
utilizzato come nome del link, invece del suo path (percorso) completo.

**Nomi locali:** non \`e necessario che i link fra metodi della stessa classe siano
assoluti. Se stai documentando ``method_one``, puoi scrivere
``:meth:`method_two```.

**Namespace globale:** se un oggetto (ad esempio ``integral``) \`e automaticamente importato
da Sage, puoi fare un link ad esso senza specificare il suo percorso completo::

    :func:`A link toward the integral function <integral>`

**Ruoli specifici di Sage:** Sage definisce parecchi specifici *ruoli* (roles):

.. list-table::
   :widths: 4 4 4
   :header-rows: 0

   * - Trac server
     - ``:trac:`17596```
     - :trac:`17596`

   * - Wikipedia
     - ``:wikipedia:`Sage_(mathematics_software)```
     - :wikipedia:`Sage_(mathematics_software)`

   * - Arxiv
     - ``:arxiv:`1202.1506```
     - :arxiv:`1202.1506`

   * - On-Line Encyclopedia of Integer Sequences
     - ``:oeis:`A000081```
     - :oeis:`A000081`

   * - Digital Object Identifier
     - ``:doi:`10.2752/175303708X390473```
     - :doi:`10.2752/175303708X390473`

   * - MathSciNet
     - ``:mathscinet:`MR0100971```
     - :mathscinet:`MR0100971`

** Link http:** puoi copiare/incollare un http link nella documentazione. Se vuoi dare al
link un nome specifico, usa ```link name <http://www.example.com>`_``

**Broken links:** Sphinx pu\`o segnalare i link non funzionanti. Vedi
:ref:`section-building-manuals`.

.. _section-add-file:

Aggiungere un nuovo file
========================

Se hai aggiunto un nuovo file a Sage (ad esempio ``sage/matroids/my_algorithm.py``) e
vuoi che il suo contenuto appaia nel manuale di riferimento, devi aggiungere il suo nome al
file ``SAGE_ROOT/src/doc/en/reference/matroids/index.rst``. Sostituisci
'matroids' con ci\`o che fa al caso tuo.

**La cartella combinat/ :** se il tuo nuovo file appartiene ad una subdirectory di combinat/
la procedura \`e differente:

* Aggiungi il tuo file all'indice memorizzato nel file ``__init__.py`` posto nella
  directory che contiene il tuo file.

* Aggiungi il tuo file all'indice contenuto in
  ``SAGE_ROOT/src/doc/en/reference/combinat/module_list.rst``.

.. _section-building-manuals:

Compilare i manuali
===================

*(Vuoi modificare la documentazione?* :ref:`Fai click qui
<section-manuals-edit>`)

Tutti i manuali di Sage sono compilati utilizzando lo script ``sage --docbuild``.
Il contenuto dello script ``sage --docbuild`` \`e definito nel file
``SAGE_ROOT/src/doc/common/builder.py``.  \`E un sottile wrapper dello script
``sphinx-build`` che fa tutto il lavoro reale. \`E stato fatto per essere una
sostituzione dei Makefile di default generati dallo script ``sphinx-quickstart``.
La forma generale del comando \`e::

    sage --docbuild <document-name> <format>

Ad esempio::

    sage --docbuild reference html

Due comandi **help** che forniscono abbondante documentazione per lo script
``sage --docbuild``::

    sage --docbuild -h # messaggio di help breve
    sage --docbuild -H # uno pi\`u completo

**Formati di output:** Tutti i formati di output supportati da Sphinx (ad esempio pdf)
possono essere usati in Sage. Vedi `<http://sphinx.pocoo.org/builders.html>`_.

**Link spezzati:** per compilare la documentazione e contemporaneamente segnalare i
link spezzati che contiene, usare il flag ``--warn-links``. Nota che Sphinx non
ricompiler\`a un documento che non \`e stato aggiornato, e quindi non riporter\`a i suoi
link spezzati::

        sage --docbuild --warn-links reference html

.. _section-manuals-names:

Nomi dei documenti
------------------

Il nome di documento ``<document-name>`` ha la forma::

    lang/name

dove ``lang`` \`e un codice di lingua di 2 lettere, e ``name`` \`e il nome
descrittivo del documento.  Se la lingua non \`e specificata, allora di default \`e
l'inglese (``en``). I seguenti 2 comandi fanno esattamente la stessa cosa::

    sage --docbuild tutorial html
    sage --docbuild en/tutorial html

Per specificare la versione francese del tutorial, ti basterebbe lanciare::

    sage --docbuild fr/tutorial html


Evidenziazione della sintassi del codice Cython
===============================================

Se vuoi scrivere codice :ref:`Cython <chapter-cython>` in un file ReST, precedi
il blocco di codice con ``.. code-block:: cython`` invece del solito ``::``. Abilita
l'evidenziazione della sintassi in un intero file con ``.. highlight:: cython``.
Ad esempio:

.. code-block:: cython

    cdef extern from "descrobject.h":
        ctypedef struct PyMethodDef:
            void *ml_meth
        ctypedef struct PyMethodDescrObject:
            PyMethodDef *d_method
        void* PyCFunction_GET_FUNCTION(object)
        bint PyCFunction_Check(object)

