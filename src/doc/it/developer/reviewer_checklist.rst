.. nodoctest

.. _chapter-review:

===========================
La check list del revisione
===========================

Tutto il codice che finisce in Sage \`e contro-verificato fra colleghi, per assicurarsi che le convenzioni presentate in questo manuale siano seguite, che ci siano sufficienti esempi nella documentazione e doctest, e per cercare di essere sicuri che il codice faccia, matematicamente, cosa si suppone che faccia.

Se qualcuno (altri che tu) ha inviato sul Trac server un ramo git per un ticket, tu puoi farne la revisione! Controlla il branch diff (l'elenco delle modifiche), cliccandoci sopra, per vedere se ha senso. Scaricalo (vedi :ref:`revisioni <section-git_trac-review>`) e compila Sage con il nuovo ramo incluso, quindi fatti delle domande come le seguenti:

* il nuovo codice sorgente ha senso?

* quando lo esegui in Sage, risolve il problema riportato nel ticket relativo?

* introduce qualche nuovo problema?

* \`e documentato a sufficienza, incluse sia le spiegazioni che i doctest? Tutto il codice in Sage deve avere dei doctest, quindi se l'autore del ticket cambia del codice che non aveva un doctest prima, la nuova versione deve includerne uno. In particolare tutto il nuovo codice deve essere provato con dei doctest, al 100%. Usa il comando ``sage -coverage <files>`` per vedere la percentuale di copertura di ``<files>``.

* in particolare, vi \`e un doctest che illustri che il baco \`e stato risolto? Se una funzione dava un risultato sbagliato e questo ticket la corregge, allora dovrebbe includere un doctest che illustri il suo successo. La doctring relativa dovrebbe includere il numero di ticket, ad esempio ``vedi :trac:'12345'``.

* se il ticket afferma di accellerare qualche calcolo, contiene anche degli esempi di codice per mostrare quanto afferma? Il ticket dovrebbe analizzare esplicitamente qual'\`e la velocit\`a prima di applicare la patch e qual'\`e dopo, e spiegare qual'\`e il guadagno di tempo.

* il manuale di riferimento compila senza errori? Puoi provare il manuale di riferimento utilizzando il comando ``sage -docbuild reference html`` per produrne la versione in HTML. Anche la versione PDF dev'essere prodotta senza errori: usa il comando ``sage -docbuild reference pdf`` per provarlo. Tale comando richiede che tu abbia installato Latex sul tuo PC.

* i doctest passano tutti senza errori? \`E difficile predire quali componenti di Sage verranno toccati da una data patch, e dovresti lanciare i test dell'intera libreria (inclusi quelli etichettati "#long") prima di segnalare esito positivo alla revisione. Poi effettuare il test della libreria Sage con "make ptestlong". Vedi (link Effettuare i doctest della Sage Library) per maggiori informazioni.
â€¢ il codice e la documentazione seguono le convenzioni documentate nelle sezioni seguenti?

* il codice e la documentazione seguono le convenzioni documentate nelle sezioni seguenti: :ref:`convenzioni di Sage <chapter-code-basics>`, :ref:`convenzioni di Python <chapter-python>`, :ref:`convenzioni di Cython <chapter-cython>`?

Se la risposta a queste ed altre domande ragionevoli simili \`e s\`i, allora puoi dar esito positivo alla revisione. Sulla pagina principale del ticket scrivi un commento nello spazio a ci\`o riservato, spiegando la tua revisione. Se ritieni di non avere abbastanza esperienza per fare ci\`o, scrivi un commento che spieghi che cosa hai verificato, e concludi chiedendo se qualcuno con pi\`u esperienza pu\`o dare un'occhiata. Se pensi che ci siano problemi con la patch, spiegali nel riquadro dei commenti e cambia lo status a ``need work`` (richiede lavorazione). Guarda altri ticket su Trac per vedere come si fa.
Se tu stesso cambi la patch, devi fare un commit sotto il tuo nome e segnarlo come patch conseguente ad una revisione. Questa va anch'essa sottoposta a revisione, per esempio dall'autore della patch originale.

Verificato quanto sopra puoi cambiare lo status del ticket (vedi
:ref:`section-trac-ticket-status`):

- **positive review** (revisione positiva): se le risposte alle domande
  suddette sono *"s\`i"*, puoi impostare il ticket a ``positive_review``.
  Aggiungi il tuo nome completo al campo"reviewer" field (vedi
  :ref:`section-trac-fields`).

- **needs_work** (richiede lavorazione): se qualcosa non \`e come
  dovrebbe, scrivi una lista di tutti i punti che vanno risolti in un
  commento e cambia lo status del ticket a ``needs_work``.

- **needs_info** (richiede informazioni): se qualcosa non ti \`e chiaro e
  ti impedisce di procedere nella revisione, fai le tue domande e poni
  lo status del ticket a ``needs_info``.

- Se non sai cosa fare, ad esempio non ti sembra di avere abbastanza
  esperienza per prendere una decisione finale, spiega quanto hai gi\`a
  fatto nel tuo commento e chiedi se qualcun altro pu\`o dare un'occhiata.

**Reviewer's commit** (commit del revisore): se puoi risolvere il problema tu stesso, puoi fare una commit nel tuo nome e segnarla come patch del revisone. Per saper come
:ref:`fai click qui <section-git_trac-push>` (git trac) oppure :ref:`qui
<section-git-push>` (git da solo). anche questo contributo deve essere sottoposto a revisione, ad esempio dall'autore della patch originale.

Per maggiori consigli sulle revisioni, vedi [WSblog]_.

.. note::

   "Il meglio \`e nemico del bene"
Lo scopo della revisione \`e assicurarsi che le lineee guida sul codice di Sage siano seguite e che l'implementazione sia matematicamente corretta. Per cortesia astieniti dalla richiesta di funzionalit\`a aggiuntive e discussioni su implementazioni alternative che non siano mirate. Se vuoi che la patch sia scritta diversamente, il tuo suggerimento dev'essere una richiesta chiara e fattibile.


REFERENCES:

.. [WSblog] William Stein, How to Referee Sage Trac Tickets,
   http://sagemath.blogspot.com/2010/10/how-to-referee-sage-trac-tickets.html
   (Caveat: mercurial was replaced with git)
