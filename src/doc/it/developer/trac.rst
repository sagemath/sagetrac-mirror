.. _chapter-sage-trac:

====================
The Sage Trac Server
====================

Tutto lo sviluppo avviene tramite il `Sage Trac server <http://trac.sagemath.org>`_. Lo scopo del server Trac \`e::

1. fornire un luogo di discussione sulle esigenze e costituire una 
   memoria permanente.

2. fornire un repository di codice sorgente e di tutti i cambiamenti 
   proposti.

3. collegare queste due funzionalit\`a.

C'\`e anche un `wiki <http://trac.sagemath.org/wiki>`_ per pagine web organizzative pi\`u generali, come i workshop di sviluppo di Sage.

Cos\`i se trovi un baco in Sage, se hai nuovo codice sorgente da proporre, vuoi verificare del nuovo codice sorgente proposto ma non ancora incluso in Sage, o se hai delle correzioni per la documentazione, dovresti effettuare un post sul server Trac. Gli elementi su tale server sono detti *ticket* (it. biglietti), e chiunque pu\`o guardare o fare ricerche fra i tickets. Per una lista dei cambiamenti recenti visita la pagina `Sage trac timeline <http://trac.sagemath.org/timeline>`_.

.. _section-trac-account:

Ottenere un account
===================

Devi innanzitutto aprire un account se vuoi *cambiare* qualcosa sul server Trac, anche se vuoi solo commentare un ticket. Parte del processo \`e fatto per provare che tu sia un essere umano, in modo da minimizzare lo spam. Per ottenere un account leggi il manuale di sviluppo (questo documento) e poi manda un'email a ``sage-trac-account@googlegroups.com`` che contenga quanto segue::

* il tuo nome completo
* il tuo username preferito
* un'email di contatto
* le ragioni per cui hai bisogno di un account Trac

Il tuo account Trac ti garantir\`a accesso anche al `wiki di Sage
<wiki.sagemath.org>`_. Accertati di conoscere il processo di revisione, e le procedure per aprire e chiudere i ticket, prima di effettuare cambiamenti. Il rimanente di questo capitolo contiene varie lineee guida su come usare il server Trac.

Autenticazione Trac attraverso SSH
==================================

Ci sono due strade per proare al server trac di essere che si pretende di essere. Il primo \`e di loggarsi in Trac usando un nome utente / password per cambiare le pagine relative ai ticket. La seconda e' la crittografia a chiave pubblica usata da Git quando si copiano nuovi files nel repository. Questa sezione ti mostrera' come configurarli entrambi.

Generare e fare l'upload con le chiavi SSH
------------------------------------------

L'installazione di Git sul server di sviluppo usa SSH per decidere se e dove ti e' permesso di fare upload di codice. Nessuna chiave SSH \`e richiesta per segnalare un baco o fare un commento su un ticket, ma non appena vuoi inviare del codice devi fornire a Trac la parte pubblica della tua chiave personale. Nelle versioni recenti di Sage puoi usare Sage stesso per generare e fare upload di una chiave SSH.

.. skip   # do not doctest

::

    sage: dev.upload_ssh_key()
    The trac git server requires your SSH public key to be able to identify you.
    Upload "/home/vbraun/.ssh/id_dsa.pub" to trac? [Yes/no] y
    Trac username: user
    Trac password:
    Your key has been uploaded.

Puoi anche generare manualmente una chiave SSH e farne l'upload su Trac. Questo e' descritto nelle due sezioni seguenti.

Generare manualmente le proprie chiavi SSH
------------------------------------------

Se non hai ancora una chiave privata, puoi crearla con lo strumento ``ssh-keygen`` ::

    [user@localhost ~]$ ssh-keygen
    Generating public/private rsa key pair.
    Enter file in which to save the key (/home/user/.ssh/id_rsa):
    Enter passphrase (empty for no passphrase):
    Enter same passphrase again:
    Your identification has been saved in /home/user/.ssh/id_rsa.
    Your public key has been saved in /home/user/.ssh/id_rsa.pub.
    The key fingerprint is:
    ce:32:b3:de:38:56:80:c9:11:f0:b3:88:f2:1c:89:0a user@localhost
    The key's randomart image is:
    +--[ RSA 2048]----+
    |  ....           |
    |   ..            |
    |   .o+           |
    | o o+o.          |
    |E + .  .S        |
    |+o .   o.        |
    |. o   +.o        |
    |      oB         |
    |     o+..        |
    +-----------------+

Questo generer\`a una nuova chiave privata casuale RSA nella cartella ``.ssh`` nella tua home directory. Di default avrai i file:

``~/.ssh/id_rsa`` 
  Qui ci sar\`a la tua chiave privata: mantienila al sicuro e non darla a **nessuno**.
``~/.ssh/id_rsa.pub``
  Qui ci sar\`a la tua chiave pubblica corrispondente: questa e solo questa puo' essere comunicata a terze parti senza problemi.

Lo strumento ``ssh-keygen`` ti permetter\`a di generare una chiave con un nome di file differente, o di proteggerla con una password. In base a quanta fiducia riponi nel tuo proprio computer o nel tuo amministratore di sistema, puoi lasciare la password vuota per essere in grado di fare login senza alcun intervento umano ulteriore.

Se hai degli account su pi\`u computer puoi usare le chiavi SSH per collegarti. Semplicemente copia la chiave **pubblica** (quella che finisce con ``.pub``) in ``~/.ssh/authorized_keys`` sulla macchina remota e verifica che solo tu abbia su tale file i permessi di lettura e scrittura. Voil\`a, la prossima volta che fai login via SSH in quella macchina non dovrai pi\`u immettere la password.


.. _section-trac-ssh-key:

Collegare manualmente la tua chiave pubblica al tuo account Trac
----------------------------------------------------------------

Il server Trac di Sage ha bisogno di conoscere una delle tue chiavi pubbliche. Puoi farne l'upload nelle preferenze, nel modo seguente::

1. Vai a http://trac.sagemath.org

2. Fai login con il tuo nome utente e password

3. Fai click su "Preferenze"

4. Vai alla scheda "Chiavi SSH"

5. Copia il contenuto del file contenente la tua chiave pubblica (ad esempio ``~/.ssh/id_rsa.pub``)

6. Clicca su "Salva cambiamenti"

Nota che questo **non** ti permette di collegarti via ssh a qualunque account su Trac, e' soltanto usato per autenticarti alla installazione gitolite su Trac. Puoi verificare di essere stato autenticato correttamente inviando qualche comando gitolite di base, ad esempio::

    [user@localhost ~]$ ssh git@trac.sagemath.org info
    hello user, this is git@trac running gitolite3 (unknown) on git 1.7.9.5

     R W      sage
    [user@localhost ~]$ ssh git@trac.sagemath.org help
    hello user, this is gitolite3 (unknown) on git 1.7.9.5

    list of remote commands available:

        desc
        help
        info
        perms
        writable

.. _trac-bug-report:

Segnalare bachi
===============

Se pensi di aver trovato un baco in Sage, dovresti innanzitutto cercare nei post dei nostri Google groups notizie relative a tale possibile baco (\`e possibile infatti che il problema che hai incontrato sia gia' stato discusso e /o risolto):

  * ``sage-devel``: `<http://groups.google.com/group/sage-devel>`_
  * ``sage-support``: `<http://groups.google.com/group/sage-support>`_

  Similmente puoi effettuare una ricerca su :ref:`chapter-sage-trac` per vedere se qualcun altro ha gia' aperto un ticket relativo a tale baco.

Se non trovi niente, e non sei sicuro di aver trovato un baco, domanda si esso su `sage-devel <http://groups.google.com/group/sage-devel>`_. Ti potrebbe essere richiesto di aprire un nuovo ticket sul server Trac (segui la sezione :ref:`section-trac-new-ticket`). Come detto sopra, hai bisogno di un account per fare ci\`o. Per segnalare un baco, fai login e clicca su "Nuovo ticket". In "Short summary" (it. riassunto breve) scrivi una riga di breve spiegazione, entrando nel dettaglio nello spazio apposito sotto. Dovresti includere almeno un esempio esplicito e **riproducibile** che dimostri il baco, con tutti i passi da seguire per causarlo. Dovresti anche includere la **versione** di Sage (ed eventuali pacchetti rilevanti) che stai usando, e informazioni sul **sistema operativo**, cercando di essere preciso il pi\`u possibile (32-bit, 64-bit, ...).

Fra "riassunto" e "descrizione completa" c'\`e un'opzione per scegliere il tipo di biglietto: "difetto", "miglioria" o "da fare" (task). Usa il buon senso: un baco dovrebbe probabilmente essere segnalato con tipo "difetto".

Inoltre scegli un componente in cui rientra il tuo baco: questo \`e spesso ovvio. Se il tuo baco ha a che fare con l'implementazione di Sage del calcolo differenziale scegli "calculus". Se non \`e ovvio, fai del tuo meglio. 

Scegli un "milestone"; se non sei sicuro su cosa scegliere, scegli semplicemente il numero di verione di Sage dal menu (ad esempio "sage-5.10").

Digita qualche parola chiave utile.

Nello riquadro etichettato "assegna a" digita "somebody (it. chiunque)se non sai cosa mettere d'altro.

Premi il bottone "anteprima" per verificare che tutto sia a posto, poi primi "invia ticket".

Se non hai un account sul sistema Trac per fare direttamente le segnalazioni, ugualmente dovresti segnalare ogni possibile baco alla mailing-list ``sage-devel`` presso ``sage-devel@googlegroups.com``. La lista \`e moderata per quanto riguarda gli utenti nuovi, e richiede di effettuare un'iscrizione. Nella segnalazione del baco su ``sage-devel`` assicurati di includere le segnenti informazioni:

* sistema operativo: sii il piu' preciso possibile ed indica l'architettura (32bit, 64bit,...)
* versione bacata: l'esatto numero di versione ed il pacchetto scaricato (sorgente, precompilato, immagine di
  macchina virtuale, o aggiornamento rispetto ad una precedente versione (quale ?))
* fornisci un esempio riproducibile e/o definisci i passaggi per riprodurre il comportamento errato.

Grazie in anticipo per la segnalazione di bachi per migliorare Sage in futuro !

.. _section-trac-new-ticket:

Linee guida sulla segnalazione di bachi
=======================================

Oltre a segnalare i bachi (vedi :ref:`trac-bug-report`), dovresti anche aprire un ticket se hai del nuovo codice che estende le capacit\`a di Sage. Se hai una richiesta di funzionalit\`a, prima inizia una discussione su ``sage-devel``, e poi, se ti sembra che tutti siano sostanzialmente daccordo che la tua sia una buona idea, apri un ticket che la descrive.

Quando pensi di aprire un nuovo ticket, **prima** per favore tieni presenti i seguenti punti::

* prima di aprire un ticket, accertati che nessun altro ha gi\`a aperto un ticket sullo stesso argomento, o su un argomento simile.

* \`e meglio aprire pi\`u ticket su questioni specifiche che uno contenente parecchie questioni. Invero un ticket contenente parecchie questioni pu\`o essere molto problematico e dovrebbe essere evitato.

* sii preciso: se la tal cosa non funzona su OSX ma funziona su Linux, tu scrivilo nel titolo. Usa l'opzione dell'immissione di parole chiave, cosicch\`e le ricerche raccolgano anche la tua richiesta.

* il problema descritto nel ticket dev'essere risolubile. Ad esempio sarebbe sciocco aprire un ticket il cui scopo fosse "rendere Sage il miglior programma per la matematica del mondo". Non c'\`e un metro di giudizio per questo e pu\`o essere molto soggettivo.

* nelle segnalazioni di bachi la descrizione del ticket deve contenere le informazioni descritte a :ref:`trac-bug-report`.

* se utile aggiungi degli URL ad altre informazioni o sezioni di mailing-list importanti per il problema che stai segnalando.

Se non sei sicuro di cosa stai facendo, lascia il campo "milestone" al suo default.
**Prima** di creare il ticket, pu\`o esserti utile leggere :ref:`section-trac-fields`.

.. _section-trac-fields:

I campi dei ticket
==================

Quando apri un nuovo ticket o cambi un ticket esistente, troverai parecchi campi da imputare. Eccone un panoramica  (per il campo 'status', vedi :ref:`section-trac-ticket-status`)::

* **Reported by** (riportato da): l'account su Trac di chi ha creato il ticket. Non pu\`o essere cambiato.

* **Owned by** (di propriet\`a di): l'account su Trac del proprietario del ticket, di default la persona incaricata della manutenzione del componente: generalmente questo campo non \`e utilizzato.

* **Type** (tipo): uno fra ``enhancement`` (miglioramento cio\`e nuova funzionalit\`a), ``defect`` (difetto cio\`e un baco), or ``task`` (obiettivo, usato raramente).

* **Priority** (priorit\`a): la priorit\`a del ticket. Tieni presente che l'etichetta "blocker" (bloccante) dovrebbe essere usata molto raramente.

* **Milestone** (pietra miliare): i milestone sono generalmente degli scopi da realizzare nel lavoro verso una nuova release del software. In Trac utilizziamo i milestone invece delle release. Ogni ticket deve essere assegnato ad una milestone. Se non sei sicuro, assegnalo alla milestone corrente.

* **Component** (componente): nella lista dei componenti di Sage, seleziona quella che pi\`u si avvicina al tuo ticket.

* **Keywords** (parole chiave): scrivi una lista di parole chiave, quelle che tu pensi possano rendere il tuo ticket pi\`u facile da trovare. I ticket su cui si \`e lavorato a qualche Sage Day hanno spesso ``sdNN`` come parola chiave, dove ``NN`` \`e il numero del Sage Day.

* **Cc** (copia carbone): lista di utenti di Trac a cui mandare emails di segnalazione di cambiamenti sul ticket. Note che gli utenti che immettono un commento sono automaticamente sottoscritti agli aggiornamenti e non hanno bisogno di essere elencati sotto Cc.

* **Merged in** (unito a): la release di Sage a cui il ticket \`e stato unito. Pu\`o essere cambiata solo dal manager di release.

* **Authors** (autori): nome reale dell'autore del ticket, o lista degli autori.

* **Reviewers** (revisori): nome reale del revisore del ticket, o lista dei revisori.

* **Report upstream** (segnala a monte): se il ticket \`e un baco in un componente a monte di Sage (ad esempio in Maxima, Pari, ecc.) questo campo \`e utilizzato per riassumere la comunicazione con gli sviluppatori di tale componente.

* **Work issues** (esigenze di lavorazione): questioni che devono essere risolte prima che il ticket possa evolvere oltre lo status "needs work".

* **Branch** (ramo): il ramo di Git che contiene il codice del ticket (vedi :ref:`section-walkthrough-branch`). \`E mostrato in colore verde, a meno che vi sia un conflitto fra il ramo e l'ultima release beta (colore rosso). In questo caso, si dovrebbe fare un merge o un rebase del ramo.

* **Dependencies** (dipendenze): il ticket dipende da un altro ticket? A volte un ticket richede che un altro venga risolto prima. Se si \`e in questa situazione, scrivere ledipendenze come una lista separata da virgole (ad esempio ``#1234, #5678``) nel campo "Dependencies".

* **Stopgaps:** See :ref:`section-trac-stopgaps`.

.. _section-trac-ticket-status:

Lo status di un ticket
======================

Lo status di un ticket appare subito vicino al suo numero, nell'angolo superiore sinistro della sua pagina. Indica che deve lavorarci sopra.

- **new** (nuovo) -- il ticket \`e solo stato creato (o l'autore ha dimenticato di cambiarne lo status a qualcos'altro).

  Se vuoi lavorarci sopra tu \`e meglio lasciare un commento per dirlo. Pu\`o evitare di avere 2 persone che lavorano sulla stessa cosa.

- **needs_review** (richiede revisione) -- il codice \`e pronto per una revisione fra pari. Se il codice non \`e tuo, allora puoi farne la revisione. Vedi :ref:`chapter-review`.

- **needs_work** (richiedere lavorazione) -- qualcosa dev'essere cambiato nel codice. La ragione dovrebbe essere visibile nei commenti.

- **needs_info** (mancano informazioni) -- qualcuno deve rispondere ad una domanda prima che qualunque altra cosa possa essere fatta. Dovrebbe essere chiaro dai commenti.

- **positive_review** (revisione positiva) -- \`e stata fatta la revisione del ticket, ed il release manager lo chiuder\`a.

Lo status di un ticket pu\`o essere cambiato usando un form in fondo alla pagina del ticket. Lascia un commento  che spieghi le tue ragioni ogni volta che fai un cambiamento.

.. _section-trac-stopgaps:

Tappabuchi
==========

Se un componente di Sage produce un errore matematico, dovresti aprire 2 ticket: il ticket principale, con tutti i dettagli, ed un ticket "tappabuchi" (ad esempio :trac:`12699`). Questo secondo ticket dovrebbe avere una patch (soluzione provvisoria) che sar\`a unita a Sage se nessuno sistema il problema principale; questa patch dovr\`a stampare un avvertimento quando qualcuno utilizza la funzionalit\`a bacata (il codice specifico). Per produrre il messaggio di avvertimento usa codice come il seguente::

    from sage.misc.stopgap import stopgap
    stopgap("This code contains bugs and may be mathematically unreliable.",
        TICKET_NUM)

Sostituisci ``TICKET_NUM`` con il numero del ticket principale. Vedi (link trac ticket #1269) per un esempio. Sul ticket principale dovresti anche immettere il numero di ticket del tappabuchi nel campo Stopgap (vedi :ref:`section-trac-fields`). I ticket tappabuchi vanno imputati come bloccanti.

.. note::
    se codice corretto matematicamente causa una segnalazione di errore in Sage o un crash allora non c'\`e
    bisogno di un tappabuchi. Essi servono per avvertire gli utenti che il codice che stanno utilizzando pu\`o
    essere difettoso: se il difetto \`e evidente perch\`e vi \`e un crash o una segnalazione di errore, non c'\`e
    quest'esigenza.

Lavorare sui ticket
===================

Se riesci a correggere un baco o a migliorare Sage, tu sei il nostro eroe. Vedi :ref:`chapter-walkthrough` per effettuare cambiamenti al codice sorgente di Sage, comunicarli al sage trac server, ed infine segnalare sul relativo ticket di Trac il nuovo ramo che hai prodotto.
Le seguenti sono altre esigenze importanti:

* il costruttore di patch automatico effettuer\`a un test sul tuo ticket. Vedi `the patchbot wiki <http://wiki.sagemath.org/buildbot>`_ per pi\`u informazioni su queste funzionalit\`a e limitazioni. Non mancare di dare un'occhiata al log, specialmente se il costruttore di patch automatico non ti d\`a semaforo verde.

* Per ogni baco risolto dev'essere prodotto un doctest.

* Questa non \`e un'esigenza con i difetti, ma ci sono molti miglioramentipossibili per Sage e troppi pochi sviluppatori per implementare tutte le buone idee. Il trac server \`e utile per tenere le ideee in un posto centralizzato perch\`e nei Google groups tendono a perdersi quando non sono pi\`u in prima pagina.

* Se sei uno sviluppatore, sii buono e cerca ogni tanto di risolvere un ticket vecchio.

* Alcune persone rivedono regolarmente le priorit\`a. In questo contesto, ci\`o significa che diamo un'occhiata ai nuovi bachi e li classifichiamo secondo quella che consideriamo esserne la priorit\`a. \`E molto probabile che altre persone possano vedere le priorit\`a in modo molto differente da noi, quindi facci sapere se vedi dei problemi con ticket specifici.

Rivedere le patch
=================

Tutto il codice che finisce in Sage \`e contro-verificato fra colleghi, per assicurarsi che le convenzioni presentate in questo manuale siano seguite, che ci siano sufficienti esempi nella documentazione e doctest, e per cercare di essere sicuri che il codice faccia, matematicamente, cosa si suppone che faccia.
Se qualcuno (altri che tu) ha inviato sul Trac server un ramo git per un ticket, tu puoi farne la revisione! Controlla il "branch diff" (l'elenco delle modifiche), cliccandoci sopra, per vedere se ha senso. Scaricalo (vedi :ref:`chapter-review`) e compila Sage con il nuovo ramo incluso, quindi fatti delle domande come le seguenti:

* il nuovo codice sorgente ha senso?

* quando lo esegui in Sage, risolve il problema riportato nel ticket relativo?

* introduce qualche nuovo problema?

* \`e documentato a sufficienza, incluse sia le spiegazioni che i doctest? Tutto il codice in Sage deve avere dei doctest, quindi se l'autore del ticket cambia del codice che non aveva un doctest prima, la nuova versione deve includerne uno. In particolare tutto il nuovo codice deve essere provato con dei doctest, al 100%. Usa il comando ``sage -coverage <files>`` per vedere la percentuale di copertura di ``<files>``.

* in particolare, vi \`e un doctest che illustri che il baco \`e stato risolto? Se una funzione dava un risultato sbagliato e questo ticket la corregge, allora dovrebbe includere un doctest che illustri il suo successo. La doctring relativa dovrebbe includere il numero di ticket, ad esempio ``vedi :trac:'12345'``.

* se il ticket afferma di accellerare qualche calcolo, contiene anche degli esempi di codice per mostrare quanto afferma? Il ticket dovrebbe analizzare esplicitamente qual'\`e la velocit\`a prima di applicare la patch e qual'\`e dopo, e spiegare qual'\`e il guadagno di tempo.

* il manuale di riferimento compila senza errori? Puoi provare il manuale di riferimento utilizzando il comando ``sage -docbuild reference html`` per produrne la versione in HTML. Anche la versione PDF dev'essere prodotta senza errori: usa il comando ``sage -docbuild reference pdf`` per provarlo. Tale comando richiede che tu abbia installato Latex sul tuo PC.

* i doctest passano tutti senza errori? \`E difficile predire quali componenti di Sage verranno toccati da una data patch, e dovresti lanciare i test dell'intera libreria (inclusi quelli etichettati ``#long``) prima di segnalare esito positivo alla revisione. Poi effettuare il test della libreria Sage con ``make ptestlong``.

* il codice e la documentazione seguono le convenzioni documentate nelle sezioni seguenti?

Se la risposta a queste ed altre domande ragionevoli simili \`e s\`i, allora puoi dar esito positivo alla revisione. Sulla pagina principale del ticket scrivi un commento nello spazio a ci\`o riservato, spiegando la tua revisione. Se ritieni di non avere abbastanza esperienza per fare ci\`o, scrivi un commento che spieghi che cosa hai verificato, e concludi chiedendo se qualcuno con pi\`u esperienza pu\`o dare un'occhiata. Se pensi che ci siano problemi con la patch, spiegali nel riquadro dei commenti e cambia lo status a "need work" (richiede lavorazione). Guarda altri ticket su Trac per vedere come si fa.
Se tu stesso cambi la patch, devi fare un commit sotto il tuo nome e segnarlo come patch conseguente ad una revisione. Questa va anch'essa sottoposta a revisione, per esempio dall'autore della patch originale.

.. note::
    "il meglio \`e nemico del bene": lo scopo della revisione \`e assicurarsi che le lineee guida sul codice
    di Sage siano seguite e che l'implementazione sia matematicamente corretta. Per cortesia astieniti dalla
    richiesta di funzionalit\`a aggiuntive e discussioni su implementazioni alternative che non siano mirate.
    Se vuoi che la patch sia scritta diversamente, il tuo suggerimento dev'essere una richiesta chiara e fattibile.

Chiusura dei ticket
===================

Solo il manager di release di Sage chiuder\`a i ticket. Molto probabilmente non sei tu ed il tuo Trac account non ha i permessi necessari. Se hai forti ragioni per ritenere che un ticket debba essere chiuso o cancellato, allora cambia il suo status  a *needs review* (richiede revisione) e cambia il "milestone" a *sage-duplicate/invalid/wontfix*. Dovresti anche aggiungere un commento, spiegando  perch\`e dovrebbe essere chiuso. Se un altro sviluppatore \`e del tuo stesso parere, cambier\`a lo status a *positive review* (revisione positiva).

Un problema simile \`e la riapertura di un ticket. Dovresti astenerti dal riaprire un ticket che \`e gi\`a stato chiuso. Apri invece un nuovo ticket e metti nella descrizione un link al vecchio ticket.

Ragioni per invalidare dei ticket
=================================

**Un problema per ticket**: un ticket deve riguardare un solo problema e non dovrebbe essere una lista della spesa di problemi scollegati fra loro. Se un ticket rigurada pi\`u di un'esigenza, non lo possiamo chiudere e sebbene alcune patch fossero state applicate in una data release, rimarr\`a in un limbo.

**No alle patch-bomb**: il codice che viene incluso in Sage \`e soggetto alla revisione fra pari. Se arrivi con 80000 linee di codice che sostituiscono un intero sottosistema con qualcos'altro, puoi immaginarti che il processo di revisione sar\`a un po' noisoso. Queste enormi patch-bomb sono problematiche per molte ragioni e preferiamo cambiamenti piccoli e graduali che possono essere rivisti ed applicati facilmente. Questo non \`e sempre possibile (ad esempio in caso di riscrittura obbligata per qualche motivo), ma \`e comunque caldamente raccomandato che eviti questo stile di sviluppo a meno che non vi siano alternative.

**Specifico per Sage**: la filosofia di Sage \`e che mettiamotutto (o quasi) in un unico archivio TAR per rendere possibile il processo di debug. Puoi immaginarti l'esplosione combinatoria che ci ritroveremmo a dover gestire se tu rimpiazzassi anche solo 10 componenti con dei pacchetti esterni. Nel momento in cui inizi a rimpiazzare i componenti pi\`u essenziali di Sage con le versioni pacchettizzate che puoi comunemente trovare (ad esempio Pari, GAP, lisp, gmp), eventuali problemi non hanno pi\`u posto sul nostro tracker. Ad esempio se utilizzi un pacchetto PARI con dei bachi, segnala il baco a loro. Di solito vogliamo e possiamo risolvere il problema, ma non garantiamo che ti aiuteremo. Dando un'occhiata al numero di tickets aperti, specifici di Sage, si spera capirai perch\`e.

**No alle discussioni di supporto**: il sistema Trac non \`e stato fatto per rispondere a difficolt\`a nell'utilizzo di Sage: i ticket devono degli evidenti bachi e non cose del tipo "a provato a fare questo e non ci sono riuscito. Come si fa?". Di solito queste cose non sono in relazione con dei bachi e verosimilmente ``sage-support`` \`e in grado di rispondere alla questione. Se poi viene fuori che ti sei imbattuto in un baco, allora qualcun aprir\`a un ticket, coincisamente e circostanziatamente.

**Le soluzioni devono essere realizzabili**: i ticket devono essere realizzabili. Spesso i ticket che ricadono in questa categoria violano qualcuna delle regole sopra elencate. Un esempio \`e il sopraddetto "rendere Sage il miglior software del mondo". Non c'\`e un criterio di misura e pu\`o essere molto soggettivo.


