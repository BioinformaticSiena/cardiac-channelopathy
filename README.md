# 🧬 Analisi di Varianti Genetiche nelle Canalopatie Cardiache

Questo progetto nasce come esame del modulo di **Fondamenti di Programmazione in Python** del Master in *Bioinformatica e Data Science* presso l'Università di Siena.

Il nostro obiettivo è sviluppare un semplice programma in Python per analizzare la sequenza del gene **SCN5A**, coinvolto nelle **canalopatie cardiache ereditarie**. Il programma utilizza come input:
- un file `.fasta` contenente la sequenza CDS del gene,
- un file `.csv` contenente varianti genetiche simulate,

e restituisce informazioni genetiche, grafici e report riassuntivi.

---

### 💻 Installazione

Per installare le dipendenze necessarie, assicurati di avere Python 3 installato, poi esegui:

```bash
pip install -r requirements.txt

## 🔄 Generazione di Varianti Simulate

Poiché non sono disponibili dati clinici reali sufficienti, questo progetto utilizza varianti genetiche simulate per testare l'analisi.
Lo script `generate_scn5a_variants.py` include le seguenti funzionalità:

1. **Lettura della sequenza proteica**  
   - Legge il file FASTA del gene SCN5A e lo traduce in proteina.

2. **Generazione casuale di varianti**  
   - Posizione casuale nella proteina  
   - Tipo di variante: `missense`, `nonsense`, `frameshift`, `inframe_deletion`  
   - Classificazione di patogenicità con probabilità configurabili

3. **Parametrizzazione da linea di comando**  
   - `--num_variants`: numero di varianti da generare (default: 20)  
   - `--output`: nome del file di output (default: variants.csv)  
   - `--random_seed`: seme per la riproducibilità (default: 42)


```bash
python src/generate_scn5a_variants.py \
    --fasta data/SCN5A.fasta \
    --output data/generated_variants.csv \
    --num_variants 50
```

---

## ▶️ Analisi delle Varianti

Una volta generate le varianti, si può eseguire lo script principale:

```bash
python src/main.py \
    --gene SCN5A \
    --sequence_file data/SCN5A.fasta \
    --variants_file data/generated_variants.csv
```

Lo script:
- Carica e traduce la sequenza
- Filtra le varianti per gene
- Calcola statistiche
- Genera grafici (distribuzione e tipi di varianti)
- Salva un report riassuntivo

---

📊 Risultati e Visualizzazioni
L'analisi produce diversi output che aiutano a comprendere la distribuzione e la natura delle varianti genetiche in SCN5A:
1. Distribuzione delle varianti lungo la proteina: Questo istogramma mostra dove si posizionano le varianti lungo la sequenza aminoacidica di SCN5A
2. Tipi di varianti: Il grafico a torta mostra la frequenza relativa dei diversi tipi di varianti
3. Report di analisi: Il report testuale fornisce statistiche dettagliate sulle varianti analizzate, inclusa la loro classificazione di patogenicità

🔬 Interpretazione Biologica
Il ruolo di SCN5A nelle canalopatie cardiache
Il gene SCN5A codifica per il canale del sodio voltaggio-dipendente Nav1.5, fondamentale per la generazione e propagazione del potenziale d'azione nei cardiomiociti. Mutazioni in questo gene sono associate a diverse patologie cardiache ereditarie, tra cui:

Sindrome di Brugada: Caratterizzata da anomalie dell'ECG e rischio di aritmie ventricolari potenzialmente fatali
Sindrome del QT lungo di tipo 3: Prolungamento dell'intervallo QT con rischio di torsioni di punta e morte cardiaca improvvisa
Malattia della conduzione cardiaca progressiva: Disturbi del sistema di conduzione che possono portare a blocchi atrioventricolari

Significato clinico dei diversi tipi di varianti

Varianti missense: Sono le più comuni e possono avere effetti molto variabili sulla funzione del canale, da lievi alterazioni della cinetica a completa perdita di funzione
Varianti nonsense: Generalmente causano la produzione di una proteina tronca non funzionale, spesso associate a fenotipi severi
Varianti frameshift: Alterano il frame di lettura, portando a proteine aberranti e generalmente a fenotipi patologici
Delezioni in-frame: L'impatto dipende dalla regione coinvolta; se interessano domini funzionali critici possono compromettere significativamente la funzionalità del canale

Domini funzionali critici
Le varianti localizzate in specifici domini funzionali di Nav1.5 spesso hanno conseguenze cliniche più severe:

Segmenti S4: Coinvolti nel sensore di voltaggio
Loop P: Determinanti per la selettività ionica
Regioni di "gating": Critiche per l'apertura e chiusura del canale
Domini di interazione: Importanti per l'assemblaggio con subunità accessorie

Correlazione genotipo-fenotipo
La correlazione tra specifiche varianti e manifestazioni cliniche è complessa e influenzata da fattori genetici ed ambientali. La stessa variante può manifestarsi diversamente in individui differenti (penetranza variabile ed espressività variabile), evidenziando l'importanza di un'analisi genetica approfondita nei pazienti con sospetta canalopatia cardiaca.

## 📁 Struttura del Progetto

```
cardiac-channelopathy/
├── data/                   # Dati di input (FASTA, CSV)
│   └── SCN5A.fasta
├── results/                # Output e grafici
├── src/                    # Codice Python
│   ├── main.py             # Script principale
│   ├── generate_scn5a_variants.py
│   ├── sequence_utils.py   # Funzioni di supporto (caricamento/traduzione sequenze)
│   └── visualize.py        # Funzioni per i grafici e report
├── README.md               # Questo file
└── requirements.txt        # (opzionale) Librerie richieste
```

---

## 📊 Tecnologie Utilizzate

- [Python 3](https://www.python.org)
- [Biopython](https://biopython.org)
- [NumPy](https://numpy.org)
- [matplotlib](https://matplotlib.org)
- Git, GitHub, GitHub Desktop
- Visual Studio Code / PyCharm
- Anaconda (facoltativo)

---

## 🔬 Obiettivi Didattici

- Elaborare sequenze proteiche da file FASTA
- Generare e gestire varianti genetiche simulate
- Analizzare la distribuzione delle varianti lungo la proteina
- Classificare e conteggiare i tipi di varianti
- Visualizzare i risultati con grafici informativi
- Generare report testuali riassuntivi
- Sviluppare codice modulare con gestione degli errori

---

🚀 Sviluppi Futuri
Questo progetto potrebbe essere esteso in diverse direzioni:

- Integrazione con dati clinici reali: Adattare il programma per utilizzare sia varianti provenienti da database pubblici (ClinVar, HGMD) che file di pazienti provenienti da coorti cliniche personali
- Analisi di conservazione evolutiva: Confrontare le posizioni delle varianti con la conservazione tra specie diverse
- Mappatura su struttura 3D: Visualizzare le varianti sulla struttura tridimensionale del canale Nav1.5
- Predizione di patogenicità: Implementare algoritmi di machine learning per predire l'impatto funzionale delle varianti
- Analisi di domini funzionali: Mappare le varianti sui domini funzionali conosciuti del canale Nav1.5

## 👥 Team

**Organizzazione GitHub:** [BioinformaticSiena](https://github.com/BioinformaticSiena)  
**Repository:** [cardiac-channelopathy](https://github.com/BioinformaticSiena/cardiac-channelopathy)

Membri del team:
- Martina (`@Modena3`)
- Marco (`@marcuzan8`)
- Francesco (`@Frampino`)
- Aldo (`@aldone1`)

---

📚 Riferimenti Bibliografici

Wilde, A. A., & Amin, A. S. (2018). Clinical Spectrum of SCN5A Mutations: Long QT Syndrome, Brugada Syndrome, and Cardiomyopathy. JACC: Clinical Electrophysiology, 4(5), 569-579.
Kroncke, B. M., Glazer, A. M., Smith, D. K., Blume, J. D., & Roden, D. M. (2018). SCN5A (NaV1. 5) Variant Functional Perturbation and Clinical Presentation: Variants of a Certain Significance. Circulation: Genomic and Precision Medicine, 11(5), e002095.
Veerman, C. C., Wilde, A. A., & Lodder, E. M. (2015). The cardiac sodium channel gene SCN5A and its gene product NaV1.5: Role in physiology and pathophysiology. Gene, 573(2), 177-187.

## 📄 Licenza

Progetto a scopo **accademico e didattico**. Nessuna licenza commerciale applicata.

---
