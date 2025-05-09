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

- Analizzare la sequenza del gene **SCN5A**
- Calcolare statistiche (es. lunghezza, %GC)
- Simulare varianti genetiche
- Produrre visualizzazioni e report riassuntivi

---

## 👥 Team

**Organizzazione GitHub:** [BioinformaticSiena](https://github.com/BioinformaticSiena)  
**Repository:** [cardiac-channelopathy](https://github.com/BioinformaticSiena/cardiac-channelopathy)

Membri del team:
- Martina (`@Modena3`)
- Marco (`@marcuzan8`)
- Francesco (`@Frampino`)
- Aldo (`@aldone1`)

---

## 📄 Licenza

Progetto a scopo **accademico e didattico**. Nessuna licenza commerciale applicata.

---
