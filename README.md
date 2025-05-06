# Analisi di Varianti Genetiche nelle Canalopatie Cardiache

Questo progetto nasce come esame del modulo di **Fondamenti di Programmazione in Python** del master: Bioinformatica e data science presso l'Università di Siena. Il nostro obiettivo è sviluppare un semplice programma in Python per analizzare la sequenza del gene **SCN5A**, coinvolto nelle canalopatie cardiache ereditarie.
Il programma è strutturato per utilizzare come input due file, uno in formato .fasta che contiene la sequenza CDS del gene e uno in formato .csv che contiene sequenze di varianti, e per poi confrontarli in modo da ottenere una serie di informazioni genetiche e generare grafici e report.
Lo script implementa 

GENERAZIONE RANDOM DI VARIANTI PER FILE .csv
Lo script generate_scn5a_variants.py include le seguenti funzionalità:
1.	Lettura della sequenza proteica: Legge direttamente il file FASTA del gene SCN5A e lo traduce in proteina.
2.	Generazione casuale di varianti: Crea varianti genetiche in base a: 
o	Posizione casuale nella proteina
o	Tipo di variante (missense, nonsense, frameshift, inserzioni/delezioni)
o	Classificazione di patogenicità con probabilità configurabili
3.	Parametrizzazione attraverso linea di comando: 
o	Numero di varianti da generare
o	Seed per la generazione casuale (per riproducibilità)
o	File di output
Comando per eseguire script
python3.8 src/generate_scn5a_variants.py --fasta data/SCN5A.fasta --output data/generated_variants.csv --num_variants 50
•	--num_variants: il numero di varianti da generare (default: 20) 
•	--output: il nome del file di output (default: variants.csv) 
•	--random_seed: un seme per la generazione casuale, utile se vuoi creare lo stesso set di varianti più volte (default: 42)
Quindi dopo aver create il files con le varianti random si esegue
python3.8 src/main.py --gene SCN5A --sequence_file data/SCN5A.fasta --variants_file data/generated_variants.csv


## 👥 Team

Organizzazione GitHub: [BioinformaticSiena](https://github.com/BioinformaticSiena)  
Repository: [cardiac-channelopathy](https://github.com/BioinformaticSiena/cardiac-channelopathy)

Membri del team:
- Martina (GitHub: `@username1`)
- Marco (GitHub: `@username2`)
- Francesco (GitHub: `@username3`)
- Aldo (GitHub: `@username4`)

## 🧬 Obiettivo

Analizzare la sequenza del gene **SCN5A** per:
- Calcolare la lunghezza e la % di GC
- Ricercare varianti genetiche note (es. sostituzioni puntiformi)
- Produrre un piccolo report testuale

## 🧰 Tecnologie usate

- **Python 3**
- [Biopython](https://biopython.org)
- [NumPy](https://numpy.org)
- Git, GitHub, GitHub Desktop
- PyCharm
- Anaconda

## 📁 Struttura del progetto

cardiac-channelopathy/
├── src/ # Codice Python
│ └── analisi_scn5a.py # Script principale
├── dati/ # File FASTA/GenBank
│ └── SCN5A.fasta
├── risultati/ # Output generati
├── README.md # Questo file
└── requirements.txt # (opzionale) Librerie richieste

cardio_variant_analysis/
├── data/                  # Dati di input
├── results/               # Risultati
└── src/
    ├── main.py            # Script principale
    ├── sequence_utils.py  # Funzioni per gestire sequenze
    └── visualize.py       # Funzioni di visualizzazione

## ▶️ Come eseguire il programma
sequence_utils.py contiene due funzioni principali:
load_sequence: Carica una sequenza da un file FASTA (o altro formato specificato) utilizzando BioPython
get_protein_sequence: Analizza la sequenza caricata e determina se è una sequenza nucleotidica o proteica. Se è nucleotidica, cerca di tradurla in proteina trovando il frame di lettura corretto

📄 Licenza
Uso accademico/didattico. Nessuna licenza commerciale applicata.
