# Analisi di Varianti Genetiche nelle Canalopatie Cardiache

Questo progetto nasce come esame del modulo di **Fondamenti di Programmazione in Python** del master: Bioinformatica e data science presso l'Università di Siena. Il nostro obiettivo è sviluppare un semplice programma in Python per analizzare la sequenza del gene **SCN5A**, coinvolto nelle canalopatie cardiache ereditarie.

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

## ▶️ Come eseguire il programma

1. Clona la repository:
   
   git clone https://github.com/BioinformaticSiena/cardiac-channelopathy.git
   Assicurati di avere Python installato. Installa Biopython:

2. Assicurati di avere Python installato. 

   Installa Biopython: pip install biopython 

3. Posiziona il file SCN5A.fasta nella cartella data/.

4. Esegui lo script: python src/analisi_scn5a.py

📄 Licenza
Uso accademico/didattico. Nessuna licenza commerciale applicata.
