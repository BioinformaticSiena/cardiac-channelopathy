from Bio import SeqIO

def load_sequence(file_path):
    """Carica una sequenza da file FASTA."""
    try:
        record = SeqIO.read(file_path, "fasta")
        return record
    except Exception as e:
        print(f"Errore nel caricamento della sequenza: {str(e)}")
        return None

def get_protein_sequence(dna_sequence_record):
    """Traduce una sequenza DNA in proteina."""
    try:
        protein_seq = dna_sequence_record.seq.translate(to_stop=True)
        return protein_seq
    except Exception as e:
        print(f"Errore nella traduzione della sequenza: {str(e)}")
        return None
