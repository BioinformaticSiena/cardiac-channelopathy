#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Modulo per la gestione delle sequenze genetiche.
"""

from Bio import SeqIO
import os

def load_sequence(file_path, format="fasta"):
    """
    Carica una sequenza da un file.
    
    Args:
        file_path (str): Percorso al file della sequenza
        format (str): Formato del file (fasta, genbank, etc.)
        
    Returns:
        SeqRecord: Oggetto contenente la sequenza
    """
    try:
        record = SeqIO.read(file_path, format)
        print(f"Sequenza caricata: {record.id}, lunghezza: {len(record.seq)} bp")
        return record
    except Exception as e:
        print(f"Errore durante il caricamento della sequenza: {str(e)}")
        return None

def get_protein_sequence(sequence_record):
    """
    Estrae la sequenza proteica (o la converte se è DNA).
    
    Args:
        sequence_record (SeqRecord): Record della sequenza
        
    Returns:
        str: Sequenza proteica
    """
    seq = sequence_record.seq
    
    # Se la sequenza è DNA, traducila in proteina
    if all(c in 'ATGCNatgcn' for c in seq):
        # Trova il frame di lettura (semplificato)
        for i in range(3):
            prot = seq[i:].translate(to_stop=True)
            if len(prot) > 100:  # Assumiamo che una proteina valida sia lunga almeno 100aa
                return prot
        print("Nessuna proteina valida trovata nella traduzione")
        return None
    else:
        # Assume che sia già una sequenza proteica
        return seq
