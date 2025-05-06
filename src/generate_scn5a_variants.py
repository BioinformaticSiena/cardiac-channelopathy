#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script per generare un file CSV di varianti genetiche di esempio per SCN5A.
Questo script crea dati simulati per testare il programma di analisi delle varianti.
"""

import random
import pandas as pd
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def parse_arguments():
    """Gestisce gli argomenti da riga di comando."""
    parser = argparse.ArgumentParser(description="Genera un file CSV di varianti genetiche di esempio per SCN5A")
    parser.add_argument("--fasta", type=str, required=True, help="File FASTA della sequenza SCN5A")
    parser.add_argument("--output", type=str, default="variants.csv", help="Nome del file CSV di output")
    parser.add_argument("--num_variants", type=int, default=20, help="Numero di varianti da generare")
    parser.add_argument("--random_seed", type=int, default=42, help="Seed per la generazione random")
    return parser.parse_args()

def load_sequence(fasta_file):
    """Carica la sequenza dal file FASTA e la traduce in proteina."""
    try:
        record = SeqIO.read(fasta_file, "fasta")
        # Assumiamo che la sequenza sia già CDS (coding sequence)
        protein_seq = record.seq.translate()
        return protein_seq
    except Exception as e:
        print(f"Errore durante il caricamento della sequenza: {str(e)}")
        return None

def generate_variants(protein_seq, num_variants, random_seed=42):
    """Genera varianti casuali basate sulla sequenza proteica."""
    random.seed(random_seed)
    
    # Definisci i possibili tipi di varianti e la loro probabilità
    variant_types = {
        "missense": 0.7,  # 70% di probabilità
        "nonsense": 0.1,  # 10% di probabilità
        "frameshift": 0.1,  # 10% di probabilità
        "inframe_deletion": 0.1  # 10% di probabilità
    }
    
    # Definisci le possibili categorie di patogenicità e la loro probabilità
    pathogenicity_types = {
        "pathogenic": 0.3,  # 30% di probabilità
        "likely_pathogenic": 0.2,  # 20% di probabilità
        "uncertain": 0.2,  # 20% di probabilità
        "likely_benign": 0.15,  # 15% di probabilità
        "benign": 0.15  # 15% di probabilità
    }
    
    # Amino acidi standard
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    
    variants = []
    protein_length = len(protein_seq)
    
    for _ in range(num_variants):
        # Scegli una posizione casuale nella proteina
        position = random.randint(1, protein_length)
        
        # Ottieni l'aminoacido originale in quella posizione
        original_aa = protein_seq[position-1]
        
        # Scegli un tipo di variante in base alla probabilità
        variant_type = random.choices(
            list(variant_types.keys()),
            weights=list(variant_types.values()),
            k=1
        )[0]
        
        # Scegli una patogenicità in base alla probabilità
        pathogenicity = random.choices(
            list(pathogenicity_types.keys()),
            weights=list(pathogenicity_types.values()),
            k=1
        )[0]
        
        # Determina il nuovo aminoacido in base al tipo di variante
        if variant_type == "missense":
            # Scegli un aminoacido diverso dall'originale
            variant_aa = random.choice(amino_acids.replace(original_aa, ""))
        elif variant_type == "nonsense":
            variant_aa = "X"  # Stop codon
        elif variant_type == "frameshift":
            variant_aa = "fs"  # Frameshift
        elif variant_type == "inframe_deletion":
            # Simuliamo una delezione di 1-3 aminoacidi
            if position + 2 <= protein_length:
                deletion_len = random.randint(1, 3)
                variant_aa = "del"
                original_aa = str(protein_seq[position-1:position-1+deletion_len])
            else:
                variant_aa = "del"
        
        # Aggiungi la variante alla lista
        variants.append({
            "gene": "SCN5A",
            "position": position,
            "variant_type": variant_type,
            "pathogenicity": pathogenicity,
            "original_aa": original_aa,
            "variant_aa": variant_aa
        })
    
    return variants

def main():
    """Funzione principale."""
    args = parse_arguments()
    
    # Carica la sequenza proteica
    protein_seq = load_sequence(args.fasta)
    if protein_seq is None:
        print("Impossibile generare varianti senza una sequenza valida.")
        return
    
    print(f"Sequenza proteica caricata, lunghezza: {len(protein_seq)} aminoacidi")
    
    # Genera le varianti
    variants = generate_variants(protein_seq, args.num_variants, args.random_seed)
    
    # Crea un DataFrame pandas
    variants_df = pd.DataFrame(variants)
    
    # Salva il DataFrame in un file CSV
    variants_df.to_csv(args.output, index=False)
    print(f"Generate {len(variants)} varianti e salvate in {args.output}")

if __name__ == "__main__":
    main()
