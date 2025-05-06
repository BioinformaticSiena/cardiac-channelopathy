#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analisi semplificata di varianti genetiche nelle canalopatie cardiache.
"""

import os
import pandas as pd
import argparse
from sequence_utils import load_sequence, get_protein_sequence
from visualize import plot_variant_distribution, plot_variant_types, generate_summary_report

# Definizione dei percorsi ai file
DATA_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")

def parse_arguments():
    """Analizza gli argomenti della riga di comando."""
    parser = argparse.ArgumentParser(description="Analisi semplificata di varianti genetiche nelle canalopatie cardiache")
    parser.add_argument("--gene", type=str, required=True, help="Gene da analizzare")
    parser.add_argument("--sequence_file", type=str, required=True, help="File FASTA della sequenza")
    parser.add_argument("--variants_file", type=str, required=True, help="File CSV delle varianti")
    return parser.parse_args()

def load_variants(file_path):
    """Carica le varianti da un file CSV."""
    try:
        variants_df = pd.read_csv(file_path)
        print(f"Caricate {len(variants_df)} varianti dal file {file_path}")
        return variants_df
    except Exception as e:
        print(f"Errore durante il caricamento delle varianti: {str(e)}")
        return None

def filter_variants_by_gene(variants_df, gene_id):
    """Filtra le varianti per un gene specifico."""
    if variants_df is None:
        return None
    
    gene_variants = variants_df[variants_df["gene"] == gene_id].copy()
    print(f"Filtrate {len(gene_variants)} varianti per il gene {gene_id}")
    
    return gene_variants

def main():
    """Funzione principale."""
    args = parse_arguments()
    
    # Verifica e crea le directory necessarie
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    gene_results_dir = os.path.join(RESULTS_DIR, args.gene)
    os.makedirs(gene_results_dir, exist_ok=True)
    
    print(f"Avvio dell'analisi delle varianti per il gene {args.gene}...")
    
    # Carica la sequenza
    sequence_file = os.path.join(DATA_DIR, args.sequence_file)
    sequence_record = load_sequence(sequence_file)
    if sequence_record is None:
        print("Impossibile procedere senza una sequenza valida.")
        return
    
    # Ottieni la sequenza proteica
    protein_seq = get_protein_sequence(sequence_record)
    if protein_seq is None:
        print("Impossibile procedere senza una sequenza proteica valida.")
        return
    
    # Carica e filtra le varianti
    variants_file = os.path.join(DATA_DIR, args.variants_file)
    variants_df = load_variants(variants_file)
    gene_variants = filter_variants_by_gene(variants_df, args.gene)
    
    if gene_variants is None or len(gene_variants) == 0:
        print(f"Nessuna variante trovata per il gene {args.gene}.")
        return
    
    # Calcola alcune statistiche di base
    print(f"\nStatistiche di base per le varianti del gene {args.gene}:")
    print(f"- Numero totale di varianti: {len(gene_variants)}")
    if "variant_type" in gene_variants.columns:
        type_counts = gene_variants["variant_type"].value_counts()
        print("- Distribuzione dei tipi di varianti:")
        for variant_type, count in type_counts.items():
            print(f"  * {variant_type}: {count} ({count/len(gene_variants)*100:.1f}%)")
    
    # Visualizza i risultati
    print("\nGenerazione delle visualizzazioni...")
    
    # Grafico della distribuzione delle varianti
    dist_plot_path = os.path.join(gene_results_dir, f"{args.gene}_variant_distribution.png")
    plot_variant_distribution(gene_variants, len(protein_seq), dist_plot_path)
    
    # Grafico dei tipi di varianti
    if "variant_type" in gene_variants.columns:
        types_plot_path = os.path.join(gene_results_dir, f"{args.gene}_variant_types.png")
        plot_variant_types(gene_variants, types_plot_path)
    
    # Genera un report di riepilogo
    report_path = os.path.join(gene_results_dir, f"{args.gene}_analysis_report.txt")
    generate_summary_report(gene_variants, args.gene, report_path)
    
    print(f"\nAnalisi completata per il gene {args.gene}!")
    print(f"I risultati sono stati salvati in: {gene_results_dir}")

if __name__ == "__main__":
    main()
