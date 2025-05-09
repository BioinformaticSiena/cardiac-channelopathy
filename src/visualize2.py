#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Modulo per la visualizzazione delle varianti genetiche.
"""

import matplotlib.pyplot as plt
import pandas as pd
import os

def plot_variant_distribution(variants_df, protein_length, output_path):
    """
    Visualizza la distribuzione delle varianti lungo la proteina.
    
    Args:
        variants_df (pd.DataFrame): DataFrame delle varianti
        protein_length (int): Lunghezza della proteina
        output_path (str): Percorso per salvare il grafico
    """
    plt.figure(figsize=(10, 6))
    
    # Crea un istogramma delle posizioni delle varianti
    plt.hist(variants_df["position"], bins=20, alpha=0.7)
    
    # Aggiungi titolo e etichette
    plt.title("Distribuzione delle varianti genetiche")
    plt.xlabel("Posizione aminoacidica")
    plt.ylabel("Numero di varianti")
    plt.grid(True, linestyle="--", alpha=0.7)
    
    # Salva la figura
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    print(f"Grafico della distribuzione delle varianti salvato in: {output_path}")

def plot_variant_types(variants_df, output_path):
    """
    Visualizza la distribuzione dei tipi di varianti.
    
    Args:
        variants_df (pd.DataFrame): DataFrame delle varianti
        output_path (str): Percorso per salvare il grafico
    """
    plt.figure(figsize=(10, 6))
    
    # Conta i tipi di varianti
    type_counts = variants_df["variant_type"].value_counts()
    
    # Crea un grafico a torta
    plt.pie(type_counts, labels=type_counts.index, autopct="%1.1f%%", 
            startangle=90, shadow=True)
    
    # Aggiungi titolo
    plt.title("Distribuzione dei tipi di varianti")
    
    # Salva la figura
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    print(f"Grafico dei tipi di varianti salvato in: {output_path}")

def generate_summary_report(variants_df, gene_id, output_path):
    """
    Genera un report di riepilogo in formato testo.
    
    Args:
        variants_df (pd.DataFrame): DataFrame delle varianti
        gene_id (str): Identificatore del gene
        output_path (str): Percorso per salvare il report
    """
    with open(output_path, "w") as f:
        f.write(f"REPORT DI ANALISI DELLE VARIANTI PER IL GENE {gene_id}\n")
        f.write("="*50 + "\n\n")
        
        f.write(f"Numero totale di varianti analizzate: {len(variants_df)}\n\n")
        
        f.write("Distribuzione dei tipi di varianti:\n")
        type_counts = variants_df["variant_type"].value_counts()
        for variant_type, count in type_counts.items():
            f.write(f"  - {variant_type}: {count} ({count/len(variants_df)*100:.1f}%)\n")
        
        f.write("\nRiassunto della patogenicità:\n")
        if "pathogenicity" in variants_df.columns:
            patho_counts = variants_df["pathogenicity"].value_counts()
            for patho, count in patho_counts.items():
                f.write(f"  - {patho}: {count} ({count/len(variants_df)*100:.1f}%)\n")
        
        f.write("\nLe 5 posizioni con il maggior numero di varianti:\n")
        pos_counts = variants_df["position"].value_counts().head(5)
        for pos, count in pos_counts.items():
            f.write(f"  - Posizione {pos}: {count} varianti\n")
    
    print(f"Report di riepilogo salvato in: {output_path}")
