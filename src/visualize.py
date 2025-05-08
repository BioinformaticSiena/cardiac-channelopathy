import matplotlib.pyplot as plt

def plot_variant_distribution(variants_df, protein_length, save_path):
    """Crea un grafico della distribuzione delle varianti sulla sequenza proteica."""
    try:
        # Assumiamo che ci sia una colonna "protein_position"
        positions = variants_df["position"].dropna()

        plt.figure(figsize=(12, 6))
        plt.hist(positions, bins=range(0, protein_length+10, 10), edgecolor='black')
        plt.title("Distribuzione delle varianti lungo la proteina")
        plt.xlabel("Residuo della proteina")
        plt.ylabel("Numero di varianti")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close()
        print(f"Salvato il grafico della distribuzione in {save_path}")
    except Exception as e:
        print(f"Errore nella creazione del grafico di distribuzione: {str(e)}")

def plot_variant_types(variants_df, save_path):
    """Crea un grafico a torta dei tipi di varianti."""
    try:
        type_counts = variants_df["variant_type"].value_counts()

        plt.figure(figsize=(8, 8))
        plt.pie(type_counts, labels=type_counts.index, autopct='%1.1f%%', startangle=140)
        plt.title("Distribuzione dei tipi di varianti")
        plt.axis('equal')  # Cerchio perfetto
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close()
        print(f"Salvato il grafico dei tipi di varianti in {save_path}")
    except Exception as e:
        print(f"Errore nella creazione del grafico dei tipi di varianti: {str(e)}")

def generate_summary_report(variants_df, gene, save_path):
    """Genera un report di riepilogo delle varianti in un file di testo."""
    try:
        with open(save_path, 'w') as report_file:
            report_file.write(f"Report analisi varianti per il gene: {gene}\n\n")
            report_file.write(f"Numero totale di varianti: {len(variants_df)}\n\n")
            
            if "variant_type" in variants_df.columns:
                report_file.write("Distribuzione dei tipi di varianti:\n")
                type_counts = variants_df["variant_type"].value_counts()
                for variant_type, count in type_counts.items():
                    report_file.write(f"- {variant_type}: {count} ({count/len(variants_df)*100:.1f}%)\n")
        print(f"Salvato il report in {save_path}")
    except Exception as e:
        print(f"Errore nella creazione del report: {str(e)}")
