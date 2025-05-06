
from Bio import SeqIO
from Bio.Align import fasta
from Bio.SeqUtils import gc_fraction
# Caricamento della sequenza dal file FASTA
record = SeqIO.read('/data/AldoSCN5A.fasta', "fasta")
seq = record.seq
print("Lunghezza della sequenza:", len(seq))
