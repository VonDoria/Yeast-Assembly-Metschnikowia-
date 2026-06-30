import os
import sys
import pandas as pd
from Bio import AlignIO
# from Bio.Align import MultipleSeqAlignment
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq


OUTPUT_DIR = sys.argv[1]
SEQ_TYPE = sys.argv[2]

aligned_folder = os.path.join(OUTPUT_DIR, "genes_aligned")

print("\n🔗 Concatenando alinhamentos (Construindo a Supermatriz)...")

samples_df = pd.read_csv(f"{OUTPUT_DIR}/gene_presence_matriz.csv")
samples = samples_df[samples_df.columns[0]]

super_alignment = {sample: "" for sample in samples}
partitions = []
current_position = 1

for gene in os.listdir(aligned_folder):
    aligned_path = os.path.join(aligned_folder, gene)
    alignment = AlignIO.read(aligned_path, "fasta")
    
    gene_len = alignment.get_alignment_length()
    
    partitions.append(f"{SEQ_TYPE.upper()}, {gene.split("_")[0]} = {current_position}-{current_position + gene_len - 1}")
    current_position += gene_len

    for sample in samples:
        found = False
        for record in alignment:
            if record.id == sample:
                super_alignment[sample] += str(record.seq)
                found = True
                break
        if not found:
            super_alignment[sample] += "-" * gene_len

supermatrix_path = os.path.join(OUTPUT_DIR, "supermatrix.fasta")
with open(supermatrix_path, "w") as f:
    for sample, sequence in super_alignment.items():
        f.write(f">{sample}\n{sequence}\n")

partitions_path = os.path.join(OUTPUT_DIR, "partitions.nex")
with open(partitions_path, "w") as f:
    f.write("#nexus\nbegin sets;\n")
    for part in partitions:
        f.write(f"    charset {part.split(', ')[1]};\n")
    f.write("end;\n")

print(f"   -> Supermatriz gerada: {supermatrix_path} ({current_position - 1} posições totais).")
print(f"   -> Arquivo de partição salvo: {partitions_path}")