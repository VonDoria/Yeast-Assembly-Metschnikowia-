#!/bin/bash
# chmod +x tree_builder.sh

# Para a execução se algum comando falhar
set -e

SCRIPTDIR="scripts" 
OUTPUT_DIR="results"

REF_ASSEMBLYS_FOLDER="../NCBI_Metschnikowias_assemblys"
SAMPLE_SC_GENES_FOLDER="../Metschnikowia_target/BUSCO_results/trimmomatic_reads_isolate/run_saccharomycetes_odb12/busco_sequences/single_copy_busco_sequences"
LINHAGEM="saccharomycetes_odb12" 
SEQ_TYPE="faa" 
THREADS=16
MIN_PRESENCE_THRESHOLD=0.8    # Limite mínimo de presença (0.8 significa que o gene deve estar presente em pelo menos 80% das amostras)

if [ -d "${OUTPUT_DIR}" ]; then
    echo "Genes de cópia única já agrupados."
else
    echo "Buscando por genes de cópia única."
    python $SCRIPTDIR/mapping_and_group_genes.py $REF_ASSEMBLYS_FOLDER $SAMPLE_SC_GENES_FOLDER $OUTPUT_DIR $SEQ_TYPE $MIN_PRESENCE_THRESHOLD
fi

if [ -d "${OUTPUT_DIR}/genes_aligned" ]; then
    echo "Alinhamento já executado."
else
    mkdir -p "${OUTPUT_DIR}/genes_aligned"
    
    echo "[MAFFT] Iniciando alinhamento."
    for unaligned_gene_file in "$OUTPUT_DIR"/genes_unaligned/*.fasta; do
        
        if [ -f "$unaligned_gene_file" ]; then

            gene_id=$(basename "$unaligned_gene_file")
            echo "Alinhando: $gene_id"
            mafft --auto --thread $THREADS "$unaligned_gene_file" > "${OUTPUT_DIR}/genes_aligned/${gene_id%.*}_aligned.fasta"

        fi
    done
fi

if [ -f "${OUTPUT_DIR}/supermatrix.fasta" ]; then
    echo "Super matrix já criada."
else
    python $SCRIPTDIR/concatenate_alignments.py $OUTPUT_DIR $SEQ_TYPE
fi

if [ -d "${OUTPUT_DIR}/arvore_metschnikowia" ]; then
    echo "Árvore filogenética já criada."
else
    echo "[IQ-TREE] Iniciando a montagem da árvore filogenética..."

    iqtree -s "${OUTPUT_DIR}/supermatrix.fasta" \
    -spp "${OUTPUT_DIR}/partitions.nex" \
    -m LG+G -bb 1000 \
    -pre "${OUTPUT_DIR}/arvore_metschnikowia" -nt $THREADS
    
    echo "Árvore Filogenética gerada com sucesso!"
    echo "-> Arquivo final da árvore (formato Newick): ${OUTPUT_DIR}/arvore_metschnikowia.treefile"
fi



