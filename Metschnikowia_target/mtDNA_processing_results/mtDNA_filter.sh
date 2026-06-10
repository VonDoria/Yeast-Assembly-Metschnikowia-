#!/bin/bash
# chmod +x mtDNA_filter.sh

# Para a execução se algum comando falhar
set -e


BASEDIR="/home/italofaria/y7002/Metschnikowia_target"
SCRIPTDIR="${BASEDIR}/mtDNA_processing_results/scripts/python_scripts"
ASSEMBLY_ORIGINAL_PATH="${BASEDIR}/SPAdes_results/trimmomatic_reads_isolate/clean_contigs.fasta"
ASSEMBLY_NAME="$(basename "$(dirname "${ASSEMBLY_ORIGINAL_PATH}")")"
ASSEMBLY_PATH="${BASEDIR}/mtDNA_processing_results/contigs/${ASSEMBLY_NAME}_original_contigs.fasta"

echo "Copiando montagem para o diretório local."

mkdir -p $(dirname "$ASSEMBLY_PATH")
cp "$ASSEMBLY_ORIGINAL_PATH" "${ASSEMBLY_PATH}"

echo "Corrigindo formatação de caracteres .fasta"
python "$SCRIPTDIR/fixContigs.1.py" "$ASSEMBLY_PATH" 

echo "Iniciando busca Mt Feature Query. (Blast 1)"

mito_features="${BASEDIR}/mtDNA_processing_results/mito_anchors_V4/select_mito_core_cds.fasta"
evalue=0.001
word_size=16
dust="yes"

blastn \
-query $mito_features \
-subject $ASSEMBLY_PATH  \
-outfmt '6 qseqid sseqid qcovs pident evalue score qstart qend sstart send' \
-out $ASSEMBLY_PATH.mt_feature_query.tab \
-evalue $evalue \
-word_size $word_size \
-dust $dust

cov_cutoff=70
echo "Filtrando resultados por cobertura >${cov_cutoff}%."

python $SCRIPTDIR/find_mito_contigs.5.py $ASSEMBLY_PATH.mt_feature_query.tab $cov_cutoff

echo "Sumarizando resultados do blastn"
python $SCRIPTDIR/summarizeMtBlast.1.2.py $ASSEMBLY_PATH $ASSEMBLY_PATH.mt_feature_query.cov_filtered.tab
echo "Sumarização completa"

echo "Busca Mt Feature Query completa. (Blast 1)"



