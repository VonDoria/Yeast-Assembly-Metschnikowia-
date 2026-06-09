#!/bin/bash
# chmod +x mtDNA_filter.sh

# Para a execução se algum comando falhar
set -e

ASSEMBLY_ORIGINAL_PATH="/home/italofaria/y7002/Metschnikowia_target/Spades_results/trimmomatic_reads_isolate/clean_contigs.fasta"
ASSEMBLY="$(basename "$(dirname "${ASSEMBLY_ORIGINAL_PATH}")")"

get_and_fix_contigs() {        
    echo "Copiando a montagem para o diretório local."
    cp "$ASSEMBLY_ORIGINAL_PATH" "./contigs/"${ASSEMBLY}"_original_contigs.fasta"

}

echo "Iniciando programa: "${ASSEMBLY}""

exit 0
get_and_fix_contigs

if [ ! -f $SPECIES.fasta ] ;  then
    cp $ASSEMBLY_ORIGINAL_PATH $SPECIES.fasta
    python $scripts_dir/python_scripts/fixContigs.1.py $SPECIES.fasta    
fi


echo "Iniciando busca Mt Feature Query."

cov_cutoff=70 #percent out of 100
mito_features="/home/italofaria/y7002/Metschnikowia_target/mtDNA_results/mito_anchors_V4/select_mito_core_cds.fasta"
evalue=0.001
word_size=16
dust="yes"

$blastn \
-query $mito_features \
-subject $SPECIES  \
-outfmt '6 qseqid sseqid qcovs pident evalue score qstart qend sstart send' \
-out $SPECIES.mt_feature_query.tab \
-evalue $evalue \
-word_size $word_size \
-dust $dust


python $scripts_dir/python_scripts/find_mito_contigs.5.py $SPECIES.mt_feature_query.tab $cov_cutoff
echo "Mt Feature Query Search Complete"

#Summarizing blast results
echo "Summarizing Blast results for mt feature query search"
python $scripts_dir/python_scripts/summarizeMtBlast.1.2.py $assembly $SPECIES.mt_feature_query.cov_filtered.tab
echo "summarizing complete"

