#!/bin/bash
# chmod +x mtDNA_filter.sh

# Para a execução se algum comando falhar
set -e

BASEDIR="/home/italofaria/y7002/Metschnikowia_target"
SCRIPTDIR="${BASEDIR}/mtDNA_processing_results/scripts/python_scripts"
CONTIGSDIR="${BASEDIR}/mtDNA_processing_results/contigs"
RESULTSDIR="${BASEDIR}/mtDNA_processing_results/results"
ASSEMBLY_ORIGINAL_PATH="${BASEDIR}/SPAdes_results/trimmomatic_reads_isolate/clean_contigs.fasta"
ASSEMBLY_NAME="$(basename "$(dirname "${ASSEMBLY_ORIGINAL_PATH}")")"
ASSEMBLY_PATH="${CONTIGSDIR}/${ASSEMBLY_NAME}.fasta"
ASSEMBLY_RESULTS_PATH="${RESULTSDIR}"/"${ASSEMBLY_NAME}"


if [ -f "${ASSEMBLY_PATH}" ]; then
    echo "Arquivo de montagem encontrado."
else
    echo "Copiando montagem para o diretório local."

    mkdir -p $(dirname "$ASSEMBLY_PATH")
    cp "$ASSEMBLY_ORIGINAL_PATH" "${ASSEMBLY_PATH}"

    echo "Corrigindo formatação de caracteres .fasta"
    python "$SCRIPTDIR/fixContigs.1.py" "$ASSEMBLY_PATH" 
fi

echo "Iniciando busca Mt Feature Query. (Blast 1)"

if [ -f "${ASSEMBLY_PATH}.blast_1.mt_feature_query.tab" ]; then
    echo "blastn já executado para ${ASSEMBLY_PATH}."
else

    mito_features="${BASEDIR}/mtDNA_processing_results/mito_anchors_V4/select_mito_core_cds.fasta"
    evalue=0.001
    word_size=16
    dust="yes"

    mkdir -p $RESULTSDIR

    blastn \
    -query $mito_features \
    -subject $ASSEMBLY_PATH  \
    -outfmt '6 qseqid sseqid qcovs pident evalue score qstart qend sstart send' \
    -out "$ASSEMBLY_RESULTS_PATH".blast_1.mt_feature_query.tab \
    -evalue $evalue \
    -word_size $word_size \
    -dust $dust
fi

if [ -f "${ASSEMBLY_RESULTS_PATH}.blast_1.mt_feature_query.cov_filtered.tab" ]; then
    echo "Resultado da filtragem encontrado."
else
    cov_cutoff=70
    echo "Filtrando resultados por cobertura >${cov_cutoff}%."

    python $SCRIPTDIR/find_mito_contigs.5.py "$ASSEMBLY_RESULTS_PATH".blast_1.mt_feature_query.tab $cov_cutoff
fi

if [ -f "${ASSEMBLY_RESULTS_PATH}.blast_1.mt_feature_query.cov_filtered.results_summary.txt" ]; then
    echo "Resultado da filtragem encontrado."
else
    echo "Sumarizando resultados do blastn"
    python $SCRIPTDIR/summarizeMtBlast.1.2.py $ASSEMBLY_PATH $ASSEMBLY_RESULTS_PATH.blast_1.mt_feature_query.cov_filtered.tab
    echo "Sumarização completa"
fi

echo "Busca Mt Feature Query completa. (Blast 1)"

echo "-------------------------------------------------------"

echo "Iniciando busca Contig Query. (Blast 2)"

if [ -f "${ASSEMBLY_PATH}.blast_2.contig_query.tab" ]; then
    echo "blastn já executado para ${ASSEMBLY_PATH}."
else

    mito_sources="${BASEDIR}/mtDNA_processing_results/mito_anchors_V4/select_sources.fasta"
    evalue=0.001
    word_size=16
    dust="yes"

    blastn \
    -query $ASSEMBLY_PATH \
    -subject $mito_sources  \
    -outfmt '6 qseqid sseqid qcovs pident evalue score qstart qend sstart send' \
    -out "$ASSEMBLY_RESULTS_PATH".blast_2.contig_query.tab \
    -evalue $evalue \
    -word_size $word_size \
    -dust $dust
fi

if [ -f "${ASSEMBLY_RESULTS_PATH}.blast_2.contig_query.cov_filtered.tab" ]; then
    echo "Resultado da filtragem encontrado."
else
    cov_cutoff=25
    echo "Filtrando resultados por cobertura >${cov_cutoff}%."

    python $SCRIPTDIR/find_mito_contigs.5.py "$ASSEMBLY_RESULTS_PATH".blast_2.contig_query.tab $cov_cutoff
fi

if [ -f "${ASSEMBLY_RESULTS_PATH}.blast_2.contig_query.cov_filtered.results_summary.txt" ]; then
    echo "Resultado da filtragem encontrado."
else
    echo "Sumarizando resultados do blastn"
    python $SCRIPTDIR/summarizeMtBlast.2.2.py $ASSEMBLY_PATH $ASSEMBLY_RESULTS_PATH.blast_2.contig_query.cov_filtered.tab
    echo "Sumarização completa"
fi

echo "Busca Contig Query completa. (Blast 2)"

echo "-------------------------------------------------------"

echo "Extraindo contigs mitocondriais."

if [ -f "${CONTIGSDIR}/found_mito_contig.fasta" ]; then
    echo "Contigs mitocondriais já extraidos."
else
    echo "Gerando lista de contigs a serem extraidos."
    awk 'NR > 1 {print $1}' "${ASSEMBLY_RESULTS_PATH}".blast_1.contig_query.cov_filtered.results_summary.txt | uniq >> "${RESULTSDIR}"/mito_contigs_id_list.txt
    awk 'NR > 1 {print $1}' "${ASSEMBLY_RESULTS_PATH}".blast_2.contig_query.cov_filtered.results_summary.txt | uniq >> "${RESULTSDIR}"/mito_contigs_id_list.txt

    seqtk subseq "$ASSEMBLY_PATH" "${RESULTSDIR}"/mito_contigs_id_list.txt > "${CONTIGSDIR}"/found_mito_contig.fasta
    echo "Contigs mitocondriais salvos em: ${CONTIGSDIR}/found_mito_contig.fasta"
fi

echo "-------------------------------------------------------"

echo "Iniciando montagem de genoma mitocondrial com NOVOPlasty."

if [ -d "${RESULTSDIR}/NOVOPlasty_results" ]; then
    echo "Montagem de genoma motocondrial já realizada."
else
    perl "${BASEDIR}/mtDNA_processing_results/scripts/NOVOPlasty_scripts/NOVOPlasty.pl" -c "${BASEDIR}/mtDNA_processing_results/scripts/NOVOPlasty_scripts/config.txt"
    echo "Montagem disponivel em: ${RESULTSDIR}/NOVOPlasty_results"
fi

echo "Montagem realizada com sucesso."












