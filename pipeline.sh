#!/bin/bash
# chmod +x pipeline.sh

# Para a execução se algum comando falhar
set -e

BASENAME="Metschnikowia"
GENOME_SIZE="15980000"                            # Tamanho aproximado do genoma em pb
R1_DIR="./reads/Y7002_S238_L003_R1_001.fastq"     # Caminho para read forward
R2_DIR="./reads/Y7002_S238_L003_R2_001.fastq"     # Caminho para read reverse
OUTPUT_DIR="."                                    # Pasta onde os resultados serão salvos
THREADS=8                                         # Número de threads/núcleos a utilizar
MEMORY=64                                         # Quantidade de RAM máxima a utilizar
TRIMMOMATIC_ADAPTERS="/opt/apps/envs/envs/trimmomatic/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa"


calc_theor_cov() {
    local fastqc_data_txt="$1"
    local origin="$2"
    local output_file="${OUTPUT_DIR}/fastqc_results/theor_cov.txt"

    if [ ! -f "$fastqc_data_txt" ]; then
        echo "[Aviso] Arquivo $fastqc_data_txt não encontrado. Pulando cobertura para $origin."
        return 0
    fi

    local total_sequences=$(grep "Total Sequences" "$fastqc_data_txt" | cut -f2)
    local raw_length=$(grep "Sequence length" "$fastqc_data_txt" | cut -f2)

    local sequence_length=$(echo "$raw_length" | awk -F'-' '{if(NF==2) print ($1+$2)/2; else print $1}')

    if [ -z "$total_sequences" ] || [ -z "$sequence_length" ] || [ -z "$GENOME_SIZE" ]; then
        return 1
    fi

    local teorical_cov=$(awk -v total="$total_sequences" -v len="$sequence_length" -v size="$GENOME_SIZE" 'BEGIN {printf "%.2f", (total * len) / size}')

    if [ ! -f "$output_file" ]; then
        echo -e "origin\ttotal_sequences\tsequence_length\tgenome_size\tteorical_cov" > "$output_file"
    fi
    
    if grep -q "^${origin}\t" "$output_file"; then
        return 0
    fi

    if ! grep -q "^${origin}\t" "$output_file" 2>/dev/null; then
        echo -e "${origin}\t${total_sequences}\t${raw_length}\t${GENOME_SIZE}\t${teorical_cov}" >> "$output_file"
    fi
}

run_fastqc() {
    local INPUT_R1=$1
    local INPUT_R2=$2
    local OUTDIR=$3

    local base_r1=$(basename "${INPUT_R1}" | sed -E 's/\.fastq(\.gz)?$//')
    local base_r2=$(basename "${INPUT_R2}" | sed -E 's/\.fastq(\.gz)?$//')

    local skip_r1=0
    local skip_r2=0

    if [ -f "${OUTDIR}/${base_r1}_fastqc.html" ]; then
        echo "[SKIP] FastQC para ${INPUT_R1} já foi executado."
        skip_r1=1
    fi
    if [ -f "${OUTDIR}/${base_r2}_fastqc.html" ]; then
        echo "[SKIP] FastQC para ${INPUT_R2} já foi executado."
        skip_r2=1
    fi
    
    if [ $skip_r1 -eq 0 ] || [ $skip_r2 -eq 0 ]; then
        mkdir -p "$OUTDIR"
        echo "Executando FastQC em ${INPUT_R1} e ${INPUT_R2}"
        fastqc --extract --threads ${THREADS} --outdir "${OUTDIR}" "${INPUT_R1}" "${INPUT_R2}"
    fi

    calc_theor_cov "${OUTDIR}/${base_r1}_fastqc/fastqc_data.txt" "${INPUT_R1}"
    calc_theor_cov "${OUTDIR}/${base_r2}_fastqc/fastqc_data.txt" "${INPUT_R2}"
}

run_fastp() {
    local INPUT_R1=$1
    local INPUT_R2=$2
    local OUTDIR="${OUTPUT_DIR}/fastp_results"
    
    if [ -f "${OUTDIR}/${BASENAME}_fastp_R1.fastq.gz" ]; then
        echo "[SKIP] Fastp para ${BASENAME} já foi executado."
        return 0
    fi

    mkdir -p "$OUTDIR"
    echo "[Fastp] Trimando ${BASENAME}..."

    fastp \
        -i "${INPUT_R1}" -I "${INPUT_R2}" \
        -o "${OUTDIR}/${BASENAME}_fastp_R1.fastq.gz" \
        -O "${OUTDIR}/${BASENAME}_fastp_R2.fastq.gz" \
        --thread ${THREADS} \
        --html "${OUTDIR}/${BASENAME}_fastp.html" \
        --json "${OUTDIR}/${BASENAME}_fastp.json"
}

run_trimmomatic() {
    local INPUT_R1=$1
    local INPUT_R2=$2
    local OUTDIR="${OUTPUT_DIR}/trimmomatic_results"

    if [ -f "${OUTDIR}/${BASENAME}_trim_R1_paired.fastq.gz" ]; then
        echo "[SKIP] Trimmomatic para ${BASENAME} já foi executado."
        return 0
    fi

    mkdir -p "$OUTDIR"
    echo "[Trimmomatic] Trimando ${BASENAME}..."

    trimmomatic PE -threads ${THREADS} \
        "${INPUT_R1}" "${INPUT_R2}" \
        "${OUTDIR}/${BASENAME}_trim_R1_paired.fastq.gz" \
        "${OUTDIR}/${BASENAME}_trim_R1_unpaired.fastq.gz" \
        "${OUTDIR}/${BASENAME}_trim_R2_paired.fastq.gz" \
        "${OUTDIR}/${BASENAME}_trim_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:"${TRIMMOMATIC_ADAPTERS}":2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
}

refining_contigs(){
    local INPUT_R1=$1
    local INPUT_R2=$2
    local SPAdes_OUT=$3

    local CONTIGS="${SPAdes_OUT}/contigs.fasta"
    local FINAL_BAM="${SPAdes_OUT}/assemble_mapping.bam"
    local REPORT="${SPAdes_OUT}/assemble_coverage.txt"

    if [ -f "$REPORT" ]; then
        echo "[SKIP] Mapeamento e cobertura para ${SPAdes_OUT} já concluídos."
    else
        bwa index "${CONTIGS}"

        echo "[BWA + Samtools] Mapeando reads e gerando BAM ordenado..."
        bwa mem -t ${THREADS} "${CONTIGS}" "${INPUT_R1}" "${INPUT_R2}" | \
            samtools view -@ ${THREADS} -b - | \
            samtools sort -@ ${THREADS} -o "${FINAL_BAM}" -

        echo "[Samtools] Indexando o BAM ordenado..."
        samtools index "${FINAL_BAM}"

        echo "[Samtools] Calculando a cobertura..."
        samtools coverage "${FINAL_BAM}" > "${REPORT}"
        echo "Relatório de cobertura salvo em: ${REPORT}"
    fi

    local CLEAN_CONTIGS="${SPAdes_OUT}/clean_contigs.fasta"
    local CLEAN_FINAL_BAM="${SPAdes_OUT}/clean_assemble_mapping.bam"
    local CLEAN_REPORT="${SPAdes_OUT}/clean_assemble_coverage.txt"

    echo "[SKIP] Filtrando contigs menores que 500 Pb."
    seqtk seq -L 500 "${CONTIGS}" > "${CLEAN_CONTIGS}"

    if [ -f "$CLEAN_REPORT" ]; then
        echo "[SKIP] Filtragem de reads para ${SPAdes_OUT} já concluídos."
    else
        bwa index "${CLEAN_CONTIGS}"

        echo "[BWA + Samtools] Mapeando reads e gerando BAM ordenado..."
        bwa mem -t ${THREADS} "${CLEAN_CONTIGS}" "${INPUT_R1}" "${INPUT_R2}" | \
            samtools view -@ ${THREADS} -b - | \
            samtools sort -@ ${THREADS} -o "${CLEAN_FINAL_BAM}" -

        echo "[Samtools] Indexando o BAM ordenado..."
        samtools index "${CLEAN_FINAL_BAM}"

        echo "[Samtools] Calculando a cobertura..."
        samtools coverage "${CLEAN_FINAL_BAM}" > "${CLEAN_REPORT}"
        echo "Relatório de cobertura salvo em: ${CLEAN_REPORT}"
    fi

}

run_assembly_and_eval() {
    local INPUT_R1=$1
    local INPUT_R2=$2
    local SPAdes_OUT=$3
    local QUAST_OUT=$4
    local PARAMS=$5

    if [ -f "${SPAdes_OUT}/contigs.fasta" ]; then
        echo "[SKIP] SPAdes em ${SPAdes_OUT} já concluído."
    else
        mkdir -p "${SPAdes_OUT}"
        echo "[SPAdes] Montando em ${SPAdes_OUT}..."
        ../programs/SPAdes-4.2.0-Linux/bin/spades.py \
            -1 "${INPUT_R1}" \
            -2 "${INPUT_R2}" \
            -o "${SPAdes_OUT}" \
            -m "${MEMORY}" \
            -t "${THREADS}" \
            ${PARAMS}
    fi
    
    refining_contigs "${INPUT_R1}" "${INPUT_R2}" "${SPAdes_OUT}"
    
    if [ -f "${QUAST_OUT}/report.txt" ]; then
        echo "[SKIP] QUAST em ${QUAST_OUT} já concluído."
    else
        mkdir -p "${QUAST_OUT}"
        echo "[QUAST] Avaliando montagem em ${QUAST_OUT}..."
        ../programs/quast-5.3.0/quast.py \
            "${SPAdes_OUT}/clean_contigs.fasta" \
            -o "${QUAST_OUT}" \
            -t "${THREADS}"
    fi
}

run_busco_assembly(){
    local INPUT=$1
    local BUSCO_OUT=$2
    local PARAMS=$3

    if [ -f "${BUSCO_OUT}/run_saccharomycetaceae_odb12/short_summary.txt" ]; then
        echo "[SKIP] BUSCO em ${BUSCO_OUT} já concluído."
    else
        mkdir -p "${BUSCO_OUT}"
        echo "[BUSCO] Avaliando em ${BUSCO_OUT}..."
        busco -i "${INPUT}" -o "${BUSCO_OUT}" \
            -m genome -l saccharomycetaceae_odb12 \
            -c "${THREADS}" -f # "${PARAMS}"
        
    fi
}

downsampling() {
    local INPUT_R1=$1
    local INPUT_R2=$2
    local THEOR_COV=$3

    local OUTDIR=$(dirname "${INPUT_R1}")
    local BASE_NAME_R1=$(basename "${INPUT_R1}" | sed -E 's/\.fastq(\.gz)?$//')
    local BASE_NAME_R2=$(basename "${INPUT_R2}" | sed -E 's/\.fastq(\.gz)?$//')

    local GOAL_COV="80"
    local PROPORTION=$(awk -v goal_cov="$GOAL_COV" -v theor_cov="$THEOR_COV" 'BEGIN {printf "%.2f", goal_cov / theor_cov}')
    local SEED="42"

    echo "Subamostrando ${INPUT_R1}"
    seqtk sample -s${SEED} ${INPUT_R1} ${PROPORTION} | gzip > "${OUTDIR}/${BASE_NAME_R1}_80x.fastq.gz"
    seqtk sample -s${SEED} ${INPUT_R2} ${PROPORTION} | gzip > "${OUTDIR}/${BASE_NAME_R2}_80x.fastq.gz"
}

echo "Iniciando QC e Trimagem..."

DIR_RAW_QC="${OUTPUT_DIR}/fastqc_results/raw_reads"
run_fastqc "${R1_DIR}" "${R2_DIR}" "${DIR_RAW_QC}"

echo "--------------------------------------------------"
echo "Processando amostra: ${BASENAME}"
echo "--------------------------------------------------"

run_fastp "${R1_DIR}" "${R2_DIR}"
run_trimmomatic "${R1_DIR}" "${R2_DIR}"

FASTP_R1="${OUTPUT_DIR}/fastp_results/${BASENAME}_fastp_R1.fastq.gz"
FASTP_R2="${OUTPUT_DIR}/fastp_results/${BASENAME}_fastp_R2.fastq.gz"
TRIM_R1="${OUTPUT_DIR}/trimmomatic_results/${BASENAME}_trim_R1_paired.fastq.gz"
TRIM_R2="${OUTPUT_DIR}/trimmomatic_results/${BASENAME}_trim_R2_paired.fastq.gz"

echo "--------------------------------------------------"
echo "Executando FastQC nas reads pós-processamento..."

DIR_FASTP_QC="${OUTPUT_DIR}/fastqc_results/fastp_reads"
run_fastqc "${FASTP_R1}" "${FASTP_R2}" "${DIR_FASTP_QC}"

DIR_TRIM_QC="${OUTPUT_DIR}/fastqc_results/trimmomatic_reads"
run_fastqc "${TRIM_R1}" "${TRIM_R2}" "${DIR_TRIM_QC}"

echo "--------------------------------------------------"
echo "Iniciando Montagens e Avaliações (SPAdes + QUAST)..."

# Cenário 1: RawReads + Isolate
run_assembly_and_eval \
    "${R1_DIR}" "${R2_DIR}" \
    "${OUTPUT_DIR}/SPAdes_results/raw_reads_isolate" \
    "${OUTPUT_DIR}/quast_results/raw_reads_isolate" \
    "--isolate"

# Cenário 2: RawReads + Careful
run_assembly_and_eval \
    "${R1_DIR}" "${R2_DIR}" \
    "${OUTPUT_DIR}/SPAdes_results/raw_reads_careful" \
    "${OUTPUT_DIR}/quast_results/raw_reads_careful" \
    "--careful"

# Cenário 3: Fastp + Isolate
run_assembly_and_eval \
    "${FASTP_R1}" "${FASTP_R2}" \
    "${OUTPUT_DIR}/SPAdes_results/fastp_reads_isolate" \
    "${OUTPUT_DIR}/quast_results/fastp_reads_isolate" \
    "--isolate"

# Cenário 4: Fastp + Careful
run_assembly_and_eval \
    "${FASTP_R1}" "${FASTP_R2}" \
    "${OUTPUT_DIR}/SPAdes_results/fastp_reads_careful" \
    "${OUTPUT_DIR}/quast_results/fastp_reads_careful" \
    "--careful"

# Cenário 5: Trimmomatic + Isolate
run_assembly_and_eval \
    "${TRIM_R1}" "${TRIM_R2}" \
    "${OUTPUT_DIR}/SPAdes_results/trimmomatic_reads_isolate" \
    "${OUTPUT_DIR}/quast_results/trimmomatic_reads_isolate" \
    "--isolate"

# Cenário 6: Trimmomatic + Careful
run_assembly_and_eval \
    "${TRIM_R1}" "${TRIM_R2}" \
    "${OUTPUT_DIR}/SPAdes_results/trimmomatic_reads_careful" \
    "${OUTPUT_DIR}/quast_results/trimmomatic_reads_careful" \
    "--careful"

# Cenário 7: Trimmomatic + Isolate + cov-cutoff
run_assembly_and_eval \
    "${TRIM_R1}" "${TRIM_R2}" \
    "${OUTPUT_DIR}/SPAdes_results/trimmomatic_reads_isolate_cutoff_80" \
    "${OUTPUT_DIR}/quast_results/trimmomatic_reads_isolate_cutoff_80" \
    "--isolate --cov-cutoff 80"

echo "--------------------------------------------------"
echo "Avaliando completude das montagens (BUSCO)"

run_busco_assembly \
    "${OUTPUT_DIR}/SPAdes_results/fastp_reads_careful/clean_contigs.fasta" \
    "${OUTPUT_DIR}/BUSCO_results/fastp_reads_careful"

run_busco_assembly \
    "${OUTPUT_DIR}/SPAdes_results/fastp_reads_isolate/clean_contigs.fasta" \
    "${OUTPUT_DIR}/BUSCO_results/fastp_reads_isolate"

run_busco_assembly \
    "${OUTPUT_DIR}/SPAdes_results/raw_reads_careful/clean_contigs.fasta" \
    "${OUTPUT_DIR}/BUSCO_results/raw_reads_careful"

run_busco_assembly \
    "${OUTPUT_DIR}/SPAdes_results/raw_reads_isolate/clean_contigs.fasta" \
    "${OUTPUT_DIR}/BUSCO_results/raw_reads_isolate"

run_busco_assembly \
    "${OUTPUT_DIR}/SPAdes_results/trimmomatic_reads_careful/clean_contigs.fasta" \
    "${OUTPUT_DIR}/BUSCO_results/trimmomatic_reads_careful"

run_busco_assembly \
    "${OUTPUT_DIR}/SPAdes_results/trimmomatic_reads_isolate/clean_contigs.fasta" \
    "${OUTPUT_DIR}/BUSCO_results/trimmomatic_reads_isolate"
    
run_busco_assembly \
    "${OUTPUT_DIR}/SPAdes_results/trimmomatic_reads_isolate_cutoff_80/clean_contigs.fasta" \
    "${OUTPUT_DIR}/BUSCO_results/trimmomatic_reads_isolate_cutoff_80"


echo "=================================================="
echo "Pipeline concluído com sucesso!"
echo "=================================================="