#!/bin/bash
# chmod +x download_NCBI_assemblys.sh

# Para a execução se algum comando falhar
set -e

THREADS=8                                         # Número de threads/núcleos a utilizar
MEMORY=64                                         # Quantidade de RAM máxima a utilizar

GENOME_LIST=(
    "GCA_900106115.1"	#Sungouiella intermedia(OUTGROUP)			
    "GCA_030556055.1"	#Metschnikowia bicuspidata var chathamia			
    "GCA_030555835.1"	#Metschnikowia bicuspidata var californica			
    "GCA_030444895.1"	#Metschnikowia sp. yHMJ9
    "GCA_030444905.1"	#Metschnikowia sp. yHKB443
    "GCA_030444905.1"	#Metschnikowia sp. yHKB443
    "GCA_003708715.3"	#Candida wancherniae			
    "GCA_030579135.1"	#Candida danieliae			
    "GCA_030561765.1"	#Candida hainanensis			
    "GCA_030582695.1"	#Candida golubevii			
    "GCA_002370145.1"	#Candida ipomoeae			
    # "GCA_939531405.1"	#Metschnikowia zobellii			
    # "GCA_003401635.1"	#Metschnikowia reukaufii			
    # "GCA_001664035.1"	#Metschnikowia bicuspidata var. bicuspidata NRRL YB-4993			
    # "GCA_030563445.1"	#Metschnikowia gruessii			
    # "GCA_002893725.1"	#Metschnikowia torresii			
    # "GCA_002370615.1"	#Metschnikowia aberdeeniae			
    # "GCA_002374725.1"	#Metschnikowia hibisci			
    # "GCA_002370515.1"	#Metschnikowia proteae			
    # "GCA_002893665.1"	#Metschnikowia orientalis			
    # "GCA_002374535.1"	#Metschnikowia kamakouana			
    # "GCA_030583255.1"	#Metschnikowia kunwiensis			
    # "GCA_002374555.1"	#Metschnikowia mauinuiana			
    # "GCA_030561945.1"	#Metschnikowia krissii			
    # "GCA_030583235.1"	#Metschnikowia lunata			
    # "GCA_002370695.1"	#Metschnikowia matae var. maris			
    # "GCA_030561745.1"	#Metschnikowia gelsemii			
    # "GCA_030582795.1"	#Metschnikowia chrysomelidarum			
    # "GCA_002374405.1"	#Metschnikowia cubensis			
    # "GCA_030569435.1"	#Metschnikowia koreensis			
    # "GCA_021272955.1"	#Metschnikowia ahupensis			
    # "GCA_030583145.1"	#Metschnikowia vanudenii			
    # "GCA_002374485.1"	#Metschnikowia santaceciliae			
    # "GCA_002370635.1"	#Metschnikowia cerradonensis			
    # "GCA_030578735.1"	#Metschnikowia noctiluminum			
    # "GCA_030565705.1"	#Metschnikowia baotianmanensis			
    # "GCA_030573015.1"	#Metschnikowia peoriensis			
    # "GCA_030556455.1"	#Metschnikowia pimensis			
    # "GCA_030572615.1"	#Metschnikowia lachancei			
    # "GCA_030556465.1"	#Metschnikowia picachoensis			
    # "GCA_030581935.1"	#Metschnikowia corniflorae			
    # "GCA_030674755.1"	#Metschnikowia henanensis			
    # "GCA_030557065.1"	#Metschnikowia rubicola			
    # "GCA_030563125.1"	#Metschnikowia laotica			
    # "GCA_030556725.1"	#Metschnikowia viticola			
    # "GCA_030674525.1"	#Metschnikowia chrysoperlae			
    # "GCA_009746055.1"	#Metschnikowia pulcherrima			
    # "GCA_019285895.1"	#Metschnikowia australis			
    # "GCA_963678735.1"	#Metschnikowia agaves			
    # "GCA_008065175.1"	#Metschnikowia caudata			
    # "GCA_002370135.1"	#Metschnikowia kipukae			
    # "GCA_002893645.1"	#Metschnikowia hawaiiana			
    # "GCA_030563105.1"	#Metschnikowia lopburiensis			
    # "GCA_002893705.1"	#Metschnikowia drosophilae			
    # "GCA_002374645.1"	#Metschnikowia shivogae			
    # "GCA_002370475.1"	#Metschnikowia drakensbergensis			
    # "GCA_030573055.1"	#Metschnikowia anglica			
    # "GCA_002370175.1"	#Metschnikowia colocasiae			
    # "GCA_009756545.1"	#Metschnikowia dekortorum			
    # "GCA_009756645.1"	#Metschnikowia lacustris			
    # "GCA_002370295.1"	#Metschnikowia bowlesiae			
    # "GCA_002370765.1"	#Metschnikowia similis			
    # "GCA_002370325.1"	#Metschnikowia hawaiiensis			
    # "GCA_002370875.1"	#Metschnikowia arizonensis			
    # "GCA_002370815.1"	#Metschnikowia hamakuensis			
    # "GCA_002370915.1"	#Metschnikowia lochheadii			
    # "GCA_002370835.1"	#Metschnikowia continentalis			
    # "GCA_002374385.1"	#Metschnikowia borealis			
    # "GCA_030564885.1"	#Metschnikowia kofuensis (nom. inval.)			
    # "GCA_008065195.1"	#Metschnikowia amazonensis		
)

ncbi_dataset_get_assembly(){
    local ACCESSION=$1
    
    species_name=$(datasets summary genome accession "${ACCESSION}" | grep -oP '"organism_name":"[^"]+"' | head -1 | cut -d'"' -f4)
    
    if [ -z "$species_name" ]; then
        species_name="Yeast_${ACCESSION}"
    fi
    
    safe_folder_name=$(echo "$species_name" | tr ' ' '_')

    
    if [ -d "${safe_folder_name}" ]; then
        echo "[SKIP] Download para ${safe_folder_name} já concluído."
    else       
        echo "📁 Criando pasta e baixando dados para: $species_name"
        mkdir -p "$safe_folder_name"

        # gff3 (anotação), genome (DNA fasta), gbff (genbank format)
        # gff3,genome,gbff 
        datasets download genome accession "${ACCESSION}" \
            --assembly-source all \
            --include genome \
            --filename "${safe_folder_name}/${ACCESSION}_dataset.zip"
            
        echo "📦 Descompactando arquivos..."
        cd "$safe_folder_name"
        unzip -q "${ACCESSION}_dataset.zip"
        
        mv ncbi_dataset/data/$ACCESSION/* .
        rm -rf ncbi_dataset README.md dataset_catalog.json md5sum.txt "${ACCESSION}_dataset.zip"
        

        cd ..
    fi
    
    run_busco_assembly ${safe_folder_name}/${ACCESSION}*_genomic.fna ${safe_folder_name}/BUSCO_results
}

run_busco_assembly(){
    local INPUT=$1
    local BUSCO_OUT=$2
    local PARAMS=$3

    if [ -f "${BUSCO_OUT}/run_saccharomycetes_odb12/short_summary.txt" ]; then
        echo "[SKIP] BUSCO em ${BUSCO_OUT} já concluído."
    else
        mkdir -p "${BUSCO_OUT}"
        echo "[BUSCO] Avaliando em ${BUSCO_OUT}..."
        busco -i "${INPUT}" -o "${BUSCO_OUT}" \
            -m genome -l ../busco_downloads/lineages/saccharomycetes_odb12 \
            -c "${THREADS}" -f # "${PARAMS}"
        
    fi
}


for accession in "${GENOME_LIST[@]}"; do
    echo "----------------------------------------"
    echo "Processando montagem: $accession"
    
    ncbi_dataset_get_assembly "${accession}"

    echo "✅ Concluído para $accession"
done

echo "🎉 Todos os downloads foram finalizados!"