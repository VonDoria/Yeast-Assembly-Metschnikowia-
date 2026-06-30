import os
import sys
import glob
import csv
from collections import defaultdict

REF_ASSEMBLYS_FOLDER = sys.argv[1]
SAMPLE_SC_GENES_FOLDER = sys.argv[2]
OUTPUT_DIR = sys.argv[3]
SEQ_TYPE = sys.argv[4]
MIN_PRESENCE_THRESHOLD = float(sys.argv[5])


def generate_assemblys_dict():
    assemblys_dict = {}
    
    assemblys_dict["Sample_Metschnikowia"] = SAMPLE_SC_GENES_FOLDER

    if os.path.exists(REF_ASSEMBLYS_FOLDER):
        for item in os.listdir(REF_ASSEMBLYS_FOLDER):
            full_path = os.path.join(REF_ASSEMBLYS_FOLDER, item)
            
            if os.path.isdir(full_path):
                target_folder = os.path.join(
                    full_path, 
                    "BUSCO_results", 
                    "run_saccharomycetes_odb12", 
                    "busco_sequences", 
                    "single_copy_busco_sequences"
                )
                assemblys_dict[item] = target_folder
    else:
        print(f"⚠️A pasta de referências '{REF_ASSEMBLYS_FOLDER}' não foi encontrada.")

    print(f"   -> {len(assemblys_dict)} amostras registradas para análise.")
    return assemblys_dict

def mapping_single_copy_genes(assemblys_dict):    
    genes_by_sample = defaultdict(dict)
    all_keys_genes = set()

    for sample, sc_genes_folder in assemblys_dict.items():
        search_path = os.path.join(sc_genes_folder, f"*.{SEQ_TYPE}")
        genes_file = glob.glob(search_path)
        
        if not genes_file:
            print(f"   ⚠️Nenhum arquivo '{SEQ_TYPE}' em {sc_genes_folder}")
            continue

        for file_path in genes_file:
            gene_id = os.path.basename(file_path).split(f".{SEQ_TYPE}")[0]
            genes_by_sample[sample][gene_id] = file_path
            all_keys_genes.add(gene_id)

        print(f"   -> Amostra '{sample}': {len(genes_file)} genes encontrados.")

    # CRIAÇÃO DO ARQUIVO CSV DE PRESENÇA
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    csv_path = os.path.join(OUTPUT_DIR, "gene_presence_matriz.csv")
    
    genes_list = sorted(list(all_keys_genes))
    
    with open(csv_path, mode="w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        
        writer.writerow(["Amostra"] + genes_list)
        
        for sample in assemblys_dict.keys():
            row = [sample]
            for gene in genes_list:
                valor = 1 if gene in genes_by_sample[sample] else 0
                row.append(valor)
            writer.writerow(row)

    print(f"\n📊 Matriz de presença gerada com sucesso em: {csv_path}")

    return genes_by_sample, all_keys_genes

def filter_common_genes(genes_by_sample, all_genes, assemblys_dict):
    num_total_sample = len(assemblys_dict)
    common_genes = []

    if num_total_sample == 0:
        return common_genes

    for gene_id in all_genes:
        presence = sum(1 for sample in assemblys_dict if gene_id in genes_by_sample[sample])
        presence_tax = presence / num_total_sample

        if presence_tax >= MIN_PRESENCE_THRESHOLD:
            common_genes.append(gene_id)

    print(f"\n🎯 Filtragem de Ortólogos Concluída:")
    print(f"   -> Genes totais analisados: {len(all_genes)}")
    print(f"   -> Genes retidos (presença >= {MIN_PRESENCE_THRESHOLD*100}%): {len(common_genes)}")
    
    return common_genes

def group_sequences(genes_by_sample, common_genes, assemblys_dict):
    unaligned_folder = os.path.join(OUTPUT_DIR, "genes_unaligned")
    os.makedirs(unaligned_folder, exist_ok=True)

    print("\n🧬 Criando arquivos FASTA agrupados (não alinhados)...")
    
    for idx, gene_id in enumerate(common_genes):
        out_fasta = os.path.join(unaligned_folder, f"{gene_id}.fasta")

        with open(out_fasta, "w", encoding="utf-8") as f_out:
            for sample in assemblys_dict.keys():
                if gene_id in genes_by_sample.get(sample, {}):
                    gene_path = genes_by_sample[sample][gene_id]
                    
                    with open(gene_path, "r") as f_in:
                        lines = f_in.readlines()
                        seq = "".join(l.strip() for l in lines if not l.startswith(">"))
                        f_out.write(f">{sample}\n{seq}\n")

        if (idx + 1) % 100 == 0 or (idx + 1) == len(common_genes):
            print(f"   -> Processados: {idx + 1}/{len(common_genes)} arquivos.")

    print(f"\n✅ Agrupamento finalizado! Arquivos salvos na pasta: {unaligned_folder}")

# -------------------------------------------------------------------------------

assemblys_dict = generate_assemblys_dict()

if not assemblys_dict:
    print("❌ Erro: Nenhuma amostra foi encontrada nos diretórios especificados.")
    exit(1)
    
genes_by_sample, all_genes = mapping_single_copy_genes(assemblys_dict)

common_genes = filter_common_genes(genes_by_sample, all_genes, assemblys_dict)

if not common_genes:
    print("❌ Erro: Nenhum gene compartilhado atingiu o limite mínimo de presença.")
    exit(1)
    
group_sequences(genes_by_sample, common_genes, assemblys_dict)