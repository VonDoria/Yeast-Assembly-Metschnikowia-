#!/usr/bin/env python3
import os
import json
import statistics

OUTPUT_DIR = "./"
DIR_FASTQC = os.path.join(OUTPUT_DIR, "fastqc_results")
DIR_QUAST = os.path.join(OUTPUT_DIR, "quast_results")
DIR_SPADES = os.path.join(OUTPUT_DIR, "SPAdes_results")
DIR_BUSCO = os.path.join(OUTPUT_DIR, "BUSCO_results")

# Função auxiliar para calcular estatísticas completas sem dependências externas (como pandas/numpy)
def calc_stats(valores):
    if not valores:
        return {"mean": 0, "std": 0, "min": 0, "q1": 0, "median": 0, "q3": 0, "max": 0}
    valores.sort()
    n = len(valores)
    return {
        "mean": statistics.mean(valores),
        "std": statistics.stdev(valores) if n > 1 else 0,
        "min": min(valores),
        "q1": statistics.quantiles(valores, n=4)[0] if n > 3 else valores[0],
        "median": statistics.median(valores),
        "q3": statistics.quantiles(valores, n=4)[2] if n > 3 else valores[-1],
        "max": max(valores)
    }

def extrair_dados_fastqc(caminho_pasta):
    dados = { "total_sequences": "-", "gc_content": "-", "base_quality": "-", "adapters": "-", "total_bases": "-", "poor_quality_sequences": "-", "sequence_length": "-" }
    fastqc_data = os.path.join(caminho_pasta, "fastqc_data.txt")
    summary = os.path.join(caminho_pasta, "summary.txt")

    if os.path.exists(fastqc_data):
        with open(fastqc_data, 'r', encoding='utf-8') as f:
            for linha in f:
                if linha.startswith("Total Sequences"): dados["total_sequences"] = f"{int(linha.strip().split('\t')[1]):,}".replace(",", ".")
                elif linha.startswith("Total Bases"): dados["total_bases"] = linha.strip().split('\t')[1]
                elif linha.startswith("Sequences flagged as poor quality"): dados["poor_quality_sequences"] = linha.strip().split('\t')[1]
                elif linha.startswith("Sequence length"): dados["sequence_length"] = linha.strip().split('\t')[1]
                elif linha.startswith("%GC"): dados["gc_content"] = linha.strip().split('\t')[1] + "%"

    if os.path.exists(summary):
        with open(summary, 'r', encoding='utf-8') as f:
            for linha in f:
                partes = linha.strip().split('\t')
                if len(partes) >= 2:
                    status, modulo = partes[0], partes[1]
                    if modulo == "Per base sequence quality": dados["base_quality"] = status
                    elif modulo == "Adapter Content": dados["adapters"] = status
    return dados

def extrair_dados_quast(caminho_tsv):
    dados = { 
        "contigs_count": "-", 
        "largest_contig": "-", 
        "total_bp": "-", 
        "GC(%)": "-", 
        "n50": "-", 
        "n90": "-", 
        "l50": "-", 
        "l90": "-", 
        "contigs_0": "-", 
        "contigs_1000": "-", 
        "contigs_5000": "-", 
        "contigs_10000": "-", 
        "contigs_25000": "-", 
        "contigs_50000": "-", 
        "total_bp_0": "-", 
        "total_bp_1000": "-", 
        "total_bp_5000": "-", 
        "total_bp_10000": "-", 
        "total_bp_25000": "-", 
        "total_bp_50000": "-", 
    }
    try:
        with open(caminho_tsv, 'r', encoding='utf-8') as f:
            linhas = f.readlines()
            if len(linhas) < 2: return dados
            headers, valores = linhas[0].strip().split('\t'), linhas[1].strip().split('\t')
            for i, h in enumerate(headers):
                if "# contigs" in h and ">=" not in h: dados["contigs_count"] = valores[i]
                elif "Largest contig" in h: dados["largest_contig"] = valores[i]
                elif "GC (%)" in h: dados["GC(%)"] = valores[i]
                elif "N50" in h: dados["n50"] = valores[i]
                elif "N90" in h: dados["n90"] = valores[i]
                elif "L50" in h: dados["l50"] = valores[i]
                elif "L90" in h: dados["l90"] = valores[i]
                elif "Total length" in h and ">=" not in h: dados["total_bp"] = valores[i]
                elif "# contigs (>= 0 bp)" in h: dados["contigs_0"] = valores[i]
                elif "# contigs (>= 1000 bp)" in h: dados["contigs_1000"] = valores[i]
                elif "# contigs (>= 5000 bp)" in h: dados["contigs_5000"] = valores[i]
                elif "# contigs (>= 10000 bp)" in h: dados["contigs_10000"] = valores[i]
                elif "# contigs (>= 25000 bp)" in h: dados["contigs_25000"] = valores[i]
                elif "# contigs (>= 50000 bp)" in h: dados["contigs_50000"] = valores[i]
                elif "Total length (>= 0 bp)" in h: dados["total_bp_0"] = valores[i]
                elif "Total length (>= 1000 bp)" in h: dados["total_bp_1000"] = valores[i]
                elif "Total length (>= 5000 bp)" in h: dados["total_bp_5000"] = valores[i]
                elif "Total length (>= 10000 bp)" in h: dados["total_bp_10000"] = valores[i]
                elif "Total length (>= 25000 bp)" in h: dados["total_bp_25000"] = valores[i]
                elif "Total length (>= 50000 bp)" in h: dados["total_bp_50000"] = valores[i]
                elif h in dados: dados[h] = valores[i]
    except Exception: pass
    return dados

def extrair_dados_spades_coverage(caminho_cov):
    parametros = ["startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"]
    dados_brutos = {p: [] for p in parametros}
    
    try:
        with open(caminho_cov, 'r', encoding='utf-8') as f:
            linhas = f.readlines()
            for linha in linhas[1:]: # Pula cabeçalho
                partes = linha.strip().split('\t')
                if len(partes) < 9: continue
                
                # Regra: descartar linhas com endpos < 500
                endpos = float(partes[2])
                if endpos < 500: continue
                
                # Índices baseados na estrutura do samtools coverage (ignorando #rname no index 0)
                dados_brutos["startpos"].append(float(partes[1]))
                dados_brutos["endpos"].append(endpos)
                dados_brutos["numreads"].append(float(partes[3]))
                dados_brutos["covbases"].append(float(partes[4]))
                dados_brutos["coverage"].append(float(partes[5]))
                dados_brutos["meandepth"].append(float(partes[6]))
                dados_brutos["meanbaseq"].append(float(partes[7]))
                dados_brutos["meanmapq"].append(float(partes[8]))
    except Exception: return {}

    # Calcula as estatísticas de cada coluna (Média, Quartis, etc)
    return {p: calc_stats(dados_brutos[p]) for p in parametros}

def mapear_busco(diretorio_base):
    resultados = {}
    if not os.path.exists(diretorio_base): return resultados
    for root, dirs, files in os.walk(diretorio_base):
        for file in files:
            if file.startswith("short_summary") and file.endswith(".json"):
                caminho_completo = os.path.join(root, file)
                nome_cenario = os.path.basename(root)
                try:
                    with open(caminho_completo, 'r', encoding='utf-8') as f:
                        conteudo = json.load(f)
                        if "results" in conteudo:
                            resultados[nome_cenario] = conteudo["results"]
                except Exception as e: print(f"Erro no BUSCO json: {e}")
    return resultados

def ler_log(caminho_log, max_linhas=500):
    try:
        with open(caminho_log, 'r', encoding='utf-8', errors='ignore') as f: return "".join(f.readlines()[-max_linhas:])
    except Exception as e: return f"Erro ao ler log: {e}"

def get_logs():
    resultados = {}
    # QUAST Logs
    if os.path.exists(DIR_QUAST):
        for root, dirs, files in os.walk(DIR_QUAST):
            if "quast.log" in files: resultados[f"QUAST_{os.path.basename(root)}_log"] = ler_log(os.path.join(root, "quast.log"), 1000)
    # SPAdes Logs
    if os.path.exists(DIR_SPADES):
        for root, dirs, files in os.walk(DIR_SPADES):
            if "spades.log" in files: resultados[f"SPAdes_{os.path.basename(root)}_log"] = ler_log(os.path.join(root, "spades.log"), 1000)
    # BUSCO Logs
    if os.path.exists(DIR_BUSCO):
        for root, dirs, files in os.walk(DIR_BUSCO):
            if os.path.basename(root) == "logs":
                cenario = os.path.basename(os.path.dirname(root))
                for file in files:
                    resultados[f"BUSCO_{cenario}_{file}"] = ler_log(os.path.join(root, file), 1000)
    return resultados

def main():
    dados_spades = {}
    if os.path.exists(DIR_SPADES):
        for root, dirs, files in os.walk(DIR_SPADES):
            if "clean_assemble_coverage.txt" in files:
                dados_spades[os.path.basename(root)] = extrair_dados_spades_coverage(os.path.join(root, "clean_assemble_coverage.txt"))

    dados_finais = {
        "fastqc": {os.path.relpath(r, DIR_FASTQC): extrair_dados_fastqc(r) for r, _, f in os.walk(DIR_FASTQC) if "fastqc_data.txt" in f},
        "quast": {os.path.basename(r): extrair_dados_quast(os.path.join(r, "transposed_report.tsv")) for r, _, f in os.walk(DIR_QUAST) if "transposed_report.tsv" in f},
        "SPAdes": dados_spades,
        "busco": mapear_busco(DIR_BUSCO),
        "logs": get_logs()
    }
    
    with open(os.path.join("./dashboard", "dashboard_data.js"), 'w', encoding='utf-8') as f:
        f.write(f"const dadosProjeto = {json.dumps(dados_finais, indent=4, ensure_ascii=False)};")
    print("Sucesso! Dados compilados salvos.")

if __name__ == "__main__":
    main()