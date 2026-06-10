# Sheetcheat Ferramentas
(This file is simply a small notepad for referencing the definition of certain parameters and results. It can be disregarded or translated into your preferred language.)

## Fastqc
Avalia a qualidade das reads (adaptadores não removidos, baixa qualidade nas extremidades das leituras...). Retorna um relatório .html para cada uma das reads analizadas. E um arquivo .zip contendo os graficos do relatório e um .txt com os dados brutos do relatório.

**Comandos:**
```bash
    fastqc 
    amostra1.fastq.gz amostra2.fastq.gz             # Caminho das amostras
    -t 4                                            # Número de threads a serem usadas
    -o ./resultados_qc                              # Pasta de destino dos resultados
    --extract                                       # Extrai automaticamente a saída .zip
    > fastqc_execucao.log 2>&1                      # Direcionamento do log para um arquivo de log
```
#### Resultados
**Per base sequence quality:**
Mostra a qualidade por base ao longo das sequências. 
- Boa qualidade >30 (zona verde). É normal as ultimas bases cairem um pouco para o zona amarela no Illumina.
- Má qualidade <30 (zona amarela ou vermelha). Quando a qualidade cai antes do meio da leitura para a zona amarela ou vermelha. Indica problemas ocorridos durante a leitura (falta de reagente ou perda de foco do laser).

**Per Tile Sequence Quality:**
A flowcell da Illumina é dividida em pequenas áreas chamadas tiles (quadrantes). Este gráfico é um heatmap: o eixo X mostra a posição na leitura e o eixo Y mostra o código de cada tile.
- Resultado Bom: O gráfico inteiro com uma coloração azul homogênea. Isso significa que a qualidade foi uniforme em toda a superfície física da flowcell.
- Resultado Ruim: O surgimento de listras ou manchas amarelas, laranjas ou vermelhas. Isso indica uma falha local e física no sequenciador durante aquela corrida. Pode ser causado por uma bolha de ar que passou pelo canal, um cisco de poeira na lente do microscópio, ou uma sobreposição de clusters (overclustering).

**Per sequence quality scores:**
Mostra a quantidade de sequências que apresentam cada qualidade. 
- O ideal é que tenha um grande pico acima de 30 ou 35. Isso indica que muitas sequências tem boa qualidade.
- Um gráfico com 2 picos (um deles em Q15 por exemplo) indica que há um grupo de reads que não foi bem sequnciada.

**Per base sequence content:**
Mostra a prorpoção de cada base ao longo da sequência.
O resultado ideal é que todas as linhas estejam retas e coerentes com a proporção do genoma do organismo em questão (Organismo => 60%GC e 40%AT, então no gráfico G e C devem estar em 30% e A e T em 20%).
É normal que no começo da leitura ocorra ocilações. Isso é um defeito da transposase do kit do Illumina, ela pode apresentar preferência por determinada base no início.

**Per sequence GC content:**
Normalmente, a distribuição de GCs nas reads segue uma distribuição normal. Caso a distribuição de GC apresente mais de um pico, ou não se assemelhe a uma curva normal, isso pode indicar que ouve contaminação no seu sequenciamento, e essa distorção da curva representa o conteudo GC do contaminante.

**Per base N content:**
Mostra quais bases não poderam ser determinadas.
É comum ver pequenos picos de "N" logo no início da leitura se o equipamento falhar nos ciclos iniciais de calibração, ou uma subida no final se a reação química perder força. Se passar de 2-5%, essas leituras devem ser limpas. Uma linha horizontal no 0 é o ideial.

**Sequence Length Distribution:**
Mostra o tamanho em pares das leituras geradas. O ideal é uma linha vertical reta no tamanho de pb contratado no sequenciador. Uma distribuição de vários tamanhos é mau sinal. Após o processamento das leituras é normal aparecer tamanhos variados e o gráfico destoar do ideal.

**Sequence Duplication Levels:**
Mostra a proporção de sequências que são unicas, estão em duplicata, triplicata... O ideal é um pico a esquerda no grau 1 e caindo rapidamente para 0 nos graus posteriores.

**Overrepresented sequences:**
Sequências idênticas que correspondem a mais de 0,1% do total de leituras. Provavelmente, adaptadores ou contaminação. No caso de RNA-seq, pode ser RNA Ribossômico não filtrado na bancada.

**Adapter Content:**
O ideal é que esteja totalmente no 0. Caso o fragmento de DNA seja menor do que o previsto na hora de programar os ciclos do sequenciador, o laser passa pelo DNA e começa a sequênciar os adaptadores.

----

## Fastp
Faz o controle de qualidade (antes e depois do processamento), filtragem por qualidade, e remove os adaptadores em uma única passagem. Retorna um arquivo .html com o relatório e um .json com os dados do relatório, além dos arquivos .fastq trimados e limpos.

**Comandos:**
```bash
    fastp -i amostra_R1.fq.gz -I amostra_R2.fq.gz   # Arquivos de entrada forward e reverse
    -o limpo_R1.fq.gz -O limpo_R2.fq.gz             # Arquivos .fastq de saída (limpos e trimados)
    -w 8 -q 20                                      # Númeor de threads (8) e corte de qualidade (Q >= 20)
    -j fastp_relatorio.json -h fastp_relatorio.html # Nome dos arquivos de saída .html e .json
    2> fastp_execucao.log                           # Direcionamento do log para um arquivo de log
```

#### Resultados

**Insert size estimation**
Exclusivo para dados Paired-End. O fastp tenta sobrepor o R1 e o R2 para descobrir o tamanho real do fragmento de DNA que estava entre os adaptadores.
- Bom Resultado: Um pico bem definido que coincide exatamente com o tamanho de fragmento que você selecionou na bancada do laboratório (geralmente entre 250bp e 400bp).
- Mau Resultado: Um pico grande próximo a zero. Isso significa que sua biblioteca está cheia de dímeros de adaptadores (adaptadores que se ligaram uns aos outros sem nenhum DNA de interesse no meio). A etapa de purificação falhou no laboratório.

**Quality:**
Mostra a qualidade das bases ao longo da sequência de forma separada para R1, R2, antes e depois. É normal a R2 antes(before) demostrar uma queda de qualidade maior no final, mas isso é corrigido pelo fastp no R2 depois(after).

**Base contents:**
O gráfico plota a porcentagem de cada base (A, T, C, G e N) no eixo Y contra a posição da base na leitura no eixo X. As linhas das bases A e T devem se sobrepor, assim como C e G.

**KMER Counting:**
A frequência com que diferentes K-mers aparecem.

---

## Trimmomatc
Focada estritamente no corte (trimming) e filtragem de leituras Illumina. Requer que você conheça as sequências dos adaptadores que deseja remover. Retorna 2 arquivos para cada read inserida, um com as reads que continuaram pareadas (P) e um com as reads que não continuaram pareadas (U).

#### Comando
```bash
    java -jar trimmomatic.jar 
    PE                                              # Define o modo (Paired-End = PE ou Single-End = SE)
    -threads 8 -phred33                             # Número de threads e especifica a codificação de qualidade
    amostra_R1.fastq.gz amostra_R2.fastq.gz         # Arquivos de entrada
    pareado_R1.fastq.gz desemparelhado_R1.fastq.gz 
    pareado_R2.fastq.gz desemparelhado_R2.fastq.gz  # Arquivos de saída
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10              # Arquivo de adaptadores para serem cortados
    LEADING:3 TRAILING:3                            # Corta bases do início (leading) ou fim (trailing) se a qualidade for menor que o valor
    SLIDINGWINDOW:4:15                              #Corta a extremidade 3' se a qualidade média da janela cair abaixo do limite
    MINLEN:36                                       # Descarta a leitura inteira se ficar menor que este tamanho após os cortes
    > trimmomatic_execucao.log 2>&1                 # Direcionamento do log para um arquivo de 
```
---

## SPAdes 
(St. Petersburg genome assembler)

Montador ideal para genomas pequenos (Bactéria, Virus, Fungos) e metagenomas.
Pipeline do SPAdes:
**1 Montagem de Grafos com Múltiplos K-mers:** 
O SPAdes quebra as reads em pedaços menores (K-mers) para montar o genoma e depois compara com uma montagem feita com K-mers maiores. Ele faz isso para vários tamanhos de k-mers diferentes (21, 33, 55, 77...). Isso garante que o genoma fique continuo e resolve regiões repetitivas.
**2 Simplificação do Grafo:**
Caminhos do grafo que não levam a nenhum outro ou ficam redundantes são podados.
**3 Mapeamento de Leituras:**
O SPAdes aplica as reads originais no grafo criado e limpo nas etapas anteriores para reforçar os caminhos certos.
**4 Scaffolding (Andaime) e Resolução de Paired-End:**
Usando o pareamento das reads, o SPAdes une os contigs em scaffolds maiores para passar por regiôes repetitivas.

#### Parâmetros:
**Modos de Operação:**
--isolate: Use isto se estiver sequenciando uma bactéria/organismo padrão cultivado em laboratório. Ele otimiza o algoritmo para genomas com cobertura alta e uniforme, desligando heurísticas feitas para single-cell que poderiam causar erros aqui.

--meta: Modo Metagenômico. Essencial se você tem uma comunidade de várias bactérias na mesma amostra. Ele impede que o SPAdes presuma que a diferença de cobertura entre as sequências seja um erro (já que em metagenomas, bactérias diferentes têm abundâncias diferentes). Nota: Não funciona junto com --careful.

--sc: Single-cell. O modo original. Útil quando você sequenciou a partir de MDA (Multiple Displacement Amplification), que gera uma cobertura altamente enviesada.

--plasmid: Tenta extrair e montar apenas plasmídeos a partir do dado genômico total, separando-os do cromossomo principal com base na topologia do grafo e na cobertura.

**Controle de Qualidade e Correção:**
--careful: Este parâmetro aciona um módulo extra (MismatchCorrector) que roda depois da montagem terminar. Ele usa a ferramenta BWA para mapear as leituras de volta nos contigs montados e corrige pequenos erros (mismatches e pequenos indels - inserções/deleções) que sobraram. É altamente recomendado para pequenos genomas (isolados), mas aumenta significativamente o tempo de processamento.

--only-assembler: Desliga a etapa inicial do SPAdes de correção de erros nas leituras originais (o BayesHammer). Use isso apenas se você já fez uma limpeza ultra-rigorosa (com fastp/Trimmomatic) e confia plenamente nos seus dados.

**Parâmetros Técnicos Básicos:**

-o: Onde os resultados (contigs e scaffolds) serão salvos.

-t: Número de núcleos de processamento. O SPAdes é muito pesado; use o máximo que tiver disponível (ex: -t 16).

-m: Limite de memória RAM em Gigabytes (ex: -m 64). O SPAdes consome muita RAM. Se ele estourar o limite da sua máquina, a montagem será abortada.

-k: Permite escolher manualmente os tamanhos dos K-mers (devem ser números ímpares e separados por vírgula, ex: -k 21,33,55,77).

--cov-cutoff: Após montar o grafo, calcula a cobertura média de cada scaffold. Se a cobertura for menor que o valor X, assume que é lixo (ruído de fundo) e delete esse caminho do resultado final.

---

## quast

---

## samtools

#### Comando
```bash
samtools coverage assemble_mapping.bam > assemble_coverage.txt
```

#### Resultados

**#rname** (Reference Name): O nome do contig

**startpos e endpos**: A posição inicial e final considerada para o cálculo. Por padrão, se você mapeou contra o contig inteiro, startpos será 1 e endpos será o tamanho total do contig em pares de bases (bp).

**numreads**: O número exato de reads que mapearam (alinharam) contra esse contig específico. Quanto maior o contig ou quanto mais expressivo for o DNA daquela região na amostra, maior será esse número.

**covbases** (Covered Bases): O número de posições (pares de bases) do contig que receberam pelo menos uma read mapeada.Exemplo: Se um contig tem 1.000 bp e covbases é 950, significa que existem 50 posições nele que ficaram completamente "descobertas" (sem nenhuma read alinhada).

**coverage** (% de Cobertura Horizontal): A proporção do contig que foi coberta por pelo menos uma read.

**meandepth** (Profundidade Média ou Cobertura Vertical): É o número médio de vezes que cada base daquele contig foi sequenciada/mapeada. Esta é a famosa "Cobertura X" (ex: 35.4x). Valores acima de 30x são excelentes para genomas bacterianos ou fúngicos pequenos, dando alta confiabilidade às bases chamadas.

**meanbaseq** (Mean Base Quality): A qualidade Phred média das bases das reads que mapearam ali. É a qualidade que veio do sequenciador (e passou pelos seus filtros do Fastp/Trimmomatic). Valores >30 indicam excelente qualidade de sequenciamento.

**meanmapq** (Mean Mapping Quality): A qualidade média de alinhamento das reads. É uma escala Phred gerada pelo BWA que mede a confiança de que a read realmente pertence àquele lugar e não foi mapeada ali por erro ou aleatoriedade (geralmente vai de 0 a 60). Se for muito baixo (ex: < 20), significa que as reads que mapearam naquele contig são muito repetitivas e poderiam pertencer a qualquer outra parte do genoma (regiões de repetição ou transposons).


#### Comando
```bash
samtools view -@ 8 -bS alinhamento.sam > alinhamento.bam
```

#### Parâmetros
-b: diz para sair em formato Binário (BAM, menor tamanho). 
-S: avisa que a entrada é SAM. 
-@ 8: usa 8 threads para ir mais rápido.

#### Resultados
Ver e Filtrar Alinhamentos (view)
É o comando principal. Ele converte formatos e filtra o que você quer ver.


#### Comando
```bash
samtools sort -@ 8 alinhamento.bam -o alinhamento_ordenado.bam
```

#### Parâmetros
-o: nome do arquivo de saida.
-@ 8: usa 8 threads para ir mais rápido.

#### Resultados
Ordenar Alinhamentos (sort)
Organiza as reads mapeadas pela posição delas no cromossomo/contig, do início ao fim. Quase todas as outras ferramentas exigem que o BAM esteja ordenado.


#### Comando
```bash
samtools flagstat alinhamento_ordenado.bam
```

#### Resultados
Obter Estatísticas Rápidas (flagstat)
Dá um resumo rápido de quantas reads mapearam, quantas não mapearam e quantas mapearam com o par correto (resultados <20% indicam que ouve erro na montagem).


#### Comando
```bash
samtools fastq alinhamento.bam > reads_recuperadas.fastq
```

#### Resultados
Converter BAM de volta para FASTQ (fastq)
Se você perdeu seus arquivos de reads originais, mas tem o arquivo BAM de um alinhamento antigo, você pode extrair as reads de volta.


#### Comando
```bash
samtools depth alinhamento_ordenado.bam > profundidade.txt
```

#### Resultados
Profundidade Base por Base (depth)
Diz quantas reads cobrem cada nucleotídeo individualmente.
Útil para gerar gráficos detalhados e procurar quebras de montagem (onde a profundidade cai repentinamente para zero).

---

## seqkt

#### Comando
```bash
seqtk seq -L 500 contigs.fasta > clean_contigs.fasta
```

#### Parâmetros
-L 500: ignora e descarte qualquer sequência que tenha menos de 500 pares de bases.

#### Resultados
Arquivo de contigs filtrados.


#### Comando
```bash
seqtk sample -s100 reads.fastq 0.5 > reads_metade.fastq
```

#### Parâmetros
-s: é a semente randômica (essencial usar a mesma para R1 e R2). 
0.5: retem 50% das reads. Também aceita número inteiro, como 10000, para extrair exatamente dez mil reads.

#### Resultados
Downsampling: Reduz o número de reads de forma aleatória para diminuir a cobertura.


#### Comando
```bash
seqtk subseq todas_as_reads.fastq lista_de_nomes.txt > apenas_algumas.fastq
```

#### Resultados
Extrair sequências específicas (por nome)
Se você tem uma lista com o ID de algumas reads (ou contigs) que deseja separar do resto do arquivo.
Muito útil quando o BUSCO ou o QUAST acusa que há um contig contaminante e você quer isolá-lo para investigar.


#### Comando
```bash
seqtk seq -r reads.fastq > reads_reversas.fastq
```

#### Resultados
Reverse Complement
Inverte a sequência de trás para frente e troca as bases pelos seus pares (A por T, C por G).

---

## BUSCO
(Benchmarking Universal Single-Copy Orthologs)

O BUSCO procura genes ortólogos de cópia única dentro dos contigs. Genes esses que todo organismo possue em cópia única.

#### Parâmetros:
busco -i contigs.fasta -o resultado_busco -m genome -l saccharomycetes_odb10 -c 8 -f

**-i**: arquivo de entrada (contigs.fasta)
**-o**: nome da pasta onde serão gerados os resultados
**-m**: tipo de dado a ser analizada (DNA=>genome, RNA=>transcriptome, Proteina=>proteins)
**-l**: banco de dados a ser usado como referencia (saccharomycetes_odb10=>leveduras)
**-c**: numero de threads a ser usado
**-f**: sobrescreve a pasta de destino caso tenha algo lá

#### Resultados:
ex: C:97.5%[S:96.0%,D:1.5%], F:1.0%, M:1.5%, n:2137

n: = número de genes essenciais do banco usado.
C: = porcentagem dos genes essenciais que foram encontrados inteiros na sua montagem.
S: = porcentagem dos genes essenciais que foram encontrados inteiros e apenas uma vez. (cenário perfeito)
D: = porcentagem dos genes essenciais que foram encontrados inteiros e mais uma vez. (valores > 5% podem idicar contaminação ou que o organismo é diploide)
F: = porcentagem dos genes essenciais que foram encontrados fragmentados. (valores muito altos podem refletir genoma muito picotado, o que leva a um n50 baixo também)
M: = porcentagem dos genes essenciais que não foram encontrados, devido a baixa cobertura ou complexidade do genoma.

## blastn
O blastn compara sequências de nucleotídeos (DNA/RNA) contra outras sequências de nucleotídeos. Se ele tentasse alinhar letra por letra de uma sequência gigante contra outra diretamente, o processo demoraria dias. Por isso, ele usa uma heurística baseada em sementes (seeds).

**Semeadura (Seeding):** Ele quebra a sua sequência de busca em pedacinhos minúsculos (words) de tamanho fixo.

**Busca por Alvos:** Ele escaneia a sequência de referência procurando por correspondências exatas dessas palavras.

**Extensão:** Ao achar uma correspondência (um hit), o algoritmo tenta estender o alinhamento para a esquerda e para a direita, dando pontos positivos para letras iguais (matches) e penalidades para letras diferentes (mismatches) ou lacunas (gaps).

**Filtragem:** Alinhamentos estatisticamente fracos são descartados, e os melhores são reportados no seu arquivo de saída.

#### Parâmetros:

**-query:** Define o arquivo FASTA contendo a sequência que você tem em mãos e quer investigar (a "query").

**-subject:** Define o arquivo FASTA contendo a sequência de referência onde a busca será feita (o "banco de dados" ou feature). Nota técnica: Usar -subject é ótimo para comparar um arquivo contra o outro diretamente sem precisar indexar um banco de dados com makeblastdb antes.

**-out:** O nome e caminho do arquivo onde o BLAST vai salvar os resultados gerados.

**-evalue:** O limite do Expect Value (Valor E). Ele é uma métrica estatística que diz quantos alinhamentos com aquela pontuação você esperaria encontrar ao acaso naquele banco de dados. Quanto menor o número (ex: 1e-5, 1e-20 ou 0.0), mais confiável é o alinhamento. Se o BLAST achar um hit com E-value maior do que o seu corte, ele simplesmente omite do resultado. (ex: -evalue 0.001)

**-word_size:** É o tamanho daquela "word" inicial da semeadura. Para o blastn padrão, o valor padrão é 11. Se você reduzir esse número (ex: para 7), o BLAST se torna muito mais sensível, encontrando sequências mais distantes evolutivamente, porém o processo fica mais lento. Se aumentar, ele fica rápido, mas só acha sequências muito idênticas.

**-dust:** Ativa/Desativa o filtro DUST para regiões de baixa complexidade genética (como repetições do tipo AAAAAAA ou ATATATAT). Passar -dust yes mascara essas regiões para evitar falsos positivos biológicos. Passar -dust no desativa o filtro, permitindo alinhar repetições (comum em genomas mitocondriais que têm regiões repetitivas nas bordas).

**-outfmt:** Expecifia o formato tabelado para os resultado e a ordem das colunas que serão geradas. (ex: -outfmt '6 qseqid sseqid qcovs pident evalue score qstart qend sstart send')

#### Resultados:

**qseqid:** Query Sequence ID: O nome da sua sequência query.
**sseqid:** Subject Sequence ID: O nome do contig/referência onde ela grudou.
**qcovs:** Query Coverage per Subject: Qual porcentagem da sua query conseguiu alinhar no alvo (vai de 0 a 100%). Excelente para ver se o gene está inteiro.
**pident:** Percentage of Identical Matches: A identidade. De todas as bases alinhadas, qual a porcentagem de letras idênticas (ex: 98.5%).
**evalue:** Expect value: A significância estatística do hit (quanto menor, melhor).
**score:** Bit Score: A pontuação bruta do alinhamento baseada na matriz de pontos. Não depende do tamanho do banco de dados.
**qstart:** Posição (número da base) onde o alinhamento começa na sua sequência Query.
**qend:** Posição onde o alinhamento termina na sua sequência Query.
**sstart:** Posição onde o alinhamento começa na sequência Subject (contig).
**send:** Posição onde o alinhamento termina na sequência Subject (contig).

## Conda

**Comando**
```bash
which python # Mostra de qual ambiente a ferramenta está sendo usada.
```

```bash
conda info --envs # Lista todos os ambientes.
```

```bash
conda activate python_env # Ativa o ambiente.
```

```bash
conda list biopython # Lista todas as versões do programa instaladas no ambiente atual.
```

```bash
conda install -c conda-forge biopython # Instala um novo programa.
```


