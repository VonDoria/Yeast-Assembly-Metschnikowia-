# Metschnikowia
---
## Assembly

O montagem foi realisada seguindo os seguintes passos:
Trimagem das reads originais usando o fastp e o trimmomatic;
Montagem dos genomas usando as reads originais, trimadas pelo fastp e trimadas pelo trimmomatic, usando o SPAdes. Essas montagens foram feitas com os parâmetros '--careful' e o '--isolate'.
Os contigs gerados pelo SPAdes foram filtrados usando o seqtk, de modo a se remover todos os contigs menores que 500 pb.
A qualidade dos contigs gerados foram avaliados usando o samtools coverage e o quast.
A completude dos genomas foi avaliada usando o BUSCO com o banco de linhagens 'saccharomycetes_odb12'.
Ao final do processo, foi repetido o processo de montagem para as reads trimadas pelo trimmomatic, usando os parâmetros '--cov-cutoff 80' e '--cov-cutoff auto'.

*Devido ao elevado tamanho de muitos dos arquivos de entrada e saida das ferramentas, apenas os arquivos expecificos dos resultados foram mapeados para subir para o repositório. Consulte o arquivo assembly_pipeline.sh para reproduzir os comandos usados e recriar todos os arquivos que foram criados nesse trabalho.
--------------------------------------------------------------------------------------------------------------------















