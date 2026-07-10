#%%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_table('../results/fastani.out', header=None, sep=r'\s+')

df.columns = ['Query', 'Reference', 'ANI', 'Mapped_Fragments', 'Total_Fragments']

df[['Query', 'Reference']] = df[['Query', 'Reference']].apply(lambda x: x.str.split('/').str[-2])

df_matriz = df.pivot(index='Query', columns='Reference', values='ANI')

df_matriz = df_matriz.fillna(70)

sns.set_theme(style="white")

g = sns.clustermap(
    df_matriz, 
    cmap='viridis',          # Paleta viridis (ótima para transições de dados contínuos)
    metric='euclidean',      # Métrica de distância para o agrupamento
    method='ward',           # Método de ligação para agrupar amostras similares
    vmin=70, vmax=100,       # Ajusta o contraste (foca na variação biológica relevante)
    figsize=(24, 22),        # Tamanho ideal para acomodar os nomes e o dendrograma
    annot=df_matriz.shape[0] <= 15, # Adiciona os números dentro dos quadrados apenas se houver poucas amostras
    fmt='.1f',               # Formato com uma casa decimal caso os números apareçam
    cbar_kws={'label': 'ANI (%)', 'orientation': 'vertical'}, # Customiza a barra de legenda
    linewidths=.5            # Linha fina separando os blocos para dar um aspecto mais limpo
)

plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=10)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)

g.figure.suptitle('Matriz de Identidade Nucleotídica Média (FastANI)', y=1.02, fontsize=14, weight='bold')

plt.show()
