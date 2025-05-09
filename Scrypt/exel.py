import pandas as pd
from sequana import DNA
from sequana import FastA
import pandas as pd
import matplotlib.pyplot as plt


f = FastA("../data/Fasta/TriTrypDB-68_LmajorFriedlin_Genome.fasta")
df = pd.read_csv("../data/Centromere_Positions/centromeres.csv")
centromeres = {
    row['Chromosome']: (row['Centromere_Start'], row['Centromere_End'])
    for _, row in df.iterrows()
}

results = []

for maseq in f:

    df_karlin = pd.read_csv(f"../data/karlin/{maseq.name}_karlin_100.csv")
    karlin = df_karlin['Karlin'].tolist()
    karlin = karlin[5000:-5000]

    mean_karlin = sum(karlin) / len(karlin)
    min_karlin = min(karlin)
    max_karlin = max(karlin)

    result = {
        'Chromosome': maseq.name,
        'Mean_Karlin': mean_karlin,
        'Min_Karlin': min_karlin,
        'Max_Karlin': max_karlin,
        'Mean_Centromere': None,
        'Min_Centromere': None,
        'Max_Centromere': None
    }

    if maseq.name in centromeres:
        cent_start, cent_end = centromeres[maseq.name]
        offset = 100 // 2
        cent_start = cent_start-5000
        cent_end = cent_end -5000
        start_idx = max(0, cent_start - offset) 
        end_idx = max(0, cent_end - offset)

        karlin_centromere = karlin[start_idx:end_idx]
        if len(karlin_centromere) > 0:
            result['Mean_Centromere'] = sum(karlin_centromere) / len(karlin_centromere)
            result['Min_Centromere'] = min(karlin_centromere)
            result['Max_Centromere'] = max(karlin_centromere)

    results.append(result)

# Création du DataFrame et export vers Excel
df_results = pd.DataFrame(results)
df_results = df_results.round(2) 
df_results.to_excel("../output/karlin_stats.xlsx", index=False)
