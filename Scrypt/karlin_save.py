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


#read = next(f)

for maseq in f:
    print(maseq.name)

    s = DNA(maseq.sequence)
    s.window = 100
    karlin = s.get_karlin_signature_difference()

    df_karlin = pd.DataFrame({'Karlin': karlin})
    df_karlin.to_csv(f"../data/karlin/{maseq.name}_karlin_100.csv", index=False)

