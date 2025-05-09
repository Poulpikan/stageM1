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

    s = DNA(maseq.sequence.upper())
    s.window = 100
    s._compute_skews()
    plt.figure(figsize=(12, 4))
    positions = range(len(s._Yn))
    plt.plot(positions, s._Yn, label=f'Karlin - {maseq.name}', color='teal')
    if maseq.name in centromeres:
        cent_start, cent_end = centromeres[maseq.name]
        plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
        plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    plt.title(f'_Yn {maseq.name}')
    plt.xlabel('Position in sequence')
    plt.ylabel('_Yn')
    plt.ylim(0, 10)
    plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"../output/dataExploration/Y/{maseq.name}.png")
    plt.close()
    #df_karlin = pd.read_csv(f"../data/karlin/{maseq.name}_karlin_100.csv")
    #karlin = df_karlin['Karlin'].tolist()
    # karlin = s.get_karlin_signature_difference()
    # positions = range(len(karlin))

    # mean_karlin = sum(karlin) / len(karlin)
    # min_karlin = min(karlin)
    # max_karlin = max(karlin)
    # print(f"{maseq.name} moy = {mean_karlin:.3f}  min = {min_karlin:.3f}  max = {max_karlin:.3f}")

    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     # Ajustement des bornes : on suppose que karlin commence à la position 0 + window // 2
    #     # donc index 0 de karlin correspond à position genome = window // 2
    #     offset = s.window // 2
    #     start_idx = max(0, cent_start - offset)
    #     end_idx = max(0, cent_end - offset)

    #     # Extraction des valeurs dans la zone centromérique
    #     karlin_centromere = karlin[start_idx:end_idx]

    #     # Calcul des stats
    #     min_k = min(karlin_centromere)
    #     max_k = max(karlin_centromere)
    #     mean_k = sum(karlin_centromere) / len(karlin_centromere)

    #     print(f"Centromere {maseq.name}: moy = {mean_k:.3f}, min = {min_k:.3f}, max = {max_k:.3f}")

    #    # === Entropy Plot - Zoom centromere ===
    #     if maseq.name in centromeres:
    #         cent_start, cent_end = centromeres[maseq.name]
    #         zoom_margin = 10000  # tu peux ajuster cette marge

    #         # Filtrer les données dans cette zone
    #         zoom_positions = []
    #         zoom_entropies = []
    #         for pos, entropy in zip(positions, karlin):
    #             if cent_start - zoom_margin <= pos <= cent_end + zoom_margin:
    #                 zoom_positions.append(pos)
    #                 zoom_entropies.append(entropy)

    #         plt.figure(figsize=(12, 4))
    #         plt.plot(zoom_positions, zoom_entropies, label='Entropy (zoom)', color='darkgreen')
    #         plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #         plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    #         plt.title(f'Entropy plot (zoom) for {maseq.name}')
    #         plt.xlabel('Position in sequence')
    #         plt.ylabel('Shannon Entropy')
    #         plt.grid(True)
    #         plt.legend()
    #         plt.tight_layout()
    #         plt.savefig(f"../output/dataExploration/karlin_zoom/{maseq.name}.png")
    #         plt.close()
    # # Tracé
    # plt.figure(figsize=(12, 4))
    # plt.plot(positions, karlin, label=f'Karlin - {maseq.name}', color='teal')
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    # plt.title(f'Karlin {maseq.name}')
    # plt.xlabel('Position in sequence')
    # plt.ylabel('Karlin')
    # plt.ylim(0, 10)
    # plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f"../output/dataExploration/karlin/{maseq.name}.png")
    # plt.close()


    # dnaFlex = s.get_dna_flexibility()
    # # Tracé
    # plt.figure(figsize=(12, 4))
    # plt.plot(positions, dnaFlex, label=f'DNA Flex - {maseq.name}', color='teal')
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    # plt.title(f'DNA Flex {maseq.name}')
    # plt.xlabel('Position in sequence')
    # plt.ylabel('DNA Flex')
    # plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f"../output/dataExploration/dna_Flex/{maseq.name}.png")
    # plt.close()

    # information_entropy = s.get_informational_entropy()
    # # Tracé
    # plt.figure(figsize=(12, 4))
    # plt.plot(positions, information_entropy, label=f'informational entropy- {maseq.name}', color='teal')
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    # plt.title(f'informational entropy {maseq.name}')
    # plt.xlabel('Position in sequence')
    # plt.ylabel('informational entropy')
    # plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f"../output/dataExploration/information_entropy/{maseq.name}.png")
    # plt.close()
