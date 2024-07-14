import seaborn as sns
import cobra
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np




if __name__ == "__main__":
    pathways = {
        "PPP": [ "zwf", "opcA", "gnd", "rpe", "rpi", "tkt_1", "tal", "tkt_2"],
        "EMP": ["pgi", "pfkA", "fbp", "fda", "tpiA", "gapA", "gapB", "eno", "pgk", "pgm", "pyk", "pps", "pdh", ],
        "ANA": [ "odx", "pyc", "mez", "ppc", "pckG"],
        "TCA": [ "gltA", "acnA", "acnB", "icd", "odhA", "sucD", "sdhCAB", "fumC", "mqo", "mdh", ],
        "GLX": ["aceB", "aceA"],
    }

    selection = set()

    labels = []
    num_fluxes = 0
    mapped_fluxes = []
    for pathway, fluxes in pathways.items():
        num_fluxes += len(fluxes)
        mapped_fluxes += fluxes
    print('num fluxes', num_fluxes)
    print(mapped_fluxes)

    order = ['PCA', 'GLC', 'CIT']
    data_set = np.zeros((num_fluxes*3, num_fluxes*3))


    PCA_GLC = pd.read_csv(f'data/wasserstein_pca-glc.csv', names=['wasserstein'], header=None)
    PCA_CIT = pd.read_csv(f'data/wasserstein_cit-pca.csv', names=['wasserstein'], header=None)
    CIT_GLC = pd.read_csv(f'data/wasserstein_glc_cit.csv', names=['wasserstein'], header=None)

    PCA_GLC['flux'] = mapped_fluxes
    PCA_CIT['flux'] = mapped_fluxes
    CIT_GLC['flux'] = mapped_fluxes
    PCA_GLC = PCA_GLC.set_index('flux')
    PCA_CIT = PCA_CIT.set_index('flux')
    CIT_GLC = CIT_GLC.set_index('flux')

    for pathway, fluxes in pathways.items():
        print('pathway', pathway)
        indices = [mapped_fluxes.index(f) for f in fluxes]
        print('PCA_GLC', PCA_GLC['wasserstein'].iloc[indices].idxmax(), PCA_GLC['wasserstein'].max())
        print('CIT_GLC', CIT_GLC['wasserstein'].iloc[indices].idxmax(), CIT_GLC['wasserstein'].max())
        print('PCA_CIT', PCA_CIT['wasserstein'].iloc[indices].idxmax(), PCA_CIT['wasserstein'].max())
        selection.add(PCA_GLC['wasserstein'].iloc[indices].idxmax())
        selection.add(CIT_GLC['wasserstein'].iloc[indices].idxmax())
        selection.add(PCA_CIT['wasserstein'].iloc[indices].idxmax())

    selection = list(selection)
    print(selection)

    pca_glc_bars = PCA_GLC.values.flatten()
    pca_cit_bars = PCA_CIT.values.flatten()
    cit_glc_bars = CIT_GLC.values.flatten()

    mapped_fluxes.reverse()
    splits = []
    remaining_fluxes = len(mapped_fluxes)
    for p, fs in pathways.items():
        if p == 'GLX':
            continue
        remaining_fluxes = remaining_fluxes - len(fs)
        splits.append(remaining_fluxes-0.5)
        print(p, len(fs), splits[-1])
    x_fluxes = np.arange(len(mapped_fluxes))
    plt.figure(figsize=(3, 6))
    plt.barh(x_fluxes - 0.2, pca_glc_bars[::-1], height=0.2, align='center', label=r'PCA vs. GLC')
    plt.barh(x_fluxes, cit_glc_bars[::-1], height=0.2, align='center', label=r'GLC vs. CIT')
    plt.barh(x_fluxes + 0.2, pca_cit_bars[::-1], height=0.2, align='center', label=r'CIT vs. PCA')
    for split in splits:
        plt.axhline(split, color='k', linestyle='dotted')

    ys = [
        36,
        26,
        13,
        10,
        0,
    ]
    for i, pathway in enumerate(pathways.keys()):
        plt.text(40, ys[i], pathway)
    plt.xlabel('Wasserstein-1-Distance')
    plt.ylabel('fluxes')
    plt.yticks(x_fluxes, mapped_fluxes)
    ticks, labels = plt.yticks()
    for i in range(len(ticks)):
        tick = ticks[i]
        label = labels[i]
        if mapped_fluxes[tick] in selection:
            label.set_color('r')
    plt.legend(framealpha=1, bbox_to_anchor=(0.35, 0.2))
    plt.tight_layout()
    plt.savefig(f'wasserstein-bars.svg', transparent=True, bbox_inches='tight')
    plt.savefig(f'wasserstein-bars.pdf', transparent=True, bbox_inches='tight')
    plt.show()