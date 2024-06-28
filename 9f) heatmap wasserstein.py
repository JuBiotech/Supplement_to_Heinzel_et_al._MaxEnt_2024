import seaborn as sns
import cobra
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np




if __name__ == "__main__":
    reactions_to_plot = set()
    model_path = os.path.join("models/iEZ481_PCA_Gluc.xml")
    cobra_model = cobra.io.read_sbml_model(model_path)
    reaction_names = [rxn.id for rxn in cobra_model.reactions]

    pathways = {
        "PPP": [ "zwf", "opcA", "gnd", "rpe", "rpi", "tkt_1", "tal", "tkt_2"],
        "EMP": ["pgi", "pfkA", "fbp", "fda", "tpiA", "gapA", "gapB", "eno", "pgk", "pgm", "pyk", "pps", "pdh", ],
        "ANA": [ "odx", "pyc", "mez", "ppc", "pckG"],
        "TCA": [ "gltA", "acnA", "acnB", "icd", "odhA", "sucD", "sdhCAB", "fumC", "mqo", "mdh", ],
        "GLX": ["aceB", "aceA"],
    }

    selection = ['tkt_2', 'pgi', 'sdhCAB', 'rpi', 'fda', 'sucD', 'gltA', 'pyc', 'pgk', 'aceA']

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


    PCA_GLC = pd.read_csv(f'PCA_Glucose_Boltzmann.txt', names=['wasserstein'])
    PCA_CIT = pd.read_csv(f'Citrat_PCA_Boltzmann.txt', names=['wasserstein'])
    CIT_GLC = pd.read_csv(f'Glucose_Citrat_Boltzmann.txt', names=['wasserstein'])
    ind_list = [reaction_names.index(c) for c in mapped_fluxes]
    PCA_GLC = PCA_GLC.iloc[ind_list]  # .set_index(central_reactions)
    PCA_CIT = PCA_CIT.iloc[ind_list]  # .set_index(central_reactions)
    CIT_GLC = CIT_GLC.iloc[ind_list]  # .set_index(central_reactions)
    # for offset, flux in enumerate(mapped_fluxes):
        # for i, flux in enumerate(fluxes):
        #     print(PCA_GLC.iloc[i])
    pca_glc_bars = PCA_GLC.values.flatten()
    pca_cit_bars = PCA_CIT.values.flatten()
    cit_glc_bars = CIT_GLC.values.flatten()

    x_fluxes = np.arange(len(mapped_fluxes))
    plt.figure(figsize=(6.4, 3))
    plt.bar(x_fluxes - 0.2, pca_glc_bars, width=0.2, align='center', label=r'PCA $\leftrightarrow$ GLC')
    plt.bar(x_fluxes, cit_glc_bars, width=0.2, align='center', label=r'GLC $\leftrightarrow$ CIT')
    plt.bar(x_fluxes + 0.2, pca_cit_bars, width=0.2, align='center', label=r'CIT $\leftrightarrow$ PCA')
    plt.xticks(x_fluxes, mapped_fluxes, rotation=90)
    plt.ylabel('Wasserstein-1-Distance')
    plt.xlabel('Central fluxes')
    ticks, labels = plt.xticks()
    for i in range(len(ticks)):
        tick = ticks[i]
        label = labels[i]
        print(tick)
        if mapped_fluxes[tick] in selection:
            label.set_color('r')
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'wasserstein-bars.svg', bbox_inches='tight')
    plt.savefig(f'wasserstein-bars.pdf', bbox_inches='tight')
    plt.show()