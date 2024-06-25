import matplotlib.pyplot as plt
import hopsy
import cobra
import os
import numpy as np
import helpers
import pandas as pd
import time



if __name__ == "__main__":
    models = [
        'iEZ481_PCA_Gluc',
        'iEZ481_Glucose-MOPS',
        'iEZ481_Citrat-MOPS',
    ]
    b_samples = {}
    for i, model in enumerate(models):
        b_samples[model] = np.load(file=os.path.join('data', f'{model}_boltzmann_samples_thinned.npz'))['samples']

    media = {
        'iEZ481_PCA_Gluc': 'PCA',
        'iEZ481_Glucose-MOPS': 'Glucose + PCA',
        'iEZ481_Citrat-MOPS': 'Citrate + PCA',
    }

    flux_names = {}
    for m in models:
        model_path = os.path.join("models", m + ".xml")
        cobra_model = cobra.io.read_sbml_model(model_path)
        # open up growth rate constraints
        flux_names[m] = []
        for i, rxn in enumerate(cobra_model.reactions):
            flux_names[m].append(rxn.id)

    assert flux_names[models[0]] == flux_names[models[1]]
    assert flux_names[models[1]] == flux_names[models[2]]
    flux_names = flux_names[models[0]]

    fluxes_to_plot = \
        [
            "GLC_t_PEP",
            "pcaGH", "acnA", "icd", "odhA", "sucD", "sdhCAB", "fumC", "mqo", "gltA", "pgi", "pyk", "pdh", "aceB", "aceA", "mdh", "odx", "pyc", "mez", "pckG", "ppc", "pfkA", "fda", "gapA", "pgk", "eno", "zwf", "opcA", "gnd", "pgm", "rpe", "rpi", "tkt_1", "tal", "tkt_2", "tpiA", "pps", "fbp", "pts", "acnB", "gapB", "actA", "pcaGF"]

    flux_samples = {}
    dfs = []
    # TODO GH durch GF in title ersetzen
    for m in models:
        for i in range(b_samples[m].shape[2]):
            if flux_names[i] not in fluxes_to_plot:
                print('skipping', flux_names[i])
                continue
            flux_samples[flux_names[i]] = [np.mean(b_samples[m][:, :, i])]
        flux_samples['medium'] = media[m]

        df = pd.DataFrame.from_dict(flux_samples)
        df.set_index('medium', inplace=True)
        dfs.append(df)
        print(df)

        df = pd.concat(dfs)
        # df.rename(columns={'pcaGF': 'pcaGH'}, inplace=True)
        df.rename(columns={'GLC_t_PEP': 'pts'}, inplace=True)
        print(df.columns)
        df.to_csv("fluxmaps.csv", index=False)
