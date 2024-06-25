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
        'iEZ481_Glucose-MOPS': 'GLC',
        'iEZ481_Citrat-MOPS': 'CIT',
    }

    fba_keys = { 'iEZ481_PCA_Gluc': 'PCA-Gluc_fluxes.csv', 'iEZ481_Glucose-MOPS': 'glucose-MOPS_fluxes.csv', 'iEZ481_Citrat-MOPS': 'citrat-MOPS_fluxes.csv'}

    fba_results = {}
    for m in models:
        fba_results[m] = pd.read_csv(f'data/fba_results/{fba_keys[m]}')[['Unnamed: 0','mean+0stds']]
        fba_results[m].columns = ['flux', 'fba result']
        fba_results[m].set_index('flux', inplace=True)
        fba_results[m] = fba_results[m].T

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

    n_bins = 20
    # fluxes_to_plot = \
        # ["pcaGH", "biomass_a", "acnA", "icd", "odhA", "sucD", "sdhCAB", "fumC", "mqo", "gltA", "pgi", "pyk", "pdh", "aceB", "aceA", "mdh", "odx", "pyc", "mez", "pckG", "ppc", "pfkA", "fda", "gapA", "pgk", "eno", "zwf", "opcA", "gnd", "pgm", "rpe", "rpi", "tkt_1", "tal", "tkt_2", "tpiA", "pps", "fbp", "pts", "acnB", "gapB", "actA", "pcaGF"]
    fluxes_to_plot = [
        # "acnB",
        # "gnd",
        # "icd",
        # "ldh",
        # "mdh",
        # "pgk",
        'tkt_2', 'pgi', 'sdhCAB', 'rpi', 'fda', 'sucD', 'gltA', 'pyc', 'pgk', 'aceA',
    ]
    # print('plotting', len( fluxes_to_plot ), 'fluxes')

    n_bins = 20

    for i, flux_name in enumerate(fluxes_to_plot):
        time.sleep(0.5)
        flux_index = flux_names.index(flux_name)
        maximum_flux = np.max([np.max(b_samples[m][:,:,flux_index].flatten()) for m in models])
        minimum_flux = np.min([np.min(b_samples[m][:,:,flux_index].flatten()) for m in models])
        minimum_flux = np.min([0, minimum_flux])
        maximum_flux = np.max([maximum_flux] + [fba_results[m][flux_name].values.flatten()[0] for m in models])
        flux_range = (maximum_flux - minimum_flux)
        maximum_flux += 0.01 * flux_range
        minimum_flux -= 0.01 * flux_range

        # densly dashed, loosely dashed, losely dotted
        style = (0, (5, 1))
        for j, m in enumerate(models):
            plt.figure(figsize=(2.5, 2.5))
            plt.title(flux_name)
            # print('range', m, np.max(b_samples[m][:, :, i].flatten())), np.min(b_samples[m][:, :, i].flatten())
            fba = fba_results[m][flux_name].values.flatten()[0]
            print(m, flux_name, fba)

            _, bins, _ = plt.hist(
                b_samples[m][:, :, flux_index].flatten(),
                bins=n_bins, # if bins is None else bins,
                alpha=0.25,
                density=True,
                # label=f'MaxEnt {media[m]}',
                color=f'C{j}',
            )
            _, bins, _ = plt.hist(
                b_samples[m][:, :, flux_index].flatten(),
                bins=n_bins,  # if bins is None else bins,
                histtype='step',
                linewidth=3,
                density=True,
                color=f'C{j}',
            )

            plt.axvline(x=fba_results[m][flux_name].values.flatten(), color=f'k', linewidth=3, linestyle=style,
                        label=f'FBA'
                        )

            plt.xlabel('flux value [mmol/gwd/h]')
            plt.ylabel('density')

            # leg = plt.legend(bbox_to_anchor=(1.0, 1))
            plt.tight_layout()
            plt.title(flux_name)
            plt.savefig(f'flux_map_components_svg/{m}_{flux_name}.svg', bbox_inches='tight',
                        # bbox_extra_artists=(leg,))
                        )
            plt.show()

