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

    n_cols = 5
    n_rows = 2
    plt.figure(figsize=(n_cols * 2.5, n_rows * 2.5))
    plt.suptitle("Selection of flux distributions")
    n_bins = 20

    for i, flux_name in enumerate(fluxes_to_plot):
        plt.subplot(n_rows, n_cols, 1+i)
        flux_index = flux_names.index(flux_name)
        maximum_flux = np.max([np.max(b_samples[m][:,:,flux_index].flatten()) for m in models])
        minimum_flux = np.min([np.min(b_samples[m][:,:,flux_index].flatten()) for m in models])
        minimum_flux = np.min([0, minimum_flux])
        maximum_flux = np.max([maximum_flux] + [fba_results[m][flux_name].values.flatten()[0] for m in models])
        flux_range = (maximum_flux - minimum_flux)
        maximum_flux += 0.01 * flux_range
        minimum_flux -= 0.01 * flux_range
        # print('flux', flux_name, 'has range', maximum_flux, minimum_flux)
        # print(flux_name)
        if flux_name == "biomass_a":
             plt.title(r"growth rate $\lambda$")
        else:
            plt.title(flux_name)

        plt.xlim(minimum_flux, maximum_flux)
        for j, m in enumerate(models):
            # print('range', m, np.max(b_samples[m][:, :, i].flatten())), np.min(b_samples[m][:, :, i].flatten())
            fba = fba_results[m][flux_name].values.flatten()[0]
            print(m, flux_name, fba)

            _, bins, _ = plt.hist(
                b_samples[m][:, :, flux_index].flatten(),
                bins=n_bins, # if bins is None else bins,
                alpha=0.25,
                density=True,
                label=f'MaxEnt {media[m]}' if i==0 else None,
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

        # densly dashed, loosely dashed, losely dotted
        styles = [(0, (5, 1)), (0, (5, 3)), (0, (5, 6))]
        for j, m in enumerate(models):
            plt.axvline(x=fba_results[m][flux_name].values.flatten(), color=f'C{j}', linewidth=3-0.5*j, linestyle=styles[j], label=f'FBA for {media[m]}' if i==0 else None)

        if flux_name == "biomass_a":
            plt.xlabel(r'growth rate [1/h]')
        else:
            plt.xlabel('flux value [mmol/gwd/h]')

        plt.ylabel('density')
        # plt.xscale('symlog')
    leg = plt.figlegend(bbox_to_anchor=(1.2, 1))
    plt.tight_layout()
    plt.savefig(f'flux_images_pdf/selection2.pdf', bbox_inches='tight', bbox_extra_artists=(leg,))
    plt.savefig(f'flux_images_svg/selection2.svg', bbox_inches='tight', bbox_extra_artists=(leg,))
    plt.savefig(f'flux_images_pdf/selection2.png', bbox_inches='tight', bbox_extra_artists=(leg,))

    plt.show()
