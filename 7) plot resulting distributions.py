import matplotlib.pyplot as plt
import hopsy
import os
import numpy as np
import helpers
import pandas as pd


class Boltzmann:
    def __init__(self, beta, index):
        self.beta = beta
        self.index = index

    def log_density(self, x):
        return self.beta * x[self.index]


if __name__ == "__main__":
    mean_growth_rates = {}
    maxent_growth_rates = {}
    std_growth_rates = {}
    all_mu = pd.read_csv('data/extracted_growth_rates_with_media.csv', index_col=0)
    # select area method (we could have used others, but area works well)"
    method = "area"
    all_mu = all_mu[all_mu['method'].str.contains("area")]
    pca_gluc_mu = all_mu[all_mu['medium'] == "PCA-Gluc"].reset_index()
    gluc_mu = all_mu[all_mu['medium'] == "Glucose-MOPS"].reset_index()
    citr_mu = all_mu[all_mu['medium'] == "Citrat-MOPS"].reset_index()

    mean_growth_rates['iEZ481_PCA_Gluc'] = float(pca_gluc_mu['mu'].mean())
    std_growth_rates['iEZ481_PCA_Gluc'] = float(pca_gluc_mu['mu'].std())

    mean_growth_rates['iEZ481_Glucose-MOPS'] = float(gluc_mu['mu'].mean())
    std_growth_rates['iEZ481_Glucose-MOPS'] = float(gluc_mu['mu'].std())

    mean_growth_rates['iEZ481_Citrat-MOPS'] = float(citr_mu['mu'].mean())
    std_growth_rates['iEZ481_Citrat-MOPS'] = float(citr_mu['mu'].std())

    models = [
        'iEZ481_PCA_Gluc',
        'iEZ481_Glucose-MOPS',
        'iEZ481_Citrat-MOPS',
    ]

    media = {
        'iEZ481_PCA_Gluc': 'PCA',
        'iEZ481_Glucose-MOPS': 'Glucose + PCA',
        'iEZ481_Citrat-MOPS': 'Citrate + PCA',
    }

    n_cols = 3
    n_rows = 1
    plt.figure(figsize=(n_cols * 3.5, n_rows * 3.5))
    plt.subplot(n_rows, n_cols, 1)
    n_bins = 20
    for i, model in enumerate(models):
        poltope = helpers.load_polytope('data', model)
        biomass_index = 300
        beta = float(pd.read_csv(f'data/{model}_beta.csv')['best_beta'][0][1:-1])
        u_samples = np.load(file=os.path.join('data', f'{model}_uniform_samples_thinned.npz'))['samples'][:, :, :]
        b_samples = np.load(file=os.path.join('data', f'{model}_boltzmann_samples_thinned.npz'))['samples'][:, :, :]
        print(f'beta is {beta}')
        print('u samples shape', u_samples.shape)
        print('b samples shape', b_samples.shape)
        maxent_growth_rates[model] = np.mean(b_samples[:, :, biomass_index])

        beta_string = "%.3f" % beta
        lambda_bar_string = "%.3f" % mean_growth_rates[model]
        mean_string = "%.3f" % maxent_growth_rates[model]
        mean_uniform = "%.3f" % np.mean(u_samples[:, :, biomass_index])
        plt.subplot(n_rows, n_cols, i + 1)
        plt.title(f'Growth medium {media[model]}')
        _, bins, _ = plt.hist(
            b_samples[:, :, biomass_index].flatten(),
            bins=n_bins,
            alpha=0.25,
            density=True,
            # label=fr"Boltzmann ($\beta$={beta_string})" if i == 2 else None,
            label=fr"Maxent (Boltzmann)" if i == 2 else None,
        )
        _ = plt.hist(
            u_samples[:, :, biomass_index].flatten(),
            bins=bins,
            alpha=0.25,
            density=True,
            label="uniform" if i == 2 else None,
        )
        _ = plt.hist(
            b_samples[:, :, biomass_index].flatten(),
            histtype='step',
            bins=bins,
            density=True,
            linewidth=2,
            color='C0'
        )
        _ = plt.hist(
            u_samples[:, :, biomass_index].flatten(),
            histtype='step',
            bins=bins,
            density=True,
            linewidth=2,
            color='C1'
        )

        f_length = lambda f: f()[1] - f()[0]
        f_offset = lambda f: f()[0]
        x_scale = 0.55
        y_scale = 0.69
        x = x_scale * f_length(plt.xlim) + f_offset(plt.xlim)
        y = y_scale * f_length(plt.ylim) + f_offset(plt.ylim)
        plt.text(x, y, rf'$\beta=${beta_string}', horizontalalignment='left', c='k')
        plt.text(x, y * .85, r'$\bar\lambda_{measured}$'f'={lambda_bar_string}', horizontalalignment='left', c='C2')
        plt.text(x, y * 0.7, r'$\bar\lambda_{maxent}$'f'={mean_string}', horizontalalignment='left', c='C0')
        plt.text(x, y * 0.55, r'$\bar\lambda_{uniform}$'f'={mean_uniform}', horizontalalignment='left', c='C1')


        styles = [(0, (5, 1)), (0, (5, 10)), (0, (5, 15))]
        plt.axvline(
            x=maxent_growth_rates[model],
            # label=r"Boltzmann mean $\bar\lambda_{maxent}$",
            color='C0',
            linewidth=3,
            linestyle=styles[0]
        )
        plt.axvline(
            x=mean_growth_rates[model],
            # label=r"measured growth rate $\bar\lambda_{exp}$",
            linewidth=2,
            color='C2',
            linestyle=styles[1]
        )
        plt.axvline(
            float(mean_uniform),
            # label=r"measured growth rate $\bar\lambda_{exp}$",
            linewidth=2,
            color='C1',
            linestyle=styles[2]
        )

        if i == 2:
            plt.legend(loc='upper right', bbox_to_anchor=(1.25, 1))

        plt.xlabel(r'growth rate $\lambda$ [1/h]')
        plt.ylabel('density')

    plt.tight_layout()
    plt.savefig('maxentResults.pdf')
    plt.savefig('maxentResults.svg')
    plt.show()
