import matplotlib.pyplot as plt
import hopsy
import os
import numpy as np
import helpers
import pandas as pd

class MarginalGauss:
    def __init__(self, mu, sigma, index):
        self.mu = mu
        self.index = index
        self.sigma = sigma

    def log_density(self, x):
        return -0.5 * (x[self.index] - self.mu) ** 2 / (self.sigma)**2


if __name__ == "__main__":
    optimal_growth_rates = {}
    mean_growth_rates = {}
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
    optimal_growth_rates['iEZ481_PCA_Gluc'] = float(pca_gluc_mu['mu'].max())

    mean_growth_rates['iEZ481_Glucose-MOPS'] = float(gluc_mu['mu'].mean())
    std_growth_rates['iEZ481_Glucose-MOPS'] = float(gluc_mu['mu'].std())
    optimal_growth_rates['iEZ481_Glucose-MOPS'] = float(gluc_mu['mu'].max())

    mean_growth_rates['iEZ481_Citrat-MOPS'] = float(citr_mu['mu'].mean())
    std_growth_rates['iEZ481_Citrat-MOPS'] = float(citr_mu['mu'].std())
    optimal_growth_rates['iEZ481_Citrat-MOPS'] = float(citr_mu['mu'].max())

    models = [
        'iEZ481_Citrat-MOPS',
        'iEZ481_Glucose-MOPS',
        'iEZ481_PCA_Gluc',
    ]
    for model in models:
        polytope = helpers.load_polytope('data', model)
        biomass_index = 300
        m = MarginalGauss(mean_growth_rates[model], std_growth_rates[model], biomass_index)
        problem = hopsy.Problem(polytope.A, polytope.b, model=m)
        problem = hopsy.add_equality_constraints(problem, A_eq=polytope.S, b_eq=polytope.h)
        problem = hopsy.round(problem)
        starting_point = hopsy.compute_chebyshev_center(problem)
        # Gleichverteilte Samples
        n_samples = 65_000
        n_procs = 4
        chains = [
            # UniformCoordinate should work well because there should not really be parameter correlations
            hopsy.MarkovChain(problem, proposal=hopsy.UniformCoordinateHitAndRunProposal, starting_point=starting_point) for
            i in range(n_procs)]
        rng = [hopsy.RandomNumberGenerator(seed=i + 1123) for i in range(n_procs)]
        print(f'start sampling {model}')
        _, samples = hopsy.sample(chains, rng, n_samples=n_samples, thinning=20, n_procs=n_procs, progress_bar=True)
        print(f'finished sampling {model}')
        print('thinned ESS ', np.min(hopsy.ess(samples[:,::1,:])))
        print('thinned rhat ', np.max(hopsy.rhat(samples[:,::1,:])))
        samples = samples[:, 1000:, :]
        print(samples.shape)
        np.savez_compressed(file=os.path.join('data', f'{model}_gauss_samples.npz'), samples=samples)
        # # Histogramm anzeigen lassen
        plt.figure()
        plt.title(f'{model} growth rate')
        plt.hist(samples[0, :, 300], density=True, bins=100, alpha=0.5)
        plt.hist(samples[1, :, 300], density=True, bins=100, alpha=0.5)
        plt.hist(samples[2, :, 300], density=True, bins=100, alpha=0.5)
        plt.hist(samples[3, :, 300], density=True, bins=100, alpha=0.5)
        plt.tight_layout()
        plt.savefig(f"gauss_{model}.png")
        plt.show(block=False)
        del samples
