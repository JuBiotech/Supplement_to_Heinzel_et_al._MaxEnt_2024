import matplotlib.pyplot as plt
import os
import numpy as np
import hopsy
from scipy.optimize import fsolve, fmin
from scipy.special import logsumexp
import pandas as pd
import helpers


class Boltzmann:
    def __init__(self, beta, index):
        self.beta = beta
        self.index = index

    def log_density(self, x):
        return self.beta * x[self.index]



def sample(weight, polytope, n_samples, rounds=5):
    biomass_index = 300
    m = Boltzmann(weight, biomass_index)
    problem = hopsy.Problem(polytope.A, polytope.b, model=m)
    problem = hopsy.add_equality_constraints(problem, A_eq=polytope.S, b_eq=polytope.h)
    problem = hopsy.round(problem)
    n_procs = 4
    starting_point = hopsy.compute_chebyshev_center(problem)
    chains = [
        # UniformCoordinate should work well because there should not really be parameter correlations
        hopsy.MarkovChain(problem, proposal=hopsy.GaussianHitAndRunProposal, starting_point=starting_point) for
        i in range(n_procs)]
    for i in range(n_procs):
        chains[i].proposal.stepsize = 10
    rng = [hopsy.RandomNumberGenerator(seed=i + 1123) for i in range(n_procs)]
    print(f'start sampling beta={weight}')
    samples = []
    for round in range(rounds):
        print('round', round)
        acc, s = hopsy.sample(chains, rng, n_samples=n_samples, thinning=5, n_procs=n_procs, progress_bar=False)
        print('acceptance rates', acc)
        # burn in
        s = s[:, int(n_samples / 2):, :]
        print(f'finished sampling beta={weight}')
        print('ess was', np.min(hopsy.ess(s[:, :, biomass_index:biomass_index + 1])))
        s = s[:, :, biomass_index]
        samples += list(s.flatten())
    print('samples mean is ', np.mean(samples))
    return samples


def create_importance_sampling_estimator(importance_beta):
    _importance_beta = importance_beta

    def estimator(_samples, _beta):
        shifted = (_beta - _importance_beta) * _samples + np.log(_samples)
        log_res1 = logsumexp(shifted)
        log_Z = logsumexp((_beta - _importance_beta) * _samples)
        log_res = log_res1 - log_Z
        return log_res

    return estimator


if __name__ == "__main__":
    biomass_index = 300
    models = [
        'iEZ481_Glucose-MOPS',
        'iEZ481_PCA_Gluc',
        'iEZ481_Citrat-MOPS',
    ]

    for model in models:
        print(f'finding beta for {model}')
        mean_growth_rates = {}
        all_mu = pd.read_csv('data/extracted_growth_rates_with_media.csv', index_col=0)
        # select area method (we could have used others, but area works well)"
        method = "area"
        all_mu = all_mu[all_mu['method'].str.contains("area")]
        pca_gluc_mu = all_mu[all_mu['medium'] == "PCA-Gluc"].reset_index()
        gluc_mu = all_mu[all_mu['medium'] == "Glucose-MOPS"].reset_index()
        citr_mu = all_mu[all_mu['medium'] == "Citrat-MOPS"].reset_index()
        mean_growth_rates['iEZ481_PCA_Gluc'] = float(pca_gluc_mu['mu'].mean())
        mean_growth_rates['iEZ481_Glucose-MOPS'] = float(gluc_mu['mu'].mean())
        mean_growth_rates['iEZ481_Citrat-MOPS'] = float(citr_mu['mu'].mean())
        bar_lambda = mean_growth_rates[model]
        print(f'\t{model} mean growth rate', bar_lambda)

        best_error = 100
        importance_beta = float(pd.read_csv(f'data/{model}_beta.csv')['best_beta'][0][1:-1])
        polytope = helpers.load_polytope('data', model)
        while np.abs(best_error) > 0.0001:
            samples = sample(importance_beta, n_samples=100000, polytope=polytope)

            lambda_estimator = create_importance_sampling_estimator(importance_beta)
            quad_integrand = lambda x: lambda_estimator(samples, x)

            eq_to_solve = lambda x: -np.log(bar_lambda) + quad_integrand(x)

            initial_guesses = [0, 1, 2, 3]
            best_beta = None
            best_error = None
            best_initial_guess = None
            for ig in initial_guesses:
                initial_guess = np.array([ig])
                result = fsolve(eq_to_solve, initial_guess, full_output=True, xtol=1e-18, maxfev=int(1e5))

                print(f'\t{model} with initial {initial_guess} returns solve info {result}')
                # Das Ergebnis enthält den gefundenen Wert für beta
                beta = np.array(result[0])
                error = np.abs(eq_to_solve(beta))
                print(f"\t\t found beta={beta} with error {error}")
                if best_error is None:
                    best_beta = beta
                    best_error = error
                    best_initial_guess = initial_guess
                elif error < best_error:
                    best_beta = beta
                    best_error = error
                    best_initial_guess = initial_guess
            importance_beta = best_beta

        print(f"\t{model} best found beta={beta} with error {error} using guess {initial_guess}")
        optim_result = pd.DataFrame(
            {"best_beta": [best_beta], "best_error": [best_error], "best_initial_guess": [best_initial_guess]},
            columns=['best_beta', 'best_error', 'best_initial_guess'])
        print(optim_result)
        optim_result.to_csv(f'data/{model}_beta.csv', index=False)
