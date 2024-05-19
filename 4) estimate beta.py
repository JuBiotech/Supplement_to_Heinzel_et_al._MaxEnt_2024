import matplotlib.pyplot as plt
import os
import numpy as np
import pandas
from scipy.optimize import fsolve, fmin
from scipy.special import logsumexp
import pandas as pd



if __name__ == "__main__":
    biomass_index = 300
    models = [
        'iEZ481_Glucose-MOPS',
        'iEZ481_PCA_Gluc',
        'iEZ481_Citrat-MOPS',
    ]

    for model in models:
        print(f'finding beta for {model}')
        model_path = os.path.join("models", f"{model}.xml")

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

        bar_lambda = mean_growth_rates[model]
        var_lambda = std_growth_rates[model]**2
        print(f'\t{model}mean growth rate', bar_lambda)
        print(f'\t{model}var growth rate', var_lambda)
        samples = np.load(os.path.join('data', f'{model}_samples.npz'))['samples']
        n_samples = samples.shape[1]
        n_procs = samples.shape[0]
        samples = samples[:, :, biomass_index].reshape(n_samples * n_procs, -1)
        res_samples = []
        for i in samples:
            res_samples.append(i[0])

        # Gleichung 11.35 aus Martino wird nach beta aufgelöst
        def lambda_function(flux_samples, _beta):
            shifted = _beta * flux_samples + np.log(flux_samples)
            log_res1 = logsumexp(shifted)
            log_res2 = logsumexp(_beta*flux_samples)
            log_res = log_res1-log_res2
            # print('log_res', log_res1, log_res2, log_res, 'for beta', _beta)
            return log_res

        # Definieren der Funktion, die das Integral berechnet - Abhängig von beta
        def quad_integrand(beta):
            return lambda_function(res_samples, beta)

        # Definieren der Funktion, die integriert wird - Abhängig von beta
        def equation_to_solve(_beta):
            return -np.log(bar_lambda) + quad_integrand(_beta)

        initial_guesses = [-1e6, -1e4, -1e2, 0, 1e2, 1e4, 1e6]
        best_beta = None
        best_error = None
        best_initial_guess = None
        for ig in initial_guesses:
            initial_guess = np.array([ig])
            result = fsolve(equation_to_solve, initial_guess, full_output=True, xtol=1e-18, maxfev=int(1e5))

            print(f'\t{model} with initial {initial_guess} returns solve info {result}')
            # Das Ergebnis enthält den gefundenen Wert für beta
            beta = np.array(result[0])
            error = np.abs(lambda_function(res_samples, beta) - np.log(bar_lambda))
            print(f"\t\t found beta={beta} with error {error}")
            if best_error is None:
                best_beta = beta
                best_error = error
                best_initial_guess = initial_guess
            elif error < best_error:
                best_beta = beta
                best_error = error
                best_initial_guess = initial_guess

        print(f"\t{model} best found beta={beta} with error {error} using guess {initial_guess}")
        optim_result = pd.DataFrame({"best_beta": [best_beta], "best_error": [best_error], "best_initial_guess": [best_initial_guess]}, columns=['best_beta', 'best_error', 'best_initial_guess'])
        print(optim_result)
        optim_result.to_csv(f'data/{model}_beta.csv', index=False)
