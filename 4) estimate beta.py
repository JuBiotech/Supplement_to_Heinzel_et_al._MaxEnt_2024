import matplotlib.pyplot as plt
import hopsy
import os
import sys
import numpy as np
from scipy.optimize import fsolve, fmin
import pandas as pd
import logging
from PolyRound.api import PolyRoundApi
import scipy as sc
import helpers

# 1) Lade die benötigten Dateien und lege die Paramter fest
# 2) Bestimme die durchschnittliche Wachstumsrate und die Varianz dieser.
# 3) Bestimme \beta aus der Boltzmann-Verteilung + Bestimme den
# Monte-Carlo-Fehler.
# 4) Stelle die Flüsse grafisch dar.
# %%
# Hier stehen alle änderbaren Variablen.
# Laden der Daten


if __name__ == "__main__":
    biomass_index = 300
    model_path = os.path.join("models", "iEZ481_Glc.xml")
    growth_rates = pd.read_csv("data/extracted_growth_rates.csv")
    # use area for growth rate estimation, it is most reliable according to johannes
    growth_rates = growth_rates[(growth_rates['method'] == "area")]
    print(growth_rates)
    bar_lambda = growth_rates['mu'].mean()
    var_lambda = growth_rates['mu'].std()**2
    print('mean growth rate', bar_lambda)
    print('var growth rate', var_lambda)
    samples = np.load(os.path.join('data', 'samples.npz'))['samples']
    print(samples.shape)
    n_samples = samples.shape[1]
    n_procs = samples.shape[0]
    samples = samples[:, :, biomass_index].reshape(n_samples * n_procs, -1)
    print(samples.shape)
    # Histogramm anzeigen lassen
    plt.hist(samples)
    plt.show()
    res_samples = []
    for i in samples:
        res_samples.append(i[0])

    # %%
    # Bestimmung der zu lösenden Gleichung
    # Gleichung 11.35 aus Martino wird nach beta aufgelöst
    def lambda_function(flux_samples, _beta):
        res = 0
        norm = 0
        for flux in flux_samples:
            log_term = flux * _beta
            norm += np.exp(log_term)
            res += flux * norm
        if norm < 1e-4:
            raise RuntimeError("dividing by 0")
        return res / (norm)

    # Definieren der Funktion, die das Integral berechnet - Abhängig von beta
    def quad_integrand(beta):
        return lambda_function(res_samples, beta)

    # Definieren der Funktion, die integriert wird - Abhängig von beta
    def equation_to_solve(_beta):
        return -bar_lambda + quad_integrand(_beta)


    # Schätze Startwert für beta
    initial_guess = np.array([1e2])
    result = fsolve(equation_to_solve, initial_guess, full_output=True, xtol=1e-13, maxfev=int(1e6))
    # result = fmin(equation_to_solve, initial_guess, args=(bar_lambda,), full_output=True, maxiter=int(1e7))

    print('solve info')
    print(result)
    # Das Ergebnis enthält den gefundenen Wert für beta
    beta = np.array(result[0])
    print(f"Der gefundene Wert für beta ist: {beta}")
    # %%
    # Wie große ist der Fehler bei der Schätzung von beta?
    print('error is around', lambda_function(res_samples, beta) - bar_lambda)
