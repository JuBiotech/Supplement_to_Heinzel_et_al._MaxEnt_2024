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
    var_lambda = growth_rates['mu'].std() ** 2
    print('mean growth rate', bar_lambda)
    print('var growth rate', var_lambda)
    samples = np.load(os.path.join('data', 'samples.npz'))['samples']
    print(samples.shape)
    n_samples = samples.shape[1]
    n_procs = samples.shape[0]
    samples = samples[:, :, biomass_index].reshape(16, 5000, 1)
    print(samples.shape)
    print("ess", hopsy.ess(samples))
    print(samples.shape)