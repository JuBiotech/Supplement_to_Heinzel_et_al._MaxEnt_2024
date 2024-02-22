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
    polytope = helpers.load_polytope('data')
    problem = hopsy.Problem(polytope.A, polytope.b)
    biomass_index = 300
    problem = hopsy.add_box_constraints(problem, upper_bound=10, lower_bound=-10)
    problem = hopsy.round(problem)
    starting_point = hopsy.compute_chebyshev_center(problem)
    # Gleichverteilte Samples
    n_samples = 16_000
    n_procs = 16
    chains = [
        hopsy.MarkovChain(problem, proposal=hopsy.UniformHitAndRunProposal, starting_point=starting_point) for
        i in range(n_procs)]
    rng = [hopsy.RandomNumberGenerator(seed=i + 1123) for i in range(n_procs)]
    print('start sampling')
    _, samples = hopsy.sample(chains, rng, n_samples=n_samples, thinning=100, n_procs=n_procs, progress_bar=False)
    samples = samples[:, 1000:, :]
    print(samples.shape)
    np.savez_compressed(file=os.path.join('data', 'samples.npz'), samples=samples)
    # # Histogramm anzeigen lassen
    plt.hist(samples[0, :, 300], density=True, alpha=0.5)
    plt.hist(samples[1, :, 300], density=True, alpha=0.5)
    plt.hist(samples[2, :, 300], density=True, alpha=0.5)
    plt.hist(samples[3, :, 300], density=True, alpha=0.5)
    plt.show()
