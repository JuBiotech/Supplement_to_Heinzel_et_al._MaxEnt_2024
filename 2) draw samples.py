import matplotlib.pyplot as plt
import hopsy
import os
import numpy as np
import helpers

if __name__ == "__main__":
    models = [
        'iEZ481_Citrat-MOPS',
        'iEZ481_Glucose-MOPS',
        'iEZ481_PCA_Gluc',
    ]
    for model in models:
        polytope = helpers.load_polytope('data', model)
        problem = hopsy.Problem(polytope.A, polytope.b)
        problem = hopsy.add_equality_constraints(problem, A_eq=polytope.S, b_eq=polytope.h)
        biomass_index = 300
        problem = hopsy.round(problem)
        starting_point = hopsy.compute_chebyshev_center(problem)
        # Gleichverteilte Samples
        n_samples = 641_000
        n_procs = 4
        chains = [
            hopsy.MarkovChain(problem, proposal=hopsy.UniformCoordinateHitAndRunProposal, starting_point=starting_point) for
            i in range(n_procs)]
        rng = [hopsy.RandomNumberGenerator(seed=i + 1123) for i in range(n_procs)]
        print(f'start sampling {model}')
        _, samples = hopsy.sample(chains, rng, n_samples=n_samples, thinning=10, n_procs=n_procs, progress_bar=False)
        print(f'finished sampling {model}')
        print('thinned ESS ', np.min(hopsy.ess(samples[:,::1000,:])))
        print('thinned rhat ', np.max(hopsy.rhat(samples[:,::1000,:])))
        samples = samples[:, 1000:, :]
        print(samples.shape)
        np.savez_compressed(file=os.path.join('data', f'{model}_samples.npz'), samples=samples)
        # # Histogramm anzeigen lassen
        plt.figure()
        plt.title(f'{model} growth rate')
        plt.hist(samples[0, :, 300], density=True, bins=1000, alpha=0.5)
        plt.hist(samples[1, :, 300], density=True, bins=1000, alpha=0.5)
        plt.hist(samples[2, :, 300], density=True, bins=1000, alpha=0.5)
        plt.hist(samples[3, :, 300], density=True, bins=1000, alpha=0.5)
        plt.show(block=False)
        del samples
