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
        print(f'model {model}')
        biomass_index = 300
        samples = np.load(file=os.path.join('data', f'{model}_boltzmann_samples.npz'))['samples']
        print('shape samples', samples.shape)
        print('ess', np.min(hopsy.ess(samples)))
        print('rhat', np.max(hopsy.rhat(samples)))
        print('mean', np.mean(samples[:, :, biomass_index]))
        print('std', np.std(samples[:, :, biomass_index]))

        # Histogramm anzeigen lassen
        # plt.figure()
        # plt.title(f'{model} growth rate')
        # plt.hist(samples[0, :, 300], density=True, bins=10, alpha=0.5)
        # plt.hist(samples[1, :, 300], density=True, bins=10, alpha=0.5)
        # plt.hist(samples[2, :, 300], density=True, bins=10, alpha=0.5)
        # plt.hist(samples[3, :, 300], density=True, bins=10, alpha=0.5)
        # plt.tight_layout()
        # plt.savefig(f"boltzmann_{model}.png")
        # plt.show(block=False)
        # del samples
