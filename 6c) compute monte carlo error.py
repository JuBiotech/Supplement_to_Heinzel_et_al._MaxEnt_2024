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
        biomass_index=300
        samples = np.load(file=os.path.join('data', f'{model}_boltzmann_samples.npz'))['samples']
        print('mean', np.mean(samples[: , : , biomass_index]))

        std = np.std(samples[: , :, biomass_index])
        print('std', std)

        ess = hopsy.ess(samples)[0][biomass_index]
        print('ess', ess)

        print('For', model, ' Monte Carlo error is', std/ess)

