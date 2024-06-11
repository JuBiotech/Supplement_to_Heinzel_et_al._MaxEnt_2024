import matplotlib.pyplot as plt
import hopsy
import os
import numpy as np
import helpers
import pandas as pd



if __name__ == "__main__":
    models = [
        'iEZ481_PCA_Gluc',
        'iEZ481_Glucose-MOPS',
        'iEZ481_Citrat-MOPS',
    ]
    for i, model in enumerate(models):
        u_samples = np.load(file=os.path.join('data', f'{model}_uniform_samples.npz'))['samples'][:, ::40, :]
        b_samples = np.load(file=os.path.join('data', f'{model}_boltzmann_samples.npz'))['samples'][:, ::35, :]

        print('shape uniform', u_samples.shape)
        print('shape boltzmann', b_samples.shape)

        np.savez_compressed(file=os.path.join('data', f'{model}_uniform_samples_thinned.npz'), samples=u_samples)
        np.savez_compressed(file=os.path.join('data', f'{model}_boltzmann_samples_thinned.npz'), samples=b_samples)
