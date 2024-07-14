import cobra
import os
import numpy as np
import scipy

if __name__ == "__main__":
    models = [
        'iEZ481_PCA_Gluc',
        'iEZ481_Glucose-MOPS',
        'iEZ481_Citrat-MOPS',
    ]
    b_samples = {}
    for i, model in enumerate(models):
        b_samples[model] = np.load(file=os.path.join('data', f'{model}_boltzmann_samples_thinned.npz'))['samples']
    reactions_to_plot = set()
    model_path = os.path.join("models/iEZ481_PCA_Gluc.xml")
    cobra_model = cobra.io.read_sbml_model(model_path)
    reaction_names = [rxn.id for rxn in cobra_model.reactions]

    pathways = {
        "PPP": ["zwf", "opcA", "gnd", "rpe", "rpi", "tkt_1", "tal", "tkt_2"],
        "EMP": ["pgi", "pfkA", "fbp", "fda", "tpiA", "gapA", "gapB", "eno", "pgk", "pgm", "pyk", "pps", "pdh", ],
        "ANA": ["odx", "pyc", "mez", "ppc", "pckG"],
        "TCA": ["gltA", "acnA", "acnB", "icd", "odhA", "sucD", "sdhCAB", "fumC", "mqo", "mdh", ],
        "GLX": ["aceB", "aceA"],
    }

    num_fluxes = 0
    mapped_fluxes = []
    for pathway, fluxes in pathways.items():
        num_fluxes += len(fluxes)
        mapped_fluxes += fluxes
    print('num fluxes', num_fluxes)
    print(mapped_fluxes)

    indices = []
    for i, r in enumerate(reaction_names):
        if r in mapped_fluxes:
            print(r)
            indices.append(i)

    print(len(indices), len(mapped_fluxes))
    print(indices)
    pca_glc = []
    cit_pca = []
    glc_cit = []
    for i, r in enumerate(mapped_fluxes):
        print('wasserstein for', r)
        print(indices[i])
        pca = b_samples[models[0]][:, :, indices[i]].flatten()
        gluc = b_samples[models[1]][:, :, indices[i]].flatten()
        cit = b_samples[models[2]][:, :, indices[i]].flatten()
        print('shapes', pca.shape, gluc.shape, cit.shape)
        pca_glc.append(scipy.stats.wasserstein_distance(pca, gluc))
        cit_pca.append(scipy.stats.wasserstein_distance(cit, pca))
        glc_cit.append(scipy.stats.wasserstein_distance(gluc, cit))

    print(mapped_fluxes)
    print(pca_glc)
    print(cit_pca)
    print(glc_cit)

    np.savetxt('data/wasserstein_pca-glc.csv', pca_glc, delimiter=',')
    np.savetxt('data/wasserstein_cit-pca.csv', cit_pca, delimiter=',')
    np.savetxt('data/wasserstein_glc_cit.csv', glc_cit, delimiter=',')
