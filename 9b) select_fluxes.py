import matplotlib.pyplot as plt
import hopsy
import cobra
import os
import numpy as np
import helpers
import pandas as pd
import time




if __name__ == "__main__":
    reactions_to_plot = set()
    model_path = os.path.join("models/iEZ481_PCA_Gluc.xml")
    cobra_model = cobra.io.read_sbml_model(model_path)
    reaction_names = [rxn.id for rxn in cobra_model.reactions]
    print(reaction_names)

    pathways = {
        "PPP": [ "zwf", "opcA", "gnd", "rpe", "rpi", "tkt_1", "tal", "tkt_2"],
        "EMP": ["pgi", "pfkA", "fbp", "fda", "tpiA", "gapA", "gapB", "eno", "pgk", "pgm", "pyk", "pps", "pdh", ],
        "ANA": [ "odx", "pyc", "mez", "ppc", "pckG"],
        "TCA": [ "gltA", "acnA", "acnB", "icd", "odhA", "sucD", "sdhCAB", "fumC", "mqo", "mdh", ],
        # "GLX": ["aceB", "aceA"],
    }

    for pathway, central_reactions in pathways.items():
        for c in central_reactions:
            if c not in reaction_names:
                print(c)
            assert c in reaction_names
        # print(len(central_reactions))

        PCA_GLC = pd.read_csv(f'PCA_Glucose_Boltzmann.txt', names=['wasserstein'])
        PCA_CIT = pd.read_csv(f'Citrat_PCA_Boltzmann.txt', names=['wasserstein'])
        CIT_GLC = pd.read_csv(f'Glucose_Citrat_Boltzmann.txt', names=['wasserstein'])

        assert len(reaction_names) == len(PCA_GLC.index)
        assert len(reaction_names) == len(PCA_CIT.index)
        assert len(reaction_names) == len(CIT_GLC.index)

        ind_list = [reaction_names.index(c) for c in central_reactions]
        PCA_GLC = PCA_GLC.iloc[ind_list] #.set_index(central_reactions)
        PCA_CIT = PCA_CIT.iloc[ind_list] #.set_index(central_reactions)
        CIT_GLC = CIT_GLC.iloc[ind_list] #.set_index(central_reactions)
        PCA_GLC['flux'] = central_reactions
        PCA_CIT['flux'] = central_reactions
        CIT_GLC['flux'] = central_reactions
        PCA_GLC = PCA_GLC.set_index('flux')
        PCA_CIT = PCA_CIT.set_index('flux')
        CIT_GLC = CIT_GLC.set_index('flux')

        print('pathway', pathway)
        print('PCA_GLC', PCA_GLC['wasserstein'].idxmax(), PCA_GLC['wasserstein'].max())
        print('CIT_GLC', CIT_GLC['wasserstein'].idxmax(), CIT_GLC['wasserstein'].max())
        print('PCA_CIT', PCA_CIT['wasserstein'].idxmax(), PCA_CIT['wasserstein'].max())
        reactions_to_plot.add(PCA_GLC['wasserstein'].idxmax())
        reactions_to_plot.add(CIT_GLC['wasserstein'].idxmax())
        reactions_to_plot.add(PCA_CIT['wasserstein'].idxmax())

    print(list(reactions_to_plot))
