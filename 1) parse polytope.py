import os
import logging
logging.getLogger('sbml').setLevel(logging.ERROR)
from PolyRound.api import PolyRoundApi
import cobra
import pandas as pd


if __name__ == "__main__":
    models = [
        'iEZ481_Glucose-MOPS',
        'iEZ481_PCA_Gluc',
        'iEZ481_Citrat-MOPS',
    ]

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

    for model in models:
        model_path = os.path.join("models", model + ".xml")
        cobra_model = cobra.io.read_sbml_model(model_path)
        # open up growth rate constraints
        for i, rxn in enumerate(cobra_model.reactions):
            if 'biomass' in rxn.id:
                rxn.upper_bound = optimal_growth_rates[model]
                rxn.lower_bound = 0
                print(rxn.id, f'({rxn.name})', 'has index', i, 'bounds', rxn.lower_bound, rxn.upper_bound)
        #     else:
        #         pass
        #         # print(rxn.id, 'has bounds', rxn.lower_bound, rxn.upper_bound)
        #         # rxn.upper_bound = optimal_growth_rates[model]
        #         # rxn.lower_bound = 0
        # print(model, 'optimal', optimal_growth_rates[model], 'mean', mean_growth_rates[model], 'std', std_growth_rates[model])
        polytope = PolyRoundApi.cobra_model_to_polytope(cobra_model)
        biomass_index = polytope.A.columns.tolist().index('biomass_a')
        print('biomass index is', biomass_index)
        polytope.A.to_csv(os.path.join('data', f'{model}_A.csv'))
        polytope.b.to_csv(os.path.join('data', f'{model}_b.csv'))
        polytope.S.to_csv(os.path.join('data', f'{model}_S.csv'))
        polytope.h.to_csv(os.path.join('data', f'{model}_h.csv'))
        polytope.shift.to_csv(os.path.join('data', f'{model}_shift.csv'))
        polytope.transformation.to_csv(os.path.join('data', f'{model}_transformation.csv'))

