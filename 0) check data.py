import pandas as pd

if __name__ == "__main__":
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

    models = [
        'iEZ481_Glucose-MOPS',
        'iEZ481_PCA_Gluc',
        'iEZ481_Citrat-MOPS',
    ]
    for model in models:
        print(f'{model} growth rates:',
              'mean',
              mean_growth_rates[model],
              'optimal',
              optimal_growth_rates[model],
              'diff',
              -mean_growth_rates[model] + optimal_growth_rates[model],
              'std',
              std_growth_rates[model],
              )
