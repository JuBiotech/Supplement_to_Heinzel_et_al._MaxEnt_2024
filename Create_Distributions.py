import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy
import numpy as np
#%%
# 0) Lade die benötigten Daten
# 1) Plotte die Wasserstein Distanzen
# 2) Plotte die Flux Distributions
#%%
# Ergebnisse FBA Fredrik
fba = pd.read_csv("C:\\Users\\carol\\OneDrive\\Desktop\\Juelich\\maxentflux\\PCA-Gluc_fluxes.csv")
df = pd.read_csv("C:\\Users\\carol\\OneDrive\\Desktop\\Juelich\\extracted_growth_rates_with_media.csv")
data = np.load("C:\\Users\\carol\\Downloads\\iEZ481_PCA_Gluc_gauss_samples.npz", allow_pickle=True)

print(df.head())
#%%
# 1. Bestimme die maximale und durchschnittliche Wachstumsrate
# Glucose
filtered_df_glucose = df[df['medium'] == 'Glucose-MOPS']
filtered_df_glucose = filtered_df_glucose[filtered_df_glucose['method'] == 'area']
lambda_glucose = filtered_df_glucose['mu']
lambda_max_clucose = max(lambda_glucose)
lambda_bar_glucose = np.mean(lambda_glucose)

# PCA
filtered_df_pca = df[df['medium'] == 'PCA-Gluc']
filtered_df_pca = filtered_df_pca[filtered_df_pca['method'] == 'area']
lambda_pca = filtered_df_pca['mu']
lambda_max_pca = max(lambda_pca)
lambda_bar_pca = np.mean(lambda_pca)
# Citrat
filtered_df_citrat = df[df['medium'] == 'Citrat-MOPS']
filtered_df_citrat = filtered_df_citrat[filtered_df_citrat['method'] == 'area']
lambda_citrat = filtered_df_citrat['mu']
lambda_max_citrat = max(lambda_citrat)
lambda_bar_citrat  = np.mean(lambda_citrat)
#%%
# Wachstumsraten plotten

data_pca = lambda_pca
num_bins = 3
categories = pd.cut(data_pca, bins=num_bins, labels=False, retbins=True)
binned_data = categories[0]
bin_edges = categories[1]
bin_mids = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in range(len(bin_edges) - 1)]
df_lambda = pd.DataFrame({'Data': data, 'Binned': binned_data})
bin_counts = df_lambda['Binned'].value_counts().sort_index()
#%%
value = fba.loc[301, 'mean+0stds']
plt.hist(data['samples'][0, :, 300], density=True, color="green", label="Normal-Distribution")
plt.axvline(value, color='r', linestyle='dashed', linewidth=1, label = "FBA Solution")
plt.bar(bin_mids, bin_counts.values, width=0.0129, align='center', edgecolor='black', label = "Empirical Growth Rate")
plt.title('Biomass, PCA')
plt.xlabel("Flux")
plt.ylabel("Probability")

plt.legend()