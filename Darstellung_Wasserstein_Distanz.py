import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy

#%%
df = pd.read_csv("C:\\Users\\carol\\OneDrive\\Desktop\\Juelich\\extracted_growth_rates_with_media.csv")
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
# 2. Bestimme die dazugehörigen Flüsse
fluss_glucose_max1 = samples['boltzmann'][0, :, biomass_index]
#%%
fluss_pca_max1 = samples['boltzmann'][0, :, biomass_index]
samples_pca_max1 = samples['boltzmann'][0, :, :]
#%%
fluss_citrat_max1 = samples['boltzmann'][0, :, biomass_index]

#%%
# 3. Bestimme die Wasserstein-Distanz
res_w1 = scipy.stats.wasserstein_distance(glucose[0, :, biomass_index],
                                          pca[0, :, biomass_index])

res_w2 = scipy.stats.wasserstein_distance(glucose[0, :, biomass_index],
                                          citrat[0, :, biomass_index])

res_w3 = scipy.stats.wasserstein_distance(pca[0, :, biomass_index],
                                          citrat[0, :, biomass_index])
#%%
# Noch für alle möglichen Flüsse


#%%
# 4. Stelle die Wasserstein-Distanz grafisch dar als Heatmap

data = {
    'Glucose': [0, res_w1, res_w2],
    'PCA': [res_w1, 0, res_w3],
    'Citrat': [res_w2, res_w3, 0]
}
index = ['GLucose', 'PCA', 'Citrat']
df = pd.DataFrame(data, index=index)

plt.figure(figsize=(8, 6))  
sns.heatmap(df, cmap='coolwarm', cbar = True, cbar_kws={'label': 'Wasserstein-Distance'})  
plt.title('Biomass, Boltzmann-Distribution') 
plt.xlabel("Growth Medium")
plt.ylabel("Growth Medium")

plt.show()  
