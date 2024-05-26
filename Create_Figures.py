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
fba = pd.read_csv("C:\\Users\\carol\\OneDrive\\Desktop\\Juelich\\maxentflux\\citrat-MOPS_fluxes.csv")
df = pd.read_csv("C:\\Users\\carol\\OneDrive\\Desktop\\Juelich\\extracted_growth_rates_with_media.csv")
print(df.head())
biomass_index = raw_polytope.A.columns.tolist().index('biomass_a')
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
lst_plot = ['EX_ac_e','EX_acgam_e','EX_ala_L_e','EX_arg_L_e','EX_asn_L_e','EX_asp_L_e',
'EX_cit_e','EX_co2_e','EX_cys_L_e','EX_etoh_e','EX_fe2_e','EX_for_e','EX_fru_e',
'EX_fum_e','EX_gal_e','EX_glc_e','EX_glcn_e','EX_gln_L_e','EX_glu_L_e','EX_gly_e',
'EX_glyc_e','EX_h2o_e','EX_his_L_e','EX_ile_L_e','EX_lac_D_e','EX_lac_L_e',
'EX_leu_L_e','EX_lys_L_e','EX_mal_L_e','EX_man_e','EX_met_L_e','EX_na1_e',
'EX_nh3_e','EX_no2_e','EX_no3_e','EX_o2_e','EX_orn_e','EX_phe_L_e','EX_ppi_e',
'EX_pnto_R_e','EX_pro_L_e','EX_pyr_e','EX_rib_D_e','EX_ser_L_e','EX_so3_e',
'EX_suc_e','EX_sucr_e','EX_thr_L_e','EX_tre_e','EX_trp_L_e','EX_tyr_L_e',
'EX_ura_e','EX_urea_e','EX_val_L_e','EX_xyl_D_e','EX_ncam_e','EX_btn_e','EX_ppa_e',
'pyk','rpe','mez']
#%%
# Samples für die unterschiedlichen Media

glucose = samples['boltzmann']
#%%
pca = samples['normal']
#%%
citrat = samples['uniform']
#%%
# Das sind die Listen mit den Wasserstein-Distanzen
w_glucose_pca = []
w_glucose_citrat = []
w_pca_citrat = []
for i in lst_plot[0:1]:
    ind = raw_polytope.A.columns.tolist().index(i)
    # 4. Stelle die Wasserstein-Distanz grafisch dar als Heatmap
    res_w1 = scipy.stats.wasserstein_distance(glucose[0, :, ind],
                                          pca[0, :, ind])

    res_w2 = scipy.stats.wasserstein_distance(glucose[0, :, ind],
                                          citrat[0, :, ind])

    res_w3 = scipy.stats.wasserstein_distance(pca[0, :, ind],
                                          citrat[0, :, ind])
    w_glucose_pca.append(res_w1)
    w_glucose_citrat.append(res_w2)
    w_pca_citrat.append(res_w3)

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

#%%
fig, axs = plt.subplots(3, 1, figsize=(8, 12), sharex=True)
plt.subplots_adjust(hspace=0.3)  # Abstand zwischen Subplots einstellen
# Subplot 1: Boltzmann-Distribution
axs[0].hist(samples['boltzmann'][0, :, biomass_index], density=True, color="green", label="Boltzmann-Distribution")
axs[0].axvline(value, color='r', linestyle='dashed', linewidth=1)
axs[0].plot(x, norm.pdf(x, mu, sigma), 'b-', linewidth=2)
axs[0].set_title('Boltzmann-Distribution')
axs[0].legend()

# Subplot 2: Uniform-Distribution
axs[1].hist(samples['uniform'][0, :, biomass_index], density=True, color="blue", label="Uniform-Distribution", alpha=0.4)
axs[1].axvline(value, color='r', linestyle='dashed', linewidth=1)
axs[1].plot(x, norm.pdf(x, mu, sigma), 'b-', linewidth=2)
axs[1].set_title('Uniform-Distribution')
axs[1].legend()

# Subplot 3: Normal-Distribution
axs[2].hist(samples['normal'][0, :, biomass_index], density=True, color="orange", label="Normal-Distribution", alpha=0.4)
axs[2].axvline(value, color='r', linestyle='dashed', linewidth=1)
axs[2].plot(x, norm.pdf(x, mu, sigma), 'b-', linewidth=2)
axs[2].set_title('Normal-Distribution')
axs[2].legend()

# Beschriftung der x-Achse hinzufügen
axs[-1].set_xlabel('Biomass Flux')

# Beschriftung der y-Achse hinzufügen
for ax in axs:
    ax.set_ylabel('Probability')
    ax.legend()

plt.tight_layout()
plt.show()
#%%

# Beispiel-Liste
data = lambda_citrat

# Anzahl der Intervalle
num_bins = 10
# Diskretisierung mit pd.cut und automatischen Bins
categories = pd.cut(data, bins=num_bins, labels=False, retbins=True)
binned_data = categories[0]
bin_edges = categories[1]
bin_mids = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in range(len(bin_edges) - 1)]

# Erstellen eines DataFrames für die Visualisierung
df = pd.DataFrame({'Data': data, 'Binned': binned_data})

# Zählen der Häufigkeiten der Bins
bin_counts = df['Binned'].value_counts().sort_index()

# Plotten des Balkendiagramms
plt.bar(bin_mids, bin_counts.values, width=0.02, align='center', edgecolor='black')
plt.xlabel('Bins')
plt.ylabel('Frequency')
plt.title('Binned Data Histogram')


plt.show()
#%%
fig, axs = plt.subplots(3, 1, figsize=(8, 12), sharex=True)
plt.subplots_adjust(hspace=0.3)  # Abstand zwischen Subplots einstellen
# Subplot 1: Boltzmann-Distribution
axs[0].hist(samples['boltzmann'][0, :, biomass_index], density=True, color="green", label="Boltzmann-Distribution")
axs[0].axvline(value, color='r', linestyle='dashed', linewidth=1, label = "FBA Solution")
axs[0].plot(x, norm.pdf(x, mu, sigma), 'b-', linewidth=2, label = "Normal Approximation Measurement")
axs[0].set_title('Boltzmann-Distribution')
axs[0].legend()

# Subplot 2: Uniform-Distribution
axs[1].hist(samples['uniform'][0, :, biomass_index], density=True, color="blue", label="Uniform-Distribution", alpha=0.4)
axs[1].axvline(value, color='r', linestyle='dashed', linewidth=1, label = "FBA Solution")
axs[1].plot(x, norm.pdf(x, mu, sigma), 'b-', linewidth=2, label = "Normal Approximation Measurement")
axs[1].set_title('Uniform-Distribution')
axs[1].legend()

# Subplot 3: Normal-Distribution
axs[2].hist(samples['normal'][0, :, biomass_index], density=True, color="orange", label="Normal-Distribution", alpha=0.4)
axs[2].axvline(value, color='r', linestyle='dashed', linewidth=1, label = "FBA Solution")
axs[2].plot(x, norm.pdf(x, mu, sigma), 'b-', linewidth=2, label = "Normal Approximation Measurement")
# Plotten des Balkendiagramms
axs[2].bar(bin_mids, bin_counts.values, width=0.02, align='center', edgecolor='black')

# Setzen der xticks auf die Mitte der Bins
#bin_labels = [f'{bin_edges[i]:.1f}-{bin_edges[i+1]:.1f}' for i in range(len(bin_edges)-1)]

axs[2].set_title('Normal-Distribution')
axs[2].legend()

# Beschriftung der x-Achse hinzufügen
axs[-1].set_xlabel('Biomass Flux')

# Beschriftung der y-Achse hinzufügen
for ax in axs:
    ax.set_ylabel('Probability')
    ax.legend()

plt.tight_layout()
plt.show()
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%
# Verteilungen untereinander


for i in lst_plot:
    ind = raw_polytope.A.columns.tolist().index(i)
    plt.title(i)
    fig, axs = plt.subplots(3, 1, figsize=(8, 12), sharex=True) 
    plt.subplots_adjust(hspace=0.1)  # Abstand zwischen Subplots einstellen
    axs[0].hist(samples['boltzmann'][0, :, ind], density=True, color="green", label="Boltzmann-Distribution")
    axs[1].hist(samples['uniform'][0, :, ind], density=True, color="blue", label="Uniform-Distribution", alpha=0.4)
    axs[2].hist(samples['normal'][0, :, ind], density=True, color="orange", label="Normal-Distribution", alpha=0.4)
    # Beschriftung der x-Achse hinzufügen
    axs[-1].set_xlabel('Biomass Flux')
    # Alle Diagramme anzeigen
    for ax in axs:
        ax.legend()
    plt.tight_layout()
    plt.show()