import matplotlib.pyplot as plt
import hopsy
import os
import numpy as np
from scipy.optimize import fsolve
import pandas as pd
from PolyRound.api import PolyRoundApi
import seaborn as sns
import scipy
import pta
#%%
# Dieser Code ist folgendermaßen gegliedert:
# 1) Lade die benötigten Dateien und lege die Paramter fest
# 2) Bestimme die durchschnittliche Wachstumsrate und die Varianz dieser.
# 3) Bestimme \beta aus der Boltzmann-Verteilung
# 4) Stelle die Flüsse grafisch dar (Verteilung und Vgl mit Wasserstein Distanz).
#%%
# Hier stehen alle änderbaren Variablen.
# Laden der Daten
model_path = os.path.join("C:/Users/carol/OneDrive/Desktop/Juelich", "iEZ481_Glc.xml")
# Datei zur Bestimmung der Wachstumsrate
df = pd.read_csv("C:\\Users\\carol\\OneDrive\\Desktop\\Juelich\\extracted_growth_rates_with_media.csv")
print(df.head())
#%%
raw_polytope = PolyRoundApi.sbml_to_polytope(model_path)
biomass_index = raw_polytope.A.columns.tolist().index('biomass_a')
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
# Abstand der Messungen des Zellwachstums in Minuten
ab_mess = 15
# Parameter für die Simulation
n_chains = 4
n_samples = 10000
thinning = 10
# Bestimmung des Monte-Carlo-Fehlers
n_int = 10 # Anzahl der simulierten Integrale
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
# Bestimmung der zu lösenden Gleichung
# Gleichung 11.35 aus Martino wird nach beta aufgelöst
def lambda_function(v1, beta):
    res = 0
    res1 = 0
    for i in v1:
        log_term = np.float64(i * beta)
        res +=  i*np.exp(log_term)
        res1 +=  np.exp(log_term)
    return res / (res1 + 0.0000001) # +0.0000001 notwendig, damit wir auf 
# keinen Fall durch 0 teilen

# Definieren der Funktion, die das Integral berechnet - Abhängig von beta
def quad_integrand(beta, res_samples):
    result = lambda_function(res_samples, beta)
    return result

# Definieren der Funktion, die integriert wird - Abhängig von beta
def equation_to_solve(beta, c, res_samples):
    return c - quad_integrand(beta, res_samples)
#%% 
# Definiere die beiden verschiedenen Modelle
class Boltzman_Modell:
    def __init__(self, beta, index):
        self.beta = beta
        self.index = index
        
    def compute_negative_log_likelihood(self, x):
        return - self.beta * x[self.index]
    
class Normal_Modell:
    def __init__(self, mu, sigma, index):
        self.mu = mu
        self.index = index
        self.sigma = sigma
        
    def compute_negative_log_likelihood(self, x):
        return 0.5 * (x[self.index]-self.mu)**2/(self.sigma)

#%%
def main(distributions, change_flux, plot, lamb, var, var_lambda):
    '''
    

    Parameters
    ----------
    ditributions : List
        Which distributions (boltzmann, uniform, normal) should be considered?
        Input ["boltzmann", "uniform", "normal"], if every distribution should be considered,
        ["boltzmann", 0, "normal"] if the uniform distribution not, order is important!
    change_flux : List
        change_flux = []: no changes for the maximal and minimal flux. 
        change_flux = [lb, ub]: lower bound and upper bound for the flux
    plot : Int
        1: Plots are printed.
    lamb : Int
        Mean biomass flux.
    var : Int
        Variance biomass flux.

    Returns
    -------
    samples : Dict
        Distributions of the flux for the considered distributions.

    '''
    
    if(len(change_flux) >0):
        model = pta.load_example_model(change_flux[0])
        model.reactions.get_by_id('biomass_a').lower_bound = change_flux[0]
        model.reactions.get_by_id('biomass_a').upper_bound = change_flux[1]
    if(distributions[0] == "Boltzmann"):
        # Bestimme den Parameter der Boltzmann-Verteilung
        polytope = PolyRoundApi.sbml_to_polytope(model_path)
        problem = hopsy.Problem(polytope.A, polytope.b)
        problem = hopsy.round(problem)
        starting_point = hopsy.compute_chebyshev_center(problem)
        # Gleichverteilte Samples
        chains = [hopsy.MarkovChain(problem, starting_point = starting_point) for i in range(1)]
        rng = [hopsy.RandomNumberGenerator(seed= i) for i in range(1)]
        accrate_e_coli, samples = hopsy.sample(chains, rng, n_samples=100, thinning=10)
        samples = samples[:,:,biomass_index][0]
        res_samples = []
        for i in samples:
            res_samples.append(i)
        initial_guess = 0
        bar_lambda = lamb
        result_beta = fsolve(equation_to_solve, initial_guess, args=(bar_lambda/60,res_samples))
        beta = np.array(result_beta[0])
        a = np.array(beta) 
        model_Boltzmann = Boltzman_Modell(a, biomass_index)
        polytope = PolyRoundApi.simplify_transform_and_round(raw_polytope)
        boltzmann = hopsy.Problem(A=polytope.A, b=polytope.b, model = model_Boltzmann, transformation=polytope.transformation, shift=polytope.shift)
        starting_point_2 = hopsy.compute_chebyshev_center(boltzmann)
    # Simuliere für jede Verteilung die entsprechenden Daten
    # Für die Uniform-Verteilung
    if(distributions[1] == "uniform"):
        uniform = hopsy.Problem(A=polytope.A, b=polytope.b, transformation=polytope.transformation, shift=polytope.shift)
        starting_point = hopsy.compute_chebyshev_center(uniform)
    if(distributions[2] == "Normal"):
    # Für die Normalverteilung
        var_lambda = 1
        sigma = np.array(var_lambda/3600) # Umrechnung von Minuten^2 in Stunden^2
        mu = np.array(lamb/60) # Umrechnung von Minuten in Stunden 
        model_normal = Normal_Modell(mu, sigma, biomass_index)
        normal = hopsy.Problem(A=polytope.A, b=polytope.b, model = model_normal, transformation=polytope.transformation, shift=polytope.shift)
        starting_point_3 = hopsy.compute_chebyshev_center(normal)
    if(distributions[0] == "Boltzmann" and distributions[1] == "uniform" and distributions[2] == "Normal"):
        problems = {
            'uniform': [uniform, starting_point],
            'boltzmann': [boltzmann, starting_point_2],
            'normal': [normal, starting_point_3]
        }
    elif(distributions[0] == 0 and distributions[1] == "uniform" and distributions[2] == "Normal"):
        problems = {
            'uniform': [uniform, starting_point],
            'normal': [normal, starting_point_3]
        }
    elif(distributions[0] == "Boltzmann" and distributions[1] == 0 and distributions[2] == "Normal"):
        problems = {
            'boltzmann': [boltzmann, starting_point_2],
            'normal': [normal, starting_point_3]
        }
    elif(distributions[0] == "Boltzmann" and distributions[1] == "uniform" and distributions[2] == 0):
        problems = {
            'uniform': [uniform, starting_point],
            'boltzmann': [boltzmann, starting_point_2],
        }
    elif(distributions[0] == "Boltzmann" and distributions[1] == 0 and distributions[2] == 0):
        problems = {
            'boltzmann': [boltzmann, starting_point_2],
        }
    elif(distributions[0] == 0 and distributions[1] == "uniform" and distributions[2] == 0):
        problems = {
            'uniform': [uniform, starting_point],
        }
    elif(distributions[0] == 0 and distributions[1] == 0 and distributions[2] == "Normal"):
        problems = {
            'normal': [normal, starting_point_3]
        }
    rhat = {}
    ess = {}
    samples = {}
    acceptance_rate = {}

    for p, liste in problems.items():
        v = liste[0]
        s = liste[1]
        print(p)
        proposal = hopsy.UniformCoordinateHitAndRunProposal(v, starting_point=s)
        mcs = [hopsy.MarkovChain(problem=v, proposal=proposal) for i in range(n_chains)]
        rngs = [hopsy.RandomNumberGenerator(i) for i in range(n_chains)]
        acceptance_rate[p], samples[p] = hopsy.sample(mcs, rngs, n_samples=n_samples, thinning=thinning, n_procs=1)
        ess[p] = hopsy.ess(samples[p])
        print('\tess', np.min(ess[p]))
        rhat[p] = hopsy.rhat(samples[p])
        print('\trhat', np.max(rhat[p]))

    if(plot == "all"):
        # Hier werden die Bilder für die einzelnen Flüsse geplottet
        # Dabei werden die einzelnen Flüsse alle in einer Abbildung dargestellt.
        # Je nachdem, welcher Fluss ausgewählt wurde, sind die Unterschiede sehr gering
        for i in lst_plot:
            ind = raw_polytope.A.columns.tolist().index(i)
            plt.title(i)
            for m in distributions:
                if(m != 0):
                    plt.hist(samples[m][0, :, ind], density = True, color = "yellow", label = "Normal-Verteilung")
                    plt.legend(loc='upper right')
                    plt.show()
    elif(plot == "biomass"):
        # Nur der Biomasse-Fluss
        # Verteilungen untereinander
        fig, axs = plt.subplots(3, 1, figsize=(8, 12), sharex=True)
        plt.subplots_adjust(hspace=0.1)  # Abstand zwischen Subplots einstellen
        axs[0].hist(samples['boltzmann'][0, :, biomass_index], density=True, color="green", label="Boltzmann-Distribution")
        axs[1].hist(samples['uniform'][0, :, biomass_index], density=True, color="blue", label="Uniform-Distribution", alpha=0.4)
        axs[2].hist(samples['normal'][0, :, biomass_index], density=True, color="orange", label="Normal-Distribution", alpha=0.4)
        # Beschriftung der x-Achse hinzufügen
        axs[-1].set_xlabel('Biomass Flux')
        # Alle Diagramme anzeigen
        for ax in axs:
            ax.legend()
        plt.tight_layout()
        plt.show()
    return samples
#%%
samples = main(["Boltzmann", 0, 0], [], 0, 1, 1, 1)
#%%
# Führe den Code für die unterschiedlichen Media aus
glucose = samples['boltzmann']
#%%
pca = samples['boltzmann']
#%%
citrat = samples['boltzmann']
#%%
# Das sind die Listen mit den Wasserstein-Distanzen
w_glucose_pca = []
w_glucose_citrat = []
w_pca_citrat = []
for i in lst_plot:
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
