import hopsy
import os
import numpy as np
from scipy.optimize import fsolve
import pandas as pd
from scipy.optimize import curve_fit
from PolyRound.api import PolyRoundApi
import scipy as sc
#%%
# To Do: ESS ist viel zu klein, da müsste man die Anzahl der Stichproben erhöhen.
# Bin mir beim relativen MC-Fehler unsicher, ob dieser richtig berechnet wurde.
# Wir wollen uns ja noch anschauen, ob die ein-dimensionalen Vereilungen bei den verschiedenen 
# Verteilungen gleich sind.
# Dazu würde ich den Kolmogorov Smirnov-Test verwenden. Im hochdimensionalen ist das aber schwierig.
# Bei den Einheiten bin ich mir auch unsicher. Es ist jetzt alles in Stunden umgerechnet.
#%%
# Dieser Code ist folgendermaßen gegliedert:
# 1) Lade die benötigten Dateien und lege die Paramter fest
# 2) Bestimme die durchschnittliche Wachstumsrate und die Varianz dieser.
# 3) Bestimme \beta aus der Boltzmann-Verteilung + Bestimme den 
# Monte-Carlo-Fehler.
# 4) Stelle die Flüsse grafisch dar.
#%%
# Hier stehen alle änderbaren Variablen.
# Laden der Daten


if __name__ == "__main__":
    model_path = os.path.join("models", "iEZ481_Glc.xml")
    # Datei zur Bestimmung der Wachstumsrate
    with open("data/growth_rate_ds.csv", newline='') as csvfile:
        # Erstellen Sie ein CSV-Leserobjekt
        csvreader = pd.read_csv(csvfile)
    # Bestimmung des  Mediums - hängt auch von der Datei ab
    Acetate = csvreader[csvreader['medium'] == 'Acetate-MOPS']
    Citrat = csvreader[csvreader['medium'] == 'Citrat-MOPS']
    Gluconate = csvreader[csvreader['medium'] == 'Gluconate-MOPS']
    CA = csvreader[csvreader['medium'] == 'PCA-Gluc']
    BHI = csvreader[csvreader['medium'] == 'BHI']
    Fructose = csvreader[csvreader['medium'] == 'Fructose-MOPS']
    Glucose = csvreader[csvreader['medium'] == 'Glucose-MOPS']
    Pyruvate = csvreader[csvreader['medium'] == 'Pyruvate-MOPS']
    # Hier wird Glucose als Medium betrachtet
    data = Glucose["cell_count"]
    laenge_teilliste = max(Glucose["frame"])
    raw_polytope = PolyRoundApi.sbml_to_polytope(model_path)
    biomass_index = raw_polytope.A.columns.tolist().index('biomass_a')
    # Zu betrachtende Flüsse
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
    n_samples = 100_000
    thinning = 100
    # Bestimmung des Monte-Carlo-Fehlers
    n_int = 10 # Anzahl der simulierten Integrale
    #%%
    # Bestimmung der Wachstumsrate durch exponentiellen Fit an die Daten
    maximaler_wert = (laenge_teilliste + 1) * ab_mess
    t = []
    for i in range(0, maximaler_wert, ab_mess):
        t.append(i)

    teillisten = []

    # Schleife, um die Ausgangsliste in Teillisten aufzuteilen
    def teile_liste(input_list, n):
        # Initialisiere eine leere Ergebnisliste
        output_list = []
        # Schleife, um die Eingabeliste in Teillisten aufzuteilen
        for i in range(0, len(input_list), n):
            teil_liste = input_list[i:i+n]
            output_list.append(teil_liste)
        return output_list

    teillisten = teile_liste(data, laenge_teilliste + 1)

    # Ansatz für die Regression
    def exponential_growth(t, r):
        return c * np.exp(r * t)

    # Epxonentieller Fit
    mw = []
    for i in teillisten:
        c = i.tolist()[0]
        N = i.tolist()
        # Schätze Parameter c und r mit Regression
        params, covariance = curve_fit(exponential_growth, t, N, maxfev = 100000)
        # Extrahieren der geschätzten Parameter
        r_estimated = params
        mw.append(r_estimated)
        print('est', r_estimated, 'covariance', covariance)
    # Mittelwert der Wachstumsrate
    bar_lambda = np.mean(mw)
    # Varianz der Wachstumsrate
    var_lambda = np.var(mw)
    print('bar_lambda', bar_lambda)
    print('var_lambda', var_lambda)
    # Plotten
    #y_fit = []
    #for i in t:
    #    y_fit.append(c * np.exp(bar_lambda * i))
    #plt.plot(t, teillisten[0], marker='o')
    #plt.plot(t, y_fit)
    #plt.show()
    # Wichtig: Hierbei ist die Einheit Minuten^2
    #print(bar_lambda, var_lambda)
    #%%
    # Bestimme den Parameter der Boltzmann-Verteilung
    polytope = PolyRoundApi.sbml_to_polytope(model_path)
    problem = hopsy.Problem(polytope.A, polytope.b)
    problem = hopsy.round(problem)
    starting_point = hopsy.compute_chebyshev_center(problem)
    # Gleichverteilte Samples
    chains = [hopsy.MarkovChain(problem, starting_point = starting_point) for i in range(4)]
    rng = [hopsy.RandomNumberGenerator(seed=i*90123) for i in range(4)]
    accrate_e_coli, samples = hopsy.sample(chains, rng, n_samples=10000, thinning=100, n_procs=4)
    samples = samples[:,:,biomass_index][0]
    # Histogramm anzeigen lassen
    plt.hist(samples)
    plt.show()
    res_samples = []
    # Ziel: Python kann die Werte wieder ausrechnen
    # ggf. notwendig i/100 o.ä reinzuschreiben
    for i in samples:
        res_samples.append(i)
    #%%
    # Bestimmung der zu lösenden Gleichung
    # Gleichung 11.35 aus Martino wird nach beta aufgelöst
    def lambda_function(v1, beta):
        res = 0
        res1 = 0
        for i in v1:
            log_term = np.float64(i * beta)
            res +=  i*sc.special.expit(log_term)
            res1 +=  sc.special.expit(log_term)
        return res / (res1 + 0.0000001) # +0.0000001 notwendig, damit wir auf
    # keinen Fall durch 0 teilen

    # Definieren der Funktion, die das Integral berechnet - Abhängig von beta
    def quad_integrand(beta):
        result = lambda_function(res_samples, beta)
        return result

    # Definieren der Funktion, die integriert wird - Abhängig von beta
    def equation_to_solve(beta, c):
        return c - quad_integrand(beta)

    # Schätze Startwert für beta
    initial_guess = 0

    # Verwende fsolve, um die Gleichung zu lösen und den Wert für beta zu finden
    # Die Einheit ist jetzt Stunden
    # Achtung: Die Ändeung von bar_lambda muss schon sehr groß sein,
    # danit sich result_beta ändert
    result_beta = fsolve(equation_to_solve, initial_guess, args=(bar_lambda/60,))

    # Das Ergebnis enthält den gefundenen Wert für beta
    print(f"Der gefundene Wert für beta ist: {result_beta[0]}")
    beta = np.array(result_beta[0])
    #%%
    # Wie große ist der Fehler bei der Schätzung von beta?
    print('fehler ist', lambda_function(res_samples, beta) - c)
    exit(0)
    #%%
    # Bestimmung MC-Fehler
    def lambda_function(v1, beta):
        res = 0
        res1 = 0
        for i in v1:
            res+= i * sc.special.expit(i*beta)
            res1 += sc.special.expit(i*beta)
        return res, res1
    # Liste, in der der relative MC-Fehler steht
    var = []
    c = [beta]
    for j in c:
        # Integral im Zähler
        result = []
        # Integral im Nenner
        result1 = []
        # Effecive Samples Size
        ess = []
        for k in range(n_int):
            starting_point = hopsy.compute_chebyshev_center(problem)
            chains = [hopsy.MarkovChain(problem, starting_point = starting_point) for i in range(1)]
            rng = [hopsy.RandomNumberGenerator(seed= k) for i in range(1)]
            accrate, samples_uniform = hopsy.sample(chains, rng, n_samples=10000, thinning=100)
            ess.append(min(min(hopsy.ess(samples_uniform))))
            res, res1 = lambda_function(samples_uniform, j)
            result.append(res)
            result1.append(res1)
        v = 0
        for i in range(n_int):
            v += (np.std(result)/ess[i])/np.mean(result) - (np.std(result1)/ess[i])/np.mean(result1)
        var.append(v)
    #%%
    print(var)
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
    # Simuliere für jede Verteilung die entsprechenden Daten
    polytope = PolyRoundApi.simplify_transform_and_round(raw_polytope)
    # Für die Uniform-Verteilung
    uniform = hopsy.Problem(A=polytope.A, b=polytope.b, transformation=polytope.transformation, shift=polytope.shift)
    starting_point = hopsy.compute_chebyshev_center(uniform)
    #%%
    # Für die Boltzmann-Verteilung
    a = np.array(beta) # Geschätzter Parameter beta
    model_Boltzmann = Boltzman_Modell(a, biomass_index)
    boltzmann = hopsy.Problem(A=polytope.A, b=polytope.b, model = model_Boltzmann, transformation=polytope.transformation, shift=polytope.shift)
    starting_point_2 = hopsy.compute_chebyshev_center(uniform)
    #%%
    # Für die Normalverteilung
    sigma = np.array(var_lambda/3600) # Umrechnung von Minuten^2 in Stunden^2
    mu = np.array(bar_lambda/60) # Umrechnung von Minuten in Stunden
    model_normal = Normal_Modell(mu, sigma, biomass_index)
    normal = hopsy.Problem(A=polytope.A, b=polytope.b, model = model_normal, transformation=polytope.transformation, shift=polytope.shift)
    starting_point_3 = hopsy.compute_chebyshev_center(normal)
    #%%
    problems = {
        'uniform': [uniform, starting_point],
         'boltzmann': [boltzmann, starting_point_2],
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
        acceptance_rate[p], samples[p] = hopsy.sample(mcs, rngs, n_samples=n_samples, thinning=thinning, n_procs=4)
        ess[p] = hopsy.ess(samples[p])
        print('\tess', np.min(ess[p]))
        rhat[p] = hopsy.rhat(samples[p])
        print('\trhat', np.max(rhat[p]))
    #%%
    # Alle Verteilungen in einer Abbildung
    # Nur der Biomasse-Fluss
    plt.title('Biomasse')
    plt.hist(samples['boltzmann'][0, :, biomass_index], density = True, color = "green", label = "Boltzmann-Distribution")
    plt.hist(samples['uniform'][0, :, biomass_index], density = True, color = "blue", label = "Uniform-Distribution", alpha = 0.4)
    plt.hist(samples['normal'][0, :, biomass_index], density = True, color = "blue", label = "Normal-Distribution", alpha = 0.4)
    plt.xlabel('Biomass Flux')
    plt.ylabel('Probability')

    plt.legend(loc='upper right')

    plt.show()
    #%%
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
    #%%
    # Hier werden die Bilder für die einzelnen Flüsse geplottet
    # Dabei werden die einzelnen Flüsse alle in einer Abbildung dargestellt.
    # Je nachdem, welcher Fluss ausgewählt wurde, sind die Unterschiede sehr gering
    for i in lst_plot:
        ind = raw_polytope.A.columns.tolist().index(i)
        plt.title(i)
        plt.hist(samples['normal'][0, :, ind], density = True, color = "yellow", label = "Normal-Verteilung")
        plt.hist(samples['boltzmann'][0, :, ind], density = True, color = "green", label = "Boltzmann-Verteilung", alpha = 0.4)
        plt.hist(samples['uniform'][0, :, ind], density = True, label = "Gleichverteilung", alpha = 0.4)
        plt.legend(loc='upper right')
        plt.show()