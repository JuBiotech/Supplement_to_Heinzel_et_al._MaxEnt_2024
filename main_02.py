import matplotlib.pyplot as plt
import hopsy
import os
import numpy as np
from scipy.optimize import fsolve
import pandas as pd
from scipy.optimize import curve_fit
from PolyRound.api import PolyRoundApi
#%%
# Dieser Code ist folgendermaßen gegliedert:
# 1. Bestimmung der durchschnittlichen Wachstumsrate mit den Daten
# 2. Bestimmung von beta mit Glg (11.35) aus dem Buch
# 3. Lösen von Ax \leq b mit hopsy und der Annahme dass v boltzmann-verteilt ist.
# 4. Ansatz zum Lösen von Ax \leq b mit hopsy und der Annahme, dass v normalverteilt ist.
#%%
model_path = os.path.join("C:/Users/carol/OneDrive/Desktop/Juelich", "e_coli_core.xml")
#%%
# Bestimme die Wachstumsrate (d.h. die Anzahl der neuen Zellen pro Minute durch
# exponentiellen Fit an die Daten)

with open("C:\\Users\\carol\OneDrive\\Desktop\\Juelich\\maxentflux\\growth_rate_ds.csv", newline='') as csvfile:
    # Erstellen Sie ein CSV-Leserobjekt
    csvreader = pd.read_csv(csvfile)
    print(csvreader)

#%%
# Aufteilen der Daten
# Hier: Beispielhaft mit einem Medium
Acetate = csvreader[csvreader['medium'] == 'Acetate-MOPS']
print(Acetate)
#%%
# Daten für Regression
data = Acetate["cell_count"]
print((data))
#%%
# Aufteilen nach Kammern
# Länge der Teillisten
laenge_teilliste = 40

# Teillisten initialisieren
teillisten = []

# Schleife, um die Ausgangsliste in Teillisten aufzuteilen
for i in range(0, len(data), laenge_teilliste):
    teilliste = data[i:i+laenge_teilliste]
    teillisten.append(teilliste)

# Die aufgeteilten Teillisten anzeigen
for i, teilliste in enumerate(teillisten):
    print(f"Teilliste {i + 1}: {teilliste}")
#%% 
# Daten für die Zeit erstellen
maximaler_wert = 39 * 15 + 1  

t = []

# Schleife, um Vielfache von 15 hinzuzufügen
for i in range(0, maximaler_wert + 1, 15):
    t.append(i)
print(t)
#%%
y_fit = []
for i in t:
    y_fit.append(5 * np.exp(0.00717829 * i))
# Plot erstellen
plt.plot(t, teillisten[0], marker='o')  # Linie mit Punkten
plt.plot(t, y_fit)
plt.show()
c = teillisten[0][0]
#%%
def exponential_growth(t, r):
    return c * np.exp(r * t)
#%%
# Epxonentieller Fit
mw = []
# Funktion für das exponentielle Wachstumsmodell
for i in teillisten:
    c = i.tolist()[0]
    N = i.tolist()
    # Schätzen der Parameter c und r mit Regression
    params, covariance = curve_fit(exponential_growth, t, N, maxfev = 10000)
    # Extrahieren der geschätzten Parameter
    r_estimated = params
    mw.append(r_estimated)
print(np.mean(mw))
#%%
###############################################################################
# Bestimme den Wert für beta
###############################################################################
polytope = PolyRoundApi.sbml_to_polytope(model_path)
problem_e_coli = hopsy.Problem(polytope.A, polytope.b)
problem_e_coli = hopsy.round(problem_e_coli)

starting_point = hopsy.compute_chebyshev_center(problem_e_coli)
chains_e_coli = [hopsy.MarkovChain(problem_e_coli, starting_point = starting_point) for i in range(1)]

rng = [hopsy.RandomNumberGenerator(seed= i) for i in range(1)]

accrate_e_coli, samples_e_coli = hopsy.sample(chains_e_coli, rng, n_samples=1000, thinning=10)

# Die Wachstumsrate ist in der 25. Spalte
samples = samples_e_coli[:,:,24][0]
print(samples.shape)
# Histogramm anzeigen lassen
#%%
res_samples = []
# teile die Werte durch 100: funktioniert durch Einheiten-Wechsel
# Ziel: Python kann die Werte wieder ausrechnen
for i in samples:
    res_samples.append(i/100)
#%%
# Bestimmung der zu lösenden Gleichung
# Wachstumsrate lambda(v) = v (hier)
def lambda_function(v1, beta):
    res = 0
    res1 = 0
    
    for i in v1:
        res+= i * np.exp(i*beta)
        res1 += np.exp(i*beta)
    return res/res1 

# Definieren der Funktion, die integriert wird - Abhängig von beta
def equation_to_solve(beta, c):
    return c - quad_integrand(beta)

# Definieren der Funktion, die das Integral berechnet - Abhängig von beta
def quad_integrand(beta):
    result = lambda_function(res_samples, beta)
    return result
# Definiere Konstante c
# Hier die Wachstumsrate von den Zellen in einem bestimmten Medium
# Kann problematisch werden bei sehr großen oder sehr kleinen Werten!
c =  0.1 

# Schätzen Sie einen Startwert für beta
# beta ist (hier) normalerweise negativ
initial_guess = 0

# Verwende fsolve, um die Gleichung zu lösen und den Wert für beta zu finden
result = fsolve(equation_to_solve, initial_guess, args=(c,))

# Das Ergebnis enthält den gefundenen Wert für beta
print(f"Der gefundene Wert für beta ist: {result[0]}")
#%%
#################################################################################
# Verwende die Boltzmann-Verteilung 
###############################################################################
# Verwende Gleichung 11.34 aus Kapitel 11 des Buchs
# x ist lambda(v) und eindimensional
class Boltzman_Modell:
    def __init__(self, a):
        self.a = a
    def compute_negative_log_likelihood(self, x):
        #a = self.a * np.ones(len(x))
        return - self.a *x[24] #(np.dot(a, x))
a = np.array(result)
model_Boltzmann = Boltzman_Modell(a)
#%%
problem = hopsy.Problem(polytope.A, polytope.b, model_Boltzmann)
problem = hopsy.round(problem)
starting_point = hopsy.compute_chebyshev_center(problem)
#%%

# number of chains we will use
n_chains = 4

# set up a few parallel markov chains using a gaussian proposal and starting in the origin
chains = [hopsy.MarkovChain(problem, starting_point=starting_point)
          for i in range(n_chains)]

# set up the random number generators
rngs = [hopsy.RandomNumberGenerator(seed=42, stream=i) for i in range(n_chains)]
#%%
accrate, states = hopsy.sample(chains, rngs, n_samples=100, thinning=10)
#%%
# Ausdünnung, damit rhat() anwendbar ist
# Anzahl zu entfernender Elemente
n = 99000

# Entferne die ersten n Elemente aus jeder inneren Liste
result_list = [outer_list[n:] for outer_list in states]

print(result_list)
#%%
y = []
for i in range(1):
    for j in range(99000,100000):
        y.append(states[3][j][:])
print(y)
        #%%
rhat = hopsy.rhat(result_list, series=10)
print(rhat)
#%%
print('Acceptance rates (per chain):', *accrate)
#%%
fig, ax = plt.subplots(1, 1, figsize=(12,7))

ax.plot(np.linspace(10, 1000, 100), rhat[:,0])
#%%
ax.plot(np.linspace(10, 1000, 100), rhat[:,1])
ax.plot(np.linspace(10, 1000, 100), rhat[:,2])


ax.plot([0, 1000], [1.05, 1.05], linestyle='dashed', color='gray',
        label='Convergence threshold')

ax.set_xlabel(r"$n_{samples}$", fontsize=14)
ax.set_ylabel(r"$\hat{R}$", fontsize=14)
ax.legend(fontsize=14)

plt.show()
#%%
print(len(states[0][0]))
#%%
y = []
for i in range(1):
    for j in range(100000):
        y.append(states[3][j][1])
#%%
x = np.linspace(90000,100000,10000)
plt.plot(x,y[90000:100000])

plt.show()
#%%
# Verwende die Gauß-Verteilung
class GaussianModel2:
    def __init__(self, mu, cov):
        self.mu = mu
        self.cov = cov

# ACHTUNG: Hierbei wurde der Name (ursprünglich log_density) geändert!
# Hier ist das auch nicht die normale Likelihood-Funktion und nicht das
# negative der Likelihood-Funktion
    def compute_negative_log_likelihood(self, x):
        return (
            0.5
            * (x.reshape(-1, 1) - self.mu).T
            @ np.linalg.inv(self.cov)
            @ (x.reshape(-1, 1) - self.mu)
        )[0, 0]

    def hessian(self, x):
        return np.linalg.inv(self.cov)

    def grad_log_density(self, x):
        return -np.linalg.inv(self.cov) @ (x - self.mu)

# Hierbei: Noch herausfinden, wie mu und cov zu bestimmen sind!
mu = np.array(res_mu)
cov = np.array(res_cov)
#%%
# Bestimme den Monte-Carlo-Fehler
# = Standardabweichung durch Effektive Stichprobengröße

data = np.zeros((4,2,3))
print(hopsy.ess(data))