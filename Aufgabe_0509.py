import numpy as np
import matplotlib.pyplot as plt
import hopsy
import os
from PolyRound.api import PolyRoundApi
import time
#%%
# Aufgabe 1 
# Darstellen der Daten unter Nebenbedingungen
# Definiere A und b
A = np.array([[1., -1.], [-1., 0.], [0., -1.], [0., 1.]])
b = np.array([0., 0., 0., 5.])

# Anzahl der zu ziehenden Punkte
num_points = 100

# Listen zur Aufbewahrung der x-Werte, die die Bedingung erfüllen, 
# und der x-Werte, die es nicht tun
x_values_satisfy = []
x_values_do_not_satisfy = []

# Ziehe 100 verschiedene Werte für x
for _ in range(num_points):
    x = np.random.rand(A.shape[1])  # Zufällige Werte für x generieren
    
    # Überprüfe, ob die Bedingung Ax < b erfüllt ist
    if np.all(np.dot(A, x) < b):
        x_values_satisfy.append(x)
    else:
        x_values_do_not_satisfy.append(x)

# Konvertiere die Listen in NumPy-Arrays
x_values_satisfy = np.array(x_values_satisfy)
x_values_do_not_satisfy = np.array(x_values_do_not_satisfy)

# Plotte die x-Werte, die die Bedingung erfüllen, in Rot
plt.scatter(x_values_satisfy[:, 0], x_values_satisfy[:, 1], c='red', label='Ax < b')

# Plotte die x-Werte, die die Bedingung nicht erfüllen, in Blau
plt.scatter(x_values_do_not_satisfy[:, 0], x_values_do_not_satisfy[:, 1], c='blue', label='Ax >= b')

# Legende hinzufügen
plt.legend()

# Achsenbeschriftungen hinzufügen
plt.xlabel('X-Achse')
plt.ylabel('Y-Achse')

# Anzeigen des Plots
plt.show()
#%%
# Aufgabe 2
problem = hopsy.Problem([[1, 1], [-1, 0], [0, -1]], [1, 0, 0])
print(problem)
#%%
chains = [hopsy.MarkovChain(problem, starting_point=[.5, .5]) for i in range(4)]
rng = [hopsy.RandomNumberGenerator(seed= i) for i in range(4)]
#%%
start_time = time.time()
# accrate: acceptance rate von jedem State
# samples: Simulierte States
accrate, samples = hopsy.sample(chains, rng, n_samples=1000, thinning=10)
end_time = time.time()
print(end_time - start_time)
#%%
print(len(samples[0]))
print(samples[0,0:4])
#%%
# samples sind die möglichen x-Werte 
res_konvergenz = hopsy.rhat(samples)
print(max(res_konvergenz[0]))
#%%
res_hopsy = hopsy.add_box_constraints(problem, lower_bound=0, upper_bound=5)
#%%
# Aufgabe 4
model_path = os.path.join("C:/Users/carol/OneDrive/Desktop/Juelich", "e_coli_core.xml")
#%%
polytope = PolyRoundApi.sbml_to_polytope(model_path)
problem_e_coli = hopsy.Problem(polytope.A, polytope.b)
#%%

#%%
print(polytope.A.columns)
#%%
starting_point = hopsy.compute_chebyshev_center(problem_e_coli)
chains_e_coli = [hopsy.MarkovChain(problem_e_coli, starting_point = starting_point) for i in range(4)]
#%%
rng = [hopsy.RandomNumberGenerator(seed= i) for i in range(4)]
#%%
accrate_e_coli, samples_e_coli = hopsy.sample(chains_e_coli, rng, n_samples=1000, thinning=10)
res_konvergenz_e_coli = hopsy.rhat(samples_e_coli)
#%%
print(max(res_konvergenz_e_coli[0])) # Sinkt mit zunehmender Stichprobengröße
#%%
# Gucke für jede Dimension, ob die x-Werte im Poytop enthalten sind
plt.figure(dpi=300)

plt.title(polytope.A.columns[0])
counts, bins, _ = plt.hist(samples_e_coli[0, :,  0], bins=25, alpha=0.5, label='chain 1, parameter 1', density=True)
plt.hist(samples_e_coli[1, :,  0], bins=bins, alpha=0.5, label='chain 2, parameter 1', density=True)
plt.hist(samples_e_coli[2, :,  0], bins=bins, alpha=0.5, label='chain 3, parameter 1', density=True)
plt.hist(samples_e_coli[3, :,  0], bins=bins, alpha=0.5, label='chain 4, parameter 1', density=True)

#%%
counts, bins, _ = plt.hist(samples_e_coli[0, :,  1], bins=25, alpha=0.5, label='chain 1, parameter 2', density=True)
plt.hist(samples_e_coli[1, :,  1], bins=bins, alpha=0.5, label='chain 2, parameter 2', density=True)
plt.hist(samples_e_coli[2, :,  1], bins=bins, alpha=0.5, label='chain 3, parameter 2', density=True)
plt.hist(samples_e_coli[3, :,  1], bins=bins, alpha=0.5, label='chain 4, parameter 2', density=True)


plt.legend()
plt.xlabel('parameter value [a. u.]')
plt.ylabel('density')
plt.show()
#%%
# Aufgabe 5
# Rundet die Werte dem Paper entsprechend
res = hopsy.round(problem)
print(res)



