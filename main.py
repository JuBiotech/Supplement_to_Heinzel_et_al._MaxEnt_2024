import hopsy
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
#%%
# 1. Bestimme A und b
# 2. Bestimme \beta aus der Boltzmann-Verteilung
# 3. Bestimme die Lösung der Gleichung Ax \leq b mit hopsy.
#%%
# To do
# Bestimme S, S_{xch}, C_{eq, net}, C_{eq, xch} bei echten Beispiel
# Bestimmem lambda(v) bei echten Beispiel 
###############################################################################
# Finden von A und b, wenn nur S bekannt
##############################################################################
#%%
# 1. S eingeben
S = np.zeros((5,3))
#%%
# 2. S_{xch}, C_{eq, net}, C_{eq, xch} eingeben
S_xch = np.zeros((5,3))
C_eq = np.zeros((5,3))
C_eq_xch = np.zeros((5,3))
#%%
# 3. p bestimmen
print(len(S_xch), S_xch.shape[1])
#%%
m_nullen = np.zeros((len(S_xch), S_xch.shape[1]))
# nebeneinander
S_neben = np.concatenate((S, m_nullen), axis=1)
S_xch_neben = np.concatenate(( m_nullen, S_xch), axis=1)
C_eq_neben = np.concatenate((C_eq, m_nullen), axis=1)
C_eq_xch_neben = np.concatenate((m_nullen, C_eq_xch), axis=1)
print(S_neben, S_xch_neben, C_eq_neben, C_eq_xch_neben)
#%%
# untereinander
A_temp =  np.concatenate((S_neben, S_xch_neben, C_eq_neben, C_eq_xch_neben), axis = 0) 
#%%
b_temp = np.ones((4*len(S_xch)))
b_temp = np.array(b_temp)
A_temp = np.array(A_temp)
print(b_temp, A_temp)
print(len(A_temp), len(b_temp))
#%%
p, residuals, rank, singular_values = np.linalg.lstsq(A_temp, b_temp, rcond=None)
#%%
# 4. A, b bestimmen
C_in = np.ones((20,6))
k = np.ones((6,4))
b_in = np.ones(20)
A = C_in @ k
b = b_in - C_in @ p
print(A, b)
#%%
# Beispielwerte
A, b = [[1, 1, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1]], [1, 0, 0, 0]
#%%
# Bestimmung des Parameters der Boltzman-Verteilung a
# Achtung: Es kann sein, dass die Lösung hier nicht konvergiert!!!

# Definieren der Funktion lambda(v) - Abhängig von v und beta
# Wachstumsrate lambda(v) = v1 + v2 (hier)
def lambda_function(v1, v2, beta):
    return (v1 + v2) * np.exp(v1*beta)* np.exp(v2*beta)  

# Definieren der Funktion, die integriert wird - Abhängig von beta
def equation_to_solve(beta, c):
    return c - quad_integrand(beta)

# Definieren der Funktion, die das Integral berechnet - Abhängig von beta
def quad_integrand(beta):
    
    # Definieren der Integrationsgrenzen (anpassen, wie benötigt)
    x_lower = 0  # Untere Grenze für v
    x_upper = 1  # Obere Grenze für v
    
    # Berechnen des Integrals unter Verwendung von lambda_function mit beta
    result, _ = quad(lambda v2: quad(lambda v1: lambda_function(v1, v2, beta), x_lower, x_upper)[0], x_lower, x_upper)

    # normierender Faktor
    res, _ = quad(lambda v2: quad(lambda v1: np.exp(v1 * beta) * np.exp(v2 * beta), x_lower, x_upper)[0], x_lower, x_upper)
    res_final = result/res
    # Rückgabewert des Integranden
    return res_final

# Definiere Konstante c
# Das ist der Mittelwert der Wachstumsraten

c =  0.00045141

# Schätze Startwert für beta
initial_guess = 1  

result = fsolve(equation_to_solve, initial_guess, args=(c,)) #-3071.7200534253125

print(f"Der gefundene Wert für beta ist: {result[0]}")

#%%
# a ist die erwartete Wachstumsrate
# Verwende Gleichung 11.34 aus Kapitel 11 des Buchs
# x ist lambda(v) und eindimensional
class Boltzman_Modell:
    def __init__(self, a):
        self.a = a
    def compute_negative_log_likelihood(self, x):
        n = len(x)
        lambda_f = 0
        for i in range(n):
            lambda_f += x[i] 
        self.a = self.a * lambda_f
        return -self.a * lambda_f
a = np.array(result)
model_Boltzmann = Boltzman_Modell(a)
#%%
# set up the problem
problem = hopsy.Problem(A, b, model_Boltzmann)
starting_point = hopsy.compute_chebyshev_center(problem)

# number of chains we will use
n_chains = 4

# set up a few parallel markov chains using a gaussian proposal and starting in the origin
chains = [hopsy.MarkovChain(problem, starting_point=starting_point)
          for i in range(n_chains)]


# set up the random number generators
rngs = [hopsy.RandomNumberGenerator(seed=42, stream=i) for i in range(n_chains)]
#%%
accrate, states = hopsy.sample(chains, rngs, n_samples=1000, thinning=10)

rhat = hopsy.rhat(states, series=10)

print('Acceptance rates (per chain):', *accrate)
#%%
fig, ax = plt.subplots(1, 1, figsize=(12,7))

ax.plot(np.linspace(10, 1000, 100), rhat[:,0])
ax.plot(np.linspace(10, 1000, 100), rhat[:,1])
ax.plot(np.linspace(10, 1000, 100), rhat[:,2])

ax.plot([0, 1000], [1.05, 1.05], linestyle='dashed', color='gray',
        label='Convergence threshold')

ax.set_xlabel(r"$n_{samples}$", fontsize=14)
ax.set_ylabel(r"$\hat{R}$", fontsize=14)
ax.legend(fontsize=14)

plt.show()