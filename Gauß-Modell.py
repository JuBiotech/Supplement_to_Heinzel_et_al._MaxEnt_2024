import hopsy
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
#%%
# Anpassen:
# lambda(v)
# Polytop
# Boltzman-Verteilung i.A.
# Zusammenhang S und A
#%%
# set up the polytope inequality for a 3-dimensional simplex
A, b = [[1, 1, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1]], [1, 0, 0, 0]
#%%
def lambda_(v, beta):
    return np.exp(v*beta)  
#%%
# Bestimmung des Parameters der Boltzman-Verteilung a
# Achtung: Es kann sein, dass die Lösung hier nicht konvergiert!!!

# Definieren der Funktion lambda(v) - Abhängig von v und beta
# Wachstumsrate lambda(v) = v (hier)
def lambda_function(v1, v2, beta):
    return v1 * np.exp(v1*beta) * v2 * np.exp(v2*beta)  # Hier verwenden wir sowohl v als auch beta

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

    # result, _ = quad(lambda v1: lambda_function(v1, beta), lower_limit, upper_limit)
    res, _ = quad(lambda v2: quad(lambda v1: np.exp(v1 * beta) * np.exp(v2 * beta), x_lower, x_upper)[0], x_lower, x_upper)
    res_final = result/res
    # Rückgabewert des Integranden
    return res_final

# Definiere Konstante c
# Das ist der Mittelwert der Wachstumsraten
c =  0.99

# Schätze Startwert für beta
initial_guess = 1  

result = fsolve(equation_to_solve, initial_guess, args=(c,))

print(f"Der gefundene Wert für beta ist: {result[0]}")

#%%
# a ist die erwartete Wachstumsrate
# Verwende Gleichung 11.34 aus Kapitel 11 des Buchs
# Noch anzupassen an die Gleichungen!!!
class Boltzman_Modell:
    def __init__(self, a):
        self.a = a
    def compute_negative_log_likelihood(self, x):
        self.a = self.a * np.ones(len(x))
        return -(np.dot(self.a, x))
a = np.array(result)
model_Boltzmann = Boltzman_Modell(a)

#%%
# Mathematiker-Version
class Boltzman_Modell:
    def __init__(self, a, temp):
        self.a = a
        self.temp = temp
# ACHTUNG: Hierbei wurde der Name (ursprünglich log_density) geändert!
# Ist aber garnicht die log likelihood-Funktion!
    def compute_negative_log_likelihood(self, x):
        for i in self.a:
            if(i > 0):
                self.temp+= 1
        if(self.temp == len(self.a)):
            return -(
                np.log(x.reshape(-1, 1))
                - x.reshape(-1, 1)**2/(2*self.a**2)
                - np.log(self.a**3) 
                )[0, 0]
        # [0,0] sorgt nur dafür, dass der Output die richtige Dimension hat
        else:
            return -100
a = np.array((1,2))
model_Boltzmann = Boltzman_Modell(a, 0)
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