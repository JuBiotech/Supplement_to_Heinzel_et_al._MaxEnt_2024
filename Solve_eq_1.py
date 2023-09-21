import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
# Mit Hilfe von ChatGPT geschrieben
#%%
# Achtung: Es kann sein, dass die Lösung hier nicht konvergiert!!!
# Außerdem: hierbie ist v immer nur eindimensional, normalerweise ist v ja 
# hochdimensional

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
    print(res_final)
    return res_final

# Definiere Konstante c
c =  0.99

# Schätzen Sie einen Startwert für beta
initial_guess = 1  # Hier können Sie einen beliebigen Startwert verwenden

# Verwenden Sie fsolve, um die Gleichung zu lösen und den Wert für beta zu finden
result = fsolve(equation_to_solve, initial_guess, args=(c,))

# Das Ergebnis enthält den gefundenen Wert für beta
print(f"Der gefundene Wert für beta ist: {result[0]}")
