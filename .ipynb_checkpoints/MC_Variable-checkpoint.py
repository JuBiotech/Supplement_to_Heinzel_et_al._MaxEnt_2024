import numpy as np
from scipy import integrate
#%%
# Zweidimensionales Integral mit MC-Methode bestimmen




def func1(x):
    return x[0] + x[1]
  
def mc_integrate(func, a, b, dim, n = 1000):
    # Monte Carlo integration of given function over domain from a to b (for each parameter)
    # dim: dimensions of function
    
    x_list = np.random.uniform(a, b, (n, dim))
    print(x_list)
    f = []
    for i in range(n):
        x = np.random.uniform(a[0], b[0])
        y = np.random.uniform(a[1], 2)
        f.append(func1([x,y]))
    
    summe = 0
    for i in f:
        summe+= i
    y_mean =  summe/len(f)
    flaeche = 2 # Muss noch angepasst werden an Integral!!!
    integ = y_mean * flaeche # * 2
    
    return integ

# Examples
print(f"Monte Carlo solution for : {mc_integrate(func1, [0,0], [2,1], 2, 1000): .3f}")
