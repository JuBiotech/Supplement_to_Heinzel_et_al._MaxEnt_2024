import numpy as np
# Testen der Dimensionen
#%%
x = np.linspace(1,10,10)
a = 2*np.linspace(1,10,10)

print(x)
print(np.dot(a,x))
mu = np.zeros((5, 1))
cov = 0.1 * np.identity(5)
x = np.linspace(1,10,5)
test =  (x.reshape(-1, 1) - mu).T @ np.linalg.inv(cov) @ (x.reshape(-1, 1) - mu)
print(test)