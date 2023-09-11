# Finden von A und b, wenn nur S bekannt
import numpy as np

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
A_temp =  np.concatenate((S_neben, S_xch_neben, C_eq_neben, C_eq_xch_neben), axis = 0) #, C_eq_neben, C_eq_xch), axis=0)
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