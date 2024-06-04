import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy

#%%

def compute_probabilities_uniform(samples, bins='auto'):
    # Berechne Histogramm und Wahrscheinlichkeitsverteilung
    hist, bin_edges = np.histogram(samples, bins=bins, density=True)
    prob = hist * np.diff(bin_edges)
    return prob, bin_edges 

def create_p(intervals, values):
    # Definiere die Intervallgrenzen
    x_min = min(values)
    x_max = max(values)
    i_min = min(intervals)
    i_max = max(intervals)
    print(i_max, x_max)
    if x_min < i_min:  
        intervals = [x_min-1] + intervals 
    if x_max > i_max:
        intervals = intervals + [x_max+1]
    bin_edges = intervals
    categories = pd.cut(values, bins=bin_edges, labels=False, right=False)
    array = categories
    counts = np.bincount(array)
    counts_list = counts.tolist()
    return counts_list, bin_edges


def create_intervals(values1, values2):
    p, bin_edges = compute_probabilities_uniform(values1, bins='auto')
    q, bin_edges_q = create_p(bin_edges.tolist(), values2)
    hist, bin_edges = np.histogram(values1, bins=bin_edges_q, density=True)
    p_final = hist * np.diff(bin_edges)
    q = [x / sum(q) for x in q]
    print(p_final, q)
    res_final = scipy.stats.wasserstein_distance(p_final, q)
    return res_final

values1 = [1, 2, 3, 4, 10, -10]  # p_final
values2 = [0.3, 0.1, 9, 100, 9]  # q

res = create_intervals(values1, values2)
print("Result:", res)
#%%
def kl_data(medium):
    
    res_diff_b = []
    res_diff_n = []
    res_diff_nb = []
    if(medium == "PCA"):
        data_u = np.load("C:\\Users\\carol\\Downloads\\iEZ481_PCA_Gluc_samples.npz", allow_pickle=True)
        data_b = np.load("C:\\Users\\carol\\Downloads\\iEZ481_PCA_Gluc_gauss_samples.npz", allow_pickle=True)
        data_n = np.load("C:\\Users\\carol\\Downloads\\iEZ481_PCA_Gluc_boltzmann_samples.npz", allow_pickle=True)

    
    elif(medium == "Glucose"):
        data_u = np.load("C:\\Users\\carol\\Downloads\\iEZ481_Glucose-MOPS_samples.npz", allow_pickle=True).astype(np.float32)
        data_b = np.load("C:\\Users\\carol\\Downloads\\iEZ481_Glucose-MOPS_gauss_samples.npz", allow_pickle=True).astype(np.float32)
        data_n = np.load("C:\\Users\\carol\\Downloads\\iEZ481_Glucose-MOPS_boltzmann_samples.npz", allow_pickle=True).astype(np.float32)
 
    for i in range(1, 10):
        u = data_u['samples'][0, :, i].astype(np.float32)
        b = data_b['samples'][0, :, i].astype(np.float32)
        n = data_n['samples'][0,:,i].astype(np.float32)
        print(i)
        res_temp_un = create_intervals(u, n)
        print(i)
        res_temp_bn = create_intervals(b, n)
        print(i)
        res_temp_ub = create_intervals(b, n)

        res_diff_b.append(res_temp_ub)
        res_diff_n.append(res_temp_un)
        res_diff_nb.append(res_temp_bn)
    
    return res_diff_n, res_diff_b, res_diff_nb
# 15:05
res = kl_data("PCA")
# 
# 0: ([0.00022801218161683283], [0.0006059922680412372], [0.0006059922680412372])