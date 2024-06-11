import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy
import seaborn as sns
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
    #print(i_max, x_max)
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
    #print(p_final, q)
    res_final = scipy.stats.wasserstein_distance(p_final, q)
    # res_final = kullback_leibler_divergence(p_final, q)
    return res_final

values1 = [1, 2, 3, 4, 10, -10]  # p_final
values2 = [0.3, 0.1, 9, 100, 9]  # q

res = create_intervals(values1, values2)
print("Result:", res)
#%%
def kl_data():
    
    res_diff_b = []
    res_diff_n = []
    res_diff_nb = []
    data_u = np.load("C:\\Users\\carol\\Downloads\\iEZ481_Citrat-MOPS_boltzmann_samples_thinned.npz", allow_pickle=True)
    data_b = np.load("C:\\Users\\carol\\Downloads\\iEZ481_Glucose-MOPS_boltzmann_samples_thinned.npz", allow_pickle=True)
    #data_n = np.load("C:\\Users\\carol\\Downloads\\iEZ481_PCA_Gluc_boltzmann_samples.npz", allow_pickle=True)
    #data_u = np.load("C:\\Users\\carol\\Downloads\\iEZ481_Glucose-MOPS_samples.npz", allow_pickle=True).astype(np.float32)
        #data_b = np.load("C:\\Users\\carol\\Downloads\\iEZ481_Glucose-MOPS_gauss_samples.npz", allow_pickle=True).astype(np.float32)
    data_n = np.load("C:\\Users\\carol\\Downloads\\iEZ481_Glucose-MOPS_boltzmann_samples_thinned.npz", allow_pickle=True)
 
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
res = kl_data()
# 
# 0: ([0.00022801218161683283], [0.0006059922680412372], [0.0006059922680412372])

#%%
data_u = np.load("C:\\Users\\carol\\Downloads\\iEZ481_Citrat-MOPS_boltzmann_samples_thinned.npz", allow_pickle=True)
data_b = np.load("C:\\Users\\carol\\Downloads\\iEZ481_PCA_Gluc_boltzmann_samples_thinned.npz", allow_pickle=True)
data_n = np.load("C:\\Users\\carol\\Downloads\\iEZ481_Glucose-MOPS_boltzmann_samples_thinned.npz", allow_pickle=True)
#%% 
res_diff_b = []
res_diff_n = []
res_diff_nb = []
for i in range(0, 558):
        u = data_u['samples'][0, :, i] # Cit.astype(np.float32)
        b = data_b['samples'][0, :, i]#.astype(np.float32)
        #n = data_n['samples'][0,:,i]# .astype(np.float32)
        print(i)
        #res_temp_un = create_intervals(u, n)
        #print(i)
        #res_temp_bn = create_intervals(b, n)
        #print(i)
        res_temp_ub = create_intervals(b, u)
        #print(res_temp_ub, res_temp_un, res_temp_bn)
        res_diff_b.append(res_temp_ub)
        #res_diff_n.append(res_temp_un)
        #res_diff_nb.append(res_temp_bn)
#%%
with open('res_diff_b.txt', 'w') as file:
    for item in res_diff_b:
        file.write(f"{item}\n")
#%%
with open('res_diff_n.txt', 'r') as file:
    read_values = file.readlines()

res_diff_n = [float(value.strip()) for value in read_values]

#%%
# Entferne das Newline-Zeichen und konvertiere die Strings zurück in Floats
res_diff_nb = [float(value.strip()) for value in read_values]

#%%
data = {
    'Glucose': [0, 0.1, 0.5],
    'PCA': [0.1, 0, 0.2],
    'Citrat': [0.5, 0.2, 0]
}
index = ['GLucose', 'PCA', 'Citrat']
df = pd.DataFrame(data, index=index)

plt.figure(figsize=(8, 6))  
#sns.heatmap(df, cmap='coolwarm', cbar = True, cbar_kws={'label': 'Wasserstein-Distance', 'fontsize':30})  
plt.title('Flux i', fontsize = 10) 
#sns.heatmap(df, cmap='coolwarm', cbar=True, cbar_kws={'label': 'Wasserstein-Distance'})  
#plt.xlabel("Growth Medium", fontsize = 30)
#plt.ylabel("Growth Medium", fontsize = 30)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
heatmap = sns.heatmap(df, cmap='coolwarm', cbar=True, cbar_kws={'label': 'Wasserstein-Distance'})  
cbar = heatmap.collections[0].colorbar
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_ylabel('Wasserstein-Distance', fontsize=20)
#plt.show()  
plt.savefig(f'plot_{i}.png')  # Speichere jede Kurve mit einem unterschiedlichen Namen
#%%
def plotten(a, b, c, titles):
    # Vorlage für das Dictionary
    template_data = {
        'Glucose': [0, 0.1, 0.5],
        'PCA': [0.1, 0, 0.2],
        'Citrat': [0.5, 0.2, 0]
        }

    
    # Erzeuge für jeden Index i ein neues Dictionary mit den Werten von a, b und c an den entsprechenden Positionen
    for i, title in enumerate(titles):
        new_data = {
            'Glucose': template_data['Glucose'][:],  # Kopiere die Werte
            'PCA': template_data['PCA'][:],
            'Citrat': template_data['Citrat'][:]
         }
        
        new_data['Glucose'][1] = a[i]  
        new_data['Glucose'][2] = b[i]  
        new_data['PCA'][0] = a[i]     
        new_data['PCA'][2] = c[i]     
        new_data['Citrat'][1] = c[i]     
        new_data['Citrat'][0] = b[i]   
        print(new_data)
        index = ['Glucose', 'PCA', 'Citrat']
        df = pd.DataFrame(new_data, index=index)
        plt.figure(figsize=(8, 6))  
        plt.title(title, fontsize=10)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        heatmap = sns.heatmap(df, cmap='coolwarm', cbar=True, cbar_kws={'label': 'Wasserstein-Distance'})  
        cbar = heatmap.collections[0].colorbar
        cbar.ax.tick_params(labelsize=10)
        cbar.ax.set_ylabel('Wasserstein-Distance', fontsize=10)
        #plt.show()  
        plt.savefig(f'WD{i}.png')  
        plt.close()# Speichere jede Kurve mit einem unterschiedlichen Namen
#%%
# Maximale Elemente
res = res_diff_n.index(max(res_diff_n))
res1 = res_diff_b.index(max(res_diff_b))
res2 = res_diff_nb.index(max(res_diff_nb))

#%%
names = fba['Unnamed: 0'].astype(str).tolist()[278:279]

plotten(res_diff_nb, res_diff_n, res_diff_b,names)