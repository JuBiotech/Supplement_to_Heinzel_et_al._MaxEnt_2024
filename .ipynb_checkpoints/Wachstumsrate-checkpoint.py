# Johanna-Daten 
# 1) Einlesen der Dateien

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#%%
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
print(Acetate[0,0])
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
print(teillisten[0])
print(teillisten[0].tolist())
#%% 
# Daten für die Zeit erstellen
maximaler_wert = 39 * 15 + 1  

t = []

# Schleife, um Vielfache von 15 hinzuzufügen
for i in range(0, maximaler_wert + 1, 15):
    t.append(i)
print(t)
#%% 
# Plot erstellen
plt.plot(t, teillisten[0], marker='o')  # Linie mit Punkten
#%%
print(teillisten[0])
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
    print(N, c)


    # Schätzen der Parameter c und r mittels Kurvenanpassung (Regression)
    params, covariance = curve_fit(exponential_growth, t, N, maxfev = 10000)

    # Extrahieren der geschätzten Parameter
    r_estimated = params
    mw.append(r_estimated)
#%%
print(np.mean(r_estimated)) # 0.00045140650004922713
#%%
print(t)
a = [5, 5, 5, 5, 6, 7, 7, 7, 7, 9, 13, 15, 15, 16, 17, 20, 21, 26, 29, 32, 33, 39, 44, 51, 54, 60, 72, 79, 88, 99, 109, 127, 149, 164, 183, 209, 240, 282, 317, 355]
params, covariance = curve_fit(exponential_growth, t, N, maxfev = 10000)
print(params)
