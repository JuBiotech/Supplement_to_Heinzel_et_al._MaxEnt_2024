import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
#%%
import csv
#%%
model_path = os.path.join("C:/Users/carol/OneDrive/Desktop/Juelich", "e_coli_core.xml")
#%%
# Bestimme die Wachstumsrate (d.h. die Anzahl der neuen Zellen pro Minute durch
# exponentiellen Fit an die Daten)

with open("C:\\Users\\carol\OneDrive\\Desktop\\Juelich\\maxentflux\\growth_rate_ds.csv", newline='') as csvfile:
    # Erstellen Sie ein CSV-Leserobjekt
    csvreader = pd.read_csv(csvfile)

#%%
# Aufteilen der Daten
# Hier: Beispielhaft mit einem Medium
Acetate = csvreader[csvreader['medium'] == 'Acetate-MOPS']
Citrat = csvreader[csvreader['medium'] == 'Citrat-MOPS']
#%%
Gluconate = csvreader[csvreader['medium'] == 'Gluconate-MOPS']
#%%
PCA = csvreader[csvreader['medium'] == 'PCA-Gluc']
BHI = csvreader[csvreader['medium'] == 'BHI']
Fructose = csvreader[csvreader['medium'] == 'Fructose-MOPS']
Glucose = csvreader[csvreader['medium'] == 'Glucose-MOPS']
#%%
Pyruvate = csvreader[csvreader['medium'] == 'Pyruvate-MOPS']
#%%
# Daten für Regression
data = Pyruvate["cell_count"]
print(data)
#%%
laenge_teilliste = max(Pyruvate["frame"])
print(laenge_teilliste)
#%%
# Daten für die Zeit erstellen
maximaler_wert = (laenge_teilliste + 1) * 15 

t = []

# Schleife, um Vielfache von 15 hinzuzufügen
for i in range(0, maximaler_wert, 15):
    t.append(i)

# Teillisten initialisieren
teillisten = []

# Schleife, um die Ausgangsliste in Teillisten aufzuteilen
def teile_liste(input_list, n):
    # Initialisiere eine leere Ergebnisliste
    output_list = []

    # Schleife, um die Eingabeliste in Teillisten aufzuteilen
    for i in range(0, len(input_list), n):
        teil_liste = input_list[i:i+n]
        output_list.append(teil_liste)

    return output_list


# Funktion aufrufen
teillisten = teile_liste(data, laenge_teilliste + 1)
#%%
# Die aufgeteilten Teillisten anzeigen
for i, teilliste in enumerate(teillisten):
    if(i == 0):
        print(f"Teilliste {i + 1}: {teilliste}")
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
    # Schätzen der Parameter c und r mit Regression
    params, covariance = curve_fit(exponential_growth, t, N, maxfev = 100000)
    # Extrahieren der geschätzten Parameter
    r_estimated = params
    mw.append(r_estimated)
print(np.mean(mw))
#%%
# Daten für die CSV-Datei
daten = [
    ("Acetate-MOPS", 0.007320106743488314),
("Citrat-MOPS", 0.011970883236652241),
("Gluconate-MOPS", 0.010259351950089916),
("PCA-Gluc",),
("BHI", 0.008654428418461922),
("Fructose-MOPS", 0.009408669046910657), 
("Glucose-MOPS", 0.0094969666209908),
("Pyruvate-MOPS", )]

# Name der CSV-Datei
dateiname = "Wachstumsrate.csv"

# CSV-Datei erstellen und schreiben
with open(dateiname, mode='w', newline='') as datei:
    schreiber = csv.writer(datei)
    
    # Spaltennamen (optional)
    schreiber.writerow(["Medium", "durchschnitlliche Wachstumsrate"])
    
    # Daten schreiben
    schreiber.writerows(daten)

print(f"Die Datei '{dateiname}' wurde erfolgreich erstellt.")
#%%
# Plotten
y_fit = []
for i in t:
    y_fit.append(5 * np.exp(0.00717829 * i))
# Plot erstellen
plt.plot(t, teillisten[0], marker='o')  # Linie mit Punkten
plt.plot(t, y_fit)
plt.show()
c = teillisten[0][0]