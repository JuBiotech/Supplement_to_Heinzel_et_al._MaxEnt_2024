# Code zum Ziehen von Stichproben aus der Verteilung p(v)
# Dabei wurde die Boltzmann-Verteilung angenommen und der Parameter willkürlich festgelegt
#%%
import numpy as np
import sympy as sp
import random
import matplotlib as plt
#%%
# Definieren Sie die Variable v1_f als ein Symbol
v1_f = sp.symbols('v1_f')

# Definieren Sie den Wert von beta
beta = -1  # Hier als Beispiel, passen Sie den Wert entsprechend an

# Definieren Sie die PDF p(v1_f)
pdf = sp.exp(beta * v1_f) / sp.integrate(sp.exp(beta * v1_f), (v1_f, 0, 1))

# Definieren Sie die kumulative Verteilungsfunktion (CDF) durch Integration der PDF
cdf = sp.integrate(pdf, (v1_f, 0, v1_f))

# Anzahl der gezogenen Stichproben
num_samples = 100  # Passen Sie die Anzahl der gewünschten Stichproben an

# Zufällige Stichproben ziehen
samples = []
for _ in range(num_samples):
    # Zufallszahl zwischen 0 und 1
    random_value = random.uniform(0, 1)

    # Umkehrung der CDF, um den entsprechenden v1_f-Wert zu erhalten
    v1_f_sample = sp.solve(cdf - random_value, v1_f)
    
    # Fügen Sie den gezogenen Wert zur Liste der Stichproben hinzu
    samples.append(v1_f_sample[0])  # Nehmen Sie den ersten Wert, da solve eine Liste zurückgibt

# Ausgabe der gezogenen Stichproben
print(samples)
#%%
samples = [0.147380036743440, 0.0540822142731394, 0.999879680594621, 
           0.107788184043952, 0.594620301742961, 0.822835022249123, 0.611848968727389,
           0.697769611778075, 0.603004511760966, 0.140424908482988, 0.251376861468225, 
           0.0799603887042078, 0.322062603327989, 0.299707313932355, 0.114075001481386, 
           0.256679633454397, 0.149228783545939, 0.914652196270490, 0.146646531384128, 0.356042221859317, 0.331433655064382, 0.890304509470818, 0.720391764932515, 0.237249655347153, 0.0211324563984073, 0.419647034043081, 0.848663361100685, 0.677919922134448, 0.480914761389795, 0.235507447692828, 0.0918144817916133, 0.120339472121649, 0.809639772996281, 0.937432836265768, 0.572591191511411, 0.601280791324250, 0.498485558470597, 0.0329803509486135, 0.612049896799789, 0.0568228877589671, 0.0339566771815310, 0.283921410427936, 0.0803892604130801, 0.841929232045997, 0.532246494161576, 0.0743371997933838, 0.0868804880739043, 0.0720900451492925, 0.722310526244219, 0.760363768364144, 0.786002557451158, 0.444159450426055, 0.421850208056275, 0.815909999097476, 0.568150711567202, 0.300538496291470, 0.458969860776317, 0.0918103835732172, 0.789603042170578, 0.677229811244146, 0.0895395614857491, 0.742126096078126, 0.00190579260457833, 0.477677626492141, 0.123822524418407, 0.192774004391013, 0.674352415338685, 0.195130819481774, 0.122329669591432, 0.969339736842660, 0.566568264370418, 0.640300339199465, 0.991580899187068, 0.100446093916373, 0.0133612717298907, 0.0379692027294137, 0.364743667393094, 0.140564107490978, 0.000550086805217004, 0.166187694484824, 0.402240677387671, 0.940447338643657, 0.875622851077420, 0.200377857771489, 0.832504144305973, 0.401876087894089, 0.473882547270131, 0.0469491290278582, 0.973067907072403, 0.955841709167890, 0.0126103031997879, 0.939265466336955, 0.194276993616378, 0.279099295625713, 0.280055464373670, 0.282242270029679, 0.200313698403094, 
           0.784043980775645, 0.339243142538404, 0.0151965372362280,
           0.514289381701970, 0.307040397902008, 0.0984662023265144, 0.0367250836644514, 0.295339351474670, 0.546426928466977, 0.281527440613403, 0.0237822516943620, 0.945289540682205, 0.904887004574822, 0.372087466906939, 0.0610199942320535, 0.882082310165517, 0.128480261164835, 0.308709746826254, 0.211613523195055, 0.107761348116962, 0.304391872285052, 0.289623101034930, 0.107572917776230, 0.618394054527986, 0.111502810053117, 0.827560558182963, 0.665282732256902, 0.921116775940903, 0.625028378379252, 0.277284618879249, 0.596736131021726, 0.115856878420043, 0.243545660625942, 0.0414384458101002, 0.281044848305175, 0.0834240243728915, 0.464858481306431, 0.444699502874308, 0.237565319747413, 0.0138916909630280, 0.490741115252824, 0.752981709544315, 0.734209771649482, 0.209700877106794, 0.728399471244575, 0.790690907963643, 0.0547779445249474, 0.369954963466248, 0.664860032511478, 0.188201094818049, 0.461504339902188, 0.492136203748982, 0.594052645349257, 0.191944573206294, 0.493720819126074, 0.475558124985732, 0.333206567592400, 0.492607979937748, 0.837348623928624, 0.995783076888123, 0.226705682757267, 0.446500912299883, 0.00963882340938300, 0.800712644961803, 0.115975799203976, 0.658044338333557, 0.0506989132707947, 0.0735804890902556, 0.283207331155465, 0.867686179865916,
           0.560130743175425, 0.832605036412597, 0.0248238482657817, 0.705580135228056, 0.431129061162673, 0.245372800562155, 0.669636967833639, 0.879245052263325, 0.805515744576687, 0.191605974153243, 0.247670279845990, 0.743793381273504, 0.0251953314710561, 0.0361220283982172, 0.594970026610553, 0.0215912826169646, 0.965241039200020, 0.654292482862188, 0.0180393981232815, 0.565580444846110, 0.677866245579459, 0.332279202527921, 0.292440194453523, 0.693560097833211, 0.142600572400268, 0.143170919650808, 0.891604128972891, 0.618247325270482, 0.502675035055534, 0.494142535585255, 0.156203436037543, 0.00140915022669043, 0.465518226222929]
# Grafische Darstellung der Stichproben als Histogramm
plt.hist(samples, bins=10, density=True, alpha=0.7, color='blue', label='Stichproben')

# Plot der PDF
v1_f_values = np.linspace(0, 1, 10)
pdf_values = [float(pdf.subs(v1_f, value)) for value in v1_f_values]
plt.plot(v1_f_values, pdf_values, color='red', label='PDF')

# Beschriftungen hinzufügen
plt.xlabel('v1_f')
plt.ylabel('Dichte')
plt.legend()

# Anzeigen des Plots
plt.show()

