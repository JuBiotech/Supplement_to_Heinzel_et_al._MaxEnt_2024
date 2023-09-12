import numpy as np
import matplotlib.pyplot as plt

# Definiere Grenzen des Integrationsintervalls [a, b]
a, b = 0, 1

# Definiere die Funktion f(x), die Sie integrieren möchten
def f(x):
    return x**2

# Anzahl der zufälligen Punkte
num_points = 10000

# Zufällige x-Werte im Integrationsintervall generieren
x_values = np.random.uniform(a, b, num_points)

# Zufällige y-Werte zwischen 0 und max(f(x)) generieren
max_function_value = np.max(f(x_values))
y_values = np.random.uniform(0, max_function_value, num_points)

# Berechnen der Funktionswerte für die zufälligen x-Werte
function_values = f(x_values)

# Punkte zählen, die zwischen 0 und f(x) oder zwischen f(x) und 0 liegen
points_under_function = ((y_values >= 0) & (y_values <= function_values)) | ((y_values <= 0) & (y_values >= function_values))

# Schätzen Sie das Integral
integral_estimate = (b - a) * max_function_value * np.sum(points_under_function) / num_points

# Das geschätzte Integral ausgeben
print("Geschätztes Integral:", integral_estimate)

# Grafische Darstellung 
plt.figure(figsize=(6, 6))
plt.scatter(x_values, function_values, s=5, alpha=0.5, label='Punkte im Bereich')
plt.scatter(x_values[points_under_function], y_values[points_under_function], s=5, alpha=0.5, color='red', label='Punkte unter/über Funktion')
plt.xlim(a, b)
plt.ylim(0, max_function_value)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Monte-Carlo-Schätzung des Integrals')
plt.legend()
plt.grid(True)
plt.show()
