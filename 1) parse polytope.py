import matplotlib.pyplot as plt
import hopsy
import os
import sys
import numpy as np
from scipy.optimize import fsolve
import pandas as pd
import logging
logging.getLogger('sbml').setLevel(logging.ERROR)
from PolyRound.api import PolyRoundApi
import pickle
# 1) Lade die benötigten Dateien und lege die Paramter fest
# 2) Bestimme die durchschnittliche Wachstumsrate und die Varianz dieser.
# 3) Bestimme \beta aus der Boltzmann-Verteilung + Bestimme den
# Monte-Carlo-Fehler.
# 4) Stelle die Flüsse grafisch dar.
# %%
# Hier stehen alle änderbaren Variablen.
# Laden der Daten


if __name__ == "__main__":
    model_path = os.path.join("models", "iEZ481_Glc.xml")
    polytope = PolyRoundApi.sbml_to_polytope(model_path)
    # TODO: simplify, transform round
    biomass_index = polytope.A.columns.tolist().index('biomass_a')
    print('biomass index is', biomass_index)
    polytope.A.to_csv(os.path.join('data', 'A.csv'))
    polytope.b.to_csv(os.path.join('data', 'b.csv'))
    polytope.S.to_csv(os.path.join('data', 'S.csv'))
    polytope.h.to_csv(os.path.join('data', 'h.csv'))
    polytope.shift.to_csv(os.path.join('data', 'shift.csv'))
    polytope.transformation.to_csv(os.path.join('data', 'transformation.csv'))

