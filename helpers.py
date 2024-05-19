from PolyRound.mutable_classes.polytope import Polytope
import pandas as pd
import os


def load_polytope(path, model):
    A = pd.read_csv(os.path.join('data', f'{model}_A.csv'), index_col=0)
    b = pd.read_csv(os.path.join('data', f'{model}_b.csv'), index_col=0)
    S = pd.read_csv(os.path.join('data', f'{model}_S.csv'), index_col=0)
    h = pd.read_csv(os.path.join('data', f'{model}_h.csv'), index_col=0)
    shift = pd.read_csv(os.path.join('data', f'{model}_shift.csv'), index_col=0)
    trafo = pd.read_csv(os.path.join('data', f'{model}_transformation.csv'), index_col=0)
    print(b.columns)
    poly = Polytope(A=A.to_numpy(), b=b['0'].to_numpy(), S=S.to_numpy(), h=h['0'].to_numpy())
    poly.shift = shift
    poly.transformation = trafo
    return poly
