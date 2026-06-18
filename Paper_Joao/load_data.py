import scipy.sparse as sp
import numpy as np
import os

def load_data(path):
    T_bc = sp.load_npz(os.path.join(path, 'T_bc.npz'))
    b = np.load(os.path.join(path, 'b.npy'))
    OP = sp.load_npz(os.path.join(path, 'OP.npz'))
    OR = sp.load_npz(os.path.join(path, 'OR.npz'))
    resp = {
        'T': T_bc,
        'b': b,
        'OP': OP,
        'OR': OR
    }
    return resp