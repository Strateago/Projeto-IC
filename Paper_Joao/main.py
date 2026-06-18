import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from pyamg.krylov import fgmres
import os
import time
import ilupp
from load_data import load_data

threshold = 1e-12
fill = 15 # ILU1
tol = 1e-4
restart = 100

data_path = './Projeto-IC/Paper_Joao/dados'
# data = []
# for dir in os.listdir(f'{data_path}'):
#     data.append(load_data(f'{data_path}/{dir}'))
#     ## continuar em breve, depois de tentar com um único

problem = "layer_1_Cr30"
data = load_data(f'{data_path}/{problem}')

M = ilupp.ILUTPreconditioner(data['T'], fill_in=fill, threshold=threshold)
M = ilupp.ILU0Preconditioner(data['T'])
residuals = []
sol, code = fgmres(data['T'], data['b'], M=M, tol=tol, residuals=residuals, restart=restart)

print(f"Code: {code}")
print(f"Residuals: {residuals}")
print(f"Solution: {sol}")
