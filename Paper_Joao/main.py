import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from pyamg.krylov import fgmres
import os
import time
import ilupp
from load_data import load_data
import matplotlib.pyplot as plt

threshold = 1e-6
fill = 2 # ILU1
tol = 1e-4
restart = 10 

result_path = './Projeto-IC/Paper_Joao/resultados'
data_path = './Projeto-IC/Paper_Joao/dados'

problem = "layer_85_Cr30"
data = load_data(f'{data_path}/{problem}')

## Resultado fina

# M = ilupp.ILUTPreconditioner(data['T'], fill_in=fill, threshold=threshold)
# # M = ilupp.ILU0Preconditioner(data['T'])
# residuals_ref = []
# sol_ref, code_ref = fgmres(data['T'], data['b'], M=M, tol=tol, residuals=residuals_ref, restart=restart)
# sol_ref = spla.spsolve(data['T'], data['b'])

## Resultado Multiescala

plt.figure()
plt.xlabel("Iterações")
plt.ylabel("Resíduo")
plt.yscale("log")
plt.grid(True)
problems = []
layer = 'layer_85'
mode = "Sem Pré-condicionador"
plt.title(f'Resíduos: Camada 85 - {mode}')

for dir in [d for d in os.listdir(f'{data_path}') if d.startswith(layer)]:
    print(dir)
    data = load_data(f'{data_path}/{dir}')

    Ac = data['OR'] @ data['T'] @ data['OP']
    b_c = data['OR'] @ data['b']
    
    if mode == "ILU0":
        M = ilupp.ILU0Preconditioner(Ac)
    elif mode == "ILU1":
        M = ilupp.ILUTPreconditioner(Ac, fill_in=fill, threshold=threshold)
    else:
        M = None

    # residuals_ms = []
    # sol, code = fgmres(Ac, b_c, M=M, tol=tol, residuals=residuals_ms, restart=restart)
    # sol = spla.spsolve(Ac, b_c)

    # sol_f = data['OP'] @ sol

    # # print(f"Solutions:\n   ref: {sol_ref}\n   ms: {sol_f}")

    # erro_absoluto = np.linalg.norm(sol_ref - sol_f)
    # erro_relativo = erro_absoluto / np.linalg.norm(sol_ref)

    # print(f"Erro Relativo: {erro_relativo:.4f} ({erro_relativo * 100:.2f}%)")


    init = time.time()

    # linOP = ilupp.ILU0Preconditioner(Ac)
    residuals = []
    def callback(xk):
        residuals.append(xk)

    # Solve the system with GMRES
    x, info = spla.gmres(
        Ac,
        b_c,
        x0=None,
        M=M,
        rtol=tol,
        maxiter=500,
        callback=callback,
        callback_type="pr_norm"
    )

    end = time.time()

    print(info)
    # sol_f = data['OP'] @ x

    # erro_absoluto = np.linalg.norm(sol_ref - sol_f)
    # erro_relativo = erro_absoluto / np.linalg.norm(sol_ref)

    # print(f"Erro Relativo (gmres): {erro_relativo:.4f} ({erro_relativo * 100:.2f}%)\n")

    # print(residuals, '\n')

    plt.plot(residuals, label=dir)
#     problems.append(dir)

plt.legend(problems, loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0.)
plt.tight_layout()
plt.savefig(f'{result_path}/{layer}_{mode}')
plt.close()

# print(f"Code (Reference): {code_ref}")
# print(f"Code (Multiscale): {code}")
# print(f"Residuals (Reference): {residuals_ref}")
# print(f"Residuals (Multiscale): {residuals_ms}")