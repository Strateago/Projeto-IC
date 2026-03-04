import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from MatrixSpectralAnalysis import MatrixSpectralAnalysis

problem = "Barreira"  # Barreira ou Joao

if problem != 'Barreira':
    data_path = './Projeto-IC/Espectro/dados'
    colsA = np.load(f'{data_path}/cols_M.npy')
    linesA = np.load(f'{data_path}/lines_M.npy')
    dataA = np.load(f'{data_path}/data_M.npy')
    colsOP = np.load(f'{data_path}/cols_op.npy')
    linesOP = np.load(f'{data_path}/lines_op.npy')
    dataOP = np.load(f'{data_path}/data_op.npy')
    colsOR = np.load(f'{data_path}/cols_or.npy')
    linesOR = np.load(f'{data_path}/lines_or.npy')
    dataOR = np.load(f'{data_path}/data_or.npy')

    b = np.load(f'{data_path}/rhs.npy')

    OR = sp.csc_matrix((dataOR, (linesOR, colsOR)))
    OP = sp.csc_matrix((dataOP, (linesOP, colsOP)))
    A = sp.csc_matrix((dataA, (linesA, colsA)))

else:
    data_path = './SPE10-PRESSURE-MATRICES/SPE10-PRESSURE-MATRICES/100_100_barreira'
    lines=np.load(f'{data_path}/OP1_fMsRSB_lines.npy')
    cols=np.load(f'{data_path}/OP1_fMsRSB_cols.npy')
    data=np.load(f'{data_path}/OP1_fMsRSB_data.npy')
    primal=np.load(f'{data_path}/primal_id_1_fMsRSB.npy')
    centroids=np.load(f'{data_path}/centroids.npy')
    l,c,d=np.load(f'{data_path}/A_lines.npy'),np.load(f'{data_path}/A_cols.npy'),np.load(f'{data_path}/A_data.npy')
    b=np.load(f"{data_path}/b_vector.npy")

    OP=sp.csc_matrix((data,(lines,cols)),shape=(lines.max()+1,cols.max()+1))
    OR=sp.csc_matrix((np.ones(len(primal)), (primal, np.arange(len(primal)))),shape=(primal.max()+1, len(primal)))
    A=sp.csc_matrix((d,(l,c)),shape=(l.max()+1,c.max()+1))

# n = A.shape[0]
# # M1_A=OP*spla.spsolve(RAP,RA)
# # I = sp.identity(n) - M1_A


# RAP = OR @ A @ OP
# RA = OR @ A
# solver = spla.factorized(RAP)
# def apply_G(x):
#     return x - OP @ solver(RA @ x)

# I = spla.LinearOperator((A.shape[0], A.shape[0]), matvec=apply_G)

# # Maior autovalor em magnitude (Largest Magnitude)
# max_val = spla.eigs(I, k=10, which='LM', return_eigenvectors=False)

# # Menor autovalor em magnitude (Smallest Magnitude)
# # O parâmetro sigma força o uso do modo shift-invert, que converge muito melhor
# min_val = spla.eigs(I, k=10, sigma=0.001, which='LM', return_eigenvectors=False)

# vals = np.concatenate((max_val, min_val))

# print(f"Maior: {max_val[0]}, Menor: {min_val[0]}")
# print(vals)

# plt.figure()
# x = np.arange(-1.0, 1.01, 0.01)
# x[-1] = 1.0
# z = np.sqrt(1 - x**2)
# y = -np.sqrt(1 - x**2)
# plt.plot(x, z, 'b')
# plt.plot(x, y, 'b')
# plt.plot(np.real(vals), np.imag(vals), '*r', label='Autovalores')
# plt.xlabel('Re(lambda)')
# plt.ylabel('Imag(lambda)')
# plt.title(f'Espectro')
# plt.axis('equal')
# plt.grid(True)
# plt.show()

sa = MatrixSpectralAnalysis(A, OP, OR, b)
# print('verify')
# sa.VerifyMatrixStructure()
# print('preconditioned matrix analysis')
# sa.PreconditionedMatrix_Analysis()
print('solve')
# methods: 'Jacobi', 'Seidel', 'ILUfac', 'Multiscale', 'MultiscaleILUfac', 'MultiscaleJacobi', 'MultiscaleSeidel'
method = 'Jacobi'
av_error, res = sa.Solve(method)
sa.PlotSpectre(method, av_error, "./Projeto-IC/Espectro/results", problem)
# sa.PlotResidues(method, res)