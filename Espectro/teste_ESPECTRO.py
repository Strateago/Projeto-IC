import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from MatrixSpectralAnalysis import MatrixSpectralAnalysis

problem = "Joao"  # Barreira ou Joao

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

sa = MatrixSpectralAnalysis(A, OP, OR, b)
# print('verify')
# sa.VerifyMatrixStructure()
# print('preconditioned matrix analysis')
# sa.PreconditionedMatrix_Analysis()
print('solve')
# methods: 'Jacobi', 'Seidel', 'ILUfac', 'Multiscale', 'MultiscaleILUfac', 'MultiscaleJacobi', 'MultiscaleSeidel'
method = 'ILUfac'
av_error, res = sa.Solve(method)
print(len(av_error))
print(av_error[:100])
sa.PlotSpectre(method, av_error, "./Projeto-IC/Espectro/results", problem)
# sa.PlotResidues(method, res)