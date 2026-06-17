import numpy as np
import scipy.sparse as sp
import scipy.io as sio
import sys

def load_problem(problem):
    if problem == 'Esdras':
        # Caso 1
        data = sio.loadmat("Projeto-IC/arquivos/Problema_heterogeneo_buck_30x30.mat")
        # Caso 2
        # data = sio.loadmat("Projeto-IC/arquivos/Problema_homogeneo_buck_40x40.mat")

        A = sp.csc_matrix(np.asarray(data["T"], dtype=np.float64))
        OP = sp.csc_matrix(np.asarray( data['OP'], dtype=np.float64))
        OR = sp.csc_matrix(np.asarray(data['OR'], dtype=np.float64))

        b = np.asarray(data["F"])

    elif problem == 'Joao':
        data_path = './Projeto-IC/Espectro/dados/Joao'
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

    elif problem == 'Barreira':
        data_path = './SPE10-PRESSURE-MATRICES/SPE10-PRESSURE-MATRICES/100_100_barreira'
        lines=np.load(f'{data_path}/OP1_MsRSB_lines.npy')
        cols=np.load(f'{data_path}/OP1_MsRSB_cols.npy')
        data=np.load(f'{data_path}/OP1_MsRSB_data.npy')
        primal=np.load(f'{data_path}/primal_id_1_MsRSB.npy')
        l,c,d=np.load(f'{data_path}/A_lines.npy'),np.load(f'{data_path}/A_cols.npy'),np.load(f'{data_path}/A_data.npy')
        b=np.load(f"{data_path}/b_vector.npy")

        OP=sp.csc_matrix((data,(lines,cols)),shape=(lines.max()+1,cols.max()+1))
        OR=sp.csc_matrix((np.ones(len(primal)), (primal, np.arange(len(primal)))),shape=(primal.max()+1, len(primal)))
        A=sp.csc_matrix((d,(l,c)),shape=(l.max()+1,c.max()+1))

    elif problem == 'Barreira_f':
        data_path = './SPE10-PRESSURE-MATRICES/SPE10-PRESSURE-MATRICES/100_100_barreira'
        lines=np.load(f'{data_path}/OP1_fMsRSB_lines.npy')
        cols=np.load(f'{data_path}/OP1_fMsRSB_cols.npy')
        data=np.load(f'{data_path}/OP1_fMsRSB_data.npy')
        primal=np.load(f'{data_path}/primal_id_1_fMsRSB.npy')
        l,c,d=np.load(f'{data_path}/A_lines.npy'),np.load(f'{data_path}/A_cols.npy'),np.load(f'{data_path}/A_data.npy')
        b=np.load(f"{data_path}/b_vector.npy")

        OP=sp.csc_matrix((data,(lines,cols)),shape=(lines.max()+1,cols.max()+1))
        OR=sp.csc_matrix((np.ones(len(primal)), (primal, np.arange(len(primal)))),shape=(primal.max()+1, len(primal)))
        A=sp.csc_matrix((d,(l,c)),shape=(l.max()+1,c.max()+1))

    elif problem == 'SPE10_0':
        data_path = './Projeto-IC/Espectro/dados/SPE10_0'
        A = (np.load(f'{data_path}/Jpp_0.npy', allow_pickle=True)).sum()  
        OR = (np.load(f'{data_path}/OR_0.npy', allow_pickle=True)).sum()
        OP = (np.load(f'{data_path}/OP_0.npy', allow_pickle=True)).sum()
        b = np.load(f'{data_path}/b_SPE-10_0.npy')

    elif problem == 'SPE10_0_true':
        data_path = './Projeto-IC/Espectro/dados/SPE10_0'
        A = (np.load(f'{data_path}/Jpp_CPR_true_IMPES_0.npy', allow_pickle=True)).sum()  
        OR = (np.load(f'{data_path}/OR_0.npy', allow_pickle=True)).sum()
        OP = (np.load(f'{data_path}/OP_0.npy', allow_pickle=True)).sum()
        b = np.load(f'{data_path}/bp_CPR_true_IMPES_0.npy')

    elif problem == 'SPE10_0_quasi':
        data_path = './Projeto-IC/Espectro/dados/SPE10_0'
        A = (np.load(f'{data_path}/Jpp_CPR_quasi_IMPES_0.npy', allow_pickle=True)).sum()  
        OR = (np.load(f'{data_path}/OR_0.npy', allow_pickle=True)).sum()
        OP = (np.load(f'{data_path}/OP_0.npy', allow_pickle=True)).sum()
        b = np.load(f'{data_path}/bp_CPR_quasi_IMPES_0.npy')

    elif problem == 'SPE10_0_ABF':
        data_path = './Projeto-IC/Espectro/dados/SPE10_0'
        A = (np.load(f'{data_path}/Jpp_CPR_ABF_0.npy', allow_pickle=True)).sum()  
        OR = (np.load(f'{data_path}/OR_0.npy', allow_pickle=True)).sum()
        OP = (np.load(f'{data_path}/OP_0.npy', allow_pickle=True)).sum()
        b = np.load(f'{data_path}/bp_CPR_ABF_0.npy')

    elif problem == 'SPE10_85':
        data_path = './Projeto-IC/Espectro/dados/SPE10_85'
        A = (np.load(f'{data_path}/Jpp_85.npy', allow_pickle=True)).sum()  
        OR = (np.load(f'{data_path}/OR_85.npy', allow_pickle=True)).sum()
        OP = (np.load(f'{data_path}/OP_85.npy', allow_pickle=True)).sum()
        b = np.load(f'{data_path}/b_SPE-10_85.npy')

    elif problem == 'SPE10_85_true':
        data_path = './Projeto-IC/Espectro/dados/SPE10_85'
        A = (np.load(f'{data_path}/Jpp_CPR_true_IMPES_85.npy', allow_pickle=True)).sum()  
        OR = (np.load(f'{data_path}/OR_85.npy', allow_pickle=True)).sum()
        OP = (np.load(f'{data_path}/OP_85.npy', allow_pickle=True)).sum()
        b = np.load(f'{data_path}/bp_CPR_true_IMPES_85.npy')

    elif problem == 'SPE10_85_quasi':
        data_path = './Projeto-IC/Espectro/dados/SPE10_85'
        A = (np.load(f'{data_path}/Jpp_CPR_quasi_IMPES_85.npy', allow_pickle=True)).sum()  
        OR = (np.load(f'{data_path}/OR_85.npy', allow_pickle=True)).sum()
        OP = (np.load(f'{data_path}/OP_85.npy', allow_pickle=True)).sum()
        b = np.load(f'{data_path}/bp_CPR_quasi_IMPES_85.npy')

    elif problem == 'SPE10_85_ABF':
        data_path = './Projeto-IC/Espectro/dados/SPE10_85'
        A = (np.load(f'{data_path}/Jpp_CPR_ABF_85.npy', allow_pickle=True)).sum()  
        OR = (np.load(f'{data_path}/OR_85.npy', allow_pickle=True)).sum()
        OP = (np.load(f'{data_path}/OP_85.npy', allow_pickle=True)).sum()
        b = np.load(f'{data_path}/bp_CPR_ABF_85.npy')

    else:
        print('problema não reconhecido.')
        sys.exit()

    return A, OP, OR, b