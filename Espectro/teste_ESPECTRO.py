import numpy as np
import scipy.sparse as sp
from MatrixSpectralAnalysisSolve import MatrixSpectralAnalysisSolve
import scipy.io as sio
import os, sys
from LoadProblem import load_problem

save_path = './Projeto-IC/Espectro/results'
operation = 'solve' # 'spectre' ou 'solve'
problems = ['Esdras', 'Barreira', 'Barreira_f', 'SPE10_0', 'SPE10_0_true', 'SPE10_0_quasi', 'SPE10_0_ABF', 'SPE10_85', 'SPE10_85_true', 'SPE10_85_quasi', 'SPE10_85_ABF']
methods = {'spectre': ['Jacobi', 'Seidel', 'ILUfac', 'Multiscale', 'MultiscaleILUfac', 'MultiscaleJacobi', 'MultiscaleSeidel'],
           'solve': ['Jacobi', 'Seidel', 'ILUfac', 'Multiscale', 'MultiscaleILUfac', 'MultiscaleJacobi', 'MultiscaleSeidel', 'AMG']
           }
iterations = [100, 100, 100, 20, 20, 20, 20, 20, 20, 20, 20]


for problem in problems:
    try:
        os.mkdir(f'{save_path}/{problem}')
    except:
        print('', end='')

    if operation == 'spectre':
        try:
            os.mkdir(f'{save_path}/{problem}/Spectre')
        except:
            print('', end='')
    else:
        try:
            os.mkdir(f'{save_path}/{problem}/Residues')
        except:
            print('', end='')

    A, OP, OR, b = load_problem(problem)
    print(f'{problem}: {A.shape[0]}')
    sa = MatrixSpectralAnalysisSolve(A, OP, OR, b)

    if operation == 'spectre':
        print('Estrutura da Matriz')
        sa.VerifyMatrixStructure(problem, save_path)
        print('ok\n')
        print('Espectro da matriz de Transmissibilidade')
        sa.InitialSpectre(problem, save_path)
        print('ok\n')

        for method in methods[operation]:
            print(method)
            try:
                av_error = sa.GetSpectre(method)
            except Exception as error:
                print(f'Analise Espectral de {method} não pôde ser realizada devido a: {error}\n')
                continue

            print(len(av_error))
            print(f'Max: {max(abs(av_error))}\nMin: {min(abs(av_error))}\n')
            sa.PlotSpectre(method, av_error, save_path, problem)

    elif operation == 'solve':
        times = open(f'{save_path}/{problem}/Residues/times.txt', 'w', encoding='utf-8')
        res = {
            'Jacobi': [],
            'Seidel': [],
            'ILUfac': [],
            'Multiscale': [], 
            'MultiscaleILUfac': [], 
            'MultiscaleJacobi': [], 
            'MultiscaleSeidel':[],
            'AMG': []
            }
        times.write(f'{problem}: Tempos, Iterações e Convergência\n\n')
        for method in methods[operation]:
            print(method)
            try:
                res[method], x, time, info = sa.Solve(method, maxiter=iterations[problems.index(problem)])
            except Exception as error:
                print(f'{method} não pôde ser resolvido devido a: {error}\n')
                continue

            print(f'Tempo: {time}, Iterações: {len(res[method])}\n')
            times.write(f'{method}: {time}, {len(res[method])}, {info}\n')
        sa.PlotResidues(methods['solve'], res, save_path, problem)
        times.close()