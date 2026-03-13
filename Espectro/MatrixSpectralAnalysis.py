import numpy as np
import scipy.io as sio
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
import time
from SpectralAnalysisMethods import SpectralAnalysisMethods

class MatrixSpectralAnalysis:
    def __init__(self, A, OP, OR, b):
        self._OR = OR
        self._OP = OP
        self._A = A
        self._b = b

        meth = SpectralAnalysisMethods(self._A, self._b, self._OP, self._OR)
        self._methods = {'Jacobi': meth.Jacobi,
                        'Seidel': meth.Seidel,
                        'ILUfac': meth.ILUfac,
                        'Multiscale': meth.Multiscale,
                        'MultiscaleILUfac': meth.Multiscale_ILUfac,
                        'MultiscaleJacobi': meth.Multiscale_Jacobi,
                        'MultiscaleSeidel': meth.Multiscale_Seidel
                        }

    def VerifyMatrixStructure(self):
        # ----- Verificando a estrutura da matriz (TPFA) -----
        plt.figure()
        plt.subplot(1, 2, 1)
        plt.spy(self._A, markersize=1)
        plt.title("Estrutura A")

        plt.show()

    def PreconditionedMatrix_Analysis(self):
        # ----- Análise espectral da matriz de transmissibilidade pré-condicionada -----
        nnod = self._A.shape[0]                # número de graus de liberdade
        av_T, VA = spla.eigs(self._A)  # espectro da matriz Tfina e espectro da matriz Tfina
        # coordxA = np.arange(1, nnod + 1)
        # coordyA = np.arange(1, nnod + 1)
        # X, Y = np.meshgrid(coordxA, coordyA)
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.plot_surface(X, Y, np.real(VA), cmap='viridis')
        # ax.set_xlabel('x')
        # ax.set_ylabel('y')
        # ax.set_title('Autovetores da matriz de transmissibilidade')
        # plt.show()
        diff = self._A - self._A.T
        print(np.abs(diff.data).max())

        plt.figure()
        # Círculo unitário
        x = np.arange(-1.0, 1.01, 0.01)
        x[-1] = 1.0
        z = np.sqrt(1 - x**2)
        y = -np.sqrt(1 - x**2)
        plt.plot(x, z, 'b')
        plt.plot(x, y, 'b')
        # autovalores
        plt.plot(np.real(av_T), np.imag(av_T), '*r')
        plt.xlabel('Re(lambda)')
        plt.ylabel('Imag(lambda)')
        plt.title('Espectro da matriz de transmissibilidade')
        plt.axis('equal')
        plt.grid(True)
        plt.show()

    def Solve(self, method):
        # ----- Análise espectral da matriz de transmissibilidade pré-condicionada -----
        nnod = self._A.shape[0]                # número de graus de liberdade
        upperA = sp.tril(self._A, k=-1)    # matriz triangular inferior
        lowerA = sp.triu(self._A, k=1)     # matriz triangular superior
        diagA  = sp.diags(self._A.diagonal(), format="csc")             # diagonal da matriz
        print(nnod)

        args = (nnod, upperA, lowerA, diagA)
        if method in self._methods:
            return self._methods[method](args)
        else:
            raise ValueError(f"Método '{method}' não reconhecido.")
        
        
    def PlotSpectre(self, method, av_error_operator, save_path=None, problem = ""):
        plt.figure()
        x = np.arange(-1.0, 1.01, 0.01)
        x[-1] = 1.0
        z = np.sqrt(1 - x**2)
        y = -np.sqrt(1 - x**2)
        plt.plot(x, z, 'b')
        plt.plot(x, y, 'b')
        plt.plot(np.real(av_error_operator), np.imag(av_error_operator), '*r', label='Autovalores')
        plt.xlabel('Re(lambda)')
        plt.ylabel('Imag(lambda)')
        plt.title(f'Espectro: pré-condicionador {method}')
        plt.axis('equal')
        plt.grid(True)
        if save_path:
            plt.savefig(f'{save_path}/Spectre_{method}_{problem}')
        # plt.show()
    
    def PlotResidues(self, method, residues):
        plt.figure()
        plt.plot(residues)
        plt.xlabel("Iteração")
        plt.ylabel("Resíduo")
        plt.title(f'Resíduo {method}')
        plt.yscale("log")
        plt.grid(True)
        plt.show()