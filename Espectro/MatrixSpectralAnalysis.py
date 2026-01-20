import numpy as np
import scipy.io as sio
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
import time
from SpectralAnalysisMethods import SpectralAnalysisMethods

class MatrixSpectralAnalysis:
    def __init__(self, arquivo_mat):
        data = sio.loadmat(arquivo_mat)
        self._Tfine = np.asarray(data["T"])       # matriz de transmissibilidade na malha fina
        self._Tcorse = np.asarray(data["Tcorse"]) # matriz de transmissibilidade na malha grossa
        self._b = np.asarray(data["F"]).ravel()   # vetor dos carregamentos externos na malha fina
        self._OP = np.asarray(data['OP'])
        self._OR = np.asarray(data['OR'])
        meth = SpectralAnalysisMethods(self._Tfine, self._Tcorse, self._b, self._OP, self._OR)
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
        plt.spy(self._Tfine, markersize=1)
        plt.title("Estrutura Tfine")

        plt.subplot(1, 2, 2)
        plt.spy(self._Tcorse, markersize=1)
        plt.title("Estrutura Tcorse")
        plt.show()

    def PreconditionedMatrix_Analysis(self):
        # ----- Análise espectral da matriz de transmissibilidade pré-condicionada -----
        nnod = self._Tfine.shape[0]                # número de graus de liberdade
        av_T, VTfine = la.eig(self._Tfine)         # espectro da matriz Tfina e espectro da matriz Tfina

        coordxTfine = np.arange(1, nnod + 1)
        coordyTfine = np.arange(1, nnod + 1)
        X, Y = np.meshgrid(coordxTfine, coordyTfine)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, np.real(VTfine), cmap='viridis')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('Autovetores da matriz de transmissibilidade')
        plt.show()

        plt.figure()
        plt.plot(np.real(av_T), np.imag(av_T), '*r')
        plt.xlabel('Re(lambda)')
        plt.ylabel('Imag(lambda)')
        plt.title('Espectro da matriz de transmissibilidade')
        # Círculo unitário
        x = np.arange(-1.0, 1.01, 0.01)
        z = np.sqrt(1 - x**2)
        y = -np.sqrt(1 - x**2)
        plt.plot(x, z, 'b')
        plt.plot(x, y, 'b')
        plt.axis('equal')
        plt.grid(True)
        plt.show()

    def Solve(self, method):
        # ----- Análise espectral da matriz de transmissibilidade pré-condicionada -----
        nnod = self._Tfine.shape[0]                # número de graus de liberdade
        upperTfine = np.tril(self._Tfine, k=-1)    # matriz triangular inferior
        lowerTfine = np.triu(self._Tfine, k=1)     # matriz triangular superior
        diagTfine  = np.diag(np.diag(self._Tfine)) # matriz diagonal

        args = (nnod, upperTfine, lowerTfine, diagTfine)
        if method in self._methods:
            return self._methods[method](args)
        else:
            raise ValueError(f"Método '{method}' não reconhecido.")
        
    def PlotSpectre(self, method, av_error_operator):
        plt.figure()
        plt.plot(np.real(av_error_operator), np.imag(av_error_operator), '*r', label='Autovalores')
        plt.xlabel('Re(lambda)')
        plt.ylabel('Imag(lambda)')
        plt.title(f'Espectro: pré-condicionador {method}')
        x = np.arange(-1.0, 1.01, 0.01)
        z = np.sqrt(1 - x**2)
        y = -np.sqrt(1 - x**2)
        plt.plot(x, z, 'b')
        plt.plot(x, y, 'b')
        plt.axis('equal')
        plt.grid(True)
        plt.show()
    
    def PlotResidues(self, method, residues):
        plt.figure()
        plt.plot(residues)
        plt.xlabel("Iteração")
        plt.ylabel("Resíduo")
        plt.title(f'Resíduo {method}')
        plt.yscale("log")
        plt.grid(True)
        plt.show()