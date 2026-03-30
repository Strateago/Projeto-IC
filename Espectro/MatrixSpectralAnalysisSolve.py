import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from SpectralAnalysisMethods import SpectralAnalysisMethods
from SolvingMethods import SolvingMethods

class MatrixSpectralAnalysisSolve:
    def __init__(self, A, OP, OR, b):
        self._OR = OR
        self._OP = OP
        self._A = A
        self._b = b

        spectral_methods = SpectralAnalysisMethods(self._A, self._b, self._OP, self._OR)
        self._SpectralMethods = {'Jacobi': spectral_methods.Jacobi,
                        'Seidel': spectral_methods.Seidel,
                        'ILUfac': spectral_methods.ILUfac,
                        'Multiscale': spectral_methods.Multiscale,
                        'MultiscaleILUfac': spectral_methods.Multiscale_ILUfac,
                        'MultiscaleJacobi': spectral_methods.Multiscale_Jacobi,
                        'MultiscaleSeidel': spectral_methods.Multiscale_Seidel
                        }
        
        solve_methods = SolvingMethods(self._A, self._b, self._OP, self._OR)
        self._SolveMethods = {'Jacobi': solve_methods.Jacobi,
                        'Seidel': solve_methods.Seidel,
                        'ILUfac': solve_methods.ILUfac,
                        'Multiscale': solve_methods.Multiscale,
                        'MultiscaleILUfac': solve_methods.Multiscale_ILUfac,
                        'MultiscaleJacobi': solve_methods.Multiscale_Jacobi,
                        'MultiscaleSeidel': solve_methods.Multiscale_Seidel
                        }

    def VerifyMatrixStructure(self):
        # ----- Verificando a estrutura da matriz (TPFA) -----
        plt.figure()
        plt.subplot(1, 2, 1)
        plt.spy(self._A, markersize=1)
        plt.title("Estrutura A")

        plt.show()

    def PreconditionedMatrix_Analysis(self, save_path, problem):
        # ----- Análise espectral da matriz de transmissibilidade pré-condicionada -----
        nnod = self._A.shape[0]                # número de graus de liberdade
        av_T = spla.eigs(self._A, which='LM', k=100, return_eigenvectors=False)  # espectro da matriz Tfina e espectro da matriz Tfina
        av_T = np.append(av_T, spla.eigs(self._A, which='SM', k=100, return_eigenvectors=False))
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
        # diff = self._A - self._A.T
        # print(np.abs(diff.data).max())

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
        if save_path:
            plt.savefig(f'{save_path}/{problem}/Spectre_Original_{problem}')

    def GetSpectre(self, method):
        # ----- Análise espectral da matriz de transmissibilidade pré-condicionada -----
        nnod = self._A.shape[0]                # número de graus de liberdade
        lowerA = sp.tril(self._A, k=-1)    # matriz triangular inferior
        upperA = sp.triu(self._A, k=1)     # matriz triangular superior
        diagA  = sp.diags(self._A.diagonal(), format="csc")             # diagonal da matriz

        args = (nnod, upperA, lowerA, diagA)
        if method in self._SpectralMethods:
            return self._SpectralMethods[method](args)
        else:
            raise ValueError(f"Método '{method}' não reconhecido.")
        
    def Solve(self, method):
        nnod = self._A.shape[0]                # número de graus de liberdade
        lowerA = sp.tril(self._A, k=-1)    # matriz triangular inferior
        upperA = sp.triu(self._A, k=1)     # matriz triangular superior
        diagA  = sp.diags(self._A.diagonal(), format="csc")             # diagonal da matriz

        args = (nnod, upperA, lowerA, diagA)
        if method in self._SolveMethods:
            return self._SolvegMethods[method](args)
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
            plt.savefig(f'{save_path}/{problem}/Spectre/Spectre_{method}_{problem}')
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