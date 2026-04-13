import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from SpectralAnalysisMethods import SpectralAnalysisMethods
from SolvingMethods import SolvingMethods
import time

class MatrixSpectralAnalysisSolve:
    def __init__(self, A, OP, OR, b = None):
        self._OR = OR
        self._OP = OP
        self._A = A
        if b is not None:
            self._b = b.flatten()
        else:
            self._b = b

        spectral_methods = SpectralAnalysisMethods(self._A, self._OP, self._OR)
        self._SpectralMethods = {'Jacobi': spectral_methods.Jacobi,
                        'Seidel': spectral_methods.Seidel,
                        'ILUfac': spectral_methods.ILUfac,
                        'Multiscale': spectral_methods.Multiscale,
                        'MultiscaleILUfac': spectral_methods.Multiscale_ILUfac,
                        'MultiscaleJacobi': spectral_methods.Multiscale_Jacobi,
                        'MultiscaleSeidel': spectral_methods.Multiscale_Seidel
                        }
        
        solve_methods = SolvingMethods(self._A, self._OP, self._OR)
        self._SolveMethods = {'Jacobi': solve_methods.Jacobi,
                        'Seidel': solve_methods.Seidel,
                        'ILUfac': solve_methods.ILUfac,
                        'Multiscale': solve_methods.Multiscale,
                        'MultiscaleILUfac': solve_methods.Multiscale_ILUfac,
                        'MultiscaleJacobi': solve_methods.Multiscale_Jacobi,
                        'MultiscaleSeidel': solve_methods.Multiscale_Seidel,
                        'AMG': solve_methods.AMG  
                        }

    def VerifyMatrixStructure(self, problem, save_path=None):
        # ----- Verificando a estrutura da matriz (TPFA) -----
        plt.figure()
        plt.spy(self._A)
        plt.title(f'Estrutura da Matriz: {problem}')
        if save_path:
            plt.savefig(f'{save_path}/{problem}/Esparsidade_{problem}')
        else:
            plt.show()
        plt.close()

    def InitialSpectre(self, problem, save_path=None):
        # ----- Análise espectral da matriz de transmissibilidade pré-condicionada -----
        nnod = self._A.shape[0]                # número de graus de liberdade
        print('Init 1')
        # 100 maior magnitude
        BM = spla.eigs(self._A, k=100, which='LM', return_eigenvectors=False)          # espectro da matriz de iteração
        print('ok\nInit 2')
        # 100 menor magnitude
        try:
            SM = spla.eigs(self._A, k=100, which='SM', return_eigenvectors=False)
        except spla.ArpackNoConvergence as error: # Captura especificamente o erro de convergência
            SM = error.eigenvalues
            print(f'SM não convergiu totalmente. Foram obtidos {len(SM)} autovalores.')
        except (Exception, KeyboardInterrupt) as error: # Captura outros erros (memória, álgebra, etc)
            SM = np.array([])
            print(f'SM falhou por outro motivo: {error}')
        print('ok')
        
        av_T = np.concatenate([BM, SM])

        print(f'Max: {max(abs(av_T))} Min: {min(abs(av_T))}\n')

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
        else:
            plt.show()
        plt.close()

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
        
    def Solve(self, method, x0=None, rtol=1e-8, maxiter=2500):
        nnod = self._A.shape[0]                # número de graus de liberdade
        lowerA = sp.tril(self._A, k=-1)    # matriz triangular inferior
        upperA = sp.triu(self._A, k=1)     # matriz triangular superior
        diagA  = sp.diags(self._A.diagonal(), format="csc")             # diagonal da matriz
        args = (nnod, upperA, lowerA, diagA)
        if method in self._SolveMethods:

            linOP = self._SolveMethods[method](args)
            residuals = []
            def callback(xk):
                r = self._b - self._A @ xk
                residuals.append(np.linalg.norm(r))
            
            init = time.time()

            # Solve the system with GMRES
            x, info = spla.gmres(
                self._A,
                self._b,
                x0=x0,
                M=linOP,
                rtol=rtol,
                maxiter=maxiter,
                callback=callback,
                callback_type="x"
            )

            end = time.time()

            print(info)
            return residuals, x, end-init
        
        
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
        else:
            plt.show()
        plt.close()
    
    def PlotResidues(self, methods, residues, save_path=None, problem=""):
        plt.figure()
        plt.xlabel("Iterações")
        plt.ylabel("Resíduo")
        plt.title(f'Resíduos: {problem}')
        plt.yscale("log")
        plt.grid(True)
        for method in methods:
            plt.plot(residues[method])
        
        plt.legend(methods, loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0.)
        plt.tight_layout()

        if save_path:
            plt.savefig(f'{save_path}/{problem}/Residues/Residues_{problem}')
        else:
            plt.show()
        plt.close()