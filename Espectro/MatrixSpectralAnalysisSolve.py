import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from SpectralAnalysisMethods import SpectralAnalysisMethods
from SolvingMethods import SolvingMethods
import time

class MatrixSpectralAnalysisSolve:
    """
    Classe para análise espectral e resolução de sistemas lineares usando diferentes métodos de pré-condicionamento.

    Métodos disponíveis:
        - Jacobi: Implementação do método de Jacobi para pré-condicionamento.
        - Seidel: Implementação do método de Gauss-Seidel para pré-condicionamento.
        - ILUfac: Implementação da fatoração ILU0 para pré-condicionamento.
        - Multiscale: Implementação do método multiscale para pré-condicionamento.
        - MultiscaleILUfac: Implementação do método multiscale em conjunto com ILU0 (Multiplicativo) para pré-condicionamento.
        - MultiscaleJacobi: Implementação do método multiscale em conjunto com Jacobi (Multiplicativo) para pré-condicionamento.
        - MultiscaleSeidel: Implementação do método multiscale em conjunto com Gauss-Seidel (Multiplicativo) para pré-condicionamento.
    
    Methods:
        VerifyMatrixStructure (str, str | None): Verifica e plota a estrutura da matriz A.
        InitialSpectre (str, str | None): Realiza a análise espectral inicial da matriz A.
        GetSpectre (str): Obtém o espectro da matriz A usando o método especificado.
        Solve (str, List[float] | None, float, int): Resolve o sistema linear usando o método especificado.
        PlotSpectre (str, List[float], str | None, str): Plota o espectro do pré-condicionador especificado de um dos métodos.
        PlotResidues (List[str], List[np.ndarray], str | None, str): Plota os resíduos das iterações para os métodos especificados.
        
    """

    def __init__(self, A, OP, OR, b = None):
        """
        Inicializa a classe com as matrizes A, OP, OR e o vetor b.

        Args:
            A (scipy.sparse.csc_matrix): Matriz de transmissibilidade.
            OP (scipy.sparse.csc_matrix): Matriz de prolongamento.
            OR (scipy.sparse.csc_matrix): Matriz de restrição.
            b (numpy.ndarray | None): Vetor do lado direito do sistema linear. Se None, assume-se que não há vetor b.
        """

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
        """
        Verifica e plota a estrutura da matriz A.
        
        Args:
            problem (str): Problema para o qual a estrutura da matriz será verificada.
            save_path (str | None): Caminho para salvar o gráfico da estrutura da matriz. Se None, o gráfico será exibido na tela.
        """

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
        """
        Gera o espectro parcial inicial da matriz A, incluindo os 100 maiores e 100 menores autovalores.
        Atenção: Menores autovalores podem não convergir, então podem ter menos de 100 menores.

        Args:
            problem (str): Problema para o qual o espectro será gerado.
            save_path (str | None): Caminho para salvar o gráfico do espectro. Se None, o gráfico será exibido na tela.
        """

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
        """
        Gera o espectro da matriz A pré-condicionada usando o método especificado.

        Args:
            method (str): Método de pré-condicionamento a ser utilizado. Pode ser 'Jacobi', 'Seidel', 'ILUfac', 'Multiscale', 'MultiscaleILUfac', 'MultiscaleJacobi' ou 'MultiscaleSeidel'.

        Returns:
            numpy.ndarray: Array contendo os 100 maiores e 100 menores autovalores da matriz A pré-condicionada.
            -- Menores autovalores podem não convergir, então podem ter menos de 100 menores.
        """

        # ----- Análise espectral da matriz de transmissibilidade pré-condicionada -----
        lowerA = sp.tril(self._A, k=-1)    # matriz triangular inferior
        diagA  = sp.diags(self._A.diagonal(), format="csc")             # diagonal da matriz

        args = (lowerA, diagA)
        if method in self._SpectralMethods:
            return self._SpectralMethods[method](args)
        else:
            raise ValueError(f"Método '{method}' não reconhecido.")
        
    def Solve(self, method, x0=None, rtol=1e-8, maxiter=2500):
        """
        Resolve o sistema linear usando o método especificado.

        Args:
            method (str): Método de pré-condicionamento a ser utilizado. Pode ser 'Jacobi', 'Seidel', 'ILUfac', 'Multiscale', 'MultiscaleILUfac', 'MultiscaleJacobi' ou 'MultiscaleSeidel'.
            x0 (numpy.ndarray | None): Vetor inicial para o método iterativo. Se None, será utilizado o vetor nulo.
            rtol (float): Tolerância relativa para a convergência do método iterativo
            maxiter (int): Número máximo de iterações para o método iterativo.

        Returns:
            Tuple[numpy.ndarray, numpy.ndarray, float, int]: Uma tupla contendo:
                - Array de resíduos das iterações.
                - Solução do sistema linear.
                - Tempo de execução do método.
                - Informação sobre a convergência do método (0 se convergiu, >0 se não convergiu).
        """

        lowerA = sp.tril(self._A, k=-1)    # matriz triangular inferior
        diagA  = sp.diags(self._A.diagonal(), format="csc")             # diagonal da matriz
        args = (lowerA, diagA)
        if method in self._SolveMethods:
            init = time.time()

            linOP = self._SolveMethods[method](args)
            residuals = []
            def callback(xk):
                r = self._b - self._A @ xk
                residuals.append(np.linalg.norm(r))

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
            return residuals, x, end-init, info
        
        
        else:
            raise ValueError(f"Método '{method}' não reconhecido.")
        
        
    def PlotSpectre(self, method, autovalores, save_path=None, problem = ""):
        """
        Plota o espectro do pré-condicionador especificado de um dos métodos.

        Args:
            method (str): Método de pré-condicionamento utilizado (Importante apenas para título e onde será armazenado).
            autovalores (numpy.ndarray): Array contendo os autovalores.
            save_path (str | None): Caminho para salvar o gráfico do espectro. Se None, o gráfico será exibido na tela.
            problem (str): Nome do problema para o qual o espectro será gerado.
        """

        plt.figure()
        x = np.arange(-1.0, 1.01, 0.01)
        x[-1] = 1.0
        z = np.sqrt(1 - x**2)
        y = -np.sqrt(1 - x**2)
        plt.plot(x, z, 'b')
        plt.plot(x, y, 'b')
        plt.plot(np.real(autovalores), np.imag(autovalores), '*r', label='Autovalores')
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
        """
        Plota os resíduos das iterações para os métodos especificados.

        Args:
            methods (List[str]): Lista de métodos utilizados.
            residues (List[numpy.ndarray]): Lista (tipo dicionário) de arrays contendo os resíduos das iterações para cada método.
            save_path (str | None): Caminho para salvar o gráfico dos resíduos. Se None, o gráfico será exibido na tela.
            problem (str): Nome do problema para o qual os resíduos serão plotados.
        """

        plt.figure()
        plt.xlabel("Iterações")
        plt.ylabel("Resíduo")
        plt.title(f'Resíduos: {problem}')
        plt.yscale("log")
        plt.grid(True)
        markers = ['s', 'D', 'h', 'v', 'X', 'p', '*', 'o']
        every = [1, 2, 3, 4, 1, 2, 3, 4]
        for method in methods:
            if len(residues[method]) < 10:
                every[methods.index(method)] = 1
            plt.plot(residues[method], marker=markers[methods.index(method)], markevery=every[methods.index(method)], label=method)

        plt.legend(methods, loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0.)
        plt.tight_layout()

        if save_path:
            plt.savefig(f'{save_path}/{problem}/Residues/Residues_{problem}')
        else:
            plt.show()
        plt.close()