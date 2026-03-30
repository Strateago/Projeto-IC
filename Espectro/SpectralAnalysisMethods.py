import numpy as np
import scipy.sparse.linalg as spla
import ilupp

class SpectralAnalysisMethods:
    def __init__(self, A, b, OP, OR):
        self._A = A
        self._b = b
        self._OP = OP
        self._OR = OR

    def Jacobi(self, args):
        # ----- MÉTODO DE JACOBI -----
        nnod, upperA, lowerA, diagA = args

        precond_jacobi = 1/(diagA.diagonal())                                  # pré-condicionador diagonal

        def apply_G(x):
            return x - (precond_jacobi * (self._A @ x))

        error_operator = spla.LinearOperator((self._A.shape[0], self._A.shape[0]), matvec=apply_G)
        print('Init 1')
        # 100 maior magnitude
        BM = spla.eigs(error_operator, k=100, which='LM', return_eigenvectors=False)          # espectro da matriz de iteração
        print('ok\nInit 2')
        # 100 menor magnitude
        try:
            SM = spla.eigs(error_operator, k=100, which='SM', return_eigenvectors=False)
        except spla.ArpackNoConvergence as error: # Captura especificamente o erro de convergência
            SM = error.eigenvalues
            print(f'SM não convergiu totalmente. Foram obtidos {len(SM)} autovalores.')
        except Exception as error: # Captura outros erros (memória, álgebra, etc)
            SM = np.array([])
            print(f'SM falhou por outro motivo: {error}')
        print('ok')

        av_error_operator = np.concatenate([BM, SM])
        return av_error_operator

    def Seidel(self, args):
        # ----- MÉTODO DE GAUSS-SEIDEL -----
        nnod, upperA, lowerA, diagA = args

        M = lowerA + diagA
        def apply_G(x):
            precond_seidel = spla.spsolve_triangular(M, self._A @ x, lower=True)
            return x - precond_seidel
        
        error_operator = spla.LinearOperator((self._A.shape[0], self._A.shape[0]), matvec=apply_G)
        # error_operator = np.eye(nnod) - (precond_seidel @ self._A)
        print('Init 1')
        # 100 maior magnitude
        BM = spla.eigs(error_operator, k=100, which='LM', return_eigenvectors=False)          # espectro da matriz de iteração
        print('ok\nInit 2')
        # 100 menor magnitude
        try:
            SM = spla.eigs(error_operator, k=100, which='SM', return_eigenvectors=False)
        except spla.ArpackNoConvergence as error: # Captura especificamente o erro de convergência
            SM = error.eigenvalues
            print(f'SM não convergiu totalmente. Foram obtidos {len(SM)} autovalores.')
        except Exception as error: # Captura outros erros (memória, álgebra, etc)
            SM = np.array([])
            print(f'SM falhou por outro motivo: {error}')
        print('ok')
        
        av_error_operator = np.concatenate([BM, SM])
        return av_error_operator
        
    def ILUfac(self, args):
        # ----- MATRIZ A COM SUAVIZADOR ILU(0) -----
        nnod, upperA, lowerA, diagA = args

        ilu = ilupp.ILU0Preconditioner(self._A)
        def apply_G(x):
            Ax = self._A @ x
            # aplicar M^{-1}Ax com ilupp
            precond_step = Ax.copy()
            ilu.apply(precond_step)

            return x - precond_step

        # error_operator = np.eye(nnod) - precond_system
        error_operator = spla.LinearOperator((self._A.shape[0], self._A.shape[0]), matvec=apply_G)
        # av_error_operator = la.eigvals(error_operator)          # espectro da matriz de iteração
        print('Init 1')
        # 100 maior magnitude
        BM = spla.eigs(error_operator, k=100, which='LM', return_eigenvectors=False)          # espectro da matriz de iteração
        print('ok\nInit 2')
        # 100 menor magnitude
        try:
            SM = spla.eigs(error_operator, k=100, which='SM', return_eigenvectors=False)
        except spla.ArpackNoConvergence as error: # Captura especificamente o erro de convergência
            SM = error.eigenvalues
            print(f'SM não convergiu totalmente. Foram obtidos {len(SM)} autovalores.')
        except Exception as error: # Captura outros erros (memória, álgebra, etc)
            SM = np.array([])
            print(f'SM falhou por outro motivo: {error}')
        print('ok')
        
        av_error_operator = np.concatenate([BM, SM])
        return av_error_operator
    
    def Multiscale(self, args):
        # ----- MATRIZ A com pré-condicionador Multiescala -----
        nnod, upperA, lowerA, diagA = args

        RAP = self._OR @ self._A @ self._OP
        RA = self._OR @ self._A
        solver = spla.factorized(RAP)
        def apply_G(x):
            return x - self._OP @ solver(RA @ x)
        # av_error_operator = la.eigvals(np.eye(nnod) - precond_system)                 # espectro da matriz de iteração
        error_operator = spla.LinearOperator((self._A.shape[0], self._A.shape[0]), matvec=apply_G)
        # espec_error_operator = np.max(np.abs(av_error_operator))                      # raio espectral do operador de propagação de erro
        print('Init 1')
        # 100 maior magnitude
        BM = spla.eigs(error_operator, k=100, which='LM', return_eigenvectors=False)          # espectro da matriz de iteração
        print('ok\nInit 2')
        # 100 menor magnitude
        try:
            SM = spla.eigs(error_operator, k=100, which='SM', return_eigenvectors=False)
        except spla.ArpackNoConvergence as error: # Captura especificamente o erro de convergência
            SM = error.eigenvalues
            print(f'SM não convergiu totalmente. Foram obtidos {len(SM)} autovalores.')
        except Exception as error: # Captura outros erros (memória, álgebra, etc)
            SM = np.array([])
            print(f'SM falhou por outro motivo: {error}')
        print('ok')
        
        av_error_operator = np.concatenate([BM, SM])
        return av_error_operator

    def Multiscale_ILUfac(self, args):
        # # ----- MATRIZ A pré-condicionador Multiescala + SUAVIZADOR ILU(0) -----
        nnod, upperA, lowerA, diagA = args

        ilu = ilupp.ILU0Preconditioner(self._A)
        def precond_ilu(x): # Aproximação de A-1 ilu
            val = x.copy()
            ilu.apply(val)
            return val

        RAP = self._OR @ self._A @ self._OP
        solver = spla.factorized(RAP)
        def precond_multiscale(x): # Aproximação de A-1 multiscale
            return self._OP @ solver(self._OR @ x)

        def apply_G(x):
            Ax = self._A @ x
            milu_Ax = (precond_multiscale(Ax) + precond_ilu(Ax) - precond_multiscale(self._A @ precond_ilu(Ax))) # precond_milu @ Ax
            return x - milu_Ax
        
        error_operator = spla.LinearOperator((self._A.shape[0], self._A.shape[0]), matvec=apply_G)
        
        # espec_error_operator = np.max(np.abs(av_error_operator))                      # raio espectral do operador de propagação de erro
        print('Init 1')
        # 100 maior magnitude
        BM = spla.eigs(error_operator, k=100, which='LM', return_eigenvectors=False)          # espectro da matriz de iteração
        print('ok\nInit 2')
        # 100 menor magnitude
        try:
            SM = spla.eigs(error_operator, k=100, which='SM', return_eigenvectors=False)
        except spla.ArpackNoConvergence as error: # Captura especificamente o erro de convergência
            SM = error.eigenvalues
            print(f'SM não convergiu totalmente. Foram obtidos {len(SM)} autovalores.')
        except Exception as error: # Captura outros erros (memória, álgebra, etc)
            SM = np.array([])
            print(f'SM falhou por outro motivo: {error}')
        print('ok')
        
        av_error_operator = np.concatenate([BM, SM])
        return av_error_operator

    def Multiscale_Jacobi(self, args):
        # ----- MATRIZ A pré-condicionador Multiescala + Jacobi -----
        nnod, upperA, lowerA, diagA = args

        def precond_jacobi(x):
            return spla.spsolve(diagA, x)
        
        RAP = self._OR @ self._A @ self._OP
        solver = spla.factorized(RAP)
        def precond_multiscale(x): # Aproximação de A-1 multiscale
            return self._OP @ solver(self._OR @ x)

        def apply_G(x):
            Ax = self._A @ x
            multiscale_jacobi_Ax = (precond_multiscale(Ax) + precond_jacobi(Ax) - precond_multiscale(self._A @ precond_jacobi(Ax))) # precond_multiscale_jacobi @ Ax
            return x - multiscale_jacobi_Ax

        error_operator = spla.LinearOperator((self._A.shape[0], self._A.shape[0]), matvec=apply_G)
        # espec_error_operator = np.max(np.abs(av_error_operator))                      # raio espectral do operador de propagação de erro
        print('Init 1')
        # 100 maior magnitude
        BM = spla.eigs(error_operator, k=100, which='LM', return_eigenvectors=False)          # espectro da matriz de iteração
        print('ok\nInit 2')
        # 100 menor magnitude
        try:
            SM = spla.eigs(error_operator, k=100, which='SM', return_eigenvectors=False)
        except spla.ArpackNoConvergence as error: # Captura especificamente o erro de convergência
            SM = error.eigenvalues
            print(f'SM não convergiu totalmente. Foram obtidos {len(SM)} autovalores.')
        except Exception as error: # Captura outros erros (memória, álgebra, etc)
            SM = np.array([])
            print(f'SM falhou por outro motivo: {error}')
        print('ok')
        
        av_error_operator = np.concatenate([BM, SM])
        return av_error_operator

    def Multiscale_Seidel(self, args):
        nnod, upperA, lowerA, diagA = args

        M = lowerA + diagA
        def precond_seidel(x):
            return spla.spsolve_triangular(M, x, lower=True)

        RAP = self._OR @ self._A @ self._OP
        solver = spla.factorized(RAP)
        def precond_multiscale(x): # Aproximação de A-1 multiscale
            return self._OP @ solver(self._OR @ x)

        def apply_G(x):
            Ax = self._A @ x
            multiscale_jacobi_Ax = (precond_multiscale(Ax) + precond_seidel(Ax) - precond_multiscale(self._A @ precond_seidel(Ax))) # precond_multiscale_jacobi @ Ax
            return x - multiscale_jacobi_Ax
        
        error_operator = spla.LinearOperator((self._A.shape[0], self._A.shape[0]), matvec=apply_G)
        # espec_error_operator = np.max(np.abs(av_error_operator))                      # raio espectral do operador de propagação de erro
        print('Init 1')
        # 100 maior magnitude
        BM = spla.eigs(error_operator, k=100, which='LM', return_eigenvectors=False)          # espectro da matriz de iteração
        print('ok\nInit 2')
        # 100 menor magnitude
        try:
            SM = spla.eigs(error_operator, k=100, which='SM', return_eigenvectors=False)
        except spla.ArpackNoConvergence as error: # Captura especificamente o erro de convergência
            SM = error.eigenvalues
            print(f'SM não convergiu totalmente. Foram obtidos {len(SM)} autovalores.')
        except Exception as error: # Captura outros erros (memória, álgebra, etc)
            SM = np.array([])
            print(f'SM falhou por outro motivo: {error}')
        print('ok')
        
        av_error_operator = np.concatenate([BM, SM])
        return av_error_operator