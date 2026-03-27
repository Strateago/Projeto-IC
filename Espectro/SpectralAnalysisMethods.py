import numpy as np
import scipy.io as sio
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
import ilupp
import time

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
        except:
            print("Aviso: SM não convergiu.")
            SM = np.array([])
        print('ok')

        av_error_operator = np.concatenate([BM, SM])

        # # Resolvendo o sistema linear
        # t0 = time.time()

        # iter = 1
        # itermax = 100
        # tolerancia = 1.0e-3

        # xold = np.zeros(nnod)
        # xnew = xold.copy()
        # resold = self._b - self._A @ xold
        # delta = spla.norm(resold)
        # deltaresold = [delta]
        # S = 1/diagA

        # while delta > tolerancia and iter < itermax:
        #     xnew = xold + S @ resold
        #     resold = self._b - self._A @ xnew
        #     delta = spla.norm(resold)

        #     xold = xnew.copy()
        #     iter += 1
        #     deltaresold.append(delta)

        # t_final = time.time() - t0
        # print(f"Tempo total (s): {t_final:.4f}")
        # print(f"Número de iterações: {iter}")
        # print(f"Resíduo final: {delta:.3e}")

        # return av_error_operator, deltaresold
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
        except:
            print("Aviso: SM não convergiu.")
            SM = np.array([])
        print('ok')
        
        av_error_operator = np.concatenate([BM, SM])

        # # Resolvendo o sistema linear
        # t0 = time.time()

        # iter = 1
        # itermax = 100
        # tolerancia = 1.0e-3

        # xold = np.zeros(nnod)
        # xnew = xold.copy()
        # resold = self._b - self._A @ xold
        # delta = np.linalg.norm(resold)
        # deltaresold = [delta]
        # S = la.inv(diagA + lowerA)

        # while delta > tolerancia and iter < itermax:
        #     xnew = xold + S @ resold
        #     resold = self._b - self._A @ xnew
        #     delta = np.linalg.norm(resold)

        #     xold = xnew.copy()
        #     iter += 1
        #     deltaresold.append(delta)

        # t_final = time.time() - t0
        # print(f"Tempo total (s): {t_final:.4f}")
        # print(f"Número de iterações: {iter}")
        # print(f"Resíduo final: {delta:.3e}")

        # return av_error_operator, deltaresold
        return av_error_operator
        
    def ILUfac(self, args):
        # ----- MATRIZ A COM SUAVIZADOR ILU(0) -----
        nnod, upperA, lowerA, diagA = args

        # ilu = spla.spilu(self._A, fill_factor=1)
        # def apply_G(x):
        #     Ax = self._A @ x
        #     precond_step = ilu.solve(Ax.astype(np.float64))
        #     return x - precond_step

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
            SM = spla.eigs(error_operator, k=100 , which='SM', return_eigenvectors=False)
        except:
            print("Aviso: SM não convergiu.")
            SM = np.array([])
        print('ok')
        
        av_error_operator = np.concatenate([BM, SM])
        # av_error_operator = BM

        # # Resolvendo o sistema linear
        # t0 = time.time()

        # xold = np.zeros(nnod)
        # xnew = xold.copy()
        # iter = 1
        # itermax = 100
        # resold = self._b - self._A @ xold
        # tolerancia = 1.0e-3
        # delta = np.linalg.norm(resold)
        # deltaresold = [delta]
        # S = precond_ilu

        # while delta > tolerancia and iter < itermax:
        #     xnew = xold + S @ resold
        #     resold = self._b - self._A @ xnew
        #     delta = np.linalg.norm(resold)
        #     xold = xnew.copy()
        #     iter += 1
        #     deltaresold.append(delta)

        # t_final = time.time() - t0
        # print(f"Tempo total (s): {t_final:.4f}")

        # return av_error_operator, deltaresold
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
            SM = spla.eigs(error_operator, k=100 , which='SM', return_eigenvectors=False)
        except:
            print("Aviso: SM não convergiu.")
            SM = np.array([])
        print('ok')
        
        av_error_operator = np.concatenate([BM, SM])


        # # Resolvendo o sistema linear
        # t0 = time.time()

        # iter = 1
        # itermax = 100
        # tolerancia = 1.0e-3

        # xold = np.zeros(nnod)
        # xnew = xold.copy()
        # resold = self._b - self._A @ xold
        # delta = np.linalg.norm(resold)
        # deltaresold = [delta]
        # S = precond_multiscale

        # while delta > tolerancia and iter < itermax:
        #     xnew = xold + S @ resold
        #     resold = self._b - self._A @ xnew
        #     delta = np.linalg.norm(resold)

        #     xold = xnew.copy()
        #     iter += 1
        #     deltaresold.append(delta)

        # t_final = time.time() - t0
        # print(f"Tempo total (s): {t_final:.4f}")
        # print(f"Número de iterações: {iter}")
        # print(f"Resíduo final: {delta:.3e}")

        # return av_error_operator, deltaresold
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
            SM = spla.eigs(error_operator, k=100 , which='SM', return_eigenvectors=False)
        except:
            print("Aviso: SM não convergiu.")
            SM = np.array([])
        print('ok')
        
        av_error_operator = np.concatenate([BM, SM])

        # # Resolvendo o sistema linear
        # t0 = time.time()

        # iter = 1
        # itermax = 1000
        # tolerancia = 1.0e-3

        # xold = np.zeros(nnod)
        # xnew = xold.copy()
        # resold1 = self._b - self._A @ xold
        # delta = np.linalg.norm(resold1)
        # deltaresold = [delta]
        # S1 = precond_multiscale
        # S2 = precond_ilu

        # while delta > tolerancia and iter < itermax:
        #     xmed = xold + S1 @ resold1
        #     resold2 = self._b - self._A @ xnew
        #     xnew = xmed + S2 @ resold2
        #     resold1 = self._b - self._A @ xnew
        #     delta = np.linalg.norm(resold1)

        #     xold = xnew.copy()
        #     iter += 1
        #     deltaresold.append(delta)

        # t_final = time.time() - t0
        # print(f"Tempo total (s): {t_final:.4f}")
        # print(f"Número de iterações: {iter}")
        # print(f"Resíduo final: {delta:.3e}")

        # return av_error_operator_milu, deltaresold
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
            SM = spla.eigs(error_operator, k=100 , which='SM', return_eigenvectors=False)
        except:
            print("Aviso: SM não convergiu.")
            SM = np.array([])
        print('ok')
        
        av_error_operator = np.concatenate([BM, SM])

        # # Resolvendo o sistema linear
        # t0 = time.time()

        # iter = 1
        # itermax = 1000
        # tolerancia = 1.0e-3

        # xold = np.zeros(nnod)
        # xnew = xold.copy()
        # resold1 = self._b - self._A @ xold
        # delta = np.linalg.norm(resold1)
        # deltaresold = [delta]
        # S1 = precond_multiscale
        # S2 = precond_jacobi

        # while delta > tolerancia and iter < itermax:
        #     xmed = xold + S1 @ resold1
        #     resold2 = self._b - self._A @ xnew
        #     xnew = xmed + S2 @ resold2
        #     resold1 = self._b - self._A @ xnew
        #     delta = np.linalg.norm(resold1)

        #     xold = xnew.copy()
        #     iter += 1
        #     deltaresold.append(delta)

        # t_final = time.time() - t0
        # print(f"Tempo total (s): {t_final:.4f}")
        # print(f"Número de iterações: {iter}")
        # print(f"Resíduo final: {delta:.3e}")

        # return av_error_operator_mjacobi, deltaresold
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
            SM = spla.eigs(error_operator, k=100 , which='SM', return_eigenvectors=False)
        except:
            print("Aviso: SM não convergiu.")
            SM = np.array([])
        print('ok')
        
        av_error_operator = np.concatenate([BM, SM])

        # # Resolvendo o sistema linear
        # t0 = time.time()

        # iter = 1
        # itermax = 1000
        # tolerancia = 1.0e-3

        # xold = np.zeros(nnod)
        # xnew = xold.copy()
        # resold1 = self._b - self._A @ xold
        # delta = np.linalg.norm(resold1)
        # deltaresold = [delta]
        # S1 = precond_multiscale
        # S2 = precond_seidel

        # while delta > tolerancia and iter < itermax:
        #     xmed = xold + S1 @ resold1
        #     resold2 = self._b - self._A @ xnew
        #     xnew = xmed + S2 @ resold2
        #     resold1 = self._b - self._A @ xnew
        #     delta = np.linalg.norm(resold1)

        #     xold = xnew.copy()
        #     iter += 1
        #     deltaresold.append(delta)

        # t_final = time.time() - t0
        # print(f"Tempo total (s): {t_final:.4f}")
        # print(f"Número de iterações: {iter}")
        # print(f"Resíduo final: {delta:.3e}")

        # return av_error_operator, deltaresold
        return av_error_operator