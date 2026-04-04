import numpy as np
import scipy.sparse.linalg as spla
import ilupp

class SolvingMethods:
    def __init__(self, A, b, OP, OR):
        self._A = A
        self._b = b
        self._OP = OP
        self._OR = OR

    def Jacobi(self, args):
        # ----- MÉTODO DE JACOBI -----
        nnod, upperA, lowerA, diagA = args

        inv_diag = 1/diagA.diagonal()
        def precond(x):
            return inv_diag * x
        
        M = spla.LinearOperator(self._A.shape, precond)
        return M

    def Seidel(self, args):
        # ----- MÉTODO DE GAUSS-SEIDEL -----
        nnod, upperA, lowerA, diagA = args

        A = lowerA + diagA
        def precond(x):
            return spla.spsolve_triangular(A, x, lower=True)
        
        M = spla.LinearOperator(self._A.shape, precond)
        return M
        
    def ILUfac(self, args):
        # ----- MATRIZ A COM SUAVIZADOR ILU(0) -----
        nnod, upperA, lowerA, diagA = args

        ilu = ilupp.ILU0Preconditioner(self._A)
        return ilu    
    
    def Multiscale(self, args):
        # ----- MATRIZ A com pré-condicionador Multiescala -----
        nnod, upperA, lowerA, diagA = args

        RAP = self._OR @ self._A @ self._OP
        solver = spla.factorized(RAP)
        def precond(x): # Aproximação de A-1 multiscale
            return self._OP @ solver(self._OR @ x)

        M = spla.LinearOperator(self._A.shape, precond)
        return M

    def Multiscale_ILUfac(self, args):
        # # ----- MATRIZ A pré-condicionador Multiescala + SUAVIZADOR ILU(0) -----
        nnod, upperA, lowerA, diagA = args

        ilu = ilupp.ILU0Preconditioner(self._A)
        def precond_ilu(x): # Aproximação de A-1 ilu
            val = x.astype(np.float64).copy() 
            ilu.apply(val)
            return val
        
        RAP = self._OR @ self._A @ self._OP
        solver = spla.factorized(RAP)
        def precond_multiscale(x): # Aproximação de A-1 multiscale
            return self._OP @ solver(self._OR @ x)

        def precond(x):
            return (precond_multiscale(x) + precond_ilu(x) - precond_multiscale(self._A @ precond_ilu(x))) # precond_milu

        M = spla.LinearOperator(self._A.shape, precond)
        return M

    def Multiscale_Jacobi(self, args):
        # ----- MATRIZ A pré-condicionador Multiescala + Jacobi -----
        nnod, upperA, lowerA, diagA = args

        inv_diag = 1/diagA.diagonal()
        def precond_jacobi(x):
            return inv_diag * x
        
        RAP = self._OR @ self._A @ self._OP
        solver = spla.factorized(RAP)
        def precond_multiscale(x): # Aproximação de A-1 multiscale
            return self._OP @ solver(self._OR @ x)

        def precond(x):
            return (precond_multiscale(x) + precond_jacobi(x) - precond_multiscale(self._A @ precond_jacobi(x))) # precond_multiscale_jacobi

        M = spla.LinearOperator(self._A.shape, precond)
        
        return M

    def Multiscale_Seidel(self, args):
        nnod, upperA, lowerA, diagA = args

        A = lowerA + diagA
        def precond_seidel(x):
            return spla.spsolve_triangular(A, x, lower=True)

        RAP = self._OR @ self._A @ self._OP
        solver = spla.factorized(RAP)
        def precond_multiscale(x): # Aproximação de A-1 multiscale
            return self._OP @ solver(self._OR @ x)

        def precond(x):
            return (precond_multiscale(x) + precond_seidel(x) - precond_multiscale(self._A @ precond_seidel(x))) # precond_multiscale_jacobi
            
        M = spla.LinearOperator(self._A.shape, precond)
        return M