import numpy as np
import scipy.sparse.linalg as spla
import ilupp
import pyamg

class SolvingMethods:
    """
    Classe que implementa diferentes métodos de pré-condicionamento para resolver sistemas lineares.
    A classe cria linear operators que representam os pré-condicionadores, que podem ser usados em métodos iterativos para resolver sistemas lineares.
    Novos métodos podem ser adicionados facilmente, basta criar uma função que retorne um objeto LinearOperator e adicioná-la ao dicionário _SolvingMethods na classe principal.

    Métodos:
        Jacobi(args): Retorna um LinearOperator que representa o pré-condicionador de Jacobi.
        Seidel(args): Retorna um LinearOperator que representa o pré-condicionador de Gauss-Seidel.
        ILUfac(args): Retorna um LinearOperator que representa o pré-condicionador ILU(0).
        Multiscale(args): Retorna um LinearOperator que representa o pré-condicionador multiescala.
        Multiscale_ILUfac(args): Retorna um LinearOperator que representa o pré-condicionador multiescala combinado com ILU(0).
        Multiscale_Jacobi(args): Retorna um LinearOperator que representa o pré-condicionador multiescala combinado com Jacobi.
        Multiscale_Seidel(args): Retorna um LinearOperator que representa o pré-condicionador multiescala combinado com Gauss-Seidel.

    """
    def __init__(self, A, OP, OR):
        self._A = A
        self._OP = OP
        self._OR = OR

    def Jacobi(self, args):
        """
        Método de pré-condicionamento de Jacobi: Diagonal.

        Args:
            args (tuple): Uma tupla contendo:
                - lowerA (scipy.sparse.csc_matrix): A matriz triangular inferior de A.
                - diagA (scipy.sparse.csc_matrix): A matriz diagonal de A.
        """

        # ----- MÉTODO DE JACOBI -----
        lowerA, diagA = args

        inv_diag = 1/diagA.diagonal()
        def precond(x):
            return inv_diag * x
        
        M = spla.LinearOperator(self._A.shape, precond)
        return M

    def Seidel(self, args):
        """
        Método de pré-condicionamento de Gauss-Seidel: Triangular inferior.

        Args:
            args (tuple): Uma tupla contendo:
                - lowerA (scipy.sparse.csc_matrix): A matriz triangular inferior de A.
                - diagA (scipy.sparse.csc_matrix): A matriz diagonal de A.
        """

        # ----- MÉTODO DE GAUSS-SEIDEL -----
        lowerA, diagA = args

        A = lowerA + diagA
        def precond(x):
            return spla.spsolve_triangular(A, x, lower=True)
        
        M = spla.LinearOperator(self._A.shape, precond)
        return M
        
    def ILUfac(self, args):
        """
        Método de pré-condicionamento ILU(0): Fatoração LU incompleta.

        Args:
            args (tuple): Uma tupla contendo:
                - lowerA (scipy.sparse.csc_matrix): A matriz triangular inferior de A.
                - diagA (scipy.sparse.csc_matrix): A matriz diagonal de A.
        """

        # ----- MATRIZ A COM SUAVIZADOR ILU(0) -----
        lowerA, diagA = args

        ilu = ilupp.ILU0Preconditioner(self._A)
        return ilu    
    
    def Multiscale(self, args):
        """
        Método de pré-condicionamento multiescala: Aproximação de A-1 usando operadores de restrição e prolongamento.

        Args:
            args (tuple): Uma tupla contendo:
                - lowerA (scipy.sparse.csc_matrix): A matriz triangular inferior de A.
                - diagA (scipy.sparse.csc_matrix): A matriz diagonal de A.
        """

        # ----- MATRIZ A com pré-condicionador Multiescala -----
        lowerA, diagA = args

        RAP = self._OR @ self._A @ self._OP
        solver = spla.factorized(RAP)
        def precond(x): # Aproximação de A-1 multiscale
            return self._OP @ solver(self._OR @ x)

        M = spla.LinearOperator(self._A.shape, precond)
        return M

    def Multiscale_ILUfac(self, args):
        """
        Método de pré-condicionamento multiescala combinado com ILU(0): Aproximação de A-1 usando operadores de restrição e prolongamento, e suavização com ILU(0).

        Args:
            args (tuple): Uma tupla contendo:
                - lowerA (scipy.sparse.csc_matrix): A matriz triangular inferior de A.
                - diagA (scipy.sparse.csc_matrix): A matriz diagonal de A.
        """

        # # ----- MATRIZ A pré-condicionador Multiescala + SUAVIZADOR ILU(0) -----
        lowerA, diagA = args

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
        """
        Método de pré-condicionamento multiescala combinado com Jacobi: Aproximação de A-1 usando operadores de restrição e prolongamento, e suavização com Jacobi.

        Args:
            args (tuple): Uma tupla contendo:
                - lowerA (scipy.sparse.csc_matrix): A matriz triangular inferior de A.
                - diagA (scipy.sparse.csc_matrix): A matriz diagonal de A.
        """

        # ----- MATRIZ A pré-condicionador Multiescala + Jacobi -----
        lowerA, diagA = args

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
        """
        Método de pré-condicionamento multiescala combinado com Gauss-Seidel: Aproximação de A-1 usando operadores de restrição e prolongamento, e suavização com Gauss-Seidel.

        Args:
            args (tuple): Uma tupla contendo:
                - lowerA (scipy.sparse.csc_matrix): A matriz triangular inferior de A.
                - diagA (scipy.sparse.csc_matrix): A matriz diagonal de A.
        """

        lowerA, diagA = args

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

    def AMG(self, args):
        """
        Método de pré-condicionamento AMG (Algebraic Multigrid): Aproximação de A-1 usando o método multigrid algébrico.
        Método do Estado da Arte para pré-condicionamento de sistemas lineares esparsos.

        Args:
            args (tuple): Uma tupla contendo:
                - lowerA (scipy.sparse.csc_matrix): A matriz triangular inferior de A.
                - diagA (scipy.sparse.csc_matrix): A matriz diagonal de A.
        """

        lowerA, diagA = args
        A = self._A.tocsr()

        ml = pyamg.smoothed_aggregation_solver(A)
        M = ml.aspreconditioner()
        return M