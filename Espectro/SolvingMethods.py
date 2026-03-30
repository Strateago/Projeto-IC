import numpy as np
import scipy.sparse.linalg as spla
import ilupp
import time

class SolvingMethods:
    def __init__(self, A, b, OP, OR):
        self._A = A
        self._b = b
        self._OP = OP
        self._OR = OR

    def Jacobi(self, args):
        # ----- MÉTODO DE JACOBI -----
        nnod, upperA, lowerA, diagA = args

        # Resolvendo o sistema linear
        t0 = time.time()

        iter = 1
        itermax = 100
        tolerancia = 1.0e-3

        xold = np.zeros(nnod)
        xnew = xold.copy()
        resold = self._b - self._A @ xold
        delta = spla.norm(resold)
        deltaresold = [delta]
        S = 1/diagA

        while delta > tolerancia and iter < itermax:
            xnew = xold + S @ resold
            resold = self._b - self._A @ xnew
            delta = spla.norm(resold)

            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")
        print(f"Número de iterações: {iter}")
        print(f"Resíduo final: {delta:.3e}")

        return deltaresold

    def Seidel(self, args):
        # ----- MÉTODO DE GAUSS-SEIDEL -----
        nnod, upperA, lowerA, diagA = args

        # Resolvendo o sistema linear
        t0 = time.time()

        iter = 1
        itermax = 100
        tolerancia = 1.0e-3

        xold = np.zeros(nnod)
        xnew = xold.copy()
        resold = self._b - self._A @ xold
        delta = np.linalg.norm(resold)
        deltaresold = [delta]
        S = la.inv(diagA + lowerA)

        while delta > tolerancia and iter < itermax:
            xnew = xold + S @ resold
            resold = self._b - self._A @ xnew
            delta = np.linalg.norm(resold)

            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")
        print(f"Número de iterações: {iter}")
        print(f"Resíduo final: {delta:.3e}")

        return deltaresold
        
    def ILUfac(self, args):
        # ----- MATRIZ A COM SUAVIZADOR ILU(0) -----
        nnod, upperA, lowerA, diagA = args

        # Resolvendo o sistema linear
        t0 = time.time()

        xold = np.zeros(nnod)
        xnew = xold.copy()
        iter = 1
        itermax = 100
        resold = self._b - self._A @ xold
        tolerancia = 1.0e-3
        delta = np.linalg.norm(resold)
        deltaresold = [delta]
        S = precond_ilu

        while delta > tolerancia and iter < itermax:
            xnew = xold + S @ resold
            resold = self._b - self._A @ xnew
            delta = np.linalg.norm(resold)
            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")

        return deltaresold
    
    def Multiscale(self, args):
        # ----- MATRIZ A com pré-condicionador Multiescala -----
        nnod, upperA, lowerA, diagA = args

        # Resolvendo o sistema linear
        t0 = time.time()

        iter = 1
        itermax = 100
        tolerancia = 1.0e-3

        xold = np.zeros(nnod)
        xnew = xold.copy()
        resold = self._b - self._A @ xold
        delta = np.linalg.norm(resold)
        deltaresold = [delta]
        S = precond_multiscale

        while delta > tolerancia and iter < itermax:
            xnew = xold + S @ resold
            resold = self._b - self._A @ xnew
            delta = np.linalg.norm(resold)

            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")
        print(f"Número de iterações: {iter}")
        print(f"Resíduo final: {delta:.3e}")

        return deltaresold

    def Multiscale_ILUfac(self, args):
        # # ----- MATRIZ A pré-condicionador Multiescala + SUAVIZADOR ILU(0) -----
        nnod, upperA, lowerA, diagA = args

        # Resolvendo o sistema linear
        t0 = time.time()

        iter = 1
        itermax = 1000
        tolerancia = 1.0e-3

        xold = np.zeros(nnod)
        xnew = xold.copy()
        resold1 = self._b - self._A @ xold
        delta = np.linalg.norm(resold1)
        deltaresold = [delta]
        S1 = precond_multiscale
        S2 = precond_ilu

        while delta > tolerancia and iter < itermax:
            xmed = xold + S1 @ resold1
            resold2 = self._b - self._A @ xnew
            xnew = xmed + S2 @ resold2
            resold1 = self._b - self._A @ xnew
            delta = np.linalg.norm(resold1)

            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")
        print(f"Número de iterações: {iter}")
        print(f"Resíduo final: {delta:.3e}")

        return deltaresold

    def Multiscale_Jacobi(self, args):
        # ----- MATRIZ A pré-condicionador Multiescala + Jacobi -----
        nnod, upperA, lowerA, diagA = args

        # Resolvendo o sistema linear
        t0 = time.time()

        iter = 1
        itermax = 1000
        tolerancia = 1.0e-3

        xold = np.zeros(nnod)
        xnew = xold.copy()
        resold1 = self._b - self._A @ xold
        delta = np.linalg.norm(resold1)
        deltaresold = [delta]
        S1 = precond_multiscale
        S2 = precond_jacobi

        while delta > tolerancia and iter < itermax:
            xmed = xold + S1 @ resold1
            resold2 = self._b - self._A @ xnew
            xnew = xmed + S2 @ resold2
            resold1 = self._b - self._A @ xnew
            delta = np.linalg.norm(resold1)

            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")
        print(f"Número de iterações: {iter}")
        print(f"Resíduo final: {delta:.3e}")

        return deltaresold

    def Multiscale_Seidel(self, args):
        nnod, upperA, lowerA, diagA = args

        # Resolvendo o sistema linear
        t0 = time.time()

        iter = 1
        itermax = 1000
        tolerancia = 1.0e-3

        xold = np.zeros(nnod)
        xnew = xold.copy()
        resold1 = self._b - self._A @ xold
        delta = np.linalg.norm(resold1)
        deltaresold = [delta]
        S1 = precond_multiscale
        S2 = precond_seidel

        while delta > tolerancia and iter < itermax:
            xmed = xold + S1 @ resold1
            resold2 = self._b - self._A @ xnew
            xnew = xmed + S2 @ resold2
            resold1 = self._b - self._A @ xnew
            delta = np.linalg.norm(resold1)

            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")
        print(f"Número de iterações: {iter}")
        print(f"Resíduo final: {delta:.3e}")

        return deltaresold