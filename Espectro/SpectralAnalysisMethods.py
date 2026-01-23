import numpy as np
import scipy.io as sio
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
import time

class SpectralAnalysisMethods:
    def __init__(self, Tfine, Tcorse, b, OP, OR):
        self._Tfine = Tfine
        self._Tcorse = Tcorse
        self._b = b
        self._OP = OP
        self._OR = OR

    def Jacobi(self, args):
        # ----- MÉTODO DE JACOBI -----
        nnod, upperTfine, lowerTfine, diagTfine = args
        precond_jacobi = la.inv(diagTfine)                      # pré-condicionador diagonal
        precond_system = precond_jacobi @ self._Tfine                 # matriz pré-condicionada
        precond_b = precond_jacobi @ self._b                          # vetor dos carregamentos externos pré-condicionado
        posto_jacobi = np.linalg.matrix_rank(precond_system)    # posto da matriz pré-condicionada
        cond_jacobi = np.linalg.cond(precond_system)            # condicionamento da matriz de Jacobi
        av_jacobi = la.eigvals(precond_system)                  # espectro da matriz de Jacobi
        esp_jacobi = np.max(np.abs(av_jacobi))                  # raio espectral da matriz de Jacobi
        error_operator = np.eye(nnod) - (precond_jacobi @ self._Tfine)
        av_error_operator = la.eigvals(error_operator)          # espectro da matriz de iteração

        # Resolvendo o sistema linear
        t0 = time.time()

        iter = 1
        itermax = 5000
        tolerancia = 1.0e-3

        xold = np.zeros(nnod)
        xnew = xold.copy()
        resold = self._b - self._Tfine @ xold
        delta = np.linalg.norm(resold)
        deltaresold = [delta]
        S = la.inv(diagTfine)

        while delta > tolerancia and iter < itermax:
            xnew = xold + S @ resold
            resold = self._b - self._Tfine @ xnew
            delta = np.linalg.norm(resold)

            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")
        print(f"Número de iterações: {iter}")
        print(f"Resíduo final: {delta:.3e}")

        return av_error_operator, deltaresold

    def Seidel(self, args):
        # ----- MÉTODO DE GAUSS-SEIDEL -----
        nnod, upperTfine, lowerTfine, diagTfine = args
        precond_seidel = la.inv(diagTfine + lowerTfine)         # pré-condicionador Gauss-Seidel
        precond_system = precond_seidel @ self._Tfine                 # matriz pré-condicionada
        precond_b = precond_seidel @ self._b                          # vetor dos carregamentos externos pré-condicionado
        posto_seidel = np.linalg.matrix_rank(precond_system)    # posto da matriz pré-condicionada
        cond_seidel = np.linalg.cond(precond_system)            # condicionamento da matriz de Seidel
        av_seidel = la.eigvals(precond_system)                  # espectro da matriz de Seidel
        esp_seidel = np.max(np.abs(av_seidel))                  # raio espectral da matriz de Seidel
        error_operator = np.eye(nnod) - (precond_seidel @ self._Tfine)
        av_error_operator = la.eigvals(error_operator)          # espectro da matriz de iteração

        # Resolvendo o sistema linear
        t0 = time.time()

        iter = 1
        itermax = 5000
        tolerancia = 1.0e-3

        xold = np.zeros(nnod)
        xnew = xold.copy()
        resold = self._b - self._Tfine @ xold
        delta = np.linalg.norm(resold)
        deltaresold = [delta]
        S = la.inv(diagTfine + lowerTfine)

        while delta > tolerancia and iter < itermax:
            xnew = xold + S @ resold
            resold = self._b - self._Tfine @ xnew
            delta = np.linalg.norm(resold)

            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")
        print(f"Número de iterações: {iter}")
        print(f"Resíduo final: {delta:.3e}")

        return av_error_operator, deltaresold

    def ILUfac(self, args):
        # ----- MATRIZ Tfine COM SUAVIZADOR ILU(0) -----
        nnod, upperTfine, lowerTfine, diagTfine = args
        Tfine_sp = sp.csc_matrix(self._Tfine)
        ilu = spla.spilu(Tfine_sp, drop_tol=0.0, fill_factor=1.0)
        precond_ilu = np.zeros((nnod, nnod))                    # pré-condicionador ILU
        for j in range(nnod):
            ej = np.zeros(nnod)   
            ej[j] = 1.0
            precond_ilu[:, j] = ilu.solve(ej)
        precond_system = precond_ilu @ self._Tfine              # matriz pré-condicionada
        precond_b = precond_ilu @ self._b                       # vetor dos carregamentos externos pré-condicionado
        posto_ilu = np.linalg.matrix_rank(precond_system)       # posto da matriz pré-condicionada
        cond_ilu = np.linalg.cond(precond_system)               # condicionamento da matriz de iteração
        av_ilu = la.eigvals(precond_system)                     # espectro da matriz de iteração
        esp_ilu = np.max(np.abs(av_ilu))                        # raio espectral da matriz de iteração
        error_operator = np.eye(nnod) - precond_system
        av_error_operator = la.eigvals(error_operator)          # espectro da matriz de iteração

        # Resolvendo o sistema linear
        t0 = time.time()

        xold = np.zeros(nnod)
        xnew = xold.copy()
        iter = 1
        itermax = 5000
        resold = self._b - self._Tfine @ xold
        tolerancia = 1.0e-3
        delta = np.linalg.norm(resold)
        deltaresold = [delta]
        S = precond_ilu

        while delta > tolerancia and iter < itermax:
            xnew = xold + S @ resold
            resold = self._b - self._Tfine @ xnew
            delta = np.linalg.norm(resold)
            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")

        return av_error_operator, deltaresold
    
    def Multiscale(self, args):
        # ----- MATRIZ Tfine com pré-condicionador Multiescala -----
        nnod, upperTfine, lowerTfine, diagTfine = args
        precond_multiscale = self._OP @ la.inv(self._OR @ self._Tfine @ self._OP) @ self._OR                  # matriz de pré-condicionamento M^-1 Multiescala
        precond_system = precond_multiscale @ self._Tfine                             # matriz pré-condicionada
        precond_b = precond_multiscale @ self._b                                      # vetor dos carregamentos externos pré-condicionado
        posto_Tfine = np.linalg.matrix_rank(self._Tfine)                              # posto da matriz de transmissibilidade
        posto_precond_multiscale = np.linalg.matrix_rank(precond_multiscale)    # posto do pré-condicionador
        posto_multiscale = np.linalg.matrix_rank(precond_system)                # posto da matriz pré-condicionada
        cond_Tfine = np.linalg.cond(self._Tfine)                                      # condicionamento da matriz de transmissibilidade
        cond_multiscale = np.linalg.cond(precond_system)                        # condicionamento da matriz pré-condicionada
        av_multiscale = la.eigvals(precond_system)                              # espectro da matriz pré-condicionada
        esp_multiscale = np.max(np.abs(av_multiscale))                          # raio espectral da matriz pré-condicionada
        av_error_operator = la.eigvals(np.eye(nnod) - precond_system)           # espectro da matriz de iteração
        espec_error_operator = np.max(np.abs(av_error_operator))                # raio espectral do operador de propagação de erro

        # Resolvendo o sistema linear
        t0 = time.time()

        iter = 1
        itermax = 5000
        tolerancia = 1.0e-3

        xold = np.zeros(nnod)
        xnew = xold.copy()
        resold = self._b - self._Tfine @ xold
        delta = np.linalg.norm(resold)
        deltaresold = [delta]
        S = precond_multiscale

        while delta > tolerancia and iter < itermax:
            xnew = xold + S @ resold
            resold = self._b - self._Tfine @ xnew
            delta = np.linalg.norm(resold)

            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")
        print(f"Número de iterações: {iter}")
        print(f"Resíduo final: {delta:.3e}")

        return av_error_operator, deltaresold

    def Multiscale_ILUfac(self, args):
        # ----- MATRIZ Tfine pré-condicionador Multiescala + SUAVIZADOR ILU(0) -----
        nnod, upperTfine, lowerTfine, diagTfine = args
        Tfine_sp = sp.csc_matrix(self._Tfine)
        ilu = spla.spilu(Tfine_sp, drop_tol=0.0, fill_factor=1.0)
        precond_ilu = np.zeros((nnod, nnod))                        # pré-condicionador ILU
        for j in range(nnod):
            ej = np.zeros(nnod)   
            ej[j] = 1.0
            precond_ilu[:, j] = ilu.solve(ej)
        precond_system = precond_ilu @ self._Tfine                  # matriz pré-condicionada
        precond_b = precond_ilu @ self._b                           # vetor dos carregamentos externos pré-condicionado
        cond_ilu = np.linalg.cond(precond_system)                   # condicionamento da matriz pré-condicionada
        av_ilu = la.eigvals(precond_system)                         # espectro da matriz pré-condicionada
        esp_ilu = np.max(np.abs(av_ilu))                            # raio espectral da matriz pré-condicionada
        error_operator = np.eye(nnod) - precond_system
        av_error_operator = la.eigvals(error_operator)              # espectro da matriz de iteração
        espec_error_operator = np.max(np.abs(av_error_operator))    # raio espectral da matriz pré-condicionada

        precond_multiscale = self._OP @ la.inv(self._OR @ self._Tfine @ self._OP) @ self._OR                  # matriz de pré-condicionamento M^-1 Multiescala
        precond_system = precond_multiscale @ self._Tfine                             # matriz pré-condicionada
        cond_multiscale = np.linalg.cond(precond_system)                        # condicionamento da matriz Mmulti
        av_multiscale = la.eigvals(precond_system)                              # espectro da matriz Mmulti
        esp_multiscale = np.max(np.abs(av_multiscale))                          # raio espectral da matriz Mmulti
        av_error_operator = la.eigvals(np.eye(nnod) - precond_system)           # espectro da matriz Tfine
        espec_error_operator = np.max(np.abs(av_error_operator))                # raio espectral do operador de propagação de erro

        precond_multiscaleilu = (precond_multiscale + precond_ilu - self._Tfine @ precond_multiscale @ precond_ilu)
        precond_system_milu = precond_multiscaleilu @ self._Tfine
        av_precond_multiscaleilu = la.eigvals(precond_system_milu)
        espec_precond_multiscaleilu = np.max(np.abs(av_precond_multiscaleilu))
        posto_multiscaleilu = np.linalg.matrix_rank(precond_system_milu)
        cond_multiscale_ilu = np.linalg.cond(precond_system_milu)
        av_error_operator_milu = la.eigvals(np.eye(nnod) - precond_system_milu)
        espec_error_operator_milu = np.max(np.abs(av_error_operator_milu))

        # Resolvendo o sistema linear
        t0 = time.time()

        iter = 1
        itermax = 1000
        tolerancia = 1.0e-3

        xold = np.zeros(nnod)
        xnew = xold.copy()
        resold1 = self._b - self._Tfine @ xold
        delta = np.linalg.norm(resold1)
        deltaresold = [delta]
        S1 = precond_multiscale
        S2 = precond_ilu

        while delta > tolerancia and iter < itermax:
            xmed = xold + S1 @ resold1
            resold2 = self._b - self._Tfine @ xnew
            xnew = xmed + S2 @ resold2
            resold1 = self._b - self._Tfine @ xnew
            delta = np.linalg.norm(resold1)

            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")
        print(f"Número de iterações: {iter}")
        print(f"Resíduo final: {delta:.3e}")

        return av_error_operator_milu, deltaresold

    def Multiscale_Jacobi(self, args):
        # ----- MATRIZ Tfine pré-condicionador Multiescala + Jacobi -----
        nnod, upperTfine, lowerTfine, diagTfine = args
        precond_jacobi = la.inv(diagTfine)                      # pré-condicionador diagonal
        precond_system = precond_jacobi @ self._Tfine                 # matriz pré-condicionada
        precond_b = precond_jacobi @ self._b                          # vetor dos carregamentos externos pré-condicionado
        posto_jacobi = np.linalg.matrix_rank(precond_system)    # posto da matriz pré-condicionada
        cond_jacobi = np.linalg.cond(precond_system)            # condicionamento da matriz de Jacobi
        av_jacobi = la.eigvals(precond_system)                  # espectro da matriz de Jacobi
        esp_jacobi = np.max(np.abs(av_jacobi))                  # raio espectral da matriz de Jacobi
        
        precond_multiscale = self._OP @ la.inv(self._OR @ self._Tfine @ self._OP) @ self._OR                  # matriz de pré-condicionamento M^-1 Multiescala
        precond_system = precond_multiscale @ self._Tfine                             # matriz pré-condicionada
        precond_b = precond_multiscale @ self._b                                      # vetor dos carregamentos externos pré-condicionado
        cond_multiscale = np.linalg.cond(precond_system)                        # condicionamento da matriz Mmulti
        av_multiscale = la.eigvals(precond_system)                              # espectro da matriz Mmulti
        esp_multiscale = np.max(np.abs(av_multiscale))                          # raio espectral da matriz Mmulti

        precond_multiscalejacobi = (precond_multiscale + precond_jacobi - self._Tfine @ precond_multiscale @ precond_jacobi)
        precond_system_mjacobi = precond_multiscalejacobi @ self._Tfine
        cond_multiscale_jacobi = np.linalg.cond(precond_system_mjacobi)
        av_error_operator_mjacobi = la.eigvals(np.eye(nnod) - precond_system_mjacobi)

        # Resolvendo o sistema linear
        t0 = time.time()

        iter = 1
        itermax = 1000
        tolerancia = 1.0e-3

        xold = np.zeros(nnod)
        xnew = xold.copy()
        resold1 = self._b - self._Tfine @ xold
        delta = np.linalg.norm(resold1)
        deltaresold = [delta]
        S1 = precond_multiscale
        S2 = precond_jacobi

        while delta > tolerancia and iter < itermax:
            xmed = xold + S1 @ resold1
            resold2 = self._b - self._Tfine @ xnew
            xnew = xmed + S2 @ resold2
            resold1 = self._b - self._Tfine @ xnew
            delta = np.linalg.norm(resold1)

            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")
        print(f"Número de iterações: {iter}")
        print(f"Resíduo final: {delta:.3e}")

        return av_error_operator_mjacobi, deltaresold

    def Multiscale_Seidel(self, args):
        nnod, upperTfine, lowerTfine, diagTfine = args
        precond_seidel = la.inv(diagTfine + lowerTfine)         # pré-condicionador Gauss-Seidel
        precond_system = precond_seidel @ self._Tfine                 # matriz pré-condicionada
        precond_b = precond_seidel @ self._b                          # vetor dos carregamentos externos pré-condicionado
        posto_seidel = np.linalg.matrix_rank(precond_system)    # posto da matriz pré-condicionada
        cond_seidel = np.linalg.cond(precond_system)            # condicionamento da matriz de Seidel
        av_seidel = la.eigvals(precond_system)                  # espectro da matriz de Seidel
        esp_seidel = np.max(np.abs(av_seidel))                  # raio espectral da matriz de Seidel

        precond_multiscale = self._OP @ la.inv(self._OR @ self._Tfine @ self._OP) @ self._OR              # matriz de pré-condicionamento M^-1 Multiescala
        precond_system = precond_multiscale @ self._Tfine                         # matriz pré-condicionada
        precond_b = precond_multiscale @ self._b                                  # vetor dos carregamentos externos pré-condicionado
        cond_multiscale = np.linalg.cond(precond_system)                    # condicionamento da matriz Mmulti
        av_multiscale = la.eigvals(precond_system)                          # espectro da matriz Mmulti
        esp_multiscale = np.max(np.abs(av_multiscale))                      # raio espectral da matriz Mmulti

        precond_multiscaleseidel = (precond_multiscale + precond_seidel - self._Tfine @ precond_multiscale @ precond_seidel)
        precond_system_mseidel = precond_multiscaleseidel @ self._Tfine
        cond_multiscale_seidel = np.linalg.cond(precond_system_mseidel)
        posto_multiscaleseidel = np.linalg.matrix_rank(precond_system_mseidel)
        av_error_operator = la.eigvals(np.eye(nnod) - precond_system_mseidel)

        # Resolvendo o sistema linear
        t0 = time.time()

        iter = 1
        itermax = 1000
        tolerancia = 1.0e-3

        xold = np.zeros(nnod)
        xnew = xold.copy()
        resold1 = self._b - self._Tfine @ xold
        delta = np.linalg.norm(resold1)
        deltaresold = [delta]
        S1 = precond_multiscale
        S2 = precond_seidel

        while delta > tolerancia and iter < itermax:
            xmed = xold + S1 @ resold1
            resold2 = self._b - self._Tfine @ xnew
            xnew = xmed + S2 @ resold2
            resold1 = self._b - self._Tfine @ xnew
            delta = np.linalg.norm(resold1)

            xold = xnew.copy()
            iter += 1
            deltaresold.append(delta)

        t_final = time.time() - t0
        print(f"Tempo total (s): {t_final:.4f}")
        print(f"Número de iterações: {iter}")
        print(f"Resíduo final: {delta:.3e}")

        return av_error_operator, deltaresold
