import numpy as np
import scipy.io as sio
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
import time

# ----- Carregando o problema -----
# Caso 1
data = sio.loadmat("Projeto-IC/arquivos/Problema_heterogeneo_buck_30x30.mat")
# Caso 2
#data = sio.loadmat("Problema_homogeneo_buck_40x40.mat")

#Tfine = np.asarray(data["Tcorse"])                         # matriz de transmissibilidade na malha grossa
#b = np.asarray(data["OR"]) @ np.asarray(data["F"]).ravel   # vetor dos carregamentos externos na malha grossa

Tfine = np.asarray(data["T"])                               # matriz de transmissibilidade na malha fina
b = np.asarray(data["F"]).ravel()                           # vetor dos carregamentos externos na malha fina
OP = np.asarray( data['OP'])
OR = np.asarray(data['OR'])

# ----- Verificando a estrutura da matriz (TPFA) -----
plt.figure()
plt.spy(Tfine, markersize=1)
plt.title("Estrutura da matriz de transmissibilidade fina")
plt.show()

Tcorse = np.asarray(data["Tcorse"])
plt.figure()
plt.spy(Tfine, markersize=1)
plt.title("Estrutura da matriz de transmissibilidade grossa")
plt.show()

# ----- Análise espectral da matriz de transmissibilidade pré-condicionada -----
nnod = Tfine.shape[0]               # número de graus de liberdade
cond_T = np.linalg.cond(Tfine)      # condicionamento da matriz Tfine
av_T, VTfine = la.eig(Tfine)        # espectro da matriz Tfina e espectro da matriz Tfina
esp_T = np.max(np.abs(av_T))        # raio espectral da matriz Tfine

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
plt.ylim([-1, 1])   ##
plt.grid(True)
plt.show()

# ----- Decomposição Tfine = D + L + U -----
upperTfine = np.tril(Tfine, k=-1)       # matriz triangular inferior
lowerTfine = np.triu(Tfine, k=1)        # matriz triangular superior
diagTfine  = np.diag(np.diag(Tfine))    # matriz diagonal
#diagTfine = Tfine - (upperTfine + lowerTfine)

# ----- Selecionando o método que vai ser utilizado -----
# method = 'Jacobi'
#method = 'Seidel'
#method = 'ILUfac'
#method = 'Multiscale'
# method = 'MultiscaleILUfac'
#method = 'MultiscaleJacobi'
method = 'MultiscaleSeidel'

if method == "Jacobi":
    # ----- MÉTODO DE JACOBI -----
    precond_jacobi = la.inv(diagTfine)                      # pré-condicionador diagonal
    precond_system = precond_jacobi @ Tfine                 # matriz pré-condicionada
    precond_b = precond_jacobi @ b                          # vetor dos carregamentos externos pré-condicionado
    posto_jacobi = np.linalg.matrix_rank(precond_system)    # posto da matriz pré-condicionada
    cond_jacobi = np.linalg.cond(precond_system)            # condicionamento da matriz de Jacobi
    av_jacobi = la.eigvals(precond_system)                  # espectro da matriz de Jacobi
    esp_jacobi = np.max(np.abs(av_jacobi))                  # raio espectral da matriz de Jacobi
    error_operator = np.eye(nnod) - precond_system
    av_error_operator = la.eigvals(error_operator)          # espectro da matriz de iteração

    # Plot do espectro
    plt.figure()
    plt.plot(np.real(av_error_operator), np.imag(av_error_operator), '*r', label='Autovalores')
    plt.xlabel('Re(lambda)')
    plt.ylabel('Imag(lambda)')
    plt.title('Espectro: pré-condicionador Diagonal (Jacobi)')
    x = np.arange(-1.0, 1.01, 0.01)
    z = np.sqrt(1 - x**2)
    y = -np.sqrt(1 - x**2)
    plt.plot(x, z, 'b')
    plt.plot(x, y, 'b')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

    # Resolvendo o sistema linear
    t0 = time.time()

    iter = 1
    itermax = 5000
    tolerancia = 1.0e-3

    xold = np.zeros(nnod)
    xnew = xold.copy()
    resold = b - Tfine @ xold
    delta = np.linalg.norm(resold)
    deltaresold = [delta]
    S = la.inv(diagTfine)

    while delta > tolerancia and iter < itermax:
        xnew = xold + S @ resold
        resold = b - Tfine @ xnew
        delta = np.linalg.norm(resold)

        xold = xnew.copy()
        iter += 1
        deltaresold.append(delta)

    t_final = time.time() - t0
    print(f"Tempo total (s): {t_final:.4f}")
    print(f"Número de iterações: {iter}")
    print(f"Resíduo final: {delta:.3e}")

    # Plot do resíduo
    plt.figure()
    plt.plot(deltaresold)
    plt.xlabel("Iteração")
    plt.ylabel("Resíduo")
    plt.title("Resíduo Jacobi")
    plt.yscale("log")
    plt.grid(True)
    plt.show()

elif method == "Seidel":
    # ----- MÉTODO DE GAUSS-SEIDEL -----
    precond_seidel = la.inv(diagTfine + lowerTfine)         # pré-condicionador Gauss-Seidel
    precond_system = precond_seidel @ Tfine                 # matriz pré-condicionada
    precond_b = precond_seidel @ b                          # vetor dos carregamentos externos pré-condicionado
    posto_seidel = np.linalg.matrix_rank(precond_system)    # posto da matriz pré-condicionada
    cond_seidel = np.linalg.cond(precond_system)            # condicionamento da matriz de Seidel
    av_seidel = la.eigvals(precond_system)                  # espectro da matriz de Seidel
    esp_seidel = np.max(np.abs(av_seidel))                  # raio espectral da matriz de Seidel
    error_operator = np.eye(nnod) - precond_system
    av_error_operator = la.eigvals(error_operator)          # espectro da matriz de iteração

    # Plot do espectro
    plt.figure()
    plt.plot(np.real(av_error_operator), np.imag(av_error_operator), '*r', label='Autovalores')
    plt.xlabel('Re(lambda)')
    plt.ylabel('Imag(lambda)')
    plt.title('Espectro: pré-condicionador Gauss-Seidel')
    x = np.arange(-1.0, 1.01, 0.01)
    z = np.sqrt(1 - x**2)
    y = -np.sqrt(1 - x**2)
    plt.plot(x, z, 'b')
    plt.plot(x, y, 'b')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

    # Resolvendo o sistema linear
    t0 = time.time()

    iter = 1
    itermax = 5000
    tolerancia = 1.0e-3

    xold = np.zeros(nnod)
    xnew = xold.copy()
    resold = b - Tfine @ xold
    delta = np.linalg.norm(resold)
    deltaresold = [delta]
    S = la.inv(diagTfine + lowerTfine)

    while delta > tolerancia and iter < itermax:
        xnew = xold + S @ resold
        resold = b - Tfine @ xnew
        delta = np.linalg.norm(resold)

        xold = xnew.copy()
        iter += 1
        deltaresold.append(delta)

    t_final = time.time() - t0
    print(f"Tempo total (s): {t_final:.4f}")
    print(f"Número de iterações: {iter}")
    print(f"Resíduo final: {delta:.3e}")

    # Plot do resíduo
    plt.figure()
    plt.plot(deltaresold)
    plt.xlabel("Iteração")
    plt.ylabel("Resíduo")
    plt.title("Resíduo Gauss-Seidel")
    plt.yscale("log")
    plt.grid(True)
    plt.show()

elif method == "ILUfac":
    # ----- MATRIZ Tfine COM SUAVIZADOR ILU(0) -----
    precond_ilu = np.zeros((nnod, nnod))                    # pré-condicionador ILU
    precond_system = precond_ilu @ Tfine                    # matriz pré-condicionada
    precond_b = precond_ilu @ b                             # vetor dos carregamentos externos pré-condicionado
    posto_ilu = np.linalg.matrix_rank(precond_system)       # posto da matriz pré-condicionada
    cond_ilu = np.linalg.cond(precond_system)               # condicionamento da matriz de iteração
    av_ilu = la.eigvals(precond_system)                     # espectro da matriz de iteração
    esp_ilu = np.max(np.abs(av_ilu))                        # raio espectral da matriz de iteração
    error_operator = np.eye(nnod) - precond_system
    av_error_operator = la.eigvals(error_operator)          # espectro da matriz de iteração

    # Plot do espectro
    plt.figure()
    plt.plot(np.real(av_error_operator), np.imag(av_error_operator), '*r', label='Autovalores')
    plt.xlabel('Re(lambda)')
    plt.ylabel('Imag(lambda)')
    plt.title('Espectro: pré-condicionador ILU(0)')
    x = np.arange(-1.0, 1.01, 0.01)
    z = np.sqrt(1 - x**2)
    y = -np.sqrt(1 - x**2)
    plt.plot(x, z, 'b')
    plt.plot(x, y, 'b')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

    # Resolvendo o sistema linear
    t0 = time.time()

    xold = np.zeros(nnod)
    xnew = xold.copy()
    iter = 1
    itermax = 5000
    resold = b - Tfine @ xold
    tolerancia = 1.0e-3
    delta = np.linalg.norm(resold)
    deltaresold = [delta]
    S = precond_ilu

    while delta > tolerancia and iter < itermax:
        xnew = xold + S @ resold
        resold = b - Tfine @ xnew
        delta = np.linalg.norm(resold)
        xold = xnew.copy()
        iter += 1
        deltaresold.append(delta)

    t_final = time.time() - t0
    print(f"Tempo total (s): {t_final:.4f}")

        # Plot do resíduo
    plt.figure()
    plt.plot(deltaresold)
    plt.xlabel("Iteração")
    plt.ylabel("Resíduo")
    plt.title("Resíduo Gauss-Seidel")
    plt.yscale("log")
    plt.grid(True)
    plt.show()

elif method == "Multiscale":
    # ----- MATRIZ Tfine com pré-condicionador Multiescala -----
    precond_multiscale = OP @ la.inv(OR @ Tfine @ OP) @ OR                  # matriz de pré-condicionamento M^-1 Multiescala
    precond_system = precond_multiscale @ Tfine                             # matriz pré-condicionada
    precond_b = precond_multiscale @ b                                      # vetor dos carregamentos externos pré-condicionado
    posto_Tfine = np.linalg.matrix_rank(Tfine)                              # posto da matriz de transmissibilidade
    posto_precond_multiscale = np.linalg.matrix_rank(precond_multiscale)    # posto do pré-condicionador
    posto_multiscale = np.linalg.matrix_rank(precond_system)                # posto da matriz pré-condicionada
    cond_Tfine = np.linalg.cond(Tfine)                                      # condicionamento da matriz de transmissibilidade
    cond_multiscale = np.linalg.cond(precond_system)                        # condicionamento da matriz pré-condicionada
    av_multiscale = la.eigvals(precond_system)                              # espectro da matriz pré-condicionada
    esp_multiscale = np.max(np.abs(av_multiscale))                          # raio espectral da matriz pré-condicionada
    av_error_operator = la.eigvals(np.eye(nnod) - precond_system)           # espectro da matriz de iteração
    espec_error_operator = np.max(np.abs(av_error_operator))                # raio espectral do operador de propagação de erro

    # Plot do espectro
    plt.figure()
    plt.plot(np.real(av_error_operator), np.imag(av_error_operator), '*r', label='Autovalores')
    plt.xlabel('Re(lambda)')
    plt.ylabel('Imag(lambda)')
    plt.title('Espectro: usando pré-condicionador M^-1 multiscale')
    x = np.arange(-1.0, 1.01, 0.01)
    z = np.sqrt(1 - x**2)
    y = -np.sqrt(1 - x**2)
    plt.plot(x, z, 'b')
    plt.plot(x, y, 'b')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

    # Resolvendo o sistema linear
    t0 = time.time()

    iter = 1
    itermax = 5000
    tolerancia = 1.0e-3

    xold = np.zeros(nnod)
    xnew = xold.copy()
    resold = b - Tfine @ xold
    delta = np.linalg.norm(resold)
    deltaresold = [delta]
    S = precond_multiscale

    while delta > tolerancia and iter < itermax:
        xnew = xold + S @ resold
        resold = b - Tfine @ xnew
        delta = np.linalg.norm(resold)

        xold = xnew.copy()
        iter += 1
        deltaresold.append(delta)

    t_final = time.time() - t0
    print(f"Tempo total (s): {t_final:.4f}")
    print(f"Número de iterações: {iter}")
    print(f"Resíduo final: {delta:.3e}")

    # Plot do resíduo
    plt.figure()
    plt.plot(deltaresold)
    plt.xlabel("Iteração")
    plt.ylabel("Resíduo")
    plt.title("Resíduo Multiescala")
    plt.yscale("log")
    plt.grid(True)
    plt.show()

    # GMRES - Generalized Minimum Residual

    itermax = 100
    tolerancia = 1.0e-6
    restart = 16

    residuals = []
    def callback(rk):
        residuals.append(np.linalg.norm(rk))
    
    t0 = time.time()
    xgmres, info = spla.gmres(
        precond_system,
        precond_b,
        restart=restart,
        rtol=tolerancia,
        maxiter=itermax,
        callback=callback
    )

    t_final = time.time() - t0
    print(f"Tempo total (s): {t_final:.4f}")
    print(f"Flag GMRES (info): {info}")
    print(f"Número de iterações: {len(residuals)}")
    print(f"Resíduo final: {residuals[-1]:.3e}")

    # Plot do resíduo - GMRES
    plt.figure()
    plt.plot(np.log10(residuals), 'r-')
    plt.xlabel("Iteração")
    plt.ylabel("Resíduo")
    plt.title("GMRES")
    plt.grid(True)
    plt.show()

elif method == "MultiscaleILUfac":
    # ----- MATRIZ Tfine pré-condicionador Multiescala + SUAVIZADOR ILU(0) -----
    precond_ilu = np.zeros((nnod, nnod))                        # pré-condicionador ILU
    precond_system = precond_ilu @ Tfine                        # matriz pré-condicionada
    precond_b = precond_ilu @ b                                 # vetor dos carregamentos externos pré-condicionado
    cond_ilu = np.linalg.cond(precond_system)                   # condicionamento da matriz pré-condicionada
    av_ilu = la.eigvals(precond_system)                         # espectro da matriz pré-condicionada
    esp_ilu = np.max(np.abs(av_ilu))                            # raio espectral da matriz pré-condicionada
    error_operator = np.eye(nnod) - precond_system
    av_error_operator = la.eigvals(error_operator)              # espectro da matriz de iteração
    espec_error_operator = np.max(np.abs(av_error_operator))    # raio espectral da matriz pré-condicionada

    # Plot do espectro
    plt.figure()
    plt.plot(np.real(av_error_operator), np.imag(av_error_operator), '*r', label='Autovalores')
    plt.xlabel('Re(lambda)')
    plt.ylabel('Imag(lambda)')
    plt.title('Espectro do Operador de Propagação de Erro: pré-condicionador ILU(0)')
    x = np.arange(-1.0, 1.01, 0.01)
    z = np.sqrt(1 - x**2)
    y = -np.sqrt(1 - x**2)
    plt.plot(x, z, 'b')
    plt.plot(x, y, 'b')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

    precond_multiscale = OP @ la.inv(OR @ Tfine @ OP) @ OR                  # matriz de pré-condicionamento M^-1 Multiescala
    precond_system = precond_multiscale @ Tfine                             # matriz pré-condicionada
    cond_multiscale = np.linalg.cond(precond_system)                        # condicionamento da matriz Mmulti
    av_multiscale = la.eigvals(precond_system)                              # espectro da matriz Mmulti
    esp_multiscale = np.max(np.abs(av_multiscale))                          # raio espectral da matriz Mmulti
    av_error_operator = la.eigvals(np.eye(nnod) - precond_system)           # espectro da matriz Tfine
    espec_error_operator = np.max(np.abs(av_error_operator))                # raio espectral do operador de propagação de erro

    # Plot do espectro
    plt.figure()
    plt.plot(np.real(av_error_operator), np.imag(av_error_operator), '*r', label='Autovalores')
    plt.xlabel('Re(lambda)')
    plt.ylabel('Imag(lambda)')
    plt.title('Espectro do Operador de Propagação de Erro: pré-condicionador Multiscale')
    x = np.arange(-1.0, 1.01, 0.01)
    z = np.sqrt(1 - x**2)
    y = -np.sqrt(1 - x**2)
    plt.plot(x, z, 'b')
    plt.plot(x, y, 'b')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

    precond_multiscaleilu = (precond_multiscale + precond_ilu - Tfine @ precond_multiscale @ precond_ilu)
    precond_system_milu = precond_multiscaleilu @ Tfine
    av_precond_multiscaleilu = la.eigvals(precond_system_milu)
    espec_precond_multiscaleilu = np.max(np.abs(av_precond_multiscaleilu))
    posto_multiscaleilu = np.linalg.matrix_rank(precond_system_milu)
    cond_multiscale_ilu = np.linalg.cond(precond_system_milu)
    av_error_operator_milu = la.eigvals(np.eye(nnod) - precond_system_milu)
    espec_error_operator_milu = np.max(np.abs(av_error_operator_milu))

    # Plot do espectro
    plt.figure()
    plt.plot(np.real(av_precond_multiscaleilu), np.imag(av_precond_multiscaleilu), '*r', label='Autovalores')
    plt.xlabel('Re(lambda)')
    plt.ylabel('Imag(lambda)')
    plt.title('Espectro do Sistema Pré-condicionado: multiscale + ILU(0)')
    x = np.arange(-1.0, 1.01, 0.01)
    z = np.sqrt(1 - x**2)
    y = -np.sqrt(1 - x**2)
    plt.plot(x, z, 'b')
    plt.plot(x, y, 'b')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

    plt.figure()
    plt.plot(np.real(av_error_operator_milu), np.imag(av_error_operator_milu), '*r', label='Autovalores')
    plt.xlabel('Re(lambda)')
    plt.ylabel('Imag(lambda)')
    plt.title('Espectro do Sistema Pré-condicionado: multiscale + ILU(0)')
    x = np.arange(-1.0, 1.01, 0.01)
    z = np.sqrt(1 - x**2)
    y = -np.sqrt(1 - x**2)
    plt.plot(x, z, 'b')
    plt.plot(x, y, 'b')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

    # Resolvendo o sistema linear
    t0 = time.time()

    iter = 1
    itermax = 1000
    tolerancia = 1.0e-3

    xold = np.zeros(nnod)
    xnew = xold.copy()
    resold1 = b - Tfine @ xold
    delta = 1
    deltaresold = [delta]
    S1 = precond_multiscale
    S2 = precond_ilu

    while delta > tolerancia and iter < itermax:
        xmed = xold + S1 @ resold1
        resold2 = b - Tfine @ xnew
        xnew = xmed + S2 @ resold2
        resold1 = b - Tfine @ xnew
        delta = np.linalg.norm(resold1)

        xold = xnew.copy()
        iter += 1
        deltaresold.append(delta)

    t_final = time.time() - t0
    print(f"Tempo total (s): {t_final:.4f}")
    print(f"Número de iterações: {iter}")
    print(f"Resíduo final: {delta:.3e}")

    # Plot do resíduo
    plt.figure()
    plt.plot(deltaresold)
    plt.xlabel("Iteração")
    plt.ylabel("Resíduo")
    plt.title("Resíduo Multiescala + ILU")
    plt.yscale("log")
    plt.grid(True)
    plt.show()

    # GMRES - Generalized Minimum Residual
    A_pc = precond_multiscaleilu @ Tfine
    b_pc = precond_multiscaleilu @ b

    itermax = 100
    tolerancia = 1.0e-3
    restart = 16

    residuals = []
    def callback(rk):
        residuals.append(np.linalg.norm(rk))
    
    t0 = time.time()
    xgmres, info = spla.gmres(
        A_pc,
        b_pc,
        restart=restart,
        rtol=tolerancia,
        maxiter=itermax,
        callback=callback
    )

    t_final = time.time() - t0
    print(f"Tempo total (s): {t_final:.4f}")
    print(f"Flag GMRES (info): {info}")
    print(f"Número de iterações: {len(residuals)}")
    print(f"Resíduo final: {residuals[-1]:.3e}")

    # Plot do resíduo - GMRES
    plt.figure()
    plt.plot(np.log10(residuals), 'r-')
    plt.xlabel("Iteração")
    plt.ylabel("Resíduo")
    plt.title("GMRES")
    plt.grid(True)
    plt.show()

elif method == "MultiscaleJacobi":
    precond_jacobi = la.inv(diagTfine)                      # pré-condicionador diagonal
    precond_system = precond_jacobi @ Tfine                 # matriz pré-condicionada
    precond_b = precond_jacobi @ b                          # vetor dos carregamentos externos pré-condicionado
    posto_jacobi = np.linalg.matrix_rank(precond_system)    # posto da matriz pré-condicionada
    cond_jacobi = np.linalg.cond(precond_system)            # condicionamento da matriz de Jacobi
    av_jacobi = la.eigvals(precond_system)                  # espectro da matriz de Jacobi
    esp_jacobi = np.max(np.abs(av_jacobi))                  # raio espectral da matriz de Jacobi
    
    precond_multiscale = OP @ la.inv(OR @ Tfine @ OP) @ OR                  # matriz de pré-condicionamento M^-1 Multiescala
    precond_system = precond_multiscale @ Tfine                             # matriz pré-condicionada
    precond_b = precond_multiscale @ b                                      # vetor dos carregamentos externos pré-condicionado
    cond_multiscale = np.linalg.cond(precond_system)                        # condicionamento da matriz Mmulti
    av_multiscale = la.eigvals(precond_system)                              # espectro da matriz Mmulti
    esp_multiscale = np.max(np.abs(av_multiscale))                          # raio espectral da matriz Mmulti

    precond_multiscalejacobi = (precond_multiscale + precond_jacobi - Tfine @ precond_multiscale @ precond_jacobi)
    precond_system_mjacobi = precond_multiscalejacobi @ Tfine
    cond_multiscale_jacobi = np.linalg.cond(precond_system_mjacobi)
    av_error_operator_mjacobi = la.eigvals(np.eye(nnod) - precond_system_mjacobi)

    plt.figure()
    plt.plot(np.real(av_error_operator_mjacobi), np.imag(av_error_operator_mjacobi), '*r', label='Autovalores')
    plt.xlabel('Re(lambda)')
    plt.ylabel('Imag(lambda)')
    plt.title('Espectro: usando pré-condicionador multiscale + Jacobi')
    x = np.arange(-1.0, 1.01, 0.01)
    z = np.sqrt(1 - x**2)
    y = -np.sqrt(1 - x**2)
    plt.plot(x, z, 'b')
    plt.plot(x, y, 'b')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

    # Resolvendo o sistema linear
    t0 = time.time()

    iter = 1
    itermax = 1000
    tolerancia = 1.0e-3

    xold = np.zeros(nnod)
    xnew = xold.copy()
    resold1 = b - Tfine @ xold
    delta = np.linalg.norm(resold1)
    deltaresold = [delta]
    S1 = precond_multiscale
    S2 = precond_jacobi

    while delta > tolerancia and iter < itermax:
        xmed = xold + S1 @ resold1
        resold2 = b - Tfine @ xnew
        xnew = xmed + S2 @ resold2
        resold1 = b - Tfine @ xnew
        delta = np.linalg.norm(resold1)

        xold = xnew.copy()
        iter += 1
        deltaresold.append(delta)

    t_final = time.time() - t0
    print(f"Tempo total (s): {t_final:.4f}")
    print(f"Número de iterações: {iter}")
    print(f"Resíduo final: {delta:.3e}")

    # Plot do resíduo
    plt.figure()
    plt.plot(deltaresold)
    plt.xlabel("Iteração")
    plt.ylabel("Resíduo")
    plt.title("Resíduo Multiescala + Jacobi")
    plt.yscale("log")
    plt.grid(True)
    plt.show()

elif method == "MultiscaleSeidel":
    precond_seidel = la.inv(diagTfine + lowerTfine)         # pré-condicionador Gauss-Seidel
    precond_system = precond_seidel @ Tfine                 # matriz pré-condicionada
    precond_b = precond_seidel @ b                          # vetor dos carregamentos externos pré-condicionado
    posto_seidel = np.linalg.matrix_rank(precond_system)    # posto da matriz pré-condicionada
    cond_seidel = np.linalg.cond(precond_system)            # condicionamento da matriz de Seidel
    av_seidel = la.eigvals(precond_system)                  # espectro da matriz de Seidel
    esp_seidel = np.max(np.abs(av_seidel))                  # raio espectral da matriz de Seidel

    precond_multiscale = OP @ la.inv(OR @ Tfine @ OP) @ OR              # matriz de pré-condicionamento M^-1 Multiescala
    precond_system = precond_multiscale @ Tfine                         # matriz pré-condicionada
    precond_b = precond_multiscale @ b                                  # vetor dos carregamentos externos pré-condicionado
    cond_multiscale = np.linalg.cond(precond_system)                    # condicionamento da matriz Mmulti
    av_multiscale = la.eigvals(precond_system)                          # espectro da matriz Mmulti
    esp_multiscale = np.max(np.abs(av_multiscale))                      # raio espectral da matriz Mmulti

    precond_multiscaleseidel = (precond_multiscale + precond_seidel - Tfine @ precond_multiscale @ precond_seidel)
    precond_system_mseidel = precond_multiscaleseidel @ Tfine
    cond_multiscale_seidel = np.linalg.cond(precond_system_mseidel)
    posto_multiscaleseidel = np.linalg.matrix_rank(precond_system_mseidel)
    av_error_operator = la.eigvals(np.eye(nnod) - precond_system_mseidel)

    plt.figure()
    plt.plot(np.real(av_error_operator), np.imag(av_error_operator), '*r', label='Autovalores')
    plt.xlabel('Re(lambda)')
    plt.ylabel('Imag(lambda)')
    plt.title('Espectro: usando pré-condicionador multiscale + Jacobi')
    x = np.arange(-1.0, 1.01, 0.01)
    z = np.sqrt(1 - x**2)
    y = -np.sqrt(1 - x**2)
    plt.plot(x, z, 'b')
    plt.plot(x, y, 'b')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

    # Resolvendo o sistema linear
    t0 = time.time()

    iter = 1
    itermax = 1000
    tolerancia = 1.0e-3

    xold = np.zeros(nnod)
    xnew = xold.copy()
    resold1 = b - Tfine @ xold
    delta = np.linalg.norm(resold1)
    deltaresold = [delta]
    S1 = precond_multiscale
    S2 = precond_seidel

    while delta > tolerancia and iter < itermax:
        xmed = xold + S1 @ resold1
        resold2 = b - Tfine @ xnew
        xnew = xmed + S2 @ resold2
        resold1 = b - Tfine @ xnew
        delta = np.linalg.norm(resold1)

        xold = xnew.copy()
        iter += 1
        deltaresold.append(delta)

    t_final = time.time() - t0
    print(f"Tempo total (s): {t_final:.4f}")
    print(f"Número de iterações: {iter}")
    print(f"Resíduo final: {delta:.3e}")

    # Plot do resíduo
    plt.figure()
    plt.plot(deltaresold)
    plt.xlabel("Iteração")
    plt.ylabel("Resíduo")
    plt.title("Resíduo Multiescala + Jacobi")
    plt.yscale("log")
    plt.grid(True)
    plt.show()