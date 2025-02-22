import numpy as np
import line_profiler

@line_profiler.profile
def teste():
    H = A = V = 9
    nX = nY = nZ = 3
    add_x_left = add_x_right = add_y_left = add_y_right = add_z_left = add_z_right = 0

    coord_grossa = np.empty((H * A * V, 3), dtype=object)

    grossa_size = (H // nX) * (V // nY) * (A // nZ)

    grossa_centroids_array = np.empty(grossa_size, dtype=object)

    grossa_centroids_array[:] = [np.empty(0, dtype=int)] * grossa_size

    dual = np.zeros(H*A*V, dtype=int)

    # Cria uma grade de valores de i, j e k
    i_vals = np.arange(H)
    j_vals = np.arange(V)
    k_vals = np.arange(A)
    i_grid, j_grid, k_grid = np.meshgrid(i_vals, j_vals, k_vals, indexing='ij')

    # Calcula indX, indY e indZ vetorizados
    indX = i_grid // nX
    indY = j_grid // nY
    indZ = k_grid // nZ

    # Calcula os nós
    nodes = i_grid * (A * V) + j_grid * A + k_grid
    nodes = nodes.flatten()

    # Atribui a coordenada grossa para cada nó
    coord_grossa[nodes, :] = np.column_stack((indX.flatten(), indY.flatten(), indZ.flatten()))

    # Prepara grossa_centroids_array usando uma lista de compreensão
    grossa_centroids_array = np.empty(grossa_size, dtype=object)
    grossa_centroids_array[:] = [np.empty(0, dtype=int)] * grossa_size

    # Atualiza grossa_centroids_array com os nós, de forma vetorizada
    idx = indZ + indY * (A // nZ) + indX * ((A // nZ) * (V // nY))
    idx = idx.flatten()
    aux = np.arange(len(idx))
    # Ordena idx para agrupar elementos iguais e reordena aux da mesma maneira
    sorted_idx = np.argsort(idx, kind='stable')
    sorted_aux = aux[sorted_idx]
    # Encontra os limites para cada grupo dentro do array ordenado
    split_indices = np.cumsum(np.bincount(idx)[:-1])
    # Divide aux nos grupos correspondentes a cada índice único de idx
    grossa_centroids_array = np.split(sorted_aux, split_indices)
    # Calcula x, y e z para cada índice
    indices = np.arange((H // nX) * (V // nY) * (A // nZ))
    y = indices // (A // nZ) % (V // nY)
    z = indices % (A // nZ)
    x = indices // ((A // nZ) * (V // nY))
    vertices = []

    # Pontas (Adiciona diretamente sem precisar buscar os valores (são únicos e conhecidos))
    vertices.extend([0, A-1, (V-1)*A, V*A-1, (H-1)*V*A, (H-1)*V*A + A-1, V*A*H - A, V*A*H - 1])
    # Intermediários bordas
    vertices.extend(grossa_centroids_array[i][(nZ + add_z_left) * (nY // 2)] for i in indices[(x == 0) & (y > 0) & (y < (V // nY - 1)) & (z == 0)])
    vertices.extend(grossa_centroids_array[i][(nZ + add_z_left) * (nY + add_y_left) * (nX // 2)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == 0) & (z == 0)])
    vertices.extend(grossa_centroids_array[i][(nZ + add_z_left) * (nY + add_y_right) * (nX // 2) + ((nZ + add_z_left) * (nY + add_y_right - 1))] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == (V // nY - 1)) & (z == 0)])
    vertices.extend(grossa_centroids_array[i][((nZ + add_z_left) * nY * (nX + add_x_right - 1) + ((nZ + add_z_left) * (nY // 2)))] for i in indices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z == 0)])
    vertices.extend(grossa_centroids_array[i][((nZ + add_z_right) * (nY // 2)) + ((nZ + add_z_right) - 1)] for i in indices[(x == 0) & (y > 0) & (y < (V // nY - 1)) & (z == (A // nZ - 1))])
    vertices.extend(grossa_centroids_array[i][((nZ + add_z_right) * (nY + add_y_left) * (nX // 2)) + (nZ + add_z_right - 1)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == 0) & (z == (A // nZ - 1))])
    vertices.extend(grossa_centroids_array[i][((nZ + add_z_right) * (nY + add_y_right) * (nX // 2)) + ((nZ + add_z_right) * (nY + add_y_right) - 1)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == (V // nY - 1)) & (z == (A // nZ - 1))])
    vertices.extend(grossa_centroids_array[i][((nZ + add_z_right) * nY * (nX + add_x_right - 1)) + ((nZ + add_z_right) * (nY // 2)) + (nZ + add_z_right - 1)] for i in indices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z == (A // nZ - 1))])
    vertices.extend(grossa_centroids_array[i][nZ // 2] for i in indices[(x == 0) & (y == 0) & (z < (A // nZ - 1)) & (z > 0)])
    vertices.extend(grossa_centroids_array[i][(nZ * (nY + add_y_right)) - (nZ // 2 + 1)] for i in indices[(x == 0) & (y == (V // nY - 1)) & (z < (A // nZ - 1)) & (z > 0)])
    vertices.extend(grossa_centroids_array[i][(nZ * (nY + add_y_left) * (nX + add_x_right - 1)) + (nZ // 2)] for i in indices[(x == (H // nX - 1)) & (y == 0) & (z < (A // nZ - 1)) & (z > 0)])
    vertices.extend(grossa_centroids_array[i][(nZ * (nY + add_y_right) * (nX + add_x_right)) - (nZ // 2 + 1)] for i in indices[(x == (H // nX - 1)) & (y == (V // nY - 1)) & (z < (A // nZ - 1)) & (z > 0)])
    # Intermediários laterais (Faces)
    vertices.extend(grossa_centroids_array[i][(nZ * nY) // 2] for i in indices[(x == 0) & (y > 0) & (y < (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))])
    vertices.extend(grossa_centroids_array[i][(nZ * nY * (nX + add_x_right - 1)) + (nZ * (nY // 2) + (nZ // 2))] for i in indices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))])
    vertices.extend(grossa_centroids_array[i][(nZ * (nY + add_y_left) * (nX // 2)) + (nZ // 2)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == 0) & (z > 0) & (z < (A // nZ - 1))])
    vertices.extend(grossa_centroids_array[i][(nZ * (nY + add_y_right) * (nX // 2 + 1)) - (nZ // 2 + 1)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))])
    vertices.extend(grossa_centroids_array[i][((nZ + add_z_left) * nY * (nX // 2)) + ((nZ + add_z_left) * (nY // 2))] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y < (V // nY - 1)) & (y > 0) & (z == 0)])
    vertices.extend(grossa_centroids_array[i][((nZ + add_z_right) * nY * (nX // 2)) + ((nZ + add_z_right) * (nY // 2)) + nZ + add_z_right - 1] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y < (V // nY - 1)) & (y > 0) & (z == (A // nZ - 1))])
    # Intermediários internos
    vertices.extend(grossa_centroids_array[i][(nZ * nY * nX) // 2] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))])
    vertices = np.array(vertices)
    dual[vertices] = 1

    # # Adiciona valores ao conjunto
    # ja_foi_mapeado = set(vertices)
    # vertices_borda = []
    # arestas = []

    # # Preciso pegar as arestas
    # # Borda x = 0 y = 0
    # for index in range(0, A):
    #     if index in ja_foi_mapeado:
    #         vertices_borda.append(index)
    #         continue
    #     ja_foi_mapeado.add(index)
    #     arestas.append(index)

    # # Borda x = 0 z = 0
    # for index in range(0, V*A, A):
    #     if index in ja_foi_mapeado:
    #         vertices_borda.append(index)
    #         continue
    #     ja_foi_mapeado.add(index)
    #     arestas.append(index)

    # # Borda x = 0 y = V
    # for index in range((V - 1)*A, V * A):
    #     if index in ja_foi_mapeado:
    #         vertices_borda.append(index)
    #         continue
    #     ja_foi_mapeado.add(index)
    #     arestas.append(index)

    # # Borda x = 0 z = A
    # for index in range(A - 1, V * A, A):
    #     if index in ja_foi_mapeado:
    #         vertices_borda.append(index)
    #         continue
    #     ja_foi_mapeado.add(index)
    #     arestas.append(index)

    # # Borda x = H y = 0
    # for index in range(V * A * (H - 1), (V * A * (H - 1)) + A):
    #     if index in ja_foi_mapeado:
    #         vertices_borda.append(index)
    #         continue
    #     ja_foi_mapeado.add(index)
    #     arestas.append(index)

    # # Borda x = H z = 0
    # for index in range(V * A * (H - 1), (V * A * (H - 1) + (A * (V - 1)) + 1), A):
    #     if index in ja_foi_mapeado:
    #         vertices_borda.append(index)
    #         continue
    #     ja_foi_mapeado.add(index)
    #     arestas.append(index)

    # # Borda x = H y = V
    # for index in range((V * A * (H - 1)) + (A * (V - 1)), V * A * H):
    #     if index in ja_foi_mapeado:
    #         vertices_borda.append(index)
    #         continue
    #     ja_foi_mapeado.add(index)
    #     arestas.append(index)

    # # Borda x = H z = A
    # for index in range((V * A * (H - 1)) + (A - 1), V * A * H, A):
    #     if index in ja_foi_mapeado:
    #         vertices_borda.append(index)
    #         continue
    #     ja_foi_mapeado.add(index)
    #     arestas.append(index)

    # # Borda y = 0 z = 0
    # for index in range(0, V * A * (H - 1) + 1, V * A):
    #     if index in ja_foi_mapeado:
    #         vertices_borda.append(index)
    #         continue
    #     ja_foi_mapeado.add(index)
    #     arestas.append(index)

    # # Borda y = 0 z = A
    # for index in range(A-1, (V * A * (H - 1)) + A, V * A):
    #     if index in ja_foi_mapeado:
    #         vertices_borda.append(index)
    #         continue
    #     ja_foi_mapeado.add(index)
    #     arestas.append(index)

    # # Borda y = V z = 0
    # for index in range((V - 1) * A, V * A * (H - 1) + ((V - 1) * A) + 1, V * A):
    #     if index in ja_foi_mapeado:
    #         vertices_borda.append(index)
    #         continue
    #     ja_foi_mapeado.add(index)
    #     arestas.append(index)

    # # Borda y = V z = A
    # for index in range((V * A) - 1, V * A * H, V * A):
    #     if index in ja_foi_mapeado:
    #         vertices_borda.append(index)
    #         continue
    #     ja_foi_mapeado.add(index)
    #     arestas.append(index)



    print(np.arange(H*V*A)[dual == 1])
    np.save('Malha_Dual.npy', dual)

teste()

# def generate_centroids_areas_and_adjacencies(mesh):
#     nb=mesh['n_blocks']
#     lb=mesh['block_size']
#     sp=mesh['starting_point']
#     ms = np.mgrid[sp[0]+0.5*lb[0]:sp[0]+(nb[0]+0.5)*lb[0]:lb[0],
#                   sp[1]+0.5*lb[1]:sp[1]+(nb[1]+0.5)*lb[1]:lb[1],
#                   sp[2]+0.5*lb[2]:sp[2]+(nb[2]+0.5)*lb[2]:lb[2]]

#     ms = ms.flatten()
#     centroids = ms.reshape(3,int(ms.size/3)).T

# mesh_generation_parameters:
# n_blocks: [60, 220, 1]
# # n_blocks: [120, 120, 1]
# block_size: [1, 1, 1]
# starting_point: [0.0, 0.0, 0.0]