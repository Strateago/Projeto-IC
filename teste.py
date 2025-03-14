import numpy as np
import line_profiler

@line_profiler.profile
def teste():
    
    H = 9
    V = 9
    A = 9
    nX = 3
    nY = 3
    nZ = 3

    centroids = generate_centroids_areas_and_adjacencies(H, V, A)

    maxs = centroids.max(axis=0) # [max x max y max z]
    mins = centroids.min(axis=0) # [min x min y min z] -> max - min = H, V, A
    cr1 = np.array([nX, nY, nZ])
    n_blocks = np.array(maxs-mins)
    n_duals = n_blocks//cr1 + 1
    second_line = (n_blocks-(n_duals-1)*cr1)//2
    xd=[]
    for i in range(3):
        xd.append(np.arange(second_line[i]+cr1[i], (n_duals[i]-1)*cr1[i],cr1[i]))
        if len(xd[i])>0:
            xd[i]=np.unique(np.concatenate([[mins[i], maxs[i]],(xd[i]+0.5)]))
        else:
            xd[i]=np.unique(np.array([mins[i], maxs[i]]))
    d=np.zeros(len(centroids))

    for i in range(3):
        d += np.isin(centroids[:,i],xd[i])

    dual=np.int16(d)

    volumes=np.arange(len(centroids))
    Grossa_Index=-np.ones_like(volumes)
    # X primal
    xp=np.array(mins[0]-1)
    xp=np.append(xp,(xd[0][:-1]+xd[0][1:])/2)
    xp=np.append(xp,maxs[0]+1)
    # Y primal
    yp=np.array(mins[1]-1)
    yp=np.append(yp,(xd[1][:-1]+xd[1][1:])/2)
    yp=np.append(yp,maxs[1]+1)
    # Z primal  
    zp=np.array(mins[2]-1)
    zp=np.append(zp,(xd[2][:-1]+xd[2][1:])/2)
    zp=np.append(zp,maxs[2]+1)
    # Gerando os indices Fina->Primal
    count=0
    for i in range(len(xp)-1):
        vx=volumes[(centroids[:,0]>=xp[i]) & (centroids[:,0]<xp[i+1])]
        for j in range(len(yp)-1):
            vy=vx[(centroids[:,1][vx]>=yp[j]) & (centroids[:,1][vx]<yp[j+1])]
            for k in range(len(zp)-1):
                vz=vy[(centroids[:,2][vy]>=zp[k]) & (centroids[:,2][vy]<zp[k+1])]
                Grossa_Index[vz]=count
                count+=1
    
    x_idx = np.digitize(centroids[:, 0], xp) - 1
    y_idx = np.digitize(centroids[:, 1], yp) - 1
    z_idx = np.digitize(centroids[:, 2], zp) - 1

    # Calcula índice linear da célula Grossa_Index diretamente
    Grossa_Index1 = x_idx * (len(yp) - 1) * (len(zp) - 1) + y_idx * (len(zp) - 1) + z_idx
    print(np.unique(Grossa_Index == Grossa_Index1))

    # add_x_left = add_x_right = add_y_left = add_y_right = add_z_left = add_z_right = 0

    # coord_grossa = np.empty((H * A * V, 3), dtype=object)

    # grossa_size = (H // nX) * (V // nY) * (A // nZ)

    # grossa_centroids_array = np.empty(grossa_size, dtype=object)

    # grossa_centroids_array[:] = [np.empty(0, dtype=int)] * grossa_size

    # dual = np.zeros(H*A*V, dtype=int)

    # # Cria uma grade de valores de i, j e k
    # i_vals = np.arange(H)
    # j_vals = np.arange(V)
    # k_vals = np.arange(A)
    # i_grid, j_grid, k_grid = np.meshgrid(i_vals, j_vals, k_vals, indexing='ij')

    # # Calcula indX, indY e indZ vetorizados
    # indX = i_grid // nX
    # indY = j_grid // nY
    # indZ = k_grid // nZ

    # # Calcula os nós
    # nodes = i_grid * (A * V) + j_grid * A + k_grid
    # nodes = nodes.flatten()

    # # Atribui a coordenada grossa para cada nó
    # coord_grossa[nodes, :] = np.column_stack((indX.flatten(), indY.flatten(), indZ.flatten()))

    # # Prepara grossa_centroids_array usando uma lista de compreensão
    # grossa_centroids_array = np.empty(grossa_size, dtype=object)
    # grossa_centroids_array[:] = [np.empty(0, dtype=int)] * grossa_size

    # # Atualiza grossa_centroids_array com os nós, de forma vetorizada
    # idx = indZ + indY * (A // nZ) + indX * ((A // nZ) * (V // nY))
    # idx = idx.flatten()
    # aux = np.arange(len(idx))
    # # Ordena idx para agrupar elementos iguais e reordena aux da mesma maneira
    # sorted_idx = np.argsort(idx, kind='stable')
    # sorted_aux = aux[sorted_idx]
    # # Encontra os limites para cada grupo dentro do array ordenado
    # initial_pointlit_indices = np.cumsum(np.bincount(idx)[:-1])
    # # Divide aux nos grupos correinitial_pointondentes a cada índice único de idx
    # grossa_centroids_array = np.initial_pointlit(sorted_aux, initial_pointlit_indices)
    # # Calcula x, y e z para cada índice
    # indices = np.arange((H // nX) * (V // nY) * (A // nZ))
    # y = indices // (A // nZ) % (V // nY)
    # z = indices % (A // nZ)
    # x = indices // ((A // nZ) * (V // nY))
    # vertices = []

    # # Pontas (Adiciona diretamente sem precisar buscar os valores (são únicos e conhecidos))
    # vertices.extend([0, A-1, (V-1)*A, V*A-1, (H-1)*V*A, (H-1)*V*A + A-1, V*A*H - A, V*A*H - 1])
    # # Intermediários bordas
    # vertices.extend(grossa_centroids_array[i][(nZ + add_z_left) * (nY // 2)] for i in indices[(x == 0) & (y > 0) & (y < (V // nY - 1)) & (z == 0)])
    # vertices.extend(grossa_centroids_array[i][(nZ + add_z_left) * (nY + add_y_left) * (nX // 2)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == 0) & (z == 0)])
    # vertices.extend(grossa_centroids_array[i][(nZ + add_z_left) * (nY + add_y_right) * (nX // 2) + ((nZ + add_z_left) * (nY + add_y_right - 1))] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == (V // nY - 1)) & (z == 0)])
    # vertices.extend(grossa_centroids_array[i][((nZ + add_z_left) * nY * (nX + add_x_right - 1) + ((nZ + add_z_left) * (nY // 2)))] for i in indices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z == 0)])
    # vertices.extend(grossa_centroids_array[i][((nZ + add_z_right) * (nY // 2)) + ((nZ + add_z_right) - 1)] for i in indices[(x == 0) & (y > 0) & (y < (V // nY - 1)) & (z == (A // nZ - 1))])
    # vertices.extend(grossa_centroids_array[i][((nZ + add_z_right) * (nY + add_y_left) * (nX // 2)) + (nZ + add_z_right - 1)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == 0) & (z == (A // nZ - 1))])
    # vertices.extend(grossa_centroids_array[i][((nZ + add_z_right) * (nY + add_y_right) * (nX // 2)) + ((nZ + add_z_right) * (nY + add_y_right) - 1)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == (V // nY - 1)) & (z == (A // nZ - 1))])
    # vertices.extend(grossa_centroids_array[i][((nZ + add_z_right) * nY * (nX + add_x_right - 1)) + ((nZ + add_z_right) * (nY // 2)) + (nZ + add_z_right - 1)] for i in indices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z == (A // nZ - 1))])
    # vertices.extend(grossa_centroids_array[i][nZ // 2] for i in indices[(x == 0) & (y == 0) & (z < (A // nZ - 1)) & (z > 0)])
    # vertices.extend(grossa_centroids_array[i][(nZ * (nY + add_y_right)) - (nZ // 2 + 1)] for i in indices[(x == 0) & (y == (V // nY - 1)) & (z < (A // nZ - 1)) & (z > 0)])
    # vertices.extend(grossa_centroids_array[i][(nZ * (nY + add_y_left) * (nX + add_x_right - 1)) + (nZ // 2)] for i in indices[(x == (H // nX - 1)) & (y == 0) & (z < (A // nZ - 1)) & (z > 0)])
    # vertices.extend(grossa_centroids_array[i][(nZ * (nY + add_y_right) * (nX + add_x_right)) - (nZ // 2 + 1)] for i in indices[(x == (H // nX - 1)) & (y == (V // nY - 1)) & (z < (A // nZ - 1)) & (z > 0)])
    # # Intermediários laterais (Faces)
    # vertices.extend(grossa_centroids_array[i][(nZ * nY) // 2] for i in indices[(x == 0) & (y > 0) & (y < (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))])
    # vertices.extend(grossa_centroids_array[i][(nZ * nY * (nX + add_x_right - 1)) + (nZ * (nY // 2) + (nZ // 2))] for i in indices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))])
    # vertices.extend(grossa_centroids_array[i][(nZ * (nY + add_y_left) * (nX // 2)) + (nZ // 2)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == 0) & (z > 0) & (z < (A // nZ - 1))])
    # vertices.extend(grossa_centroids_array[i][(nZ * (nY + add_y_right) * (nX // 2 + 1)) - (nZ // 2 + 1)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))])
    # vertices.extend(grossa_centroids_array[i][((nZ + add_z_left) * nY * (nX // 2)) + ((nZ + add_z_left) * (nY // 2))] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y < (V // nY - 1)) & (y > 0) & (z == 0)])
    # vertices.extend(grossa_centroids_array[i][((nZ + add_z_right) * nY * (nX // 2)) + ((nZ + add_z_right) * (nY // 2)) + nZ + add_z_right - 1] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y < (V // nY - 1)) & (y > 0) & (z == (A // nZ - 1))])
    # # Intermediários internos
    # vertices.extend(grossa_centroids_array[i][(nZ * nY * nX) // 2] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))])
    # vertices = np.array(vertices)
    # dual[vertices] = 1

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



    # print(len(dual) == H*A*V)
    # print(centroids)
    np.save('Malha_Dual.npy', dual)

def generate_centroids_areas_and_adjacencies(H, V, A):
    block_num=[H, V, A] # numero de blocos
    block_size=[1, 1, 1] # considerado constante aqui (tamanho do bloco)
    initial_point=[0.0, 0.0, 0.0] # considerado constante aqui (ponto inicial)
    ms = np.mgrid[initial_point[0]+0.5*block_size[0]:initial_point[0]+(block_num[0]+0.5)*block_size[0]:block_size[0],
                  initial_point[1]+0.5*block_size[1]:initial_point[1]+(block_num[1]+0.5)*block_size[1]:block_size[1],
                  initial_point[2]+0.5*block_size[2]:initial_point[2]+(block_num[2]+0.5)*block_size[2]:block_size[2]]

    ms = ms.flatten()
    return ms.reshape(3,int(ms.size/3)).T

teste()


# mesh_generation_parameters:
# n_blocks: [60, 220, 1]
# # n_blocks: [120, 120, 1]
# block_size: [1, 1, 1]
# starting_point: [0.0, 0.0, 0.0]

# def get_dual_and_primal_1(centroids):
#     maxs=centroids.max(axis=0)
#     mins=centroids.min(axis=0)
#     engrossamento=np.array(inputs.multiscale_and_multilevel_inputs['multiscale']['coarsening_ratio_1'])
#     n_blocks=np.array(inputs.finescale_inputs['mesh_generation_parameters']['n_blocks'])
#     block_size=np.array(inputs.finescale_inputs['mesh_generation_parameters']['block_size'])
#     starting_point=np.array(inputs.finescale_inputs['mesh_generation_parameters']['starting_point'])
#     n_duals=n_blocks//engrossamento+1
#     second_line = (n_blocks-(n_duals-1)*engrossamento)//2
#     xd=[]
#     for i in range(3):
#         xd.append(np.arange(starting_point[i]+second_line[i]+engrossamento[i],starting_point[i]+(n_duals[i]-1)*engrossamento[i],engrossamento[i]))
#         if len(xd[i])>0:
#             xd[i]=np.unique(np.concatenate([[mins[i], maxs[i]],(xd[i]+0.5)*block_size[i]]))
#         else:
#             xd[i]=np.unique(np.array([mins[i], maxs[i]]))
#     d=np.zeros(len(centroids))

#     for i in range(3):
#         for x in xd[i]:
#             d+=centroids[:,i]==x
#     DUAL_1=d

    # volumes=np.arange(len(centroids))
    # GID_1=-np.ones_like(volumes)
    # xp=np.array(mins[0]-1)
    # xp=np.append(xp,(xd[0][:-1]+xd[0][1:])/2)
    # xp=np.append(xp,maxs[0]+1)

    # yp=np.array(mins[1]-1)
    # yp=np.append(yp,(xd[1][:-1]+xd[1][1:])/2)
    # yp=np.append(yp,maxs[1]+1)

    # zp=np.array(mins[2]-1)
    # zp=np.append(zp,(xd[2][:-1]+xd[2][1:])/2)
    # zp=np.append(zp,maxs[2]+1)
    # count=0
    # for i in range(len(xp)-1):
    #     vx=volumes[(centroids[:,0]>=xp[i]) & (centroids[:,0]<xp[i+1])]
    #     for j in range(len(yp)-1):
    #         vy=vx[(centroids[:,1][vx]>=yp[j]) & (centroids[:,1][vx]<yp[j+1])]
    #         for k in range(len(zp)-1):
    #             vz=vy[(centroids[:,2][vy]>=zp[k]) & (centroids[:,2][vy]<zp[k+1])]
    #             GID_1[vz]=count
    #             count+=1
    # return GID_1, DUAL_1


    # nodes, elements = Read_Mesh('Malhas_teste_3D/semi_9_9_9.msh', 3)
# print('Processando a malha fina')
# teste = MeshDataStructure('','', 3)
# teste.process(nodes, elements)
# teste.divide_mesh(teste.centroids, teste.nHorizontal, teste.nVertical, 3, 3, teste.nAltura, 3)
# print(teste.nHorizontal, teste.nVertical, teste.nAltura)

# print('Salvando malha dual 1')
# dual = np.zeros(len(teste.connectivities), dtype = int)
# for vertice in teste.vertices:
#   dual[vertice] = 3
# for aresta in teste.arestas:
#   for a in aresta:
#     dual[a] = 2
# for face in teste.faces:
#   dual[face] = 1
# np.save('Malha_Dual.npy', dual)
# print('Ok')