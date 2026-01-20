import numpy as np
import time
# import matplotlib.pyplot as plt
from MultiLevel_NUADM import MultiLevel_NUADM

def generate_adj_cent(H, V, A):
    nb=[H, V, A]
    lb=[1, 1, 1]
    sp=[0.0, 0.0, 0.0]
    ms = np.mgrid[sp[0]+0.5*lb[0]:sp[0]+(nb[0]+0.5)*lb[0]:lb[0],
                  sp[1]+0.5*lb[1]:sp[1]+(nb[1]+0.5)*lb[1]:lb[1],
                  sp[2]+0.5*lb[2]:sp[2]+(nb[2]+0.5)*lb[2]:lb[2]]

    ms = ms.flatten()
    centroids = ms.reshape(3,int(ms.size/3)).T
    GID_0 = np.arange(len(centroids))
    adjs0 = np.tile(GID_0,(np.array(nb)>1).sum())
    adjs1 = np.array([])
    sz=1
    sy=1
    if nb[2]>1:
        adjs1 = np.concatenate([adjs1, GID_0+1])
        sz=nb[2]
    if nb[1]>1:
        adjs1 = np.concatenate([adjs1, GID_0+sz])
        sy=nb[1]
    if nb[0]>1:
        adjs1 = np.concatenate([adjs1, GID_0+sz*sy])
    adjs1=adjs1.astype(int)
    adjs=np.vstack([adjs0,adjs1]).T
    adjs=adjs[adjs.max(axis=1)<=GID_0.max()]
    dif=abs(centroids[adjs[:,0]]-centroids[adjs[:,1]])
    adjacencies=adjs[((dif>0).sum(axis=1)==1) & ((dif<=lb).sum(axis=1)==3)]
    areas=np.array([lb[1]*lb[2], lb[0]*lb[2], lb[0]*lb[1]])
    areas=np.tile(areas,len(adjacencies)).reshape(len(adjacencies),3)[abs(centroids[adjacencies[:,0]]-centroids[adjacencies[:,1]])>0]
    return centroids, adjacencies

def label_levels_from_fines(adjs, fines, max_level):
    gids = np.arange(adjs.max() + 1)
    levels = np.full_like(gids, max_level + 1)  # valores default maiores que o nível máximo
    levels[fines] = 0  # nível inicial (fines)

    current_level = 0
    current_front = set(fines)

    while current_level < max_level:
        # Encontra adjacências onde pelo menos um vértice está no nível atual
        mask = np.isin(adjs[:, 0], list(current_front)) | np.isin(adjs[:, 1], list(current_front))
        neighbors = adjs[mask].flatten()

        # Remove os que já foram rotulados com nível menor ou igual
        next_front = np.setdiff1d(neighbors, np.where(levels <= current_level)[0], assume_unique=False)

        # Atribui novo nível
        levels[next_front] = current_level + 1

        # Prepara próxima iteração
        current_front = set(next_front)
        current_level += 1

        if len(current_front) == 0:
            break  # Não há mais vizinhos para expandir

    return levels

H = 81
V = 81
A = 81

# for i in range(4):
print(f'Generating level_vector')
fines = np.random.choice(H*V*A, size=round(H*V*A * 0.3), replace=False) # 30% dos volumes são finos
centroids, adjs = generate_adj_cent(H, V, A)
gids = np.arange(H*V*A)
vector = label_levels_from_fines(adjs, gids[fines],2) # 5 níveis
print('Done')

# Processamento
print(f'Start Processing')
nX = [3, 3, 3, 3]
nY = [3, 3, 3, 3]
nZ = [3, 3, 3, 3]
sla = MultiLevel_NUADM(nX, nY, nZ, levels=3) # 5 níveis -> 4 engrossamentos

init = time.time()
sla.run(centroids, vector)
end = time.time()
# print(f'Tempo de processamento: {end - init:.4f} segundos')
sla.validate()
duals = sla.get_duals()
print(f'Duals generated: {len(duals)}')
for i in range(len(duals)):
    print(f'Level {i} dual shape: {duals[i].shape}')

