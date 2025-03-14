import time
import numpy as np
import line_profiler
import random

class MeshDataStructure:
  def __init__(self):
    self.dual = None
    self.Grossa_Index = None
    self.Grossa_Size = None
    self.nextH = None
    self.nextV = None
    self.nextA = None

  @line_profiler.profile
  def divide_mesh(self, centroids, nX, nY, nZ):

    maxs = centroids.max(axis=0) # [max x, max y, max z]
    mins = centroids.min(axis=0) # [min x, min y, min z] -> max - min = [H, V, A]
    cr1 = np.array([nX, nY, nZ])
    n_blocks = np.array(maxs - mins)
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
    # Gerando os indices Fina->Primal (Vetorizando o código de Cícero)
    # Mapeia cada coordenada dos centroides para um índice de célula usando digitize
    x_idx = np.digitize(centroids[:, 0], xp) - 1
    y_idx = np.digitize(centroids[:, 1], yp) - 1
    z_idx = np.digitize(centroids[:, 2], zp) - 1

    # Calcula índice linear da célula Grossa_Index diretamente
    Grossa_Index = x_idx * (len(yp) - 1) * (len(zp) - 1) + y_idx * (len(zp) - 1) + z_idx

    self.nextH = (len(xp)-1)
    self.nextV = (len(yp)-1)
    self.nextA = (len(zp)-1)
    self.Grossa_Size = self.nextH * self.nextV * self.nextA
    self.dual = dual
    self.Grossa_Index = Grossa_Index

class NUADM_ID:
    def __init__(self, nX_all, nY_all, nZ_all = None, niveis = 2):
        self.nX = nX_all
        self.nY = nY_all
        self.nZ = nZ_all
        self.levels = niveis
        self.mesh_level = []
        self.mapping = []
        self.index_atual = 0
        self.ids = None
        self.id_NUADM = None
        self.get_index_vector = None
        self.len_centroids = None
        
    @line_profiler.profile
    def generate_levels(self, centroids):
        already_ones = False
        self.len_centroids = len(centroids)
        self.get_index_vector = np.arange(self.len_centroids)

        # Caso o primeiro nível seja 2D (nZ[0] = 1), o nZ deve sempre ser 1, e se for diferente de 1, seta para 1
        if self.nZ is None or self.nZ[0] == 1:
            self.nZ = np.ones(self.levels - 1)
            already_ones = True

        # Inicializa o nivel 0
        first = MeshDataStructure()
        first.divide_mesh(centroids, self.nX[0], self.nY[0], self.nZ[0])
        self.mesh_level.append(first)
        # Inicializa o mapeamento de indices de nivel 0: fina -> primal 0 (já existe o mapeamento para o nivel 0)
        self.mapping.append(self.mesh_level[0].Grossa_Index)

        # Inicializa os próximos níveis
        for i in range(self.levels - 2):
            nextH = self.mesh_level[i].nextH
            nextV = self.mesh_level[i].nextV
            nextA = self.mesh_level[i].nextA
            # Garantir que caso o próximo nível seja 2D, o nZ seja 1
            if nextA == 1 and not already_ones:
                self.nZ = np.ones(self.levels - 1)
                already_ones = True

            centroids_grossa = generate_centroids_areas_and_adjacencies(nextH, nextV, nextA)
            new = MeshDataStructure()
            new.divide_mesh(centroids_grossa, self.nX[i+1], self.nY[i+1], self.nZ[i+1])
            self.mesh_level.append(new)

            # criando o mapeamento de indices de níveis > 0: fina -> primal 0 -> primal 1...
            index = self.mapping[0]
            for j in range(i+1):
                index = self.mesh_level[j+1].Grossa_Index[index]
            self.mapping.append(index)

        print('done generating levels')

    @line_profiler.profile
    def generate_ids(self, level_vector):
        # inicializando vetor de ids NU-ADM
        self.id_NUADM = np.full(self.len_centroids, -1)
        self.ids = [np.full(lvl.Grossa_Size, -1) for lvl in self.mesh_level]
        print('ids initialized')
        
        index_0 = self.get_index_vector[level_vector == 0]
        self.id_NUADM[index_0] = np.arange(len(index_0))
        self.index_atual += len(index_0)

        for level in range(1, self.levels):
            index_level = self.get_index_vector[level_vector == level]
            grossa_index = np.unique(self.mapping[level - 1][index_level])
            self.ids[level - 1][grossa_index] = np.arange(self.index_atual, self.index_atual + len(grossa_index))
            self.id_NUADM[index_level] = self.ids[level - 1][self.mapping[level - 1][index_level]]
            self.index_atual += len(grossa_index)

    def validate(self):
        print(f'levels: {self.levels}')
        print(f'level 0: {len(self.get_index_vector)}')
        for i in range(self.levels - 1):
          print(f'level {i+1}: {self.mesh_level[i].Grossa_Size}')


# Usar como exemplo de geração de centroides (Identico ao vetor que vou receber de cícero)
def generate_centroids_areas_and_adjacencies(H, V, A):
    block_num=[H, V, A] # numero de blocos
    block_size=[1, 1, 1] # considerado constante aqui (tamanho do bloco)
    initial_point=[0.0, 0.0, 0.0] # considerado constante aqui (ponto inicial)
    ms = np.mgrid[initial_point[0]+0.5*block_size[0]:initial_point[0]+(block_num[0]+0.5)*block_size[0]:block_size[0],
                  initial_point[1]+0.5*block_size[1]:initial_point[1]+(block_num[1]+0.5)*block_size[1]:block_size[1],
                  initial_point[2]+0.5*block_size[2]:initial_point[2]+(block_num[2]+0.5)*block_size[2]:block_size[2]]

    ms = ms.flatten()
    return ms.reshape(3,int(ms.size/3)).T

# Setando para gerar os centroides, mas vou pegar a partir dos centroides que vou receber de cícero
A = 500
V = 500
H = 500

# Recebo pronto
print('Generating level_vector')
vector = []
for i in range(H * V * A):
    vector.append(random.randint(0, 2))
vector = np.array(vector)
print('Done')

# Recebo pronto
print('Generating centroids')
centroids = generate_centroids_areas_and_adjacencies(H, V, A)
print('Done')

# Recebo pronto (Acho)
nX = [5, 5, 5, 5]
nY = [5, 5, 5, 5]
nZ = [5, 5, 5, 5]

# Processamento
print('Start Processing')
init = time.time()
sla = NUADM_ID(nX, nY, nZ, niveis=5)
sla.generate_levels(centroids)
sla.generate_ids(vector)
fim = time.time()
print(f'Done in {fim-init :.4f} seconds')
sla.validate()
# np.save('nuadm_id.npy', sla.id_NUADM)
# np.save('level.npy', vector)
# np.save('Malha_Dual.npy', sla.mesh_level[0].dual)