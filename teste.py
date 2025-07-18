import numpy as np
import time
import matplotlib.pyplot as plt

class MeshDataStructure:
    """
    This class is an Data Structure for generate the primal and dual meshes and store the resulting data.
    It contains a method to divide the mesh into primal and dual cells based on the centroids.

    Attributes:
        dual (np.ndarray): An array representing the dual mesh.
        primal_index (np.ndarray): An array representing the indexes of the primal cells.
        primal_size (int): The size of the primal mesh.

        nextH (int): The number of primal cells in the horizontal direction.
        nextV (int): The number of primal cells in the vertical direction.
        nextA (int): The number of primal cells in the axial direction.

    Methods:
        divide_mesh(centroids, nX, nY, nZ): Divides the mesh into primal and dual cells based on the centroids and coarsening rates.
    """
    def __init__(self):
        self.dual = None
        self.primal_index = None
        self.primal_size = None
        self.H = None
        self.V = None
        self.A = None
        self.nextH = None
        self.nextV = None
        self.nextA = None

    def divide_mesh(self, centroids, nX, nY, nZ):
        """
        This method divides the mesh into primal and dual meshes based on the centroids and the coarsening rate for each dimension.
        Args:
            centroids (np.ndarray): An array of shape (N, 3) where N is the number of centroids.
                                    Each row represents a centroid in the format [x, y, z].
            nX (int): Coarsening rate in the horizontal direction.
            nY (int): Coarsening rate in the vertical direction.
            nZ (int): Coarsening rate in the axial direction.
        """
    # Compute the dual mesh
        maxs = centroids.max(axis=0) # [max x, max y, max z]
        mins = centroids.min(axis=0) # [min x, min y, min z]
        cr1 = np.array([nX, nY, nZ])
        n_blocks = np.array(maxs - mins) # max - min = [H, V, A]
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

    # Compute the primal mesh
        primal_index=-np.ones(len(centroids))

        # X primal
        self.H = n_blocks[0] + 1
        if self.H % nX != 0:
            div_x = np.int64(self.H // nX) - 1
            rest_x = self.H - (div_x * nX)
            rest_x_left = np.int64(rest_x // 2)
            rest_x_right = rest_x - rest_x_left
            
            xp = np.array([0, rest_x_left])
            xp = np.append(xp, np.arange(rest_x_left + nX, H + 1, nX))
            xp = np.append(xp, xp[-1] + rest_x_right)
        else:
            xp = np.arange(0, self.H+1, nX)

        # Y primal
        self.V = n_blocks[1] + 1
        if self.V % nY != 0:
            div_y = np.int64(self.V // nY) - 1
            rest_y = self.V - (div_y * nY)
            rest_y_left = np.int64(rest_y // 2)
            rest_y_right = rest_y - rest_y_left
            
            yp = np.array([0, rest_y_left])
            yp = np.append(yp, np.arange(rest_y_left + nY, V + 1, nY))
            yp = np.append(yp, yp[-1] + rest_y_right)
        else:
            yp = np.arange(0, self.V+1, nY)

        # Z primal  
        self.A = n_blocks[2] + 1
        if self.A % nZ != 0:
            div_a = np.int64(self.A // nZ) - 1
            rest_a = self.A - (div_a * nZ)
            rest_a_left = np.int64(rest_a // 2)
            rest_a_right = rest_a - rest_a_left
            
            zp = np.array([0, rest_a_left])
            zp = np.append(zp, np.arange(rest_a_left + nZ, A + 1, nZ))
            zp = np.append(zp, zp[-1] + rest_a_right)
        else:
            zp = np.arange(0, self.A+1, nZ)

        # Maps each fine cell to its corresponding primal cell
        x_idx = np.digitize(centroids[:, 0], xp) - 1
        y_idx = np.digitize(centroids[:, 1], yp) - 1
        z_idx = np.digitize(centroids[:, 2], zp) - 1

        # Compute the linear primal index for each centroid
        primal_index = x_idx * (len(yp) - 1) * (len(zp) - 1) + y_idx * (len(zp) - 1) + z_idx

        self.nextH = (len(xp)-1)
        self.nextV = (len(yp)-1)
        self.nextA = (len(zp)-1)
        self.primal_size = self.nextH * self.nextV * self.nextA
        self.dual = dual
        self.primal_index = primal_index


class MultiLevel_NUADM:
    """
    This class is used to generate and store the multi-level NU-ADM IDs.
    The NUADM IDs are generated based on the primal meshes of each coarsening.

    Attributes:
        nX (list): List of coarsening rates in the horizontal direction for each level.
        nY (list): List of coarsening rates in the vertical direction for each level.
        nZ (list or None): List of coarsening rates in the axial direction for each level, or None if using a 2D mesh.
        levels (int): Number of levels in the multi-level NU-ADM (Ex 3 levels means two coarsenings: Fine (level 0) -> Primal 0 -> Primal 1).
        mesh_level (list): List of MeshDataStructure objects representing the primal and dual meshes at each level.
        mapping (list): List of mappings from fine mesh indexes to its primal equivalent in each level.
        current_index (int): Current index for generating NU-ADM IDs.
        id_NUADM (np.ndarray): NU-ADM IDs for all fine cells.
        fine_mesh (np.ndarray): Array representing the fine mesh.
        len_centroids (int): Length of the centroids array.
    
    Methods:
        generate_centroids(H, V, A): Generates centroids based on the dimensions.
        generate_levels(centroids): Generates the primal and dual meshes for each level based on the centroids.
        generate_ids(level_vector): Generates NU-ADM IDs based on the level vector and the levels generated in generate_levels.
    """
    def __init__(self, nX_all, nY_all, nZ_all = None, niveis = 2):
        self.nX = nX_all
        self.nY = nY_all
        self.nZ = nZ_all
        self.levels = niveis
        self.mesh_level = []
        self.mapping = []
        self.current_index = 0
        self.id_NUADM = None
        self.fine_mesh = None
        self.len_centroids = None

    def generate_centroids(self, H, V, A):
        """
        Generate the centroids based on the dimensions.
        Args:
            H (int): Number of blocks in the horizontal direction.
            V (int): Number of blocks in the vertical direction.
            A (int): Number of blocks in the axial direction.
        Returns:
            np.ndarray: An array of shape (N, 3) where N is the total number of centroids.
                        Each row represents a centroid in the format [x, y, z].
        """
        block_num=[H, V, A] # Number of blocks in each direction
        block_size=[1, 1, 1] # Size of each block in each direction (in this case considered as 1 unit)
        initial_point=[0, 0, 0] # Initial point of the mesh (origin)
        ms = np.mgrid[initial_point[0]+0.5*block_size[0]:initial_point[0]+(block_num[0]+0.5)*block_size[0]:block_size[0],
                    initial_point[1]+0.5*block_size[1]:initial_point[1]+(block_num[1]+0.5)*block_size[1]:block_size[1],
                    initial_point[2]+0.5*block_size[2]:initial_point[2]+(block_num[2]+0.5)*block_size[2]:block_size[2]]

        ms = ms.flatten()
        return ms.reshape(3,int(ms.size/3)).T
    
    # @line_profiler.profile
    def generate_levels(self, centroids):
        """
        Generates the primal and dual meshes for each level based on the centroids.
        Args:
            centroids (np.ndarray): An array of shape (N, 3) where N is the number of centroids.
                                    Each row represents a centroid in the format [x, y, z].
        """
        self.len_centroids = len(centroids)
        self.fine_mesh = np.arange(self.len_centroids)

        # If the first level is 2D (nZ = 1) then nZ must be 1 for all levels
        if self.nZ is None or self.nZ[0] == 1:
            self.nZ = np.ones(self.levels - 1)

        # Initialize first level
        first = MeshDataStructure()
        first.divide_mesh(centroids, self.nX[0], self.nY[0], self.nZ[0])
        self.mesh_level.append(first)
        # Initialize mapping with the level 0 mapping: fine -> primal 0 (already generated in divide_mesh)
        self.mapping.append(self.mesh_level[0].primal_index)

        # Initialize next levels
        for i in range(self.levels - 2):
            nextH = self.mesh_level[i].nextH
            nextV = self.mesh_level[i].nextV
            nextA = self.mesh_level[i].nextA

            # Guarantee that if next level has any dimension equals 1, it coarsening rate vector become 1
            if nextA == 1:
                self.nZ = np.ones(self.levels - 1)
            
            if nextH == 1:
                self.nX = np.ones(self.levels - 1)

            if nextV == 1:
                self.nY = np.ones(self.levels - 1)

            centroids_grossa = self.generate_centroids(nextH, nextV, nextA)
            new = MeshDataStructure()
            new.divide_mesh(centroids_grossa, self.nX[i+1], self.nY[i+1], self.nZ[i+1])
            self.mesh_level.append(new)

            # Saving the mapping of each level: fine -> primal 1, fine -> primal 2, ...
            index = self.mapping[0]
            for j in range(i+1):
                index = self.mesh_level[j+1].primal_index[index]
            self.mapping.append(index)

        print('done generating levels')

    # @line_profiler.profile
    def generate_ids(self, level_vector):
        """
        Generates NU-ADM IDs based on the level vector and the levels generated in generate_levels.
        Args:
            level_vector (np.ndarray): An array where each element represents the level of the corresponding fine mesh cell.
                                       The values should be in the range [0, levels-1].
        """
        # Initializing the NU-ADM IDs vector
        self.id_NUADM = np.full(self.len_centroids, -1)
        ids = [np.full(lvl.primal_size, -1) for lvl in self.mesh_level]
        print('ids initialized')
        
        index_0 = self.fine_mesh[level_vector == 0]
        self.id_NUADM[index_0] = np.arange(len(index_0))
        self.current_index += len(index_0)

        for level in range(1, self.levels):
            index_level = self.fine_mesh[level_vector == level]
            primal_index = np.unique(self.mapping[level - 1][index_level])
            ids[level - 1][primal_index] = np.arange(self.current_index, self.current_index + len(primal_index))
            self.id_NUADM[index_level] = ids[level - 1][self.mapping[level - 1][index_level]]
            self.current_index += len(primal_index)

    def get_ids(self):
        """
        Returns the NU-ADM IDs for all fine cells.
        """
        if self.id_NUADM is None:
            raise ValueError("NU-ADM IDs have not been generated yet. Call run() first.")
        return self.id_NUADM
    
    def get_mapping(self):
        """
        Returns the mapping from fine mesh indexes to primal mesh indexes of each level.
        """
        if not self.mapping:
            raise ValueError("Mappings have not been generated yet. Call run() first.")
        return self.mapping
    
    def get_dual(self):
        """
        Returns a dual mesh vector for each level
        """
        if not self.mesh_level:
            raise ValueError("Mesh levels have not been generated yet. Call run() first.")
        duals = []
        duals.append(self.mesh_level[0].dual)
        for i in range(1, self.levels - 1):
            duals.append(self.mesh_level[i].dual[self.mapping[i-1]])
        return duals

    def get_fine_mesh(self):
        """
        Returns the fine mesh vector.
        """
        if self.fine_mesh is None:
            raise ValueError("Fine mesh has not been generated yet. Call run() first.")
        return self.fine_mesh

    def run(self, centroids, level_vector):
        """
        Runs the generation of levels and IDs.
        Args:
            centroids (np.ndarray): An array of shape (N, 3) where N is the number of centroids.
                                    Each row represents a centroid in the format [x, y, z].
            level_vector (np.ndarray): An array where each element represents the level of the corresponding fine mesh cell.
                                       The values should be in the range [0, levels-1].
        Returns:
            np.ndarray: NU-ADM IDs for all fine cells.
        """
        self.generate_levels(centroids)
        self.generate_ids(level_vector)
        return self.id_NUADM

####################################################################################################
# Daqui pra baixo é código de teste, verificar se preciso manter isso no código final
####################################################################################################

def generate_adj_cent(H, V, A):
    nb=[H, V, A]
    lb=[1, 1, 1]
    sp=[0.0, 0.0, 0.0]
    ms = np.mgrid[sp[0]+0.5*lb[0]:sp[0]+(nb[0]+0.5)*lb[0]:lb[0],
                  sp[1]+0.5*lb[1]:sp[1]+(nb[1]+0.5)*lb[1]:lb[1],
                  sp[2]+0.5*lb[2]:sp[2]+(nb[2]+0.5)*lb[2]:lb[2]]

    ms = ms.flatten()
    centroids = ms.reshape(3,int(ms.size/3)).T
    ############
    # ijk0 = np.array([centroids[:, 0]//lb[0], centroids[:, 1]//lb[1], centroids[:, 2]//lb[2]])
    # ee = ijk0[0] + ijk0[1]*nb[0] + ijk0[2]*nb[0]*nb[1] #60 3600
    # # import pdb; pdb.set_trace()
    # centroids=centroids[ee.astype(int)]
    ############
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

H = 10
V = 10
A = 1
tempos_lvl = []
tempos_ids = []
tamanhos = []
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
sla = MultiLevel_NUADM(nX, nY, nZ, niveis=3) # 5 níveis -> 4 engrossamentos

# @line_profiler.profile
def func_main(centroids, levels, sla):
    init = time.time()
    sla.generate_levels(centroids)
    meio = time.time()
    sla.generate_ids(levels)
    end = time.time()
    return init, meio, end

tempo = func_main(centroids, vector, sla)
print(f'Tempo de processamento: {tempo[2]-tempo[0]:.4f} segundos')
tempos_lvl.append(tempo[1] - tempo[0])
tempos_ids.append(tempo[2] - tempo[1])
tamanhos.append(H*V*A)
    # A = (A/2).__int__()

# Plotando os tempos dos ids e das gerações dos níveis lado a lado de acordo com o tamanho do mesh
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(tamanhos, tempos_lvl, marker='o', label='Geração de Níveis')
plt.xlabel('Tamanho do Mesh')
plt.ylabel('Tempo (s)')
plt.title('Tempo de Geração de Níveis')
plt.grid()
plt.legend()
plt.subplot(1, 2, 2)
plt.plot(tamanhos, tempos_ids, marker='o', label='Geração de IDs')
plt.xlabel('Tamanho do Mesh')
plt.ylabel('Tempo (s)')
plt.title('Tempo de Geração de IDs NU-ADM')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()