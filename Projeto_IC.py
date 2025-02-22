import time
import numpy as np
import re
import h5py
import meshio
import sys
import line_profiler
import random

def Read_Mesh(filePath, dimension):
    inicio = time.time()
    mesh = meshio.read(filePath)
    fim = time.time()
    print(f"Tempo de leitura da malha: {fim - inicio:.4f} segundos")

    nodes = np.round(mesh.points, 6)
    elem = 2**dimension
    # Inicializa um array para conectividades
    elements = np.empty((0, elem), dtype=int)

    # Itera sobre os tipos de células no arquivo
    for cell in mesh.cells:
        if cell.data.shape[1] == elem: 
            elements = np.concatenate((elements, cell.data)) + 1
    
    return nodes, elements

class MeshDataStructure:
  def __init__(self, nodes_file: str, connectivities_file: str, dimension = None):
    self.dimension = dimension
    self.nodes_file = nodes_file
    self.connectivities_file = connectivities_file
    self.nodes = None
    self.connectivities = None
    self.centroids = None
    self.coords_centroids = None
    self.nHorizontal = None
    self.nVertical = None
    self.nAltura = None
    self.coord_grossa = None
    self.grossa_centroids_array = None
    self.vertices = None
    self.arestas = None
    self.faces = None
    self.internos = None
    self.nos_na_proxima_hierarquia = None
    self.conexoes_na_proxima_hierarquia = None
    self.nX = None
    self.nY = None
    self.nZ = None
    self.Grossa_Index = None

  def set_nodes(self, nodes=None):
    if nodes is None:
      nos = np.array([[1, 2, 3]])
      with h5py.File(self.nodes_file, 'r') as file:
          grupos = []
          for group in file.keys():
              if("entity" in group):
                grupos.append(group)

          for grupo_nome in grupos:
            grupo = file[grupo_nome]
            for dataset in grupo.keys():
                dados = []
                try:
                  dados = grupo[dataset][:]
                  dados = np.round(dados, 6)
                  nos = np.concatenate((nos, dados))
                except:
                  dados = grupo[dataset]
          nos = nos[1:]
          self.nodes = nos
    else:
      self.nodes = nodes

  def set_connectivities(self, connections = None):
    if connections is None:
      conexao = np.array([[1,2,3,4]])
      with h5py.File(self.connectivities_file, 'r') as file:
          grupos = []
          for group in file.keys():
              if("entity" in group):
                grupos.append(group)

          for grupo_nome in grupos:
            grupo = file[grupo_nome]
            for dataset in grupo.keys():
                dados = []
                try:
                  dados = grupo[dataset][:]
                  if(isinstance(dados[0], np.ndarray)):
                    if(len(dados[0]) == 4):
                      conexao = np.concatenate((conexao, dados))
                except:
                  dados = grupo[dataset]
          conexao = conexao[1:]
          self.connectivities = conexao
    else:
      self.connectivities = connections
  
  def get_centroids(self):
    connectivities = np.array(self.connectivities) - 1
    nodes = np.array(self.nodes)

    coords = nodes[connectivities]
    sum_coords = np.sum(coords, axis=1)
    centroids = sum_coords / 4.0

    centroids_rounded = np.round(centroids, 6)
    indices = np.arange(len(centroids)).reshape(-1, 1)
    centroids_with_indices = np.hstack((centroids_rounded, indices))

    self.centroids = centroids_with_indices
    self.coords_centroids = coords

  def sort_centroids(self, centroids):
    centroids = np.asarray(centroids)
    indices_ordenados = np.lexsort((centroids[:, 2], centroids[:, 1], centroids[:, 0]))
    return centroids[indices_ordenados]


  def get_dimensions(self):
    error = 0.0005
    firstX = self.centroids[0, 0]
    firstY = self.centroids[0, 1]
    if self.dimension == 3:
      firstZ = self.centroids[0, 2]

      nVertical = np.sum((self.centroids[:, 0] >= firstX - error) & (self.centroids[:, 0] <= firstX + error) & (self.centroids[:, 2] >= firstZ - error) & (self.centroids[:, 2] <= firstZ + error))
      nHorizontal = np.sum((self.centroids[:, 1] >= firstY - error) & (self.centroids[:, 1] <= firstY + error) & (self.centroids[:, 2] >= firstZ - error) & (self.centroids[:, 2] <= firstZ + error))
      nAltura = np.sum((self.centroids[:, 1] >= firstY - error) & (self.centroids[:, 1] <= firstY + error) & (self.centroids[:, 0] >= firstX - error) & (self.centroids[:, 0] <= firstX + error))
      self.nAltura = nAltura

    else:
      nVertical = np.sum((self.centroids[:, 1] >= firstX - error) & (self.centroids[:, 1] <= firstX + error))
      nHorizontal = np.sum((self.centroids[:, 1] >= firstY - error) & (self.centroids[:, 1] <= firstY + error))

    self.nHorizontal = nHorizontal
    self.nVertical = nVertical

  # NX E NY DEVEM SER ÍMPARES!
  @line_profiler.profile
  def divide_mesh(self, centroids, H, V, nX, nY, A = None, nZ = None):
    
    if (nX % 2 == 0) or (nY % 2 == 0):
      print("Os valores nX ou nY são pares")
      return

    if self.dimension == 3:
      if A is None or nZ is None:
        print("Insira um valor de Altura (A) ou da redução nZ")
        return
      if nZ % 2 == 0:
        print("O valor A não é múltiplo de nZ, ou nZ não é ímpar")
        return

    coord_grossa = np.empty((centroids.shape[0], 3), dtype=object)
    ja_foi_mapeado = set()

    # calcula quantos nós serão adicionados em cada lado caso H % nX != 0, e o mesmo para as outras dimensões
    rest_x = (H % nX)
    add_x_left = rest_x // 2
    add_x_right = rest_x - add_x_left

    rest_y = (V % nY)
    add_y_left = rest_y // 2
    add_y_right = rest_y - add_y_left

    if self.dimension == 3:
      rest_z = (A % nZ)
      add_z_left = rest_z // 2
      add_z_right = rest_z - add_z_left

      grossa_size = (H // nX) * (V // nY) * (A // nZ)

      grossa_centroids_array = np.empty(grossa_size, dtype=object)

      grossa_centroids_array[:] = [np.empty(0, dtype=int)] * grossa_size

      vertices = np.zeros(grossa_size, dtype=int)
      arestas = np.empty(grossa_size, dtype=object)
      faces = np.empty(grossa_size, dtype=object)
      internos = np.empty(grossa_size, dtype=object)

      arestas[:] = [np.empty(0, dtype=int)] * grossa_size

      # Cria uma grade de valores de i, j e k
      i_vals = np.arange(H)
      j_vals = np.arange(V)
      k_vals = np.arange(A)
      i_grid, j_grid, k_grid = np.meshgrid(i_vals, j_vals, k_vals, indexing='ij')

      # Calcula indX, indY e indZ vetorizados
      indX = i_grid // nX
      indY = j_grid // nY
      indZ = k_grid // nZ
      
      # Se o tamanho não for divisível pelo fator de redução, temos que alterar a distribuição dos nós 
      # (Só são alterados se tiver que adicionar na esquerda ou direita)
      # Alterando indX
      if rest_x != 0:
          # Para a parte da esquerda
          indX[nX:nX+add_x_left] = 0
          # Para a parte da direita
          indX[-add_x_right:] = indX[-nX]
          # Para a parte do meio
          indices = np.arange(nX + add_x_left, H - add_x_right, nX)
          for index in indices:
              indX[index + nX - add_x_left : index + nX] = indX[index] 
          
      # Alterando indY
      if rest_y != 0:
          # Para a parte de cima
          indY[:, nY:nY + add_y_left] = 0   
          # Para a parte de baixo
          indY[:, -add_y_right:] = indY[:, [-nY]]      
          # Para a parte do meio
          indices = np.arange(nY + add_y_left, V - add_y_right, nY)
          for index in indices:
              indY[:, index + nY - add_y_left : index + nY] = indY[:, [index]] 

      # Alterando indZ
      if rest_z != 0:
          # Para a parte da frente
          indZ[:, :, nZ:nZ + add_z_left] = 0 
          # Para a parte de trás
          indZ[:, :, -add_z_right:] = indZ[:, :, [-nZ]]      
          # Para a parte do meio
          indices = np.arange(nZ + add_z_left, A - add_z_right, nZ)
          for index in indices:
              indZ[:, :, index + nZ - add_z_left : index + nZ] = indZ[:, :, [index]]

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
      self.Grossa_Index = idx
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

      # Calcula os vértices para cada condição
      # Pontas
      vertices[(x == 0) & (y == 0) & (z == 0)] = [grossa_centroids_array[i][0] for i in indices[(x == 0) & (y == 0) & (z == 0)]]
      vertices[(x == (H // nX - 1)) & (y == 0) & (z == 0)] = [grossa_centroids_array[i][((nZ + add_z_left) * (nY + add_y_left) * (nX + add_x_right - 1))] for i in indices[(x == (H // nX - 1)) & (y == 0) & (z == 0)]]
      vertices[(x == 0) & (y == (V // nY - 1)) & (z == 0)] = [grossa_centroids_array[i][((nZ + add_z_left) * (nY + add_y_right - 1))] for i in indices[(x == 0) & (y == (V // nY - 1)) & (z == 0)]]
      vertices[(x == (H // nX - 1)) & (y == (V // nY - 1)) & (z == 0)] = [grossa_centroids_array[i][((nZ + add_z_left) * (nY + add_y_right) * (nX + add_x_right - 1)) + ((nZ + add_z_left) * (nY + add_y_right - 1))] for i in indices[(x == (H // nX - 1)) & (y == (V // nY - 1)) & (z == 0)]]
      vertices[(x == 0) & (y == 0) & (z == (A // nZ - 1))] = [grossa_centroids_array[i][(nZ + add_z_right) - 1] for i in indices[(x == 0) & (y == 0) & (z == (A // nZ - 1))]]
      vertices[(x == (H // nX - 1)) & (y == 0) & (z == (A // nZ - 1))] = [grossa_centroids_array[i][((nY + add_y_left) * (nZ + add_z_right) * (nX + add_x_right - 1)) + (nZ + add_z_right - 1)] for i in indices[(x == (H // nX - 1)) & (y == 0) & (z == (A // nZ - 1))]]
      vertices[(x == 0) & (y == (V // nY - 1)) & (z == (A // nZ - 1))] = [grossa_centroids_array[i][(nZ + add_z_right) * (nY + add_y_right) - 1] for i in indices[(x == 0) & (y == (V // nY - 1)) & (z == (A // nZ - 1))]]
      vertices[(x == (H // nX - 1)) & (y == (V // nY - 1)) & (z == (A // nZ - 1))] = [grossa_centroids_array[i][-1] for i in indices[(x == (H // nX - 1)) & (y == (V // nY - 1)) & (z == (A // nZ - 1))]]
      # Intermediários bordas
      vertices[(x == 0) & (y > 0) & (y < (V // nY - 1)) & (z == 0)] = [grossa_centroids_array[i][(nZ + add_z_left) * (nY // 2)] for i in indices[(x == 0) & (y > 0) & (y < (V // nY - 1)) & (z == 0)]]
      vertices[(x > 0) & (x < (H // nX - 1)) & (y == 0) & (z == 0)] = [grossa_centroids_array[i][(nZ + add_z_left) * (nY + add_y_left) * (nX // 2)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == 0) & (z == 0)]]
      vertices[(x > 0) & (x < (H // nX - 1)) & (y == (V // nY - 1)) & (z == 0)] = [grossa_centroids_array[i][(nZ + add_z_left) * (nY + add_y_right) * (nX // 2) + ((nZ + add_z_left) * (nY + add_y_right - 1))] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == (V // nY - 1)) & (z == 0)]]
      vertices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z == 0)] = [grossa_centroids_array[i][((nZ + add_z_left) * nY * (nX + add_x_right - 1) + ((nZ + add_z_left) * (nY // 2)))] for i in indices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z == 0)]]
      vertices[(x == 0) & (y > 0) & (y < (V // nY - 1)) & (z == (A // nZ - 1))] = [grossa_centroids_array[i][((nZ + add_z_right) * (nY // 2)) + ((nZ + add_z_right) - 1)] for i in indices[(x == 0) & (y > 0) & (y < (V // nY - 1)) & (z == (A // nZ - 1))]]
      vertices[(x > 0) & (x < (H // nX - 1)) & (y == 0) & (z == (A // nZ - 1))] = [grossa_centroids_array[i][((nZ + add_z_right) * (nY + add_y_left) * (nX // 2)) + (nZ + add_z_right - 1)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == 0) & (z == (A // nZ - 1))]]
      vertices[(x > 0) & (x < (H // nX - 1)) & (y == (V // nY - 1)) & (z == (A // nZ - 1))] = [grossa_centroids_array[i][((nZ + add_z_right) * (nY + add_y_right) * (nX // 2)) + ((nZ + add_z_right) * (nY + add_y_right) - 1)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == (V // nY - 1)) & (z == (A // nZ - 1))]]
      vertices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z == (A // nZ - 1))] = [grossa_centroids_array[i][((nZ + add_z_right) * nY * (nX + add_x_right - 1)) + ((nZ + add_z_right) * (nY // 2)) + (nZ + add_z_right - 1)] for i in indices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z == (A // nZ - 1))]]
      vertices[(x == 0) & (y == 0) & (z < (A // nZ - 1)) & (z > 0)] = [grossa_centroids_array[i][nZ // 2] for i in indices[(x == 0) & (y == 0) & (z < (A // nZ - 1)) & (z > 0)]]
      vertices[(x == 0) & (y == (V // nY - 1)) & (z < (A // nZ - 1)) & (z > 0)] = [grossa_centroids_array[i][(nZ * (nY + add_y_right)) - (nZ // 2 + 1)] for i in indices[(x == 0) & (y == (V // nY - 1)) & (z < (A // nZ - 1)) & (z > 0)]]
      vertices[(x == (H // nX - 1)) & (y == 0) & (z < (A // nZ - 1)) & (z > 0)] = [grossa_centroids_array[i][(nZ * (nY + add_y_left) * (nX + add_x_right - 1)) + (nZ // 2)] for i in indices[(x == (H // nX - 1)) & (y == 0) & (z < (A // nZ - 1)) & (z > 0)]]
      vertices[(x == (H // nX - 1)) & (y == (V // nY - 1)) & (z < (A // nZ - 1)) & (z > 0)] = [grossa_centroids_array[i][(nZ * (nY + add_y_right) * (nX + add_x_right)) - (nZ // 2 + 1)] for i in indices[(x == (H // nX - 1)) & (y == (V // nY - 1)) & (z < (A // nZ - 1)) & (z > 0)]]
      # Intermediários laterais (Faces)
      vertices[(x == 0) & (y > 0) & (y < (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))] = [grossa_centroids_array[i][(nZ * nY) // 2] for i in indices[(x == 0) & (y > 0) & (y < (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))]]
      vertices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))] = [grossa_centroids_array[i][(nZ * nY * (nX + add_x_right - 1)) + (nZ * (nY // 2) + (nZ // 2))] for i in indices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))]]
      vertices[(x > 0) & (x < (H // nX - 1)) & (y == 0) & (z > 0) & (z < (A // nZ - 1))] = [grossa_centroids_array[i][(nZ * (nY + add_y_left) * (nX // 2)) + (nZ // 2)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == 0) & (z > 0) & (z < (A // nZ - 1))]]
      vertices[(x > 0) & (x < (H // nX - 1)) & (y == (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))] = [grossa_centroids_array[i][(nZ * (nY + add_y_right) * (nX // 2 + 1)) - (nZ // 2 + 1)] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y == (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))]]
      vertices[(x > 0) & (x < (H // nX - 1)) & (y < (V // nY - 1)) & (y > 0) & (z == 0)] = [grossa_centroids_array[i][((nZ + add_z_left) * nY * (nX // 2)) + ((nZ + add_z_left) * (nY // 2))] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y < (V // nY - 1)) & (y > 0) & (z == 0)]]
      vertices[(x > 0) & (x < (H // nX - 1)) & (y < (V // nY - 1)) & (y > 0) & (z == (A // nZ - 1))] = [grossa_centroids_array[i][((nZ + add_z_right) * nY * (nX // 2)) + ((nZ + add_z_right) * (nY // 2)) + nZ + add_z_right - 1] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y < (V // nY - 1)) & (y > 0) & (z == (A // nZ - 1))]]
      # Intermediários internos
      vertices[(x > 0) & (x < (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))] = [grossa_centroids_array[i][(nZ * nY * nX) // 2] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y > 0) & (y < (V // nY - 1)) & (z > 0) & (z < (A // nZ - 1))]]
      # Adiciona valores ao conjunto
      ja_foi_mapeado = set(vertices)
      vertices_borda = []

      # Preciso pegar as arestas
      # Borda x = 0 y = 0
      for index in range(0, A):
        if index in ja_foi_mapeado:
          vertices_borda.append(index)
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))

      # Borda x = 0 z = 0
      for index in range(0, V*A, A):
        if index in ja_foi_mapeado:
          vertices_borda.append(index)
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))

      # Borda x = 0 y = V
      for index in range((V - 1)*A, V * A):
        if index in ja_foi_mapeado:
          vertices_borda.append(index)
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))

      # Borda x = 0 z = A
      for index in range(A - 1, V * A, A):
        if index in ja_foi_mapeado:
          vertices_borda.append(index)
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))

      # Borda x = H y = 0
      for index in range(V * A * (H - 1), (V * A * (H - 1)) + A):
        if index in ja_foi_mapeado:
          vertices_borda.append(index)
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))

      # Borda x = H z = 0
      for index in range(V * A * (H - 1), (V * A * (H - 1) + (A * (V - 1)) + 1), A):
        if index in ja_foi_mapeado:
          vertices_borda.append(index)
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))

      # Borda x = H y = V
      for index in range((V * A * (H - 1)) + (A * (V - 1)), V * A * H):
        if index in ja_foi_mapeado:
          vertices_borda.append(index)
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))

      # Borda x = H z = A
      for index in range((V * A * (H - 1)) + (A - 1), V * A * H, A):
        if index in ja_foi_mapeado:
          vertices_borda.append(index)
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))

      # Borda y = 0 z = 0
      for index in range(0, V * A * (H - 1) + 1, V * A):
        if index in ja_foi_mapeado:
          vertices_borda.append(index)
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))

      # Borda y = 0 z = A
      for index in range(A-1, (V * A * (H - 1)) + A, V * A):
        if index in ja_foi_mapeado:
          vertices_borda.append(index)
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))

      # Borda y = V z = 0
      for index in range((V - 1) * A, V * A * (H - 1) + ((V - 1) * A) + 1, V * A):
        if index in ja_foi_mapeado:
          vertices_borda.append(index)
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))

      # Borda y = V z = A
      for index in range((V * A) - 1, V * A * H, V * A):
        if index in ja_foi_mapeado:
          vertices_borda.append(index)
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))

      # Para cada vértice não borda, andamos para os 6 lados até encontrar um vértice, ou chegar nos limites do volume
      todos_vertices_nao_borda = []
      for index in vertices:
        if index not in vertices_borda:
          todos_vertices_nao_borda.append(index)
      todos_vertices_nao_borda = np.array(todos_vertices_nao_borda)

      for index in todos_vertices_nao_borda:
        y = index // A % V
        z = index % A
        x = index // (A * V)

        # Ando para cima até encontrar um vertice.
        if y < V - 1:
          daVezU = index + A
          while (daVezU not in ja_foi_mapeado):
            pertence_na_grossa = coord_grossa[daVezU]
            index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
            arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([daVezU])))
            ja_foi_mapeado.add(daVezU)
            daVezU += A

        # Ando para baixo
        if y > 0:
          daVezD = index - A
          while (daVezD not in ja_foi_mapeado):
            pertence_na_grossa = coord_grossa[daVezD]
            index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
            arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([daVezD])))
            ja_foi_mapeado.add(daVezD)
            daVezD -= A

        # Ando para a direita
        if x < H - 1:
          daVezR = index + V*A
          while (daVezR not in ja_foi_mapeado):
            pertence_na_grossa = coord_grossa[daVezR]
            index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
            arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([daVezR])))
            ja_foi_mapeado.add(daVezR)
            daVezR += V*A

        # Ando para a esquerda
        if x > 0:
          daVezL = index - V*A
          while (daVezL not in ja_foi_mapeado):
            pertence_na_grossa = coord_grossa[daVezL]
            index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
            arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([daVezL])))
            ja_foi_mapeado.add(daVezL)
            daVezL -= V*A

        # Ando para dentro
        if z < A - 1:
          daVezI = index + 1
          while (daVezI not in ja_foi_mapeado):
            pertence_na_grossa = coord_grossa[daVezI]
            index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
            arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([daVezI])))
            ja_foi_mapeado.add(daVezI)
            daVezI += 1

        # Ando para fora
        if z > 0:
          daVezO = index - 1
          while (daVezO not in ja_foi_mapeado):
            pertence_na_grossa = coord_grossa[daVezO]
            index_na_grossa = pertence_na_grossa[0] * ((A // nZ)) * (V // nY) + pertence_na_grossa[1] * (A // nZ) + pertence_na_grossa[2]
            arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([daVezO])))
            ja_foi_mapeado.add(daVezO)
            daVezO -= 1

      # Definir as faces e os internos
      index = np.arange(V * A * H)
      mascara = ~np.isin(index, list(ja_foi_mapeado))
      index = index[mascara]
      y = index // A % V
      z = index % A
      x = index // (A * V)
      is_face =  (x == 0) | (x == H - 1) | (y == 0) | (y == V - 1) | (z == 0) | (z == A - 1)
      face_index = index[is_face]
      internos_index = index[~is_face]
      faces = np.array(face_index)
      internos = np.array(internos_index)
      
      self.coord_grossa = coord_grossa
      self.grossa_centroids_array = grossa_centroids_array
      self.vertices = vertices
      self.arestas = arestas
      self.faces = faces
      self.internos = internos
      self.nX = nX
      self.nY = nY
      self.nZ = nZ

    else:
      grossa_centroids_array = np.empty((H // nX) * (V // nY), dtype=object)

      grossa_centroids_array[:] = [np.array([], dtype=int) for _ in range((H // nX) * (V // nY))]

      vertices = np.zeros((H // nX) * (V // nY), dtype=int)

      arestas = np.empty((H // nX) * (V // nY), dtype=object)
      faces = np.empty((H // nX) * (V // nY), dtype=object)

      arestas[:] = [np.array([], dtype=int) for _ in range((H // nX) * (V // nY))]
      faces[:] = [np.array([], dtype=int) for _ in range((H // nX) * (V // nY))]

      # Cria uma grade de valores de i e j
      i_vals = np.arange(H)
      j_vals = np.arange(V)
      i_grid, j_grid = np.meshgrid(i_vals, j_vals, indexing='ij')

      # Calcula indX e indY vetorizados
      indX = i_grid // nX
      indY = j_grid // nY

      # Alterando indX
      if rest_x != 0:
        for i in range(add_x_left):
          indX[nX + i] = 0
        for i in range(nX + add_x_left, H - add_x_right, nX):
          for j in range(0, add_x_left):
            indX[i + nX - (j + 1)] = indX[i]
        for i in range(add_x_right):
          indX[-(i+1)] = indX[-nX]

      # Alterando indY
      if rest_y != 0:
        for x in range(H):
          for i in range(add_y_left):
            indY[x][nY + i] = 0
          for i in range(nY + add_y_left, V - add_y_right, nY):
            for j in range(add_y_left):
              indY[x][i + nY - (j + 1)] = indY[x][i]
          for i in range(add_y_right):
            indY[x][-(i+1)] = indY[x][-nY]

      # Calcula os nós
      nodes = i_grid * V + j_grid

      # Atribui a coordenada grossa para cada nó
      coord_grossa[nodes.flatten()] = list(zip(indX.flatten(), indY.flatten()))

      # Prepara grossa_centroids_array usando uma lista de compreensão
      grossa_centroids_array = np.empty((H // nX) * (V // nY), dtype=object)
      grossa_centroids_array[:] = [np.array([], dtype=int) for _ in range((H // nX) * (V // nY))]

      self.Grossa_Index = np.arange(H*V)
      # Atualiza grossa_centroids_array com os nós, de forma vetorizada
      for indX_val, indY_val, node in zip(indX.flatten(), indY.flatten(), nodes.flatten()):
          idx = indX_val * (V // nY) + indY_val
          self.Grossa_Index[node] = idx
          grossa_centroids_array[idx] = np.append(grossa_centroids_array[idx], node)

      # Calcula x e y para cada índice
      indices = np.arange((H // nX) * (V // nY))
      x = indices // (V // nY)
      y = indices % (V // nY)

      # Calcula os vértices para cada condição
      vertices[(x == 0) & (y == 0)] = [grossa_centroids_array[i][0] for i in indices[(x == 0) & (y == 0)]]
      vertices[(x == (H // nX - 1)) & (y == 0)] = [grossa_centroids_array[i][-(nY + add_y_left)] for i in indices[(x == (H // nX - 1)) & (y == 0)]]
      vertices[(x == 0) & (y == (V // nY - 1))] = [grossa_centroids_array[i][(nY + add_y_right - 1)] for i in indices[(x == 0) & (y == (V // nY - 1))]]
      vertices[(x == (H // nX - 1)) & (y == (V // nY - 1))] = [grossa_centroids_array[i][-1] for i in indices[(x == (H // nX - 1)) & (y == (V // nY - 1))]]
      vertices[(x == 0) & (y > 0) & (y < (V // nY - 1))] = [grossa_centroids_array[i][nY // 2] for i in indices[(x == 0) & (y > 0) & (y < (V // nY - 1))]]
      vertices[(y == (V // nY - 1)) & (x > 0) & (x < (H // nX - 1))] = [grossa_centroids_array[i][(nY + add_y_right) * (nX // 2 + 1) - 1] for i in indices[(y == (V // nY - 1)) & (x > 0) & (x < (H // nX - 1))]]
      vertices[(y == 0) & (x > 0) & (x < (H // nX - 1))] = [grossa_centroids_array[i][(nY + add_y_left) * (nX // 2)] for i in indices[(y == 0) & (x > 0) & (x < (H // nX - 1))]]
      vertices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1))] = [grossa_centroids_array[i][-(nY // 2 + 1)] for i in indices[(x == (H // nX - 1)) & (y > 0) & (y < (V // nY - 1))]]
      vertices[(x > 0) & (x < (H // nX - 1)) & (y > 0) & (y < (V // nY - 1))] = [grossa_centroids_array[i][(nX * nY) // 2] for i in indices[(x > 0) & (x < (H // nX - 1)) & (y > 0) & (y < (V // nY - 1))]]

      # Adiciona valores ao conjunto
      ja_foi_mapeado = set(vertices)

      # Preciso pegar as arestas
      # Borda lateral esquerda
      for index in range(0, V):
        if index in ja_foi_mapeado:
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0]*((V // nY)) + pertence_na_grossa[1]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))

      # Borda inferior
      for index in range(0, H*V, V):
        if index in ja_foi_mapeado:
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0]*((V // nY)) + pertence_na_grossa[1]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))


      # Borda superior
      for index in range(V-1, H*V, V):
        if index in ja_foi_mapeado:
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0]*((V // nY)) + pertence_na_grossa[1]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))


      # Borda lateral direita
      for index in range(V*(H-1), H*V):
        if index in ja_foi_mapeado:
          continue
        ja_foi_mapeado.add(index)
        pertence_na_grossa = coord_grossa[index]
        index_na_grossa = pertence_na_grossa[0]*((V // nY)) + pertence_na_grossa[1]
        arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([index])))


      #Para cada vertice não-borda, só andar pros 4 lados até atingir um vértice.

      todos_vertices = vertices
      todos_vertices_nao_borda = []

      for index in range(0, (H // nX) * (V // nY)):
        i = index // ((V // nY))
        j = index % ((V // nY))
        if(i == 0 or i == (H//nX - 1) or j == 0 or j == (V // nY - 1)):
          continue
        todos_vertices_nao_borda.append(vertices[index])

      todos_vertices_nao_borda = np.array(todos_vertices_nao_borda)

      for index in todos_vertices_nao_borda:
        # Ando para cima até encontrar um vertice.
        daVezU = index + 1

        while (daVezU not in ja_foi_mapeado):
          pertence_na_grossa = coord_grossa[daVezU]
          index_na_grossa = pertence_na_grossa[0]*((V // nY)) + pertence_na_grossa[1]
          arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([daVezU])))
          ja_foi_mapeado.add(daVezU)
          daVezU += 1

        # Ando para baixo
        daVezD = index - 1
        while (daVezD not in ja_foi_mapeado):
          pertence_na_grossa = coord_grossa[daVezD]
          index_na_grossa = pertence_na_grossa[0]*((V // nY)) + pertence_na_grossa[1]
          arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([daVezD])))
          ja_foi_mapeado.add(daVezD)
          daVezD -= 1

        # Ando para a direita
        daVezR = index + V
        while (daVezR not in ja_foi_mapeado):
          pertence_na_grossa = coord_grossa[daVezR]
          index_na_grossa = pertence_na_grossa[0]*((V // nY)) + pertence_na_grossa[1]
          arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([daVezR])))
          ja_foi_mapeado.add(daVezR)
          daVezR += V

        # Ando para a esquerda
        daVezL = index - V
        while (daVezL not in ja_foi_mapeado):
          pertence_na_grossa = coord_grossa[daVezL]
          index_na_grossa = pertence_na_grossa[0]*((V // nY)) + pertence_na_grossa[1]
          arestas[index_na_grossa] = np.concatenate((arestas[index_na_grossa], np.array([daVezL])))
          ja_foi_mapeado.add(daVezL)
          daVezL -= V

      for index in range(0, V*H):
        if(index not in ja_foi_mapeado):
          ja_foi_mapeado.add(index)
          pertence_na_grossa = coord_grossa[index]
          index_na_grossa = pertence_na_grossa[0]*((V // nY)) + pertence_na_grossa[1]
          faces[index_na_grossa] = np.concatenate((faces[index_na_grossa], np.array([index])))

      self.coord_grossa = coord_grossa
      self.grossa_centroids_array = grossa_centroids_array
      self.vertices = vertices
      self.arestas = arestas
      self.faces = faces
      self.nX = nX
      self.nY = nY
      
  def process(self, nodes = None, connections = None):
    self.set_nodes(nodes)
    self.set_connectivities(connections)
    self.get_centroids()
    self.centroids = self.sort_centroids(self.centroids)
    self.get_dimensions()
    # self.divide_mesh(self.centroids, self.nHorizontal, self.nVertical, nX, nY, self.nAltura, nZ)
    print("Process  Ended")

class NUADM_ID:
    def __init__(self, nX_all, nY_all, nZ_all, niveis = 2, dim = 3):
        self.levels = niveis
        self.mesh_level = []
        self.dimension = dim
        self.nX = nX_all
        self.nY = nY_all
        self.nZ = nZ_all
        self.index_atual = 0
        self.ids = []
        self.id_NUADM = []
        self.mapping = []
        self.get_index_vector = None
        self.len_centroids = 0
        
    @line_profiler.profile
    def generate_levels(self, centroids, H, V, A = None):
        self.len_centroids = len(centroids)
        self.get_index_vector = np.arange(self.len_centroids)
        if self.dimension == 3:
            first = MeshDataStructure('', '', self.dimension)
            first.divide_mesh(centroids, H, V, self.nX[0], self.nY[0], A, self.nZ[0])
            self.mesh_level.append(first)
            nextH = H
            nextV = V
            nextA = A
            # Inicializa o mapeamento de indices de nivel 0: fina -> grossa 0 (já existe o mapeamento para o nivel 0)
            self.mapping.append(self.mesh_level[0].Grossa_Index)
            
            for i in range(self.levels - 2):
                nextH = nextH // self.nX[i]
                nextV = nextV // self.nY[i]
                nextA = nextA // self.nZ[i]
                centroids_grossa = np.arange(nextV * nextA * nextH) # n é o centroide de verdade, mas tem o mesmo numero de centroids
                new = MeshDataStructure('', '', self.dimension)
                new.divide_mesh(centroids_grossa, nextH, nextV, self.nX[i+1], self.nY[i+1], nextA, self.nZ[i+1])
                self.mesh_level.append(new)

                # criando o mapeamento de indices de níveis > 0: fina -> grossa 0 -> grossa 1...
                index = self.mapping[0]
                for j in range(i+1):
                    index = self.mesh_level[j+1].Grossa_Index[index]
                self.mapping.append(index)
                
        if self.dimension == 2:
            first = MeshDataStructure('', '', self.dimension)
            first.divide_mesh(centroids, H, V, self.nX[0], self.nY[0])
            self.mesh_level.append(first)
            nextH = H
            nextV = V
            # Inicializa o mapeamento de indices de nivel 0: fina -> grossa 0 (já existe o mapeamento para o nivel 0)
            self.mapping.append(self.mesh_level[0].Grossa_Index)
            for i in range(self.levels - 2):
                nextH = nextH // self.nX[i]
                nextV = nextV // self.nY[i]
                centroids_grossa = np.arange(nextV * nextH)
                new = MeshDataStructure('', '', self.dimension)
                new.divide_mesh(centroids_grossa, nextH, nextV, self.nX[i+1], self.nY[i+1])
                self.mesh_level.append(new)

                # criando o mapeamento de indices de níveis > 0: fina -> grossa 0 -> grossa 1...
                index = self.mapping[0]
                for j in range(i+1):
                    index = self.mesh_level[j+1].Grossa_Index[index]
                self.mapping.append(index)

        print('done generating levels')

    @line_profiler.profile
    def generate_ids(self, level_vector):
        # inicializando vetor de ids NU-ADM
        self.id_NUADM = np.full(self.len_centroids, -1)
        self.ids = [np.full(len(lvl.grossa_centroids_array), -1) for lvl in self.mesh_level]
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
          print(f'level {i+1}: {len(self.mesh_level[i].grossa_centroids_array)}')



# nodes, elements = Read_Mesh('Malhas_teste_3D/semi_81_81_81.msh', 3)
start = MeshDataStructure('', '', 3)
# start.process(nodes, elements)

A = 100
V = 100
H = 2000

vector = []
for i in range(H * V * A):
    vector.append(random.randint(0, 2))
vector = np.array(vector)

centroids = np.arange(H * V * A)

nX = [5, 5, 5]
nY = [5, 5, 5]
nZ = [5, 5, 5]
sla = NUADM_ID(nX, nY, nZ, niveis=3, dim=3)
sla.generate_levels(centroids, H, V, A)
init = time.time()
sla.generate_ids(vector)
fim = time.time()
print(fim - init)
sla.validate()
np.save('nuadm_id.npy', sla.id_NUADM)
np.save('level.npy', vector)

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