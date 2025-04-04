import scipy.sparse as sp
import numpy as np
import matplotlib.pyplot as plt

# Dual, operador, ids -> Entradas

# Criando a matriz do operador da apresentação P^

lines_operator = np.array([1, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 9, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 17, 18, 19, 20, 20, 21, 21, 22, 22, 23, 24, 24, 25, 25, 26, 26, 27]) -1
columns_operator = np.array([0, 0, 1, 0, 1, 0, 1, 1, 1, 2, 1, 2, 1, 2, 2, 0, 0, 1, 0, 1, 0, 1, 1, 1, 2, 1, 2, 1, 2, 2, 0, 0, 1, 0, 1, 0, 1, 1, 1, 2, 1, 2, 1, 2, 2])
values_operator = np.array([1, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 1, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 1, 1, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 1, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 1, 1, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 1, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 1])
operator = sp.csr_matrix((values_operator, (lines_operator, columns_operator)))

vector_ids = np.array([1, 2, 3, 10, 10, 10, 11, 12, 13, 4, 5, 6, 10, 10, 10, 14, 15, 16, 7, 8, 9, 10, 10, 10, 17, 18, 19], dtype=int) -1

vector_level = np.array([0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0], dtype=int)

# vertices de cada grossa
primal_vertex = np.array([9, 13, 17])

# ids_0 = np.where(vector_level == 0)[0]
# ids_1 = np.where(vector_level == 1)[0]

new_lines = []
new_columns = []
new_values = []

for i in range(len(vector_level)):
    if vector_level[i] == 0:
        new_lines.append(i)
        new_columns.append(vector_ids[i])
        new_values.append(1)
    else:
        for j in range(len(operator[i].indices)):
            new_lines.append(i)
            new_columns.append(vector_ids[primal_vertex[operator[i].indices[j]]])
            new_values.append(operator[i].data[j])

# for i in operator:
#     print(i.indices, i.data, i.shape)

P = sp.csr_matrix((new_values, (new_lines, new_columns)))
# x = a.todense()

# LEMBRAR DISSO
x = P.toarray()
plt.matshow(x)
print(P)
plt.show()