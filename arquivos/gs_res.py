import numpy as np
import gstools as gs
import pyvista as pv
import scipy.sparse as sp

def get_grid():
    ordem=[2,1,0]
    nb=np.array([nx, ny, nz])[ordem]+1
    lb=np.array([hx,hy,hz])[ordem]
    sp=np.array([0,0,0])
    grid = pv.ImageData()
    grid.dimensions = nb
    grid.origin = sp  # The bottom left corner of the data set
    grid.spacing = lb  # These are the cell sizes along each axis
    return grid
# Dimensões SPE10
# nx, ny, nz = 60, 220, 85
# Lx, Ly, Lz = 1200.0, 2200.0, 850.0
nx, ny, nz = 60, 220, 10
Lx, Ly, Lz = 1200.0, 2200.0, 850.0
hx, hy, hz = Lx/nx, Ly/ny, Lz/nz

x = np.linspace(hx/2, Lx - hx/2, nx)
y = np.linspace(hy/2, Ly - hy/2, ny)
z = np.linspace(hz/2, Lz - hz/2, nz)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Modelo geoestatístico
model = gs.Exponential(
    dim=3,
    var=1.0,
    len_scale=[300, 300, 3],
    anis=[1.0, 1.0, 0.05],
)

srf = gs.SRF(model, seed=42)
log_k = srf((X, Y, Z))
log_k = 6 * (log_k - np.mean(log_k)) / np.std(log_k)
k = np.exp(log_k)  # mD

# Salvar como .txt
k_flat = k.flatten(order='F')
grid=get_grid()
centroids=np.array([X.flatten(),Y.flatten(),Z.flatten()]).T

# grid.cell_data["ks"] = np.log10(k_flat)
l,c,d=np.load('op_nu_adm.npy')
d1=np.load('DUAL_1.npy')
nid=np.load('NU_ADM_ID.npy')
lv=np.load('levels.npy')
l0,c0,d0=np.load('op_132000_576.npy')
v=0
op_nu_adm=sp.csc_matrix((d, (l.astype(int), c.astype(int))), shape = (int(l.max()+1), int(c.max()+1)), dtype=np.float32)
op_ms=sp.csc_matrix((d0, (l0.astype(int), c0.astype(int))), shape = (int(l0.max()+1), int(c0.max()+1)), dtype=np.float32)
FB=op_nu_adm[:,v].toarray().transpose()[0]
FBms=op_ms[:,v].toarray().transpose()[0]

grid.cell_data["OP_"+str(v)] = FB
grid.cell_data["OPms_"+str(v)] = FBms
grid.cell_data["DUAL"] = d1
grid.cell_data["levels"] = lv

l,c,d=np.load('new_op.npy')
new_op=sp.csc_matrix((d, (l.astype(int), c.astype(int))), shape = (int(l.max()+1), int(c.max()+1)), dtype=np.float32)
FB=new_op[:,v].toarray().transpose()[0]
grid.cell_data["NEW_OP_"+str(v)] = FB


# grid.cell_data["ks"] = np.log10(k_flat)



grid.save('results/gs_spe10.vtk') #ativar
# import pdb; pdb.set_trace()
# np.savetxt("permeabilidade_spe10_geoestatistica.txt", k_flat, fmt="%.6e")
# import pdb; pdb.set_trace()
