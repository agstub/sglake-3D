import numpy as np
from dolfin import *
from params import *
from mpi4py import MPI
from plotting import plot_fields

def get_coords(mesh):

    M = mesh.coordinates()
    xh = M[:,0][np.abs(M[:,2]-Hght)<0.25*DZ]
    yh = M[:,1][np.abs(M[:,2]-Hght)<0.25*DZ]
    xs = M[:,0][np.abs(M[:,2])<0.25*DZ]
    ys = M[:,1][np.abs(M[:,2])<0.25*DZ]

    return [xh,yh,xs,ys]

def get_fields(w,mesh):

    w_vv = w.sub(0).sub(2).compute_vertex_values(mesh)

    M = mesh.coordinates()

    h = M[:,2][np.abs(M[:,2]-Hght)<0.25*DZ]
    s = M[:,2][np.abs(M[:,2])<0.25*DZ]
    wb = w_vv[np.abs(M[:,2])<0.25*DZ]*3.154e7   # save in meters per year

    return [h,s,wb]

for i in range(nt):
    hdf5 = HDF5File(MPI.COMM_WORLD, "./results/h5files/sol"+str(i)+".h5", "r")
    mesh = Mesh()
    hdf5.read(mesh, "mesh",False)

    P1 = FiniteElement('P',mesh.ufl_cell(),1)     # pressure
    P2 = FiniteElement('P',mesh.ufl_cell(),2)     # velocity
    R  = FiniteElement("R", mesh.ufl_cell(),0)    # mean water pressure
    element = MixedElement([[P2,P2,P2],P1,R])
    W = FunctionSpace(mesh,element)
    w = Function(W)
    hdf5.read(w, "solution")

    xh,yh,xs,ys = get_coords(mesh)

    h,s,wb = get_fields(w,mesh)

    plot_fields(h,s,wb,xh,yh,xs,ys,i)
