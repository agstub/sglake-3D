import numpy as np
from dolfin import *
from params import dim
from mpi4py import MPI
from plotting import plot_fields
from mesh_fcns import get_fields
import sys
from realtime_process import realtime_plot

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

for i in range(nt):
    hdf5 = HDF5File(MPI.COMM_WORLD, "./results/h5files/sol"+str(i)+".h5", "r")
    mesh = Mesh()
    hdf5.read(mesh, "mesh",False)

    P1 = FiniteElement('P',mesh.ufl_cell(),1)     # pressure
    P2 = FiniteElement('P',mesh.ufl_cell(),2)     # velocity
    R  = FiniteElement("R", mesh.ufl_cell(),0)    # mean water pressure

    if dim != '2D':
        element = MixedElement([[P2,P2,P2],P1,R])
    else:
        element = MixedElement([[P2,P2],P1,R])

    W = FunctionSpace(mesh,element)
    w = Function(W)
    hdf5.read(w, "solution")

    realtime_plot(w,mesh,i)
