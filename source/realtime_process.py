import numpy as np
from dolfin import *
from params import dim
from mpi4py import MPI
from plotting import save_and_plot
from mesh_fcns import get_fields
import sys

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def realtime_proc(w,mesh,i):
    if dim != '2D':
        h,s,wb,xh,yh,xs,ys = get_fields(w,mesh)
    else:
        h,s,wb,xh,xs = get_fields(w,mesh)


    h = comm.gather(h,root=0)
    s = comm.gather(s,root=0)
    wb = comm.gather(wb,root=0)

    xh = comm.gather(xh,root=0)
    xs = comm.gather(xs,root=0)


    if dim != '2D':
        yh = comm.gather(yh,root=0)
        ys = comm.gather(ys,root=0)

    if rank ==0:
        h = np.concatenate(h).ravel()
        s = np.concatenate(s).ravel()
        wb = np.concatenate(wb).ravel()

        xh = np.concatenate(xh).ravel()
        xs = np.concatenate(xs).ravel()

        if dim != '2D':
            yh = np.concatenate(yh).ravel()
            ys = np.concatenate(ys).ravel()
            save_and_plot(h,s,wb,xh,yh,xs,ys,i)
            sys.stdout.flush()
        else:
            save_and_plot(h,s,wb,xh,0*xh,xs,0*xs,i)
            sys.stdout.flush()
