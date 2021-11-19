#-------------------------------------------------------------------------------
# This FEniCS program simulates ice flow over a subglacial lake undergoing water
# volume changes over time in three spatial dimensons
#
# This is the main file that calls the stoke solver and free surface evolution
# functions at each timestep, and saves the results.
#-------------------------------------------------------------------------------
import os
import sys
from dolfin import *
from stokes import stokes_solve
from geometry import bed,bed_2D
from mesh_fcns import move_mesh
from boundary_conds import mark_boundary
from params import tol,Lngth,Hght,nt,dt,Nx,Ny,Nz,save_vtk,plot_now,dim
from realtime_process import realtime_plot
from mpi4py import MPI
#-------------------------------------------------------------------------------

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# 1. make a directory for the results
resultsname = 'results'

if os.path.isdir(resultsname)==False and rank == 0:
    os.mkdir(resultsname)

set_log_level(40)    # suppress Newton convergence information if desired.

# create VTK files for the results for use in paraview, if desired
if save_vtk == 'on':
    vtkfile_u = File(resultsname+'/vtkfiles/u.pvd')
    vtkfile_p = File(resultsname+'/vtkfiles/p.pvd')
    vtkfile_marks = File(resultsname+'/vtkfiles/marks.pvd')

#-------------------------------------------------------------------------------
# 2. create mesh
if dim != '2D':
    p0 = Point((0.0,0.0,0.0))
    p1 = Point((Lngth,Lngth,Hght))
    mesh = BoxMesh(p0,p1, Nx, Ny,Nz)
    # get mesh coordinates
    M = mesh.coordinates()
    # make sure all vertices are bounded below by the bed elevation
    M[:,2][M[:,2]<bed(M[:,0],M[:,1])] = bed(M[:,0],M[:,1])[M[:,2]<bed(M[:,0],M[:,1])]
else:
    p0 = Point((0.0,0.0))
    p1 = Point((Lngth,Hght))
    mesh = RectangleMesh(p0,p1, Nx, Nz)
    M = mesh.coordinates()
    # make sure all vertices are bounded below by the bed elevation
    M[:,1][M[:,1]<bed_2D(M[:,0])] = bed_2D(M[:,0])[M[:,1]<bed_2D(M[:,0])]



#-------------------------------------------------------------------------------
# 3. solve the problem
t = 0     # time

# begin time stepping
for i in range(nt):

    if rank == 0:
        print('-----------------------------------------------')
        print('Timestep '+str(i+1)+' out of '+str(nt))
        sys.stdout.flush()

    # sovle stokes problem for solution w = (velocity,ice pressure,mean water pressure)
    w = stokes_solve(mesh,t)

    # save hdf5 files of solution (w) and mesh
    hdf5 = HDF5File(mesh.mpi_comm(), "results/h5files/sol"+str(i)+".h5", "w")
    hdf5.write(mesh, "mesh")
    hdf5.write(w, "solution")
    hdf5.close()


    # plot in real time if desired
    if plot_now == 'on':
        realtime_plot(w,mesh,i)

    # move the mesh according by solving the surface kinematic equations
    # and displacemnt the interior mesh vertices smoothly
    mesh = move_mesh(w,mesh)

    # # save Stokes solution if desired
    if save_vtk == 'on':
        marks = mark_boundary(mesh)  # boundary markers for different subdomains
        _u, _p,_pw = w.split()
        _u.rename("vel", "U")
        _p.rename("press","P")
        vtkfile_u << (_u,t)
        vtkfile_p << (_p,t)
        vtkfile_marks << (marks,t)

    # update time
    t += dt
