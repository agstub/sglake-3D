#------------------------------------------------------------------------------
# This module is used to update the mesh at each timestep by solving the
#  surface kinematic equations
#------------------------------------------------------------------------------

from params import dt,DZ
from dolfin import *
import numpy as np
from geometry import bed
from boundary_conds import mark_boundary

#------------------------------------------------------------------------------
def move_mesh(w,mesh):
    # this function computes the surface displacements and moves the mesh
    # by solving Laplace's equation for a smooth displacement function
    # defined for all mesh vertices

    V = FunctionSpace(mesh, 'CG', 1)

    # define surface elevation functions
    z_expr = Expression('x[2]',degree=1)
    z_fcn = Function(V)
    z_fcn.assign(interpolate(z_expr,V))

    z_x = Function(V)
    z_x.assign(project(Dx(z_fcn,0),V))

    z_y = Function(V)
    z_y.assign(project(Dx(z_fcn,1),V))

    # displacement at upper and lower boundaries
    disp_expr = w.sub(0).sub(2) - w.sub(0).sub(0)*z_x - w.sub(0).sub(1)*z_y

    disp_bdry = project(dt*disp_expr,V)

    boundary_markers = mark_boundary(mesh)

    # define displacement boundary conditions on upper and lower surfaces
    bc1 = DirichletBC(V,disp_bdry,boundary_markers,0)            # top boundary
    bc2 = DirichletBC(V,disp_bdry,boundary_markers,4)            # ice-water boundary
    bc3 = DirichletBC(V,disp_bdry,boundary_markers,3)            # ice-bed boundary

    bcs = [bc1,bc2,bc3]

    # solve Laplace's equation for a smooth displacement field on all vertices,
    # given the boundary displacement disp_bdry
    disp = Function(V)
    v = TestFunction(V)
    F = inner(grad(disp), grad(v))*dx
    solve(F == 0, disp, bcs=bcs)

    # get the displacement at the nodes
    disp_vv = disp.compute_vertex_values(mesh)

    # get mesh coordinates
    M = mesh.coordinates()

    # displacement the mesh vertices with the displacement function
    M[:,2] += disp_vv
    M[:,2][M[:,2]<bed(M[:,0],M[:,1])] = bed(M[:,0],M[:,1])[M[:,2]<bed(M[:,0],M[:,1])]

    return mesh
