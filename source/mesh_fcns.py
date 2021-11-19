#------------------------------------------------------------------------------
# This module is used to update the mesh at each timestep by solving the
#  surface kinematic equations
#------------------------------------------------------------------------------

from params import dt,DZ,Hght,dim
from dolfin import *
import numpy as np
from geometry import bed,bed_2D
from boundary_conds import mark_boundary

#------------------------------------------------------------------------------
def move_mesh(w,mesh):
    # this function computes the surface displacements and moves the mesh
    # by solving Laplace's equation for a smooth displacement function
    # defined for all mesh vertices

    V = FunctionSpace(mesh, 'CG', 1)

    # define surface elevation functions
    if dim != '2D':
        z_expr = Expression('x[2]',degree=1)
        z_fcn = Function(V)
        z_fcn.assign(interpolate(z_expr,V))

        z_x = Function(V)
        z_x.assign(project(Dx(z_fcn,0),V))

        z_y = Function(V)
        z_y.assign(project(Dx(z_fcn,1),V))

        # displacement at upper and lower boundaries
        disp_expr = w.sub(0).sub(2) - w.sub(0).sub(0)*z_x - w.sub(0).sub(1)*z_y
    else:
        z_expr = Expression('x[1]',degree=1)
        z_fcn = Function(V)
        z_fcn.assign(interpolate(z_expr,V))

        z_x = Function(V)
        z_x.assign(project(Dx(z_fcn,0),V))

        # displacement at upper and lower boundaries
        disp_expr = w.sub(0).sub(1) - w.sub(0).sub(0)*z_x


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
    if dim != '2D':
        M[:,2] += disp_vv
        M[:,2][M[:,2]<bed(M[:,0],M[:,1])] = bed(M[:,0],M[:,1])[M[:,2]<bed(M[:,0],M[:,1])]
    else:
        M[:,1] += disp_vv
        M[:,1][M[:,1]<bed_2D(M[:,0])] = bed_2D(M[:,0])[M[:,1]<bed_2D(M[:,0])]

    return mesh

#-------------------------------------------------------------------------------
# get upper surface elevation, lower surface elevation, basal vertical velocity,
# and mesh coordinates
def get_fields(w,mesh):

    M = mesh.coordinates()

    if dim != '2D':
        w_vv = w.sub(0).sub(2).compute_vertex_values(mesh)
        xh = M[:,0][np.abs(M[:,2]-Hght)<0.25*DZ]
        yh = M[:,1][np.abs(M[:,2]-Hght)<0.25*DZ]
        xs = M[:,0][np.abs(M[:,2])<0.25*DZ]
        ys = M[:,1][np.abs(M[:,2])<0.25*DZ]

        h = M[:,2][np.abs(M[:,2]-Hght)<0.25*DZ]
        s = M[:,2][np.abs(M[:,2])<0.25*DZ]
        wb = w_vv[np.abs(M[:,2])<0.25*DZ]*3.154e7   # save in meters per year
        fields = [h,s,wb,xh,yh,xs,ys]
    else:
        w_vv = w.sub(0).sub(1).compute_vertex_values(mesh)
        xh = M[:,0][np.abs(M[:,1]-Hght)<0.25*DZ]
        xs = M[:,0][np.abs(M[:,1])<0.25*DZ]

        h = M[:,1][np.abs(M[:,1]-Hght)<0.25*DZ]
        s = M[:,1][np.abs(M[:,1])<0.25*DZ]
        wb = w_vv[np.abs(M[:,1])<0.25*DZ]*3.154e7   # save in meters per year
        fields = [h,s,wb,xh,xs]

    return fields
