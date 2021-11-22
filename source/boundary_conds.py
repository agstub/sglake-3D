# -------------------------------------------------------------------------------
# This file contains functions that:
# (1) define the boundaries (ice-air,ice-water,ice-bed) of the mesh, AND...
# (2) mark the boundaries of the mesh
# (3) apply Dirichlet boundary conditions on the side-walls of the domain
#-------------------------------------------------------------------------------
from params import tol,Lngth,Hght,dim
from geometry import bed,bed_2D
import numpy as np
from dolfin import *

#-------------------------------------------------------------------------------
# Define SubDomains for ice-water boundary, ice-bed boundary and side-walls
# of the domain

class WaterBoundary(SubDomain):
    # Ice-water boundary.
    # This boundary is marked first and all of the irrelevant portions are
    # overwritten by the other boundary markers.
    def inside(self, x, on_boundary):
        if dim != '2D':
            return (on_boundary and (x[2]<0.5*Hght))
        else:
            return (on_boundary and (x[1]<0.5*Hght))

class BedBoundary(SubDomain):
    # Ice-bed boundary
    def inside(self, x, on_boundary):
        if dim != '2D':
            return (on_boundary and ((x[2]-bed(x[0],x[1]))<=tol))
        else:
            return (on_boundary and (x[1]-bed_2D(x[0])<=tol))

class WestBoundary(SubDomain):
    # West boundary (x=0)
    def inside(self, x, on_boundary):
        return (on_boundary and np.abs(x[0])<tol)

class EastBoundary(SubDomain):
    # East boundary (x=L)
    def inside(self, x, on_boundary):
        return (on_boundary and np.abs(x[0]-Lngth)<tol)

class SouthBoundary(SubDomain):
    # South boundary (y=0)
    def inside(self, x, on_boundary):
        return (on_boundary and np.abs(x[1])<tol)

class NorthBoundary(SubDomain):
    # North boundary (y=L)
    def inside(self, x, on_boundary):
        return (on_boundary and np.abs(x[1]-Lngth)<tol)

#-------------------------------------------------------------------------------

def mark_boundary(mesh):
    #
    # Assign markers to each boundary segment (except the upper surface).
    # This is used at each time step to update the markers.
    #
    # Boundary marker numbering convention:
    # 0 - "Top" boundary
    # 1 - "East" boundary
    # 2 - "West" boundary
    # 3 - Ice-bed boundary
    # 4 - Ice-water boundary
    # 5 - "North" boundary
    # 6 - "South" boundary
    # This function returns these markers, which are used to define the
    # boundary integrals and dirichlet conditions.
    #
    # *These markers can be saved as pvd file if "save_vtk" is turned on in params.py
    #

    if dim != '2D':
        boundary_markers = MeshFunction('size_t', mesh,dim=2)
    else:
        boundary_markers = MeshFunction('size_t', mesh,dim=1)

    boundary_markers.set_all(0)

    # Mark ice-water boundary
    bdryWater = WaterBoundary()
    bdryWater.mark(boundary_markers, 4)

    # Mark ice-bed boundary away from lake
    bdryBed = BedBoundary()
    bdryBed.mark(boundary_markers, 3)

    # Mark inflow boundary
    bdryWest = WestBoundary()
    bdryWest.mark(boundary_markers, 1)

    # Mark outflow boundary
    bdryEast = EastBoundary()
    bdryEast.mark(boundary_markers, 2)


    if dim != '2D':
        # North boundary
        bdryNorth = NorthBoundary()
        bdryNorth.mark(boundary_markers, 5)

        # South boundary
        bdrySouth = SouthBoundary()
        bdrySouth.mark(boundary_markers, 6)

    return boundary_markers



def apply_bcs(W,boundary_markers):
    # Apply Dirichlet conditions on the side-walls of the domain.
    # here we assume zero vertical velocity, which is consistent with the
    # cryostatic normal stress condition

    if dim != '2D':
        bc1 = DirichletBC(W.sub(0).sub(2), Constant(0), boundary_markers,1)
        bc2 = DirichletBC(W.sub(0).sub(2), Constant(0), boundary_markers,2)
        bc3 = DirichletBC(W.sub(0).sub(2), Constant(0), boundary_markers,5)
        bc4 = DirichletBC(W.sub(0).sub(2), Constant(0), boundary_markers,6)
        bcs = [bc1,bc2,bc3,bc4]
    else:
        bcs = []

    return bcs
