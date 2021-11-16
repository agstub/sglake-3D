# This file contains the functions needed for solving the Stokes system.

from params import rho_i,g,tol,B,rm2,rho_w,C,eps_p,eps_v,dt,quad_degree,Lngth,sigma_0,Hght
from boundary_conds import mark_boundary,apply_bcs
from geometry import bed,s_mean0,lake_vol_0
from hydrology import Vdot
import numpy as np
from dolfin import *

def Pi(u,nu):
        # penalty functional for enforcing the impenetrability constraint on the ice-bed boundary.
        # NOTE: not used unless you want to assemble(Pi(u,nu)*ds(3))
        #       to see how much the constraint is being violated
        un = dot(u,nu)
        return 0.5*(un**2.0+un*abs(un))

def dPi(u,nu):
        # derivative of penalty functional for enforcing impenetrability
        # on the ice-bed boundary.
        un = dot(u,nu)
        return un+abs(un)

def eta(u):
        # nonlinear (Glen's law) is viscosity
        return 0.5*B*((inner(sym(grad(u)),sym(grad(u)))+Constant(eps_v))**(rm2/2.0))

def sigma(u,p):
        return -p*Identity(3) + 2*eta(u)*sym(grad(u))

def weak_form(u,p,pw,v,q,qw,f,g_lake,g_cryo,ds,nu,T,lake_vol_0,t):
    # define weak form of the subglacial lake problem

    # measures of the entire lower boundary (L0) and ice-water boundary (L1)
    L0 = Constant(assemble(1*ds(4))+assemble(1*ds(3)))
    L1 = Constant(assemble(1*ds(4)))

    # Nonlinear residual
    Fw =  2*eta(u)*inner(sym(grad(u)),sym(grad(v)))*dx + (- div(v)*p + q*div(u))*dx - inner(f, v)*dx\
         + (g_lake+pw+Constant(rho_w*g*dt)*(dot(u,nu)+Constant(Vdot(lake_vol_0,t)/L1)))*inner(nu, v)*ds(4)\
         + qw*(inner(u,nu)+Constant(Vdot(lake_vol_0,t))/L0)*ds(4)\
         + (g_lake+pw-Constant(sigma_0)+Constant(rho_w*g*dt)*(dot(u,nu)+Constant(Vdot(lake_vol_0,t)/L1)))*inner(nu, v)*ds(3)\
         + qw*(inner(u,nu)+Constant(Vdot(lake_vol_0,t))/L0)*ds(3)\
         + Constant(1/eps_p)*dPi(u,nu)*dot(v,nu)*ds(3)\
         + Constant(C)*inner(dot(T,u),dot(T,v))*ds(3) \
         + g_cryo*inner(nu,v)*ds(2) + g_cryo*inner(nu,v)*ds(1) \
         + g_cryo*inner(nu,v)*ds(5) + g_cryo*inner(nu,v)*ds(6)
    return Fw


def stokes_solve(mesh,t):
        # stokes solver using Taylor-Hood elements and a Lagrange multiplier
        # for the water pressure.

        # define function space
        P1 = FiniteElement('P',mesh.ufl_cell(),1)     # pressure
        P2 = FiniteElement('P',mesh.ufl_cell(),2)     # velocity
        R  = FiniteElement("R", mesh.ufl_cell(),0)    # mean water pressure
        element = MixedElement([[P2,P2,P2],P1,R])
        W = FunctionSpace(mesh,element)

        #---------------------define variational problem------------------------
        w = Function(W)
        (u,p,pw) = split(w)             # (velocity,pressure,mean water pressure)
        (v,q,qw) = TestFunctions(W)     # test functions corresponding to (u,p,pw)

        # Gravitational body force
        f = Constant((0,0,-rho_i*g))

        nu = FacetNormal(mesh)            # Outward-pointing unit normal to the boundary
        I = Identity(3)                   # Identity tensor
        T = I - outer(nu,nu)              # Orthogonal projection (onto boundary)

        # mark the boundary and define a measure for integration
        boundary_markers = mark_boundary(mesh)
        ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)

        z_expr = Expression('x[2]',degree=1)

        s_mean = Constant(assemble(z_expr*ds(4))/assemble(1*ds(4)))

        # Define Neumann (water pressure) condition at ice-water interface
        g_lake = Expression('rho_w*g*(s_mean-x[2])',rho_w=rho_w,g=g,s_mean=s_mean,degree=1)

        # Define Neumann (cryostatic) condition at side-walls
        g_cryo = Expression('rho_i*g*(Hght-x[2])',rho_i=rho_i,g=g,Hght=Hght,degree=1)

        #Apply Dirichlet BC on side walls
        bcs =  apply_bcs(W,boundary_markers)

        # define weak form
        Fw = weak_form(u,p,pw,v,q,qw,f,g_lake,g_cryo,ds,nu,T,lake_vol_0,t)

        # solve for (u,p,pw).
        solve(Fw == 0, w, bcs=bcs,
        solver_parameters={"newton_solver":{"relative_tolerance": 1e-14,"linear_solver":"mumps","maximum_iterations":50}},
        form_compiler_parameters={"quadrature_degree":quad_degree,"optimize":True,"eliminate_zeros":False})

        # return solution w
        return w
