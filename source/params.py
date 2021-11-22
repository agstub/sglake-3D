# All model parameters and options are set here.

import numpy as np
#-------------------------------------------------------------------------------
#-----------------------------MODEL OPTIONS-------------------------------------

# save vtk files for stokes solution if 'on':
save_vtk = 'off'


# real-time plotting
plot_now = 'on'

dim = '2D'

#-------------------------------------------------------------------------------
#-----------------------------MODEL PARAMETERS----------------------------------
# physical units:
# time - seconds
# space - meters
# pressure - pascals
# mass - kg

# material parameters
A0 = 3.1689e-24                    # Glen's law coefficient (ice softness, Pa^{-n}/s)
n = 3.0                            # Glen's law exponent
rm2 = 1 + 1.0/n - 2.0              # exponent in variational forms: r-2
B0 = A0**(-1/n)                    # ice hardness (Pa s^{1/n})
B = (2**((n-1.0)/(2*n)))*B0        # coefficient in weak form (Pa s^{1/n})
rho_i = 917.0                      # density of ice (kg/m^3)
rho_w = 1000.0                     # density of water (kg/m^3)
g = 9.81                           # gravitational acceleration (m/s^2)
C = 5.0e9                          # sliding law friction coefficient (Pa s/m)

# numerical parameters
eps_p = 1.0e-13                    # penalty method parameter for unilateral condition
eps_v = (2*1e13/B)**(1/(rm2/2.0))

sigma_0 = 1e3                      # minimum separation stress:
                                   # water pressure must exceed the normal stress
                                   # by this amount for ice-bed separation to occur
                                   # (has a regularizing effect that leads to better convergence)

quad_degree = 16                   # quadrature degree for weak forms

tol = 1.0e-2                       # numerical tolerance for boundary geometry:
                                   # s(x,t) - b(x) > tol on ice-water boundary,
                                   # s(x,t) - b(x) <= tol on ice-bed boundary.

# geometry/mesh parameters
Hght = 1000.0                      # (initial) height of the domain (m)
Lngth = 80*1000.0                  # length of the domain (m)


Ny = int(Lngth/100.0)             # number of elements in vertical direction
Nx = int(Lngth/100.0)             # number of elements in horizontal direction
Nz = int(Hght/100.0)

DX = Lngth/Nx
DY = Lngth/Ny
DZ = Hght/Nz

# time-stepping parameters
t_period = 5*3.154e7                     # oscillation period (sec)
t_final = 2*t_period                     # final time
nt_per_cycle = 2000                      # number of timesteps per oscillation
nt = int(t_final/t_period*nt_per_cycle)  # number of time steps
dt = t_final/nt                          # timestep size
t_arr = np.linspace(0,t_final,nt)        # time array (mainly for plotting)

# spatial coordinates for plotting and interpolation on fine grid
nx = 100                           # number of grid points for plotting
                                   # (larger than true number of elements Nx)
ny = nx


X_fine = np.linspace(0,Lngth,nx)   # horizontal x coordinate

Y_fine = np.linspace(0,Lngth,ny)   # horizontal y coordinate

X0,Y0 = np.meshgrid(X_fine,Y_fine,indexing='ij')
