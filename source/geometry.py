#-------------------------------------------------------------------------------
# This file defines the bed topography, initial ice-water interface, and inital
# lake volume.
# Note: Bed and ice-water interface should be equal on margins of the domain!
#-------------------------------------------------------------------------------

import numpy as np
from params import Lngth,X0,Y0,tol,X_fine,Y_fine,dim,tol
import scipy.integrate as scpint

def bed(x,y):
    # generate bed topography
    ## An interesting pattern:
    #B = -8*(np.exp(-((x-Lngth/2.0)**2+(y-Lngth/2.0)**2+np.abs(x+y-Lngth)**2)/(8000**2) ))*(0.8+0.1*np.cos((y-0.5*Lngth)/2000)+0.1*np.sin((x-0.5*Lngth)/2000)+0.05*np.sin((x-y)/1000))+5

    # # Default Gaussian analogous to 2D example
    B = -8*np.exp(-((x-Lngth/2.0)**2 + (y-Lngth/2.0)**2 )/(8000**2) ) + 4

    return B

def bed_2D(x):
    # generate bed topography
    return -8*np.exp(-((x-Lngth/2.0)**2)/(8000**2) )+4

def interface(x,y):
    # generate initial ice-water/ice-bed interface
    Int = np.maximum(0*x,bed(x,y))
    return Int

def interface_2D(x):
    # generate initial ice-water/ice-bed interface
    Int = np.maximum(0*x,bed_2D(x))
    return Int

s_mean0 = 0.0            # initial mean ice-water elevation for default geometry

if dim != '2D':
    f0 = lambda x,y: interface(x,y)-bed(x,y)
    lake_vol_0 = scpint.nquad(f0,[[0,Lngth],[0,Lngth]],opts={'limit':100})[0]
else:
    f0 = lambda x: interface_2D(x)-bed_2D(x)
    lake_vol_0 = scpint.nquad(f0,[[0,Lngth]],opts={'limit':100})[0]



#-------------------------------------------------------------------------------
# # # sanity check plotting
# import matplotlib.pyplot as plt
# plt.figure(figsize=(8,6))
# #p1 = plt.contourf(X0/1e3,Y0/1e3,interface(X0,Y0)-bed(X0,Y0),cmap='coolwarm')
# #plt.contour(X0/1e3,Y0/1e3,interface(X0,Y0)-bed(X0,Y0),levels=[tol],colors='k',linewidths=3)
# p1=plt.contourf(X0/1e3,Y0/1e3,bed(X0,Y0),cmap='coolwarm')
# plt.xlabel(r'$x$ (km)',fontsize=20)
# plt.ylabel(r'$y$ (km)',fontsize=20)
# plt.colorbar(p1)
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.tight_layout()
# plt.show()

# #
# import matplotlib.pyplot as plt
# x = np.linspace(0,Lngth,1000)
# plt.plot(x,bed_2D(x),color='k',linewidth=2)
# plt.plot(x[interface_2D(x)-bed_2D(x)>tol],interface_2D(x)[interface_2D(x)-bed_2D(x)>tol],color='crimson',linewidth=2)
# plt.show()
