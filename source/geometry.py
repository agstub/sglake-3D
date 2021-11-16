#-------------------------------------------------------------------------------
# This file defines the bed topography, initial ice-water interface, and inital
# lake volume.
# Note: Bed and ice-water interface should be equal on margins of the domain!
#-------------------------------------------------------------------------------

import numpy as np
from params import Lngth,X0,Y0,tol,X_fine,Y_fine
import scipy.integrate as scpint

def bed(x,y):
    # generate bed topography
    B = -8*(np.exp(-((x-Lngth/2.0)**2+(y-Lngth/2.0)**2+np.abs(x+y-Lngth)**2)/(10000**2) ))*(0.8+0.1*np.cos((y-0.5*Lngth)/2000)+0.1*np.sin((x-0.5*Lngth)/2000)+0.05*np.sin((x-y)/1000))+7
    return B

def interface(x,y):
    # generate initial ice-water/ice-bed interface
    Int = np.maximum(0*x,bed(x,y))
    return Int

s_mean0 = 0.0               # initial mean ice-water elevation for default geometry
#s_mean0 = #np.mean(interface(X_fine,Y_fine)[interface(X_fine,Y_fine)-bed(X_fine,Y_fine)>tol])

lake_vol_0 = 2644661.0      # initial lake volume for default geometry
# if geometry is changed, calculate initial lake volume with:
# f0 = lambda x,y: interface(x,y)-bed(x,y)
# lake_vol_0 = scpint.nquad(f0,[[0,Lngth],[0,Lngth]],opts={'limit':100})[0]


# # sanity check plotting
# import matplotlib.pyplot as plt
# plt.figure(figsize=(8,6))
# p1 = plt.contourf(X0/1e3,Y0/1e3,interface(X0,Y0)-bed(X0,Y0),cmap='coolwarm')
# plt.contour(X0/1e3,Y0/1e3,interface(X0,Y0)-bed(X0,Y0),levels=[tol],colors='k',linewidths=3)
# #p1=plt.contourf(X0/1e3,Y0/1e3,bed(X0,Y0),cmap='coolwarm')
# plt.xlabel(r'$x$ (km)',fontsize=20)
# plt.ylabel(r'$y$ (km)',fontsize=20)
# plt.colorbar(p1)
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.tight_layout()
# plt.show()
