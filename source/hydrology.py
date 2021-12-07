#-------------------------------------------------------------------------------
# This function defines the rate of water volume change in the subglacial lake.
# *Default volume change = smoothed sawtooth timeseries from:
# https://mathematica.stackexchange.com/questions/38293/make-a-differentiable-smooth-sawtooth-waveform
# (of course you can calculate dVol/dt analytically for the default example)
#-------------------------------------------------------------------------------

import numpy as np
import scipy.misc as scm
from params import nt,t_period,t_final

d0 = 0.1            # Smoothing parameter

# Smoothed triangle wave
def trg(t):
    return 1 - 2*np.arccos((1 - d0)*np.sin(2*np.pi*t))/np.pi


# Smooth square wave
def sqr(t):
    return 2*np.arctan(np.sin(2*np.pi*t)/d0)/np.pi

# Smoothed sawtooth wave
def swt(t):
    return (1 + trg((2*t - 1)/4)*sqr(t/2))/2

# Sawtooth volume change
def Vol(t,lake_vol_0):
    V = 1.5*lake_vol_0*swt((t-0.11*t_period)/t_period)
    return V


def Vdot(lake_vol_0,t):
    # compute rate of subglacial lake volume change
    dt_fine = 3.154e7/5000.0       # timestep for computing derivative (1/5000 yr)
    Vd = scm.derivative(Vol,t,dx=dt_fine,args=(lake_vol_0,))
    return Vd


#-------------------------------------------------------------------------------
# #sanity check plotting
import matplotlib.pyplot as plt
ind = np.arange(0,4000,20)
t = np.linspace(0,t_final,nt)
plt.plot(t/3.154e7,Vol(t,1),linewidth=3,color='royalblue')
plt.plot(t[ind]/3.154e7,Vol(t,1)[ind],'ko',markersize=8)
plt.plot(t[ind][47]/3.154e7,Vol(t,1)[ind][47],'ro',markersize=16)
plt.plot(t[ind][100]/3.154e7,Vol(t,1)[ind][100],'ro',markersize=16)
plt.plot(t[ind][113]/3.154e7,Vol(t,1)[ind][113],'ro',markersize=16)
plt.show()
