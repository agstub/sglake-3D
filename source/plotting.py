from params import Lngth,tol,Hght,nt,X0,Y0,t_arr
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
import os
import numpy as np
import copy
from geometry import bed

def plot_fields(h_i,s_i,wb_i,xh,yh,xs,ys,i):
    if os.path.isdir('results/pngs')==False:
        os.mkdir('results/pngs')    # make a directory for the results.

    levels = np.linspace(-1,1,9)
    levels0 = np.linspace(0,4,9)
    levels0[0] = tol

    h_int = LinearNDInterpolator(list(zip(xh, yh)),h_i)
    s_int = LinearNDInterpolator(list(zip(xs, ys)),s_i)
    wb_int = LinearNDInterpolator(list(zip(xs, ys)),wb_i)

    print(r'solution properties at t='+"{:.2f}".format(t_arr[i]/3.154e7)+' yr:')
    print('max dh = '+str(np.max(np.abs(h_int(X0,Y0)-Hght))))
    print('max ds = '+str(np.max(np.abs(s_int(X0,Y0)-bed(X0,Y0)))))
    print('max wb = '+str(np.max(np.abs(wb_int(X0,Y0)))))
    print('\n')

    levels1 = np.linspace(-5,5,9)

    cmap1 = copy.copy(mpl.cm.get_cmap("Blues"))
    cmap1.set_under('burlywood')

    plt.figure(figsize=(8,10))
    plt.subplot(311)
    plt.title(r'$t=$'+"{:.2f}".format(t_arr[i]/3.154e7)+' yr',loc='left',fontsize=22)
    plt.contourf((X0-Lngth/2)/1e3,(Y0-Lngth/2)/1e3,h_int(X0,Y0)-Hght,levels=levels,vmin=levels[0],vmax=levels[-1],extend='both',cmap='coolwarm')
    plt.yticks(fontsize=16)
    plt.ylabel(r'$y$ (km)',fontsize=20)
    plt.gca().xaxis.set_ticklabels([])
    cbar = plt.colorbar(ticks=levels)
    cbar.set_label(label=r'elevation anomaly (m)',size=18)
    cbar.ax.tick_params(labelsize=16)

    plt.subplot(312)
    p1 = plt.contourf((X0-Lngth/2)/1e3,(Y0-Lngth/2)/1e3,s_int(X0,Y0)-bed(X0,Y0),levels=levels0,vmin=levels0[0],vmax=levels0[-1],extend='both',cmap=cmap1)
    l1 = plt.contour((X0-Lngth/2)/1e3,(Y0-Lngth/2)/1e3,s_int(X0,Y0)-bed(X0,Y0),levels=[tol],colors='k',linewidths=3)
    plt.yticks(fontsize=16)
    plt.gca().xaxis.set_ticklabels([])
    plt.ylabel(r'$y$ (km)',fontsize=20)
    cbar = plt.colorbar(p1,ticks=levels0)
    cbar.set_label(label=r'water thickness (m)',size=18)
    cbar.add_lines(l1)
    cbar.ax.tick_params(labelsize=16)

    plt.subplot(313)
    p1 = plt.contourf((X0-Lngth/2)/1e3,(Y0-Lngth/2)/1e3,wb_int(X0,Y0),levels=levels1,vmin=levels1[0],vmax=levels1[-1],extend='both',cmap='coolwarm')
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    plt.xlabel(r'$x$ (km)',fontsize=20)
    plt.ylabel(r'$y$ (km)',fontsize=20)
    cbar = plt.colorbar(p1,ticks=levels1)
    cbar.set_label(label=r'$w_b$ (m/yr)',size=18)
    cbar.ax.tick_params(labelsize=16)

    plt.tight_layout()
    plt.savefig('results/pngs/surfs_'+str(i))
    plt.close()
