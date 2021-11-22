from params import Lngth,tol,Hght,nt,X0,Y0,t_arr,dim,X_fine
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator, interp1d
import os
import numpy as np
import copy
from geometry import bed,bed_2D

def plot_fields(h_i,s_i,wb_i,xh,yh,xs,ys,i):
    if dim == '2D':
        plot_2D(h_i,s_i,wb_i,xh,xs,i)
    else:
        if os.path.isdir('results/pngs')==False:
            os.mkdir('results/pngs')    # make a directory for the results.

        levels = np.linspace(-1,1,9)            # contour levels for elevation
        levels0 = np.linspace(0,4,9)            # contour levels for water thickness
        levels0[0] = tol
        levels1 = np.linspace(-2,2,9)           # contour levels for basal vertical vel.

        h_int = LinearNDInterpolator(list(zip(xh, yh)),h_i)
        s_int = LinearNDInterpolator(list(zip(xs, ys)),s_i)
        wb_int = LinearNDInterpolator(list(zip(xs, ys)),wb_i)

        print(r'solution properties at t='+"{:.2f}".format(t_arr[i]/3.154e7)+' yr:')
        print('max elevation anom. = '+str(np.max(np.abs(h_int(X0,Y0)-Hght)))+' m')
        print('max water thick. = '+str(np.max(np.abs(s_int(X0,Y0)-bed(X0,Y0))))+' m')
        print('max basal vertical vel. = '+str(np.max(np.abs(wb_int(X0,Y0))))+" m/yr")
        print('\n')

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





def plot_2D(h_i,s_i,wb_i,xh,xs,i):
    if os.path.isdir('results/pngs')==False:
        os.mkdir('results/pngs')    # make a directory for the results.
    if os.path.isdir('results/arrays')==False:
        os.mkdir('results/arrays')    # make a directory for the results.

    h_int = interp1d(xh,h_i,kind='linear',fill_value='extrapolate',bounds_error=False)
    s_int = interp1d(xs,s_i,kind='linear',fill_value='extrapolate',bounds_error=False)
    wb_int = interp1d(xs,wb_i,kind='linear',fill_value='extrapolate',bounds_error=False)

    X = X_fine
    X_plt = X/1000-0.5*Lngth/1000

    dh = h_int(X)-Hght                      # elevation anaomaly
    ds = s_int(X) - bed_2D(X)               # water layer thickness
    wb = wb_int(X)                  # basal vertical velocity

    print('max |dh| = '+str(np.max(np.abs(dh))))
    print('max |ds| = '+str(np.max(np.abs(ds))))
    print('max |wb| = '+str(np.max(np.abs(wb))))


    np.savetxt('results/arrays/dh_'+str(i),dh)
    np.savetxt('results/arrays/ds_'+str(i),ds)
    np.savetxt('results/arrays/wb_'+str(i),wb)


    plt.figure(figsize=(8,10))
    plt.subplot(311)
    plt.title(r'$t=$'+"{:.2f}".format(t_arr[i]/3.154e7)+' yr',loc='left',fontsize=22)

    # Plot upper surface
    plt.plot(X_plt,dh,color='royalblue',linewidth=3)
    plt.ylabel(r'elevation anomaly (m)',fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylim(-2,1)
    plt.gca().xaxis.set_ticklabels([])
    plt.gca().yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))


    plt.subplot(312)
    plt.plot(X_plt,ds,color='royalblue',linewidth=3)
    plt.ylabel(r'water layer thickness (m)',fontsize=16)
    plt.yticks(fontsize=16)
    plt.gca().xaxis.set_ticklabels([])
    plt.gca().yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
    plt.ylim(-0.1,5)


    plt.subplot(313)
    plt.plot(X_plt,wb,color='royalblue',linewidth=3)
    plt.ylabel(r'basal vertical vel. (m/yr)',fontsize=16)
    plt.ylim(-10,4)
    plt.gca().yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))


    # Label axes and save png:
    plt.xlabel(r'$x$ (km)',fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.savefig('results/pngs/surfs_'+str(i))
    plt.close()
