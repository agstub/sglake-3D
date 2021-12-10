#1.-------------------------process solution------------------------------------
import numpy as np

ind = np.arange(0,4000,20)

dh = np.zeros((200,101))
wb = np.zeros((200,101))

j = 0

# loop over time steps
for i in ind:

    wb_i = np.loadtxt('./results/arrays/wb_'+str(i))  # basal vertical velocity (units = m/yr)
    dh_i = np.loadtxt('./results/arrays/dh_'+str(i))  # elevation anomaly (units = m)

    wb[j,:] = wb_i
    dh[j,:] = dh_i

    j += 1

p = np.ones(101)

wb_ext = np.outer(wb,p).reshape((200,101,101))
h_ext = np.outer(dh,p).reshape((200,101,101))


# save numpy files for use in inversion
np.save('wb_true.npy',wb_ext)
np.save('h_true.npy',h_ext)


# 2. --------------------- make figure of solution -----------------------------
import numpy as np
import scipy.misc as scm

t_period = 5.0*3.154e7
t_final = 2*t_period
nt = 4000
Lngth = 80*1000.0
Hght = 1000
tol = 1e-2

d0 = 0.1            # Smoothing parameter

def bed_2D(x):
    # generate bed topography
    return -8*np.exp(-((x-Lngth/2.0)**2)/(8000**2) )+4

# Smoothed triangle wave
def trg(t):
    return 1 - 2*np.arccos((1 - d0)*np.sin(2*np.pi*t))/np.pi


# Smooth square wave
def sqr(t):
    return 2*np.arctan(np.sin(2*np.pi*t)/d0)/np.pi

# Smoothed sawtooth wave
def swt(t):
    return (1 + trg((2*t - 1)/4)*sqr(t/2))/2

# Sawtooth volume change time series
def Vol(t,lake_vol_0):
    V = 1.5*lake_vol_0*swt((t-0.11*t_period)/t_period)
    return V

t = np.linspace(0,t_final,nt)
X = np.linspace(0,Lngth,100)
bed = bed_2D(X)
dH = 10
wb_true = wb_ext
wb_inf = np.max(np.abs(wb_true))

V = Vol(t,1)/Vol(t,1)[0]

import matplotlib.pyplot as plt
plt.figure(figsize=(8,8))
ind = np.arange(0,4000,20)
plt.subplot(311)
plt.annotate(r'(a)',xy=(-0.042,1.075),fontsize=18,bbox=dict(facecolor='w',alpha=1))
plt.plot(t/t_final,V,linewidth=3,color='k')
plt.annotate(r'$t_1$',xy=(t[ind][47]/t_final-0.025,V[ind][47]+0.12),fontsize=24)
plt.annotate(r'$t_2$',xy=(t[ind][100]/t_final-0.025,V[ind][100]-0.2),fontsize=24)
plt.annotate(r'$t_3$',xy=(t[ind][113]/t_final+0.008,V[ind][113]+0.08),fontsize=24)
plt.plot(t[ind][47]/t_final,V[ind][47],'o',color='crimson',markersize=12)
plt.plot(t[ind][100]/t_final,V[ind][100],'o',color='crimson',markersize=12)
plt.plot(t[ind][113]/t_final,V[ind][113],'o',color='crimson',markersize=12)
plt.gca().xaxis.tick_top()
plt.gca().xaxis.set_label_position('top')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.ylim(0.25,1.2)
plt.xlim(-0.05,1)
plt.xlabel(r'$t\,/\,T$', fontsize=20)
plt.ylabel(r'$V\,/\,V_0$', fontsize=20)
plt.tight_layout()

plt.subplot(334)
i = ind[47]
plt.title(r'$t_1$',fontsize=24)
plt.annotate(r'(b)',xy=(-38,13.85),fontsize=18,bbox=dict(facecolor='w',alpha=1))

plt.annotate(r'air',xy=(-35,10.5),fontsize=16)
plt.annotate(r'ice',xy=(-35,6),fontsize=16)
plt.annotate(r'bed',xy=(-35,0),fontsize=16)
plt.annotate(r'water',xy=(0,-3),xytext=(13,-4),fontsize=16,arrowprops=dict(facecolor='w', shrink=0.0,headwidth=10,headlength=8,width=3))


dh = np.loadtxt('./results/arrays/dh_'+str(i))  # elevation anomaly (units = m)
ds = np.loadtxt('./results/arrays/ds_'+str(i))+bed  # elevation anomaly (units = m)
plt.plot(X/1000-0.5*Lngth/1000,dh+dH,color='forestgreen',linewidth=2)
p1 = plt.fill_between(X/1000-0.5*Lngth/1000,y1=ds, y2=dh+dH,facecolor='aliceblue',alpha=1.0)
p2 = plt.fill_between(X/1000-0.5*Lngth/1000,bed,ds,facecolor='slateblue',alpha=0.5)
p3 = plt.fill_between(X/1000-0.5*Lngth/1000,-18*np.ones(np.size(X)),bed,facecolor='burlywood',alpha=1.0)
plt.plot(X/1000-0.5*Lngth/1000,bed,color='k',linewidth=2)
plt.plot(X[ds-bed>tol]/1000-0.5*Lngth/1000,ds[ds-bed>tol],'-',color='royalblue',linewidth=2)
plt.gca().xaxis.set_ticklabels([])
plt.yticks([4,dH],[r'$s$',r'$h$'],fontsize=20)
plt.ylim(np.min(bed)-2.0,dH+7,8)
plt.xlim(-0.5*Lngth/1000.0,0.5*Lngth/1000.0)
plt.ylabel('free surfaces\n',fontsize=18)
plt.tight_layout()

plt.subplot(335)
plt.title(r'$t_2$',fontsize=24)
plt.annotate(r'(c)',xy=(-38,13.85),fontsize=18,bbox=dict(facecolor='w',alpha=1))

i = ind[100]
dh = np.loadtxt('./results/arrays/dh_'+str(i))  # elevation anomaly (units = m)
ds = np.loadtxt('./results/arrays/ds_'+str(i))+bed  # elevation anomaly (units = m)
plt.plot(X/1000-0.5*Lngth/1000,dh+dH,color='forestgreen',linewidth=2)
p1 = plt.fill_between(X/1000-0.5*Lngth/1000,y1=ds, y2=dh+dH,facecolor='aliceblue',alpha=1.0)
p2 = plt.fill_between(X/1000-0.5*Lngth/1000,bed,ds,facecolor='slateblue',alpha=0.5)
p3 = plt.fill_between(X/1000-0.5*Lngth/1000,-18*np.ones(np.size(X)),bed,facecolor='burlywood',alpha=1.0)
plt.plot(X/1000-0.5*Lngth/1000,bed,color='k',linewidth=2)
plt.plot(X[ds-bed>tol]/1000-0.5*Lngth/1000,ds[ds-bed>tol],'-',color='royalblue',linewidth=2)
plt.yticks([4,dH],[r'',r''],fontsize=16)
plt.gca().xaxis.set_ticklabels([])
plt.ylim(np.min(bed)-2.0,dH+7,8)
plt.xlim(-0.5*Lngth/1000.0,0.5*Lngth/1000.0)
plt.tight_layout()

plt.subplot(336)
i = ind[113]
plt.title(r'$t_3$',fontsize=24)
plt.annotate(r'(d)',xy=(-38,13.85),fontsize=18,bbox=dict(facecolor='w',alpha=1))
dh = np.loadtxt('./results/arrays/dh_'+str(i))  # elevation anomaly (units = m)
ds = np.loadtxt('./results/arrays/ds_'+str(i))+bed  # elevation anomaly (units = m)
plt.plot(X/1000-0.5*Lngth/1000,dh+dH,color='forestgreen',linewidth=2)
p1 = plt.fill_between(X/1000-0.5*Lngth/1000,y1=ds, y2=dh+dH,facecolor='aliceblue',alpha=1.0)
p2 = plt.fill_between(X/1000-0.5*Lngth/1000,bed,ds,facecolor='slateblue',alpha=0.5)
p3 = plt.fill_between(X/1000-0.5*Lngth/1000,-18*np.ones(np.size(X)),bed,facecolor='burlywood',alpha=1.0)
plt.plot(X/1000-0.5*Lngth/1000,bed,color='k',linewidth=2)
plt.plot(X[ds-bed>tol]/1000-0.5*Lngth/1000,ds[ds-bed>tol],'-',color='royalblue',linewidth=2)
plt.yticks([4,dH],[r'',r''],fontsize=16)
plt.gca().xaxis.set_ticklabels([])
plt.ylim(np.min(bed)-2.0,dH+7,8)
plt.xlim(-0.5*Lngth/1000.0,0.5*Lngth/1000.0)
plt.tight_layout()

plt.subplot(337)
plt.annotate(r'(e)',xy=(-38,0.74),fontsize=18,bbox=dict(facecolor='w',alpha=1))
i = ind[47]
wb = np.loadtxt('./results/arrays/wb_'+str(i))/wb_inf
plt.plot(X/1000-0.5*Lngth/1000,wb,color='k',linewidth=3)
plt.xlabel(r'$x$',fontsize=20)
plt.ylabel(r'$w_b\,/\, \Vert w_b\Vert_\infty$',fontsize=20)
plt.xticks(fontsize=16)
plt.ylim(-1,1)
plt.yticks(fontsize=16)
plt.xlim(-0.5*Lngth/1000.0,0.5*Lngth/1000.0)
plt.tight_layout()

plt.subplot(338)
plt.annotate(r'(f)',xy=(-38,0.74),fontsize=18,bbox=dict(facecolor='w',alpha=1))
i = ind[100]
wb = np.loadtxt('./results/arrays/wb_'+str(i))/wb_inf
plt.plot(X/1000-0.5*Lngth/1000,wb,color='k',linewidth=3)
plt.xlabel(r'$x$',fontsize=20)
plt.xticks(fontsize=16)
plt.ylim(-1,1)
plt.xlim(-0.5*Lngth/1000.0,0.5*Lngth/1000.0)
plt.gca().yaxis.set_ticklabels([])
plt.tight_layout()

plt.subplot(339)
plt.annotate(r'(g)',xy=(-38,0.74),fontsize=18,bbox=dict(facecolor='w',alpha=1))
i = ind[113]
wb = np.loadtxt('./results/arrays/wb_'+str(i))/wb_inf
plt.plot(X/1000-0.5*Lngth/1000,wb,color='k',linewidth=3)
plt.xlabel(r'$x$',fontsize=20)
plt.xticks(fontsize=16)
plt.xlim(-0.5*Lngth/1000.0,0.5*Lngth/1000.0)
plt.ylim(-1,1)
plt.gca().yaxis.set_ticklabels([])
plt.tight_layout()
plt.savefig('fig8',bbox_inches='tight')
plt.close()


#3.---------------------make a movie of the results-----------------------------
import os
if os.path.isdir('movie')==False:
    os.mkdir('movie')    # make a directory for the results.

j=0
for i in ind:
    plt.figure(figsize=(6,6))

    plt.subplot(211)
    dh = np.loadtxt('./results/arrays/dh_'+str(i))  # elevation anomaly (units = m)
    ds = np.loadtxt('./results/arrays/ds_'+str(i))+bed  # elevation anomaly (units = m)
    plt.plot(X/1000-0.5*Lngth/1000,dh+dH,color='forestgreen',linewidth=2)

    plt.annotate(r'air',xy=(-35,10.5),fontsize=16)
    plt.annotate(r'ice',xy=(-35,6),fontsize=16)
    plt.annotate(r'bed',xy=(-35,0),fontsize=16)
    plt.annotate(r'water',xy=(0,-3),xytext=(13,-4),fontsize=16,arrowprops=dict(facecolor='w', shrink=0.0,headwidth=14,headlength=12,width=3))

    p1 = plt.fill_between(X/1000-0.5*Lngth/1000,y1=ds, y2=dh+dH,facecolor='aliceblue',alpha=1.0)
    p2 = plt.fill_between(X/1000-0.5*Lngth/1000,bed,ds,facecolor='slateblue',alpha=0.5)
    p3 = plt.fill_between(X/1000-0.5*Lngth/1000,-18*np.ones(np.size(X)),bed,facecolor='burlywood',alpha=1.0)
    plt.plot(X/1000-0.5*Lngth/1000,bed,color='k',linewidth=2)
    plt.plot(X[ds-bed>tol]/1000-0.5*Lngth/1000,ds[ds-bed>tol],'-',color='royalblue',linewidth=2)
    plt.gca().xaxis.set_ticklabels([])
    plt.yticks([4,dH],[r'$s$',r'$h$'],fontsize=20)
    plt.ylim(np.min(bed)-2.0,dH+7,8)
    plt.xlim(-0.5*Lngth/1000.0,0.5*Lngth/1000.0)
    plt.ylabel('free surfaces\n',fontsize=18)
    plt.tight_layout()

    plt.subplot(212)
    wb = np.loadtxt('./results/arrays/wb_'+str(i))/wb_inf
    plt.plot(X/1000-0.5*Lngth/1000,wb,color='k',linewidth=3)
    plt.xlabel(r'$x$',fontsize=20)
    plt.ylabel(r'$w_b\,/\, \Vert w_b\Vert_\infty$',fontsize=20)
    plt.xticks(fontsize=16)
    plt.ylim(-1,1)
    plt.yticks(fontsize=16)
    plt.xlim(-0.5*Lngth/1000.0,0.5*Lngth/1000.0)
    plt.tight_layout()
    plt.savefig('movie/'+str(j),bbox_inches='tight')
    plt.close()
    j+=1
